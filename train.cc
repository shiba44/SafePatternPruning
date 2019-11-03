
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <cassert>
#include <cmath>

#include "cmdline.h"
#include "datautil.h"
#include "miner.h"
#include "solver.h"
#include "testutil.h"
#include "type.h"

using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::string;
using std::vector;

int main(int argc, char** argv)
{
  const auto p = [argv, argc]() {
    cmdline::parser ret;

    ret.add<string>("train", 't', "file path for train data", true);
    ret.add<string>("test", 'T', "file path for test data", true);
    ret.add<string>("out", 'o', "file path for output", true);
    ret.add<string>("grid", 'g', "grid", true);
    ret.add<int>("maxpat", 'L', "max pattern", true);
    ret.add<double>("lambda_max", 'a', "lambda_max", false);
    ret.add<int>("num_lambda", 'm', "num of lambda in reg path", false, 100);
    ret.add<int>("max_iter", 'i', "max iteration", false, 10000000);
    ret.add<int>("freq", 'f', "frequency of screening", false, 50);
    ret.add<double>("eps",
                    'e',
                    "(primal - dual) / primal < eps => convergence",
                    false,
                    1e-6);
    ret.add<double>("ratio", 'r', "ratio", false, 2.0);
    ret.add<double>("l2co", 'z', "coefficient of l2 norm", false, 0.0);
    ret.parse_check(argc, argv);

    return ret;
  }();

  const string train_file_path = p.get<string>("train");
  const string test_file_path  = p.get<string>("test");
  const string out_file_path   = p.get<string>("out");

  std::ofstream fout(out_file_path, std::ios::app);
  if (!fout.is_open()) {
    cerr << "RUNTIME ERROR! Can Not Open: " << out_file_path << endl;
    exit(1);
  }

  const string grid_str = p.get<string>("grid");

  const int maxpat     = p.get<int>("maxpat");
  const int max_iter   = p.get<int>("max_iter");
  const int freq       = p.get<int>("freq");
  const double eps     = p.get<double>("eps");
  const double ratio   = p.get<double>("ratio");
  const int T          = p.get<int>("num_lambda");
  const double coef_l2 = p.get<double>("l2co");

  const auto grid = [&grid_str]() {
    auto split = [](string str, char del) -> vector<string> {
      vector<string> ret;
      std::stringstream ss(str);
      string item;
      while (getline(ss, item, del)) {
        if (!item.empty()) {
          ret.push_back(item);
        }
      }
      return ret;
    };

    vector<vector<double>> ret;
    for (const auto sub_str : split(grid_str, ';')) {
      vector<double> g;
      for (const auto s_sub_str : split(sub_str, ',')) {
        g.emplace_back(std::stod(s_sub_str));
      }
      ret.emplace_back(g);
    }

    std::sort(ret.begin(), ret.end());

    return ret;
  }();

  const int event_size = static_cast<int>(grid.size());
  /* proceasure the following **********************************************/
  auto [label, raw_transaction] =
      PPM::DataUtil::read(train_file_path, event_size);

  auto [test_label, test_raw_transaction] =
      PPM::DataUtil::read(test_file_path, event_size);

  PPM::DataUtil::binary_label_normalize(label, test_label);

  const vector<double> y(label.begin(), label.end());

  const PPM::Transaction<PPM::Range> ts =
      PPM::DataUtil::discretize(raw_transaction, grid, event_size);

  const int n = ts.size();

  const PPM::Transaction<PPM::Range> test_ts =
      PPM::DataUtil::discretize(test_raw_transaction, grid, event_size);

  const auto tree = std::make_shared<std::list<PPM::Node>>();

  const auto num_grid = [&grid]() {
    vector<int> ret;
    for (const vector<double>& g : grid) {
      ret.emplace_back(g.size());
    }
    return ret;
  }();

  PPM::Miner miner(tree, ts, event_size, maxpat, num_grid);
  miner.init_tree();

  const auto lambda_max = [&p, n, &y, &miner]() -> double {
    if (p.exist("lambda_max")) {
      return p.get<double>("lambda_max");
    } else {
      double ret = 0.0;

      double bias = 0.0;
      for (int i = 0; i < n; ++i) {
        bias += y[i];
      }
      bias /= n;

      vector<double> r(n, 0.0);
      for (int i = 0; i < n; ++i) {
        r[i] = 1 - y[i] * bias;
      }

      miner.m_checker =
          [&r, &y, &ret](PPM::Node& node) -> std::tuple<bool, bool> {
        bool can_safe_pruning = false;

        node.val = 0.0;
        double p = 0.0, m = 0.0;
        for (int i = 0; i < node.support; ++i) {
          int id = node.x_id[i];
          if (r[id] > 0) {
            double temp = r[id] * y[id];
            node.val += temp;
            (temp > 0) ? p += temp : m += temp;
          }
        }

        node.val = std::fabs(node.val);

        if (std::max(p, -m) < ret) {
          can_safe_pruning = true;
        } else {
          can_safe_pruning = false;
          if (node.val > ret) {
            ret = node.val;
          }
        }

        return {can_safe_pruning, false};
      };

      miner.traverse();
      return ret;
    }
  }();

  // regularization path
  const auto lambda = [T, lambda_max, ratio]() {
    vector<double> ret(T);
    for (int t = 0; t < T; ++t) {
      ret[t] = lambda_max * std::pow(10, -ratio * t / (T - 1));
    }
    return ret;
  }();

  PPM::Solver solver(max_iter, freq, eps, ratio, T, coef_l2);
  const vector<PPM::Result> result = solver.solve(miner, n, y, lambda);

  for (int t = 1; t < T; t++) {
    assert(result[t].model.size() == result[t].weight.size());

    const auto predict = [&test_ts, &a_result = result[t]]() {
      vector<int> ret;
      for (const auto& tts : test_ts) {
        const auto value = [&a_result, &tts]() {
          double ret = a_result.bias;
          for (size_t j = 0, j_size = a_result.model.size(); j < j_size; j++) {
            if (PPM::is_specific_and_abstract(tts,
                                              a_result.model[j]->pattern)) {
              ret += a_result.model[j]->w;
            }
          }
          return ret;
        }();

        if (value == 0.0) {
          ret.emplace_back(PPM::Label::UNDECIDABLE);
        } else if (value > 0.0) {
          ret.emplace_back(PPM::Label::POSITIVE);
        } else {
          ret.emplace_back(PPM::Label::NEGATIVE);
        }
      }
      return ret;
    }();

    assert(predict.size() == test_label.size());

    const double measure = PPM::TestUtil::macro_f1_score(predict, test_label);

    auto join = [](const vector<int>& v, const char* del) {
      std::stringstream ss;
      for (size_t i = 0, i_size = v.size(); i < i_size; i++) {
        if (i != 0) {
          ss << del;
        }
        ss << v[i];
      }
      return ss.str();
    };

    std::stringstream ss;
    ss << maxpat << "\t"  // maxpat
       << t << "\t"       // t番目の正則化係数
       << join(num_grid, ",") << "\t" << measure << endl;

    fout << ss.str();
  }

  /* cerr << "tree size = " << tree->size() << endl; */
  /* cerr << "lambda max = " << lambda_max << endl; */
  /* std::copy(label.begin(), label.end(), std::ostream_iterator<int>(std::cerr,
   * ",")); */

  /* for(auto a_ts : ts){ */
  /*   cerr << "a_ts " << PPM::to_string(a_ts); */
  /* } */

  /* for_each(tree->begin(),tree->end(),[](PPM::Node node){ */
  /*     cerr << node.pattern_str << "," << node.support << endl; */
  /* }); */
}
