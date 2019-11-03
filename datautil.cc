#include <algorithm>
#include <fstream>
#include <limits>
#include <random>

#include <cassert>

#include "datautil.h"

namespace PPM::DataUtil
{
std::tuple<vector<int>, Transaction<vector<double>>> read(
    const string file_name, const int event_size)
{
  Transaction<vector<double>> transaction;
  vector<int> label;

  auto split = [](string str, char del) {
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

  int BUFSIZE = 1300000;

  std::ifstream ifs(file_name);
  if (ifs.fail()) {
    std::cerr << "cant open file: file = " << file_name << std::endl;
    exit(1);
  }

  string line(BUFSIZE, '\0');
  while (getline(ifs, line)) {
    /* cerr << "line = " << line << endl; */
    /* int ladfa; */
    /* std::cin >> ladfa; */

    vector<vector<double>> a_ts;

    vector<string> sp_line = split(line, '\t');

    assert(!sp_line.empty());

    for (size_t i = 1, i_size = sp_line.size(); i < i_size; i++) {
      vector<string> sp_sp_line = split(sp_line[i], ':');
      vector<double> event;

      std::transform(sp_sp_line.begin(),
                     sp_sp_line.end(),
                     std::back_inserter(event),
                     [](string str) { return std::stod(str); });
      if (event_size != static_cast<int>(event.size())) {
        cerr << "DATA FORMAT ERROR! Not Consistant Dimension" << endl
             << "data dimension number = " << event.size() << endl
             << "user designed grid size = " << event_size << endl
             << event.size() << endl;
        exit(1);
      }
      a_ts.emplace_back(event);
    }

    transaction.emplace_back(a_ts);
    label.emplace_back(stoi(sp_line.front()));
  }

  if (transaction.empty()) {
    cerr << "INPUT ERROR! transaction is empty!" << endl;
    exit(1);
  }

  if (label.empty()) {
    cerr << "INPUT ERROR! label is empty!" << endl;
    exit(1);
  }

  assert(transaction.size() == label.size());

  return {label, transaction};
}

void binary_label_normalize(vector<int> &alice, vector<int> &bob)
{
  std::set<int> label_set;
  std::for_each(alice.begin(), alice.end(), [&label_set](int a) {
    label_set.emplace(a);
  });
  std::for_each(
      bob.begin(), bob.end(), [&label_set](int a) { label_set.emplace(a); });
  binary_label_normalize(label_set, alice);
  binary_label_normalize(label_set, bob);
}

void binary_label_normalize(vector<int> &label)
{
  std::set<int> label_set;
  std::for_each(label.begin(), label.end(), [&label_set](int a) {
    label_set.emplace(a);
  });
  binary_label_normalize(label_set, label);
}

void binary_label_normalize(std::set<int> label_set, vector<int> &label)
{
  if (int class_number = label_set.size(); class_number > 2) {
    cerr << "DATA FORMAT ERROR! class number > 2" << endl
         << "class number = " << class_number << endl;
    exit(1);
  } else if (class_number < 2) {
    cerr << "DATA FORMAT ERROR! class number < 2" << endl
         << "class number = " << class_number << endl;
    exit(1);
  }

  auto binary_label_normalize = [label_set,
                                 NEGATIVE = Label::NEGATIVE,
                                 POSITIVE = Label::POSITIVE](int &label) {
    if (label == *(label_set.begin())) {
      label = NEGATIVE;
    } else {
      label = POSITIVE;
    }
  };

  for_each(label.begin(), label.end(), binary_label_normalize);
}

Transaction<Range> discretize(const Transaction<vector<double>> &raw,
                              const vector<vector<double>> &grid,
                              const int event_size)
{
  auto discretization = [grid, event_size](
                            const vector<vector<double>> &raw_sequence) {
    auto discretization = [grid, event_size](const vector<double> &raw_event) {
      Range ret(event_size);
      for (int k = 0; k < event_size; ++k) {
        vector<double> bp(grid[k]);
        //先頭にデータの下界を追加
        bp.insert(bp.begin(), -std::numeric_limits<long double>::infinity());
        //末尾にデータの上界を追加
        bp.emplace_back(std::numeric_limits<long double>::infinity());

        auto first    = bp.begin();
        auto last     = bp.end();
        int upper     = distance(first, upper_bound(first, last, raw_event[k]));
        int lower     = upper - 1;
        ret[k].first  = lower;
        ret[k].second = upper;
      }
      return ret;
    };

    vector<Range> ret;
    std::transform(raw_sequence.begin(),
                   raw_sequence.end(),
                   std::back_inserter(ret),
                   discretization);
    return ret;
  };

  Transaction<Range> ret;
  std::transform(
      raw.begin(), raw.end(), std::back_inserter(ret), discretization);

  return ret;
}

}  // namespace PPM::DataUtil
