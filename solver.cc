
#include <algorithm>
#include <numeric>
#include <set>
#include <tuple>

#include <cassert>
#include <cmath>

#include "solver.h"

namespace PPM
{
Solver::Solver(
    int max_iter, int freq, double eps, double ratio, int T, double coef_l2)
    : m_max_iter(max_iter)
    , m_freq(freq)
    , m_eps(eps)
    , m_ratio(ratio)
    , m_T(T)
    , m_has_l2_norm(coef_l2 != 0.0)
    , m_eta(coef_l2){};

vector<Result> Solver::solve(Miner &miner,
                             const int n,
                             const vector<double> &y,
                             const vector<double> &lambda)
{
  if (lambda.empty()) {
    cerr << "RUNTIME ERROR!; lambda parameter empty" << endl;
    exit(1);
  }

  vector<Result> result(m_T);

  double bias = 0.0;
  for (int i = 0; i < n; ++i) {
    bias += y[i];
  }
  bias /= n;

  vector<double> r(n, 0.0);
  for (int i = 0; i < n; ++i) {
    r[i] = 1 - y[i] * bias;
  }

  size_t active_size = 0;
  double l1_norm     = 0.0;
  double l2_norm     = 0.0;
  double maxval      = lambda.front();

  vector<TreeIter> model;
  vector<TreeIter> violate;

  for (int t = 1; t < m_T; ++t) {
    double lam = lambda[t];
    cerr << "------------------------------" << endl
         << "[" << t << "] lambda: " << lam << endl;

    vector<int> index;
    vector<int> norm;
    vector<double> grad;

    /* TRAVERSE *****************************************/
    {
      /* UPDATE LOSS */
      double loss = 0;
      double oTr  = 0;
      {
        for (int i = 0; i < n; ++i) {
          if (r[i] > 0) {
            loss += r[i] * r[i];
            oTr += r[i];
          }
        }
      }
      double primal, alpha, dual;
      if (m_has_l2_norm) {
        primal = 0.5 * loss + lam * (l1_norm + 0.5 * m_eta * l2_norm);
        alpha  = 1 / lam;
        miner.m_checker =
            [&r, &alpha, &y](Node &node) -> std::tuple<bool, bool> {
          bool can_safe_pruning = false;
          bool active           = false;

          node.val = 0.0;
          double p = 0.0, m = 0.0;
          for (int i = 0; i < node.support; ++i) {
            int id = node.x_id[i];
            if (r[id] > 0) {
              double temp = alpha * r[id] * y[id];
              node.val += temp;
              (temp > 0) ? p += temp : m += temp;
            }
          }

          node.val = std::fabs(node.val);

          if (std::max(p, -m) < 1.0) {
            can_safe_pruning = true;
            active           = false;
          } else {
            if (node.val >= 1.0) {
              can_safe_pruning = false;
              active           = true;
            } else {
              can_safe_pruning = false;
              active           = false;
            }
          }

          return {can_safe_pruning, active};
        };

        violate = miner.traverse();

        double violate_term =
            0.5 * lam * (1.0 / m_eta) *
            std::accumulate(violate.begin(),
                            violate.end(),
                            0.0,
                            [](double acc, const auto iter) {
                              return acc +
                                     ((iter->val) - 1.0) * ((iter->val) - 1.0);
                            });
        dual = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr -
               violate_term;
      } else {
        primal = 0.5 * loss + lam * l1_norm;
        alpha  = std::min(std::max(oTr / (lam * loss), 0.0), 1 / maxval);
        dual   = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr;
      }

      double gap    = primal - dual;
      double radius = std::sqrt(2.0 * gap) / lam;

      // for traverse
      miner.m_checker =
          [&r, &alpha, &y, &radius, &n](Node &node) -> std::tuple<bool, bool> {
        bool can_safe_pruning = false;
        bool active           = false;

        node.val = 0.0;
        double p = 0.0, m = 0.0;
        for (int i = 0; i < node.support; ++i) {
          int id = node.x_id[i];
          if (r[id] > 0) {
            double temp = alpha * r[id] * y[id];
            (temp > 0) ? p += temp : m += temp;
            node.val += temp;
          }
        }
        node.val = std::fabs(node.val);

        double s = (static_cast<double>(node.support) / n);
        double u = (s >= 0.5) ? sqrt(n * s * (1 - s))
                              : sqrt(static_cast<double>(n)) / 2;
        if (std::max(p, -m) + radius * u < 1) {
          can_safe_pruning = true;
          active           = false;
        } else {
          can_safe_pruning = false;
          double score     = node.val + radius * sqrt(node.support * (1 - s));
          active           = (score >= 1);
        }

        return {can_safe_pruning, active};
      };

      model = miner.traverse();

      active_size = model.size();

      index.resize(active_size);
      norm.resize(active_size);
      grad.resize(active_size);

      for (size_t j = 0; j < active_size; ++j) {
        index[j] = j;
        norm[j]  = model[j]->x_id.size();
        grad[j]  = 0;
      }
    }
    /*end traverse***************************************/

    /* OPTIMIZE */
    for (int iter = 0; iter <= m_max_iter; iter++) {
      /* SHUFFLE */
      for (size_t j = 0; j < active_size; ++j) {
        size_t i = j + rand() % (active_size - j);
        std::swap(index[i], index[j]);
      }

      /* UPDATE LOSS */
      double loss = 0;
      double oTr  = 0;
      for (int i = 0; i < n; ++i) {
        if (r[i] > 0) {
          loss += r[i] * r[i];
          oTr += r[i];
        }
      }

      /* UPDATE GRADIENT */
      maxval = 0;
      for (size_t s = 0; s < active_size; ++s) {
        size_t j = index[s];
        grad[j]  = 0;

        for (int id : model[j]->x_id) {
          if (r[id] > 0) {
            grad[j] += y[id] * r[id];
          }
        }

        if (fabs(grad[j]) > maxval) maxval = std::fabs(grad[j]);
      }

      /* UPDATE PENALTY */
      l1_norm = 0;
      l2_norm = 0;
      for (size_t s = 0; s < active_size; ++s) {
        size_t j = index[s];
        l1_norm += std::fabs(model[j]->w);
        if (m_has_l2_norm) {
          l2_norm += (model[j]->w) * (model[j]->w);
        }
      }

      /* CHECK CONVERGENCE */
      double primal, alpha, dual;
      if (m_has_l2_norm) {
        primal = 0.5 * loss + lam * (l1_norm + 0.5 * m_eta * l2_norm);
        alpha  = 1.0 / lam;
        // violate
        miner.m_checker =
            [&r, &alpha, &y](Node &node) -> std::tuple<bool, bool> {
          bool can_safe_pruning = false;
          bool active           = false;

          node.val = 0.0;
          double p = 0.0, m = 0.0;
          for (int i = 0; i < node.support; ++i) {
            int id = node.x_id[i];
            if (r[id] > 0) {
              double temp = alpha * r[id] * y[id];
              (temp > 0) ? p += temp : m += temp;
              node.val += temp;
            }
          }
          node.val = std::fabs(node.val);

          if (std::max(p, -m) <= 1.0) {
            can_safe_pruning = true;
            active           = false;
          } else {
            if (node.val >= 1.0) {
              can_safe_pruning = false;
              active           = true;
            } else {
              can_safe_pruning = false;
              active           = false;
            }
          }

          return {can_safe_pruning, active};
        };

        violate = miner.traverse();
        double violate_term =
            0.5 * lam * (1.0 / m_eta) *
            std::accumulate(violate.begin(),
                            violate.end(),
                            0.0,
                            [](double acc, const auto iter) {
                              return acc +
                                     ((iter->val) - 1.0) * ((iter->val) - 1.0);
                            });
        dual = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr -
               violate_term;
      } else {
        primal = 0.5 * loss + lam * l1_norm;
        alpha  = std::min(std::max(oTr / (lam * loss), 0.0), 1 / maxval);
        dual   = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr;
      }

      double gap    = primal - dual;
      double radius = sqrt(2 * gap) / lam;

      /* cerr << "[iter " << iter << "] primal: " << primal << ", dual: " <<
       * dual */
      /*      << ", gap: " << gap / primal << endl; */
      if (gap / primal < m_eps) {
        result[t].bias = bias;
        int active     = 0;
        for (size_t s = 0; s < active_size; ++s) {
          size_t j = index[s];
          if (model[j]->w != 0.0) {
            active++;
            result[t].model.emplace_back(model[j]);
            result[t].weight.emplace_back(model[j]->w);
          }
        }

        cerr << "[iter " << iter << "] primal: " << primal << ", dual: " << dual
             << ", gap: " << gap / primal << ", activeset: " << active << endl;
        cerr << "m_active size = " << active_size << endl;

        break;
      }

      /* UPDATE weight */
      {
        for (size_t s = 0; s < active_size; ++s) {
          int j         = index[s];
          double weight = model[j]->w;

          double loss_old = lam * fabs(weight);
          if (m_has_l2_norm) {
            loss_old += 0.5 * lam * m_eta * (weight) * (weight);
          }

          for (int k = 0; k < norm[j]; ++k) {
            int idx = model[j]->x_id[k];
            if (r[idx] > 0) {
              loss_old += 0.5 * r[idx] * r[idx];
            }
          }

          double G = 0.0;
          double H = 0.0;
          for (int k = 0; k < norm[j]; ++k) {
            int idx = model[j]->x_id[k];
            if (r[idx] > 0) {
              G -= y[idx] * r[idx];
              H += 1.0;
            }
          }

          double Gp = G + lam;
          double Gn = G - lam;
          double d  = 0.0;
          if (m_has_l2_norm) {
            if (Gp <= H * weight)
              d = -(Gp + lam * m_eta * weight) / (H + m_eta * lam);
            else if (Gn >= H * weight)
              d = -(Gn + lam * m_eta * weight) / (H + m_eta * lam);
            else
              d = -weight;
          } else {
            if (Gp <= H * weight)
              d = -Gp / H;
            else if (Gn >= H * weight)
              d = -Gn / H;
            else
              d = -weight;
          }

          double d_old = 0;
          double bound =
              0.5 * (G * d + lam * fabs(weight + d) - lam * fabs(weight));
          if (m_has_l2_norm) {
            bound += 0.5 * lam * weight * d * m_eta;
          }

          for (int linesearch = 0; linesearch < 10; ++linesearch) {
            double loss_new = lam * fabs(weight + d);
            for (int k = 0; k < norm[j]; ++k) {
              int idx = model[j]->x_id[k];
              r[idx] += y[idx] * (d_old - d);
              if (r[idx] > 0) {
                loss_new += 0.5 * r[idx] * r[idx];
              }
            }
            if (loss_new - loss_old <= bound) {
              break;
            } else {
              d_old = d;
              d *= 0.5;
              bound *= 0.5;
            }
          }
          model[j]->w += d;
        }
      }  // end update

      /* UPDATE BIAS */
      {
        int nn     = 0;
        double tmp = 0;
        for (int i = 0; i < n; i++) {
          if (r[i] > 0) {
            tmp += y[i] * r[i] + bias;
            nn++;
          }
        }
        double bias_old = bias;
        bias            = tmp / nn;
        for (int i = 0; i < n; i++) {
          r[i] += y[i] * (bias_old - bias);
        }
      }  // update bias

      // DYNAMIC SAFE SCREENING
      if (iter % m_freq == 0) {
        for (size_t s = 0; s < active_size; ++s) {
          size_t j = index[s];
          double score =
              fabs(alpha * grad[j]) +
              radius * sqrt(norm[j] * (1.0 - static_cast<double>(norm[j]) / n));
          if (score < 1.0) {
            if (model[j]->w != 0.0) {
              for (int k = 0; k < norm[j]; ++k) {
                int id = model[j]->x_id[k];
                r[id] += y[id] * model[j]->w;
              }
            }
            model[j]->w = 0;
            active_size--;
            std::swap(index[s], index[active_size]);
            s--;
          }
        }
      }
    }

  }
  return result;
}

}  // namespace PPM
