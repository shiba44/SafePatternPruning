#ifndef PPM_SOLVER_H
#define PPM_SOLVER_H

#include <list>
#include "miner.h"

namespace PPM
{
class Solver
{
public:
  Solver(int max_iter, int freq, double eps, double ratio, int T, double reg);

  vector<Result> solve(Miner& miner,
                       const int n,
                       const vector<double>& y,
                       const vector<double>& lambda);

private:
  const int m_max_iter;
  const int m_freq;
  const double m_eps;
  const double m_ratio;
  const int m_T;
  const bool m_has_l2_norm;
  const double m_eta;
};
}  // namespace PPM

#endif
