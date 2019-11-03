#ifndef PPM_MINER_H
#define PPM_MINER_H

#include <functional>
#include <list>
#include <memory>
#include <unordered_map>
#include "type.h"

template <class Key, class Value>
using hash_map = std::unordered_map<Key, Value>;

namespace PPM
{
class Miner
{
public:
  void init_tree();
  vector<TreeIter> traverse();
  std::function<std::tuple<bool, bool>(Node&)> m_checker;

private:
  vector<Range> m_pattern;
  vector<TreeIter> m_active_set;

  const Transaction<Range> m_ts;
  const int m_event_size;
  const int m_maxpat;
  const vector<int> m_num_grid;
  const std::shared_ptr<std::list<Node>> m_tree;
  const int m_min_sup;
  const int MAXVAL = std::numeric_limits<int>::max();

  void project(const vector<pair<int, int>> &pdb);
  Node make_node(const vector<pair<int, int>>& pdb);

public:
  explicit Miner(const std::shared_ptr<std::list<Node>> tree,
                 const Transaction<Range>& ts,
                 const int event_size,
                 const int maxpat,
                 const vector<int>& num_grid,
                 const int min_sup = 1);
};

}  // namespace PPM

#endif
