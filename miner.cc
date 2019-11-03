#include "miner.h"
#include <algorithm>
#include <cassert>
#include <list>
#include <map>
#include <set>

using std::map;

namespace PPM
{
Miner::Miner(const std::shared_ptr<std::list<Node>> tree,
             const Transaction<Range> &ts,
             const int event_size,
             const int maxpat,
             const vector<int> &num_grid,
             const int min_sup)
    : m_ts(ts)
    , m_event_size(event_size)
    , m_maxpat(maxpat)
    , m_num_grid(num_grid)
    , m_tree(tree)
    , m_min_sup(min_sup){};

void Miner::init_tree()
{
  m_tree->clear();

  map<Range, vector<pair<int, int>>> root_pdb;
  for (size_t i = 0, i_size = m_ts.size(); i < i_size; ++i) {
    for (size_t j = 0, j_size = m_ts[i].size(); j < j_size; ++j) {
      root_pdb[m_ts[i][j]].emplace_back(pair<int, int>(i, j));
    }
  }

  for (auto it = root_pdb.begin(), end = root_pdb.end(); it != end; it++) {
    m_pattern.emplace_back(it->first);

    Node node{m_pattern,             // pattern
              to_string(m_pattern),  // pattern string
              it->second,            // pdb
              vector<int>(),         // x_id
              0.0,                   // weight
              0,                     // support
              0.0,                   // val
              true};                 // is_leaf

    auto current = m_tree->insert(m_tree->end(), node);

    int oid = MAXVAL;
    for (size_t i = 0, i_size = current->pdb.size(); i < i_size; ++i) {
      int id = current->pdb[i].first;
      if (id != oid) {
        current->support++;
        current->x_id.emplace_back(id);
      }
      oid = id;
    }

    m_pattern.pop_back();
  }
}

Node Miner::make_node(const vector<pair<int, int>> &pdb)
{
  Node node{
      m_pattern,             // pattern
      to_string(m_pattern),  // pattern_string
      pdb,                   // pdb
      vector<int>(),         // x_id
      0.0,                   // weight
      0,                     // support
      0.0,                   // val
      false                  // is_leaf
  };

  int oid = MAXVAL;
  for (int i = 0, size = node.pdb.size(); i < size; ++i) {
    int id = node.pdb[i].first;
    if (oid != id) {
      node.support++;
      node.x_id.emplace_back(id);
    }
    oid = id;
  }

  return node;
}

vector<TreeIter> Miner::traverse()
{
  m_active_set.clear();

  vector<TreeIter> leaf;
  for (TreeIter it = m_tree->begin(), end = m_tree->end(); it != end; ++it) {
    auto [can_safe_pruning, append] = m_checker(*it);
    if (!can_safe_pruning) {
      if (append) {
        m_active_set.emplace_back(it);
      }
      if (it->is_leaf) {
        it->is_leaf = false;
        leaf.emplace_back(it);
      }
    }
  }

  for (TreeIter iter : leaf) {
    m_pattern = iter->pattern;
    project(iter->pdb);
  }

  return m_active_set;
}

void Miner::project(const vector<pair<int, int>> &pdb)
{
  if (static_cast<int>(m_pattern.size()) < m_maxpat) {
    map<Range, vector<pair<int, int>>> counter;

    for (const auto [id, j] : pdb) {
      int trsize = m_ts[id].size();
      if (int next_j = j + 1; next_j < trsize) {
        counter[m_ts[id][next_j]].emplace_back(pair<int, int>(id, j));
      }
    }

    for (auto it = counter.begin(), end = counter.end(); it != end; it++) {
      m_pattern.emplace_back(it->first);

      Node node = make_node(it->second);
      if (node.support >= m_min_sup) {
        auto [can_safe_pruning, append] = m_checker(node);

        auto current = m_tree->insert(m_tree->end(), node);
        if (can_safe_pruning) {
          current->is_leaf = true;
        } else {
          if (append) {
            m_active_set.emplace_back(current);
          }
          project(it->second);
        }
      }

      m_pattern.pop_back();
    }
  }
}

}  // namespace PPM
