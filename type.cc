
#include "type.h"
#include <algorithm>

namespace PPM
{
bool is_specific_and_abstract(const vector<Range>& alice,
                              const vector<Range>& bob)
{
  if (alice.size() < bob.size()) {
    return false;
  }

  auto result =
      std::search(alice.begin(),
                  alice.end(),
                  bob.begin(),
                  bob.end(),
                  [](const Range& a, const Range& b) { return include(b, a); });

  if (result != alice.end()) {
    return true;
  }

  return false;
}

string to_string(const Range& range)
{
  auto pair2str = [](const pair<int, int>& a) {
    std::stringstream ss;
    ss << "[" << a.first << "," << a.second << "]";
    return ss.str();
  };

  std::stringstream ss;
  for (size_t k = 0, k_size = range.size(); k < k_size; ++k) {
    ss << pair2str(range[k]);
  }
  return ss.str();
}

string to_string(const vector<Range>& pat)
{
  std::stringstream ss;
  for (size_t i = 0, patsize = pat.size(); i < patsize; ++i) {
    if (i != 0) {
      ss << ",";
    }
    ss << to_string(pat[i]);
  }
  return ss.str();
}

string to_string(const vector<double>& value)
{
  std::stringstream ss;
  for (size_t i = 0, i_size = value.size(); i < i_size; ++i) {
    if (i != 0) {
      ss << ",";
    }
    ss << std::to_string(value[i]);
  }
  return ss.str();
}

string to_string(const Transaction<vector<double>>& ts)
{
  std::stringstream ss;
  {
    for (const vector<vector<double>>& a_ts : ts) {
      std::for_each(
          a_ts.begin(), a_ts.end(), [&ss](const vector<double>& value) {
            ss << to_string(value);
          });
    }
  }
  return ss.str();
}

string to_string(const Transaction<Range>& ts)
{
  std::stringstream ss;
  {
    for (const vector<Range>& a_ts : ts) {
      std::for_each(a_ts.begin(), a_ts.end(), [&ss](const Range& range) {
        ss << to_string(range);
      });
    }
  }
  return ss.str();
}

string to_string(const vector<vector<double>>& sequence)
{
  std::stringstream ss;
  {
    for (size_t i = 0, i_size = sequence.size(); i < i_size; ++i) {
      if (i != 0) {
        ss << ",";
      }
      ss << to_string(sequence[i]);
    }
  }
  return ss.str();
}

bool include(const Range& alice, const Range& bob)
{
  auto inc_range = [](const pair<int, int>& a, const pair<int, int>& b) {
    return a.first <= b.first && b.second <= a.second;
  };

  const int event_size = alice.size();
  if (static_cast<int>(bob.size()) != event_size) {
    cerr << "BUG ERROR! twin range is not same size" << endl
         << "alice size = " << alice.size() << endl
         << "bob size = " << bob.size() << endl;
    exit(1);
  }

  for (int k = 0; k < event_size; ++k) {
    if (false == inc_range(alice[k], bob[k])) {
      return false;
    }
  }

  return true;
}

}  // namespace PPM
