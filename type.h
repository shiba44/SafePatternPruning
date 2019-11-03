
#ifndef PPM_TYPE_H
#define PPM_TYPE_H

#include <iostream>

#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using std::pair;
using std::string;
using std::vector;

using std::cerr;
using std::endl;

namespace PPM
{
using Range = vector<pair<int, int>>;

bool include(const Range& alice, const Range& bob);

/* FIXME: this function should be Pattern method */
bool is_specific_and_abstract(const vector<Range>& alice,
                              const vector<Range>& bob);

bool include(const vector<Range>& alice, const vector<Range>& bob);

struct Node {
    vector<Range> pattern;
    string pattern_str;
    vector<pair<int, int>> pdb;
    vector<int> x_id;
    double w;
    int support;
    double val;
    bool is_leaf;
};

using TreeIter = std::list<Node>::iterator;

/* type of pattern */
enum TypePattern {
    UNIQ_SYMBOL_FREQ = 1,
    UNIQ_SYMBOL_CLOS = 2,
    UNIQ_SYMBOL_GENE = 3,
    ADAP_SYMBOL_FREQ = 11,
    ADAP_SYMBOL_CLOS = 12,
    ADAP_SYMBOL_GENE = 13,
    ADAP_RANGE_FREQ  = 21,
    ADAP_RANGE_CLOS  = 22,
    ADAP_RANGE_GENE  = 23
};

template <class T>
using Transaction = vector<vector<T>>;

struct Result {
    double bias;
    vector<TreeIter> model;
    vector<double> weight;
};

// label
namespace Label
{
const int POSITIVE    = 1;
const int NEGATIVE    = -1;
const int UNDECIDABLE = 0;
}  // namespace Label

string to_string(const Range& range);
string to_string(const vector<Range>& pat);
string to_string(const Transaction<Range>& ts);
string to_string(const vector<double>& value);
string to_string(const Transaction<vector<double>>& ts);
string to_string(const vector<vector<double>>& ts);






}  // namespace PPM

#endif
