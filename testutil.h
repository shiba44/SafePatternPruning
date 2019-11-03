
#ifndef PPM_TESTUTIL_H
#define PPM_TESTUTIL_H

#include <vector>

namespace PPM::TestUtil
{
double macro_f1_score(const std::vector<int>& predict,
                      const std::vector<int>& label);

}

#endif
