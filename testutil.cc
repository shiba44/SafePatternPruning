
#include "testutil.h"
#include <algorithm>
#include "type.h"

using std::vector;

namespace PPM::TestUtil
{
double macro_f1_score(const std::vector<int> &predict,
                      const std::vector<int> &correct)
{
  auto calc_f_value = [](const vector<int> &predict,
                         const vector<int> &correct,
                         const int POSITIVE) -> double {
    int num_predict_positive =
        std::count(predict.begin(), predict.end(), POSITIVE);
    int num_correct_positive =
        std::count(correct.begin(), correct.end(), POSITIVE);

    int TP = 0;
    for (int i = 0, size = correct.size(); i < size; i++) {
      if (predict[i] == POSITIVE && correct[i] == POSITIVE) TP++;
    }

    return 2.0 * TP / (num_predict_positive + num_correct_positive);
  };

  double f_value         = calc_f_value(predict, correct, Label::POSITIVE);
  double reverse_f_value = calc_f_value(predict, correct, Label::NEGATIVE);

  return (f_value + reverse_f_value) / 2.0;
}

}  // namespace PPM::TestUtil
