#include "Histogram.h"

namespace coulomb {

void histinfo_fixbar(const std::vector<double>& values,
                     std::vector<int>& counts, double minimum,
                     double maximum) {
  const int value_count = static_cast<int>(values.size());
  const int bin_count = static_cast<int>(counts.size());
  const double width = (maximum - minimum) / bin_count;

  for (int bin = 0; bin < bin_count; ++bin) counts[bin] = 0;

  for (int index = 0; index < value_count; ++index) {
    int bin = static_cast<int>((values[index] - minimum) / width);
    if (bin < 0) bin = 0;
    if (bin >= bin_count) bin = bin_count - 1;
    ++counts[bin];
  }
}

}  // namespace coulomb
