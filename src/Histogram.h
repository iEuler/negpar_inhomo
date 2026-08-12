#pragma once

#include <vector>

namespace coulomb {

class Histogram {
public:
  static void fixed_bins(const std::vector<double> &values,
                         std::vector<int> &counts, double minimum,
                         double maximum);
};

} // namespace coulomb
