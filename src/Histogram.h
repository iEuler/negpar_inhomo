#pragma once

#include <vector>

namespace coulomb {

class Histogram {
public:
  void fixed_bins(const std::vector<double> &values, std::vector<int> &counts,
                  double minimum, double maximum) const;
};

} // namespace coulomb
