#pragma once

#include <vector>

namespace coulomb {

class Histogram {
  public:
	void fixedBins(const std::vector<double>& values, std::vector<int>& counts,
				   double minimum, double maximum) const;
};

} // namespace coulomb
