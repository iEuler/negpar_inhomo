#pragma once

#include <vector>

namespace coulomb {

void histinfo_fixbar(const std::vector<double>& values,
                     std::vector<int>& counts, double minimum,
                     double maximum);

}  // namespace coulomb
