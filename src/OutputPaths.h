#pragma once

#include <string>

namespace coulomb {

struct SimulationState;

class OutputPaths {
public:
  static std::string resolve(const std::string &filename,
                             const SimulationState &state);
  static std::string format_index(int value, int digits = 3);
};

} // namespace coulomb
