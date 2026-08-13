#pragma once

#include <string>

namespace coulomb {

struct SimulationState;

class OutputPaths {
  public:
	explicit OutputPaths(const SimulationState& state) : state_(state) {}

	std::string resolve(const std::string& filename) const;
	std::string format_index(int value, int digits = 3) const;

  private:
	const SimulationState& state_;
};

} // namespace coulomb
