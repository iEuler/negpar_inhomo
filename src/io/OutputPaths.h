#pragma once

#include <string>

namespace coulomb {

struct SimulationState;

class OutputPaths {
  public:
	explicit OutputPaths(const SimulationState& state) : stateRef(state) {}

	std::string resolve(const std::string& filename) const;
	std::string formatIndex(int value, int digits = 3) const;

  private:
	const SimulationState& stateRef;
};

} // namespace coulomb
