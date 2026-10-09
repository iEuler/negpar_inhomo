#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct SimulationState;

class Initialization {
  public:
	void initialize(NumericGridClass& grid,
					std::vector<NeParticleGroup>& groups,
					SimulationState& state, double landauAmplitude = 0.4);
};

} // namespace coulomb
