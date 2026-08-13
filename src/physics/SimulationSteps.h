#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct SimulationState;

class SimulationSteps {
  public:
	SimulationSteps(NumericGridClass& grid, ParaClass& parameters,
					SimulationState& state)
		: gridRef(grid), parametersRef(parameters), stateRef(state) {}

	void advanceHdp(std::vector<NeParticleGroup>& groups);
	void advancePic(std::vector<NeParticleGroup>& groups);

  private:
	NumericGridClass& gridRef;
	ParaClass& parametersRef;
	SimulationState& stateRef;
};

} // namespace coulomb
