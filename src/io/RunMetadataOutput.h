#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;
class ParaClass;
struct SimulationState;

class RunMetadataOutput {
  public:
	void saveGrid(const NumericGridClass& grid, const SimulationState& state);
	void saveParameters(const ParaClass& parameters,
						const NumericGridClass& grid,
						const SimulationState& state);
	void saveInitialConditions(IniValClass& initialData,
							   const SimulationState& state);
};

} // namespace coulomb
