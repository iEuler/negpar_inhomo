#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;
class ParaClass;
struct SimulationState;

class RunMetadataOutput {
public:
  static void save_grid(const NumericGridClass &grid,
                        const SimulationState &state);
  static void save_parameters(const ParaClass &parameters,
                              const NumericGridClass &grid,
                              const SimulationState &state);
  static void save_initial_conditions(IniValClass &initial_data,
                                      const SimulationState &state);
};

} // namespace coulomb
