#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct SimulationState;

class SimulationSteps {
public:
  static void advance_hdp(std::vector<NeParticleGroup> &groups,
                          NumericGridClass &grid, ParaClass &parameters,
                          SimulationState &state);
  static void advance_pic(std::vector<NeParticleGroup> &groups,
                          NumericGridClass &grid, ParaClass &parameters,
                          SimulationState &state);
};

} // namespace coulomb
