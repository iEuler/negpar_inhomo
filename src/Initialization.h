#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct SimulationState;

class Initialization {
public:
  static void initialize(NumericGridClass &grid,
                         std::vector<NeParticleGroup> &groups,
                         SimulationState &state);
};

} // namespace coulomb
