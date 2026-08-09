#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct SimulationState;

void initialize_distri_Negpar(NumericGridClass& grid,
                              std::vector<NeParticleGroup>& groups,
                              SimulationState& state);

}  // namespace coulomb
