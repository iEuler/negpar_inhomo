#pragma once

#include <vector>

namespace coulomb {

class IniValClass;
class NeParticleGroup;
class NumericGridClass;
struct RandomContext;
struct SimulationState;

void initialize_Negpar(NeParticleGroup& groups, const IniValClass& initial_data,
                       double effective_particles, double effective_full_particles,
                       double dx, RandomContext& random);
void initialize_TwoStreamInstab(IniValClass& initial_data);
void initialize_distri_Negpar(NumericGridClass& grid,
                              std::vector<NeParticleGroup>& groups,
                              SimulationState& state);
void initialize_distri_Negpar_test(NumericGridClass& grid,
                                   std::vector<NeParticleGroup>& groups,
                                   SimulationState& state);

}  // namespace coulomb
