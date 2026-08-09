#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;

bool particleresample_homo(NeParticleGroup& groups,
                           const ParaClass& parameters,
                           SimulationState& state);
void particleresample_inhomo(std::vector<NeParticleGroup>& groups,
                             NumericGridClass& grid, ParaClass& parameters,
                             SimulationState& state);
void resampleF_homo(NeParticleGroup& groups, double new_effective_particles,
                    double effective_particles, int frequency,
                    double dx_space, RandomContext& random);
void resampleF_inhomo(std::vector<NeParticleGroup>& groups,
                      double new_effective_particles, NumericGridClass& grid,
                      int frequency, SimulationState& state);
void resampleF_keeptotalmass(std::vector<NeParticleGroup>& groups,
                             NumericGridClass& grid, int old_count,
                             RandomContext& random);
void sync_coarse(std::vector<NeParticleGroup>& groups, NumericGridClass& grid,
                 ParaClass& parameters, SimulationState& state);

}  // namespace coulomb
