#pragma once

#include "Classes.h"
#include "_global_variables.h"

namespace coulomb {

Particle1d3d moveparticle(const Particle1d3d& particle, double elecfield,
                          const NumericGridClass& grid, SimulationState& state);
int findparticlegroup(Particle1d3d& particle, const NumericGridClass& grid);

void relocateparticle(std::vector<ParticleGroup>& groups, int group_before,
                      int particle_index, int group_after);
void reset_flag_moved(std::vector<ParticleGroup>& groups, int grid_size);
void particleadvection(std::vector<ParticleGroup>& groups,
                       const NumericGridClass& grid, SimulationState& state);

void relocateparticle(std::vector<NeParticleGroup>& groups, char partype,
                      int group_before, int particle_index, int group_after);
void reset_flag_moved(std::vector<NeParticleGroup>& groups, char partype,
                      int grid_size);

void particleadvection(std::vector<NeParticleGroup>& groups, char partype,
                       const NumericGridClass& grid, SimulationState& state);
void particleadvection(std::vector<NeParticleGroup>& groups,
                       ParticleKind kind, const NumericGridClass& grid,
                       SimulationState& state);
void particleadvection(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid, SimulationState& state);

}  // namespace coulomb
