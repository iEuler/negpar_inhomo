#pragma once

#include <vector>

#include "Particle.h"

namespace coulomb {

class NumericGridClass;
class ParticleGroup;
class NeParticleGroup;
struct SimulationState;

void save_dist(std::vector<ParticleGroup>& groups,
               const NumericGridClass& grid, const SimulationState& state);
void save_dist(std::vector<NeParticleGroup>& groups,
               const NumericGridClass& grid, ParticleKind kind,
               const SimulationState& state);
void save_dist(std::vector<NeParticleGroup>& groups,
               const NumericGridClass& grid, const SimulationState& state);
void save_particle1d1d(std::vector<ParticleGroup>& groups,
                       const NumericGridClass& grid,
                       const SimulationState& state);
void save_particle1d1d(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid, ParticleKind kind,
                       int quantity, const SimulationState& state);
void save_particle1d1d(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid,
                       const SimulationState& state);
void save_particleenergy(std::vector<NeParticleGroup>& groups,
                         const NumericGridClass& grid,
                         const SimulationState& state);
void save_homo_rdist(const SimulationState& state);
void save_homo_rdist(int bin_count, const SimulationState& state);
void save_homo_dist(const NeParticleGroup& group, int bin_count,
                    int case_index, const SimulationState& state);

}  // namespace coulomb
