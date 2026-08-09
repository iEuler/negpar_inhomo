#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
enum class ParticleKind;

// Diagnostics are pure reductions over particle and grid data. They do not
// mutate run counters, output policy, or timing, so they intentionally remain
// state-free while stateful callers pass SimulationState to the APIs that need
// it.

double compute_elec_energy(const std::vector<NeParticleGroup>& S_x,
                           const NumericGridClass& grid);
double compute_elec_energy_F(const std::vector<NeParticleGroup>& S_x,
                             const NumericGridClass& grid);
double compute_total_energy(const std::vector<NeParticleGroup>& S_x,
                            const NumericGridClass& grid);
double compute_total_energy_F(const std::vector<NeParticleGroup>& S_x,
                              const NumericGridClass& grid);
int count_particle_number(const std::vector<NeParticleGroup>& groups, int size,
                          ParticleKind kind);

}  // namespace coulomb
