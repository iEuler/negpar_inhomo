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

class Diagnostics {
public:
  static double electric_energy(const std::vector<NeParticleGroup> &groups,
                                const NumericGridClass &grid);
  static double full_electric_energy(const std::vector<NeParticleGroup> &groups,
                                     const NumericGridClass &grid);
  static double total_energy(const std::vector<NeParticleGroup> &groups,
                             const NumericGridClass &grid);
  static double full_total_energy(const std::vector<NeParticleGroup> &groups,
                                  const NumericGridClass &grid);
  static int particle_count(const std::vector<NeParticleGroup> &groups,
                            int size, ParticleKind kind);
};

} // namespace coulomb
