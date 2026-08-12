#pragma once

#include <vector>

#include "Particle.h"

namespace coulomb {

class NumericGridClass;
class ParticleGroup;
class NeParticleGroup;
struct SimulationState;

class ParticleOutput {
public:
  void save_distribution(std::vector<ParticleGroup> &groups,
                         const NumericGridClass &grid,
                         const SimulationState &state);
  void save_distribution(std::vector<NeParticleGroup> &groups,
                         const NumericGridClass &grid, ParticleKind kind,
                         const SimulationState &state);
  void save_distribution(std::vector<NeParticleGroup> &groups,
                         const NumericGridClass &grid,
                         const SimulationState &state);
  void save_phase_space(std::vector<ParticleGroup> &groups,
                        const NumericGridClass &grid,
                        const SimulationState &state);
  void save_phase_space(std::vector<NeParticleGroup> &groups,
                        const NumericGridClass &grid, ParticleKind kind,
                        int quantity, const SimulationState &state);
  void save_phase_space(std::vector<NeParticleGroup> &groups,
                        const NumericGridClass &grid,
                        const SimulationState &state);
  void save_energy(std::vector<NeParticleGroup> &groups,
                   const NumericGridClass &grid, const SimulationState &state);
  void save_homogeneous_radial_distribution(const SimulationState &state);
  void save_homogeneous_radial_distribution(int bin_count,
                                            const SimulationState &state);
  void save_homogeneous_distribution(const NeParticleGroup &group,
                                     int bin_count, int case_index,
                                     const SimulationState &state);
};

} // namespace coulomb
