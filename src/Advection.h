#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class Particle1d3d;
class ParticleGroup;
enum class ParticleKind;

class Advection {
public:
  static Particle1d3d move_particle(const Particle1d3d &particle,
                                    double electric_field,
                                    const NumericGridClass &grid,
                                    SimulationState &state);
  static int find_particle_group(Particle1d3d &particle,
                                 const NumericGridClass &grid);

  static void relocate_particle(std::vector<ParticleGroup> &groups,
                                int group_before, int particle_index,
                                int group_after);
  static void reset_moved_flags(std::vector<ParticleGroup> &groups,
                                int grid_size);
  static void advance(std::vector<ParticleGroup> &groups,
                      const NumericGridClass &grid, SimulationState &state);

  static void relocate_particle(std::vector<NeParticleGroup> &groups,
                                ParticleKind kind, int group_before,
                                int particle_index, int group_after);
  static void reset_moved_flags(std::vector<NeParticleGroup> &groups,
                                ParticleKind kind, int grid_size);

  static void advance(std::vector<NeParticleGroup> &groups, ParticleKind kind,
                      const NumericGridClass &grid, SimulationState &state);
  static void advance(std::vector<NeParticleGroup> &groups,
                      const NumericGridClass &grid, SimulationState &state);
};

} // namespace coulomb
