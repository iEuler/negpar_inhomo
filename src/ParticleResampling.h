#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;

class ParticleResampling {
public:
  static bool resample_homogeneous(NeParticleGroup &groups,
                                   const ParaClass &parameters,
                                   SimulationState &state);
  static void resample(std::vector<NeParticleGroup> &groups,
                       NumericGridClass &grid, ParaClass &parameters,
                       SimulationState &state);
  static void resample_full_homogeneous(NeParticleGroup &groups,
                                        double new_effective_particles,
                                        double effective_particles,
                                        int frequency, double dx_space,
                                        RandomContext &random);
  static void resample_full(std::vector<NeParticleGroup> &groups,
                            double new_effective_particles,
                            NumericGridClass &grid, int frequency,
                            SimulationState &state);
  static void
  resample_full_preserving_mass(std::vector<NeParticleGroup> &groups,
                                NumericGridClass &grid, int old_count,
                                RandomContext &random);
  static void synchronize_coarse(std::vector<NeParticleGroup> &groups,
                                 NumericGridClass &grid, ParaClass &parameters,
                                 SimulationState &state);
};

} // namespace coulomb
