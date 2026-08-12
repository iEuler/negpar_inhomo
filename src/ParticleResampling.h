#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;

class ParticleResampling {
public:
  ParticleResampling(NumericGridClass &grid, ParaClass &parameters,
                     SimulationState &state)
      : grid_(grid), parameters_(parameters), state_(state) {}

  bool resample_homogeneous(NeParticleGroup &groups);
  void resample(std::vector<NeParticleGroup> &groups);
  void resample_full_homogeneous(NeParticleGroup &groups,
                                 double new_effective_particles,
                                 double effective_particles, int frequency,
                                 double dx_space);
  void resample_full(std::vector<NeParticleGroup> &groups,
                     double new_effective_particles, int frequency);
  void resample_full_preserving_mass(std::vector<NeParticleGroup> &groups,
                                     int old_count);
  void synchronize_coarse(std::vector<NeParticleGroup> &groups);

private:
  NumericGridClass &grid_;
  ParaClass &parameters_;
  SimulationState &state_;
};

} // namespace coulomb
