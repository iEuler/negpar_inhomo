#pragma once

#include <vector>

namespace coulomb {
class NeParticleGroup;
class Particle1d3d;

namespace resampling {

class ResamplingVelocity {
public:
  NeParticleGroup normalize_signed(const NeParticleGroup &particles);
  void restore(std::vector<Particle1d3d> &particles,
               const std::vector<double> &velocity_bounds);
};

} // namespace resampling
} // namespace coulomb
