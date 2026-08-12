#pragma once

// Stochastic full-particle reconstruction and sampling.

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

class FullParticleSampling {
public:
  NeParticleGroup resample(NeParticleGroup &groups, int frequency,
                           double effective_particles,
                           double effective_full_particles, double dx_space,
                           RandomContext &random);
};

} // namespace coulomb
