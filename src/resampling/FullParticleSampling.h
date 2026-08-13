#pragma once

// Stochastic full-particle reconstruction and sampling.

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

class FullParticleSampling {
  public:
	NeParticleGroup resample(NeParticleGroup& groups, int frequency,
							 double effectiveParticles,
							 double effectiveFullParticles, double dxSpace,
							 RandomContext& random);
};

} // namespace coulomb
