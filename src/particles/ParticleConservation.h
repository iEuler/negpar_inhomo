#pragma once

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

class ParticleConservation {
  public:
	void enforce(double m0, double m11, double m12, double m13, double m21,
				 double m22, double m23, NeParticleGroup& groups,
				 double effectiveParticles, bool conserveEnergyVector,
				 RandomContext& random);
	void enforceZero(NeParticleGroup& groups, double effectiveParticles,
					 RandomContext& random);
};

} // namespace coulomb
