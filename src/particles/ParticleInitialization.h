#pragma once

namespace coulomb {

class IniValClass;
class NeParticleGroup;
struct RandomContext;

class ParticleInitialization {
  public:
	void initialize(NeParticleGroup& groups, const IniValClass& initialData,
					double effectiveParticles, double effectiveFullParticles,
					double dx, RandomContext& random);
};

} // namespace coulomb
