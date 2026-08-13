#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

class NegativeParticleSampling {
  public:
	double evaluateMaxwellian(const std::vector<double>& velocity,
							  const NeParticleGroup& groups);
	double evaluateSource(const std::vector<double>& velocity,
						  const std::vector<double>& sourceVelocity,
						  const NeParticleGroup& groups,
						  const ParaClass& parameters, int mode = 0);
	void updateBounds(NeParticleGroup& groups, const ParaClass& parameters);
	void updateBounds(std::vector<NeParticleGroup>& groups,
					  const NumericGridClass& grid,
					  const ParaClass& parameters);
	int estimateVirtualCount(const NeParticleGroup& groups,
							 double effectiveParticles, RandomContext& random);
	void sampleDelta(NeParticleGroup& groups, NeParticleGroup& newGroups,
					 const ParaClass& parameters, double effectiveParticles,
					 RandomContext& random);
};

} // namespace coulomb
