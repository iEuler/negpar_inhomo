#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
enum class ParticleKind;
struct RandomContext;

class ParticleGroupOperations {
  public:
	void mergeSigned(NeParticleGroup& groups, const NeParticleGroup& newGroups);
	void mergeFull(NeParticleGroup& groups, const NeParticleGroup& newGroups);
	void merge(NeParticleGroup& groups, const NeParticleGroup& newGroups,
			   const std::vector<ParticleKind>& particleTypes);
	void assignPositions(NeParticleGroup& groups, double xmin, double xmax,
						 RandomContext& random);
};

} // namespace coulomb
