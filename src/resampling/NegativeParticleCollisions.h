#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

class NegativeParticleCollisions {
  public:
	NegativeParticleCollisions(const NumericGridClass& grid,
							   ParaClass& parameters, RandomContext& random)
		: gridRef(grid), parametersRef(parameters), randomContext(random) {}

	void collideWithFull(NeParticleGroup& groups);
	void collideHomogeneous(NeParticleGroup& groups);
	void collide(std::vector<NeParticleGroup>& groups);
	void collideParallel(std::vector<NeParticleGroup>& groups);
	void collideBgkHomogeneous(NeParticleGroup& groups);
	void collideBgk(std::vector<NeParticleGroup>& groups);

  private:
	const NumericGridClass& gridRef;
	ParaClass& parametersRef;
	RandomContext& randomContext;
};

} // namespace coulomb
