#pragma once

#include <utility>
#include <vector>

namespace coulomb {

class ParaClass;
class Particle1D3D;
struct RandomContext;

class CollisionOperator {
  public:
	CollisionOperator(const ParaClass& parameters, RandomContext& random)
		: parametersRef(parameters), randomContext(random) {}

	std::pair<std::vector<double>, std::vector<double>>
	collidePair(const std::vector<double>& velocity1,
				const std::vector<double>& velocity2);

	void collideHomogeneous(std::vector<Particle1D3D>& particles,
							int particleCount);

  private:
	const ParaClass& parametersRef;
	RandomContext& randomContext;
};

} // namespace coulomb
