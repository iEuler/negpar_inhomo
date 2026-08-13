#pragma once

#include <utility>
#include <vector>

namespace coulomb {

class ParaClass;
class Particle1d3d;
struct RandomContext;

class CollisionOperator {
  public:
	CollisionOperator(const ParaClass& parameters, RandomContext& random)
		: parameters_(parameters), random_(random) {}

	std::pair<std::vector<double>, std::vector<double>>
	collide_pair(const std::vector<double>& velocity1,
				 const std::vector<double>& velocity2);

	void collide_homogeneous(std::vector<Particle1d3d>& particles,
							 int particle_count);

  private:
	const ParaClass& parameters_;
	RandomContext& random_;
};

} // namespace coulomb
