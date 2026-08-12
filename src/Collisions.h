#pragma once

#include <utility>
#include <vector>

namespace coulomb {

class ParaClass;
class Particle1d3d;
struct RandomContext;

class CollisionOperator {
public:
  static std::pair<std::vector<double>, std::vector<double>>
  collide_pair(const std::vector<double> &velocity1,
               const std::vector<double> &velocity2,
               const ParaClass &parameters, RandomContext &random);

  static void collide_homogeneous(std::vector<Particle1d3d> &particles,
                                  int particle_count,
                                  const ParaClass &parameters,
                                  RandomContext &random);
};

} // namespace coulomb
