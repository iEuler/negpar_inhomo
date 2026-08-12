#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

class NegativeParticleCollisions {
public:
  static void collide_with_full(NeParticleGroup &groups,
                                const ParaClass &parameters,
                                RandomContext &random);
  static void collide_homogeneous(NeParticleGroup &groups,
                                  const ParaClass &parameters,
                                  double effective_particles,
                                  RandomContext &random);
  static void collide(std::vector<NeParticleGroup> &groups,
                      const NumericGridClass &grid, const ParaClass &parameters,
                      RandomContext &random);
  static void collide_parallel(std::vector<NeParticleGroup> &groups,
                               const NumericGridClass &grid,
                               const ParaClass &parameters,
                               RandomContext &random);
  static void collide_bgk_homogeneous(NeParticleGroup &groups,
                                      ParaClass &parameters,
                                      RandomContext &random);
  static void collide_bgk(std::vector<NeParticleGroup> &groups,
                          NumericGridClass &grid, ParaClass &parameters,
                          RandomContext &random);
};

} // namespace coulomb
