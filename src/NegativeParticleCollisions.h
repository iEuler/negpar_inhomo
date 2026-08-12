#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

class NegativeParticleCollisions {
public:
  NegativeParticleCollisions(const NumericGridClass &grid,
                             ParaClass &parameters, RandomContext &random)
      : grid_(grid), parameters_(parameters), random_(random) {}

  void collide_with_full(NeParticleGroup &groups);
  void collide_homogeneous(NeParticleGroup &groups);
  void collide(std::vector<NeParticleGroup> &groups);
  void collide_parallel(std::vector<NeParticleGroup> &groups);
  void collide_bgk_homogeneous(NeParticleGroup &groups);
  void collide_bgk(std::vector<NeParticleGroup> &groups);

private:
  const NumericGridClass &grid_;
  ParaClass &parameters_;
  RandomContext &random_;
};

} // namespace coulomb
