#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct RandomContext;

class ProjectionSampling {
public:
  static void sample_homogeneous(NeParticleGroup &groups,
                                 const NumericGridClass &grid,
                                 RandomContext &random);
  static void sample(std::vector<NeParticleGroup> &groups,
                     const NumericGridClass &grid, RandomContext &random);
};

} // namespace coulomb
