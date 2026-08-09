#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct RandomContext;

void sample_from_MMprojection_homo(NeParticleGroup& groups,
                                   const NumericGridClass& grid,
                                   RandomContext& random);
void sample_from_MMprojection(std::vector<NeParticleGroup>& groups,
                              const NumericGridClass& grid,
                              RandomContext& random);

}  // namespace coulomb
