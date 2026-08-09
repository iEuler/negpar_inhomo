#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

void coulomb_collision_homo_PFNF(NeParticleGroup& groups,
                                  const ParaClass& parameters,
                                  RandomContext& random);
void NegPar_collision_homo(NeParticleGroup& groups,
                           const ParaClass& parameters,
                           double effective_particles,
                           RandomContext& random);
void NegPar_collision(std::vector<NeParticleGroup>& groups,
                      const NumericGridClass& grid,
                      const ParaClass& parameters, RandomContext& random);
void NegPar_collision_openmp(std::vector<NeParticleGroup>& groups,
                             const NumericGridClass& grid,
                             const ParaClass& parameters,
                             RandomContext& random);
void NegPar_BGK_collision_homo(NeParticleGroup& groups, ParaClass& parameters,
                               RandomContext& random);
void NegPar_BGK_collision(std::vector<NeParticleGroup>& groups,
                          NumericGridClass& grid, ParaClass& parameters,
                          RandomContext& random);

}  // namespace coulomb
