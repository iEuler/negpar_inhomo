#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

double evaluateM(const std::vector<double>& velocity,
                 const NeParticleGroup& groups);
double evaluateH(const std::vector<double>& velocity,
                 const std::vector<double>& source_velocity,
                 const NeParticleGroup& groups, const ParaClass& parameters,
                 int mode = 0);
void finddeltambound(NeParticleGroup& groups, const ParaClass& parameters);
void finddeltambound_inhomo(std::vector<NeParticleGroup>& groups,
                            const NumericGridClass& grid,
                            const ParaClass& parameters);
int samplefromDeltamp_Npv(const NeParticleGroup& groups,
                          double effective_particles,
                          RandomContext& random);
void samplefromDeltam(NeParticleGroup& groups, NeParticleGroup& new_groups,
                      const ParaClass& parameters, double effective_particles,
                      RandomContext& random);

}  // namespace coulomb
