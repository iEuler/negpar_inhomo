#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
enum class ParticleKind;
struct RandomContext;

void merge_NeParticleGroup(NeParticleGroup& groups,
                           const NeParticleGroup& new_groups);
void mergeF_NeParticleGroup(NeParticleGroup& groups,
                            const NeParticleGroup& new_groups);
void mergeNeParticleGroup(NeParticleGroup& groups,
                          const NeParticleGroup& new_groups,
                          const std::vector<ParticleKind>& particle_types);
void assign_positions(NeParticleGroup& groups, double xmin, double xmax,
                      RandomContext& random);

}  // namespace coulomb
