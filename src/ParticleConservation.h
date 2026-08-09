#pragma once

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

void enforce_conservation(double m0, double m11, double m12, double m13,
                          double m21, double m22, double m23,
                          NeParticleGroup& groups, double effective_particles,
                          bool conserve_energy_vector, RandomContext& random);
void enforce_conservation_zero(NeParticleGroup& groups,
                               double effective_particles,
                               RandomContext& random);

}  // namespace coulomb
