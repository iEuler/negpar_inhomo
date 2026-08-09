#pragma once

// Stochastic full-particle reconstruction orchestration.

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

NeParticleGroup resample_F_from_MPN(NeParticleGroup& groups, int frequency,
                                    double effective_particles,
                                    double effective_full_particles,
                                    double dx_space, RandomContext& random);

}  // namespace coulomb
