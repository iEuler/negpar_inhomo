#pragma once

namespace coulomb {

class IniValClass;
class NeParticleGroup;
struct RandomContext;

class ParticleInitialization {
public:
  void initialize(NeParticleGroup &groups, const IniValClass &initial_data,
                  double effective_particles, double effective_full_particles,
                  double dx, RandomContext &random);
};

} // namespace coulomb
