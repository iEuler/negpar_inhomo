#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct RandomContext;

class NegativeParticleSampling {
public:
  double evaluate_maxwellian(const std::vector<double> &velocity,
                             const NeParticleGroup &groups);
  double evaluate_source(const std::vector<double> &velocity,
                         const std::vector<double> &source_velocity,
                         const NeParticleGroup &groups,
                         const ParaClass &parameters, int mode = 0);
  void update_bounds(NeParticleGroup &groups, const ParaClass &parameters);
  void update_bounds(std::vector<NeParticleGroup> &groups,
                     const NumericGridClass &grid, const ParaClass &parameters);
  int estimate_virtual_count(const NeParticleGroup &groups,
                             double effective_particles, RandomContext &random);
  void sample_delta(NeParticleGroup &groups, NeParticleGroup &new_groups,
                    const ParaClass &parameters, double effective_particles,
                    RandomContext &random);
};

} // namespace coulomb
