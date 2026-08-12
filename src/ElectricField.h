#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParticleGroup;

class ElectricFieldSolver {
public:
  static std::vector<double> solve_poisson(const std::vector<double> &rho,
                                           int grid_size, double domain_size,
                                           double lambda);
  static void update(std::vector<ParticleGroup> &groups,
                     const NumericGridClass &grid);
  static void update(std::vector<NeParticleGroup> &groups,
                     const NumericGridClass &grid);
  static void update_pic(std::vector<NeParticleGroup> &groups,
                         const NumericGridClass &grid);
  static void update_from_coarse(std::vector<NeParticleGroup> &groups,
                                 const NumericGridClass &grid);
  static void clear(std::vector<NeParticleGroup> &groups,
                    const NumericGridClass &grid);
  static void update_from_density(std::vector<NeParticleGroup> &groups,
                                  const NumericGridClass &grid);
};

} // namespace coulomb
