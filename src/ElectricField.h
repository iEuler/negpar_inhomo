#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParticleGroup;

class ElectricFieldSolver {
public:
  explicit ElectricFieldSolver(const NumericGridClass &grid) : grid_(grid) {}

  std::vector<double> solve_poisson(const std::vector<double> &rho,
                                    double lambda);
  void update(std::vector<ParticleGroup> &groups);
  void update(std::vector<NeParticleGroup> &groups);
  void update_pic(std::vector<NeParticleGroup> &groups);
  void update_from_coarse(std::vector<NeParticleGroup> &groups);
  void clear(std::vector<NeParticleGroup> &groups);
  void update_from_density(std::vector<NeParticleGroup> &groups);

private:
  const NumericGridClass &grid_;
};

} // namespace coulomb
