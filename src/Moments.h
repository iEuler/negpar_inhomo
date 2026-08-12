#pragma once

#include <tuple>
#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParticleGroup;
struct SimulationState;

class MomentOperations {
public:
  static void primitive_to_conserved(const double &rho, const double *velocity,
                                     const double &temperature,
                                     double *momentum, double &energy);
  static void conserved_to_primitive(const double &rho, const double *momentum,
                                     const double &energy, double *velocity,
                                     double &temperature);

  static void primitive_to_conserved(int grid_size,
                                     const std::vector<double> &rho,
                                     const std::vector<double> &velocity,
                                     const std::vector<double> &temperature,
                                     std::vector<double> &momentum,
                                     std::vector<double> &energy);
  static void conserved_to_primitive(int grid_size,
                                     const std::vector<double> &rho,
                                     const std::vector<double> &momentum,
                                     const std::vector<double> &energy,
                                     std::vector<double> &velocity,
                                     std::vector<double> &temperature);

  static void particle_to_conserved(const ParticleGroup &group,
                                    double effective_particles, double &rho,
                                    double *momentum, double &energy);
  static void
  compute_primitive(int grid_size, const std::vector<ParticleGroup> &groups,
                    double effective_particles, std::vector<double> &rho,
                    std::vector<double> &u1, std::vector<double> &u2,
                    std::vector<double> &u3, std::vector<double> &temperature);

  static void update_primitive(std::vector<NeParticleGroup> &groups,
                               const NumericGridClass &grid);
  static void update_full_primitive(std::vector<NeParticleGroup> &groups,
                                    const NumericGridClass &grid);
  static void
  update_maxwellian_derivatives(std::vector<NeParticleGroup> &groups,
                                const NumericGridClass &grid);
  static void update_macro(std::vector<NeParticleGroup> &groups,
                           const NumericGridClass &grid);
  static void update_maxwellian(std::vector<NeParticleGroup> &groups,
                                const NumericGridClass &grid);
  static void compute_macro_change(std::vector<NeParticleGroup> &groups,
                                   const NumericGridClass &grid,
                                   SimulationState &state);
  static void compute_kinetic_macro_change(std::vector<NeParticleGroup> &groups,
                                           const NumericGridClass &grid);

  static std::tuple<std::vector<double>, std::vector<double>,
                    std::vector<double>>
  moment_change(const NeParticleGroup *groups, const NumericGridClass &grid);
};

} // namespace coulomb
