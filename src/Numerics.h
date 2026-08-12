#pragma once

#include <vector>

namespace coulomb {

enum class BoundaryCondition;

class Numerics {
public:
  static double error_function(double x);
  static std::vector<double> circular_shift(const std::vector<double> &values,
                                            int size, int shift_distance);
  static std::vector<double> boundary_shift(const std::vector<double> &values,
                                            int size, int shift_distance,
                                            double boundary);
  static std::vector<double>
  central_difference(const std::vector<double> &values, int size,
                     BoundaryCondition boundary);
  static std::vector<double>
  limited_flux(const std::vector<double> &values, int size, double dx,
               double dt, int velocity_sign,
               const std::vector<double> &boundary_values,
               BoundaryCondition boundary);
  static void advance_kinetic_euler(
      const std::vector<double> &density, const std::vector<double> &velocity,
      const std::vector<double> &temperature, int size, double dx, double dt,
      BoundaryCondition boundary, std::vector<double> &density_change,
      std::vector<double> &momentum_change, std::vector<double> &energy_change);
};

} // namespace coulomb
