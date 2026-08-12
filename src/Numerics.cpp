#include "Numerics.h"

#include "SimulationTypes.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace coulomb {

double Numerics::error_function(double x) const {
  // Abramowitz and Stegun formula 7.1.26, retained for numerical continuity.
  constexpr double a1 = 0.254829592;
  constexpr double a2 = -0.284496736;
  constexpr double a3 = 1.421413741;
  constexpr double a4 = -1.453152027;
  constexpr double a5 = 1.061405429;
  constexpr double p = 0.3275911;

  const int sign = x < 0 ? -1 : 1;
  x = std::abs(x);
  const double t = 1.0 / (1.0 + p * x);
  const double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t *
                             std::exp(-x * x);
  return sign * y;
}

std::vector<double> Numerics::circular_shift(const std::vector<double> &values,
                                             int size,
                                             int shift_distance) const {
  if (size < 0 || static_cast<std::size_t>(size) != values.size())
    throw std::invalid_argument(
        "Numerics::circular_shift size does not match values");
  if (size == 0)
    return {};
  std::vector<double> shifted(size);
  shift_distance %= size;
  for (int index = 0; index < size; ++index) {
    int old_index = index + shift_distance;
    if (old_index >= size)
      old_index -= size;
    if (old_index < 0)
      old_index += size;
    shifted[index] = values[old_index];
  }
  return shifted;
}

std::vector<double> Numerics::boundary_shift(const std::vector<double> &values,
                                             int size, int shift_distance,
                                             double boundary) const {
  if (size < 0 || static_cast<std::size_t>(size) != values.size())
    throw std::invalid_argument(
        "Numerics::boundary_shift size does not match values");
  if (size == 0)
    return {};
  std::vector<double> shifted(size);
  for (int index = 0; index < size; ++index) {
    const int old_index = index + shift_distance;
    shifted[index] =
        (old_index >= size || old_index < 0) ? boundary : values[old_index];
  }
  return shifted;
}

std::vector<double>
Numerics::central_difference(const std::vector<double> &values, int size,
                             BoundaryCondition boundary) const {
  if (size < 0 || static_cast<std::size_t>(size) != values.size())
    throw std::invalid_argument(
        "Numerics::central_difference size does not match values");
  if (boundary != BoundaryCondition::Periodic &&
      boundary != BoundaryCondition::Reflective)
    throw std::invalid_argument(
        "Numerics::central_difference boundary is invalid");
  if (size == 0)
    return {};
  const auto left = boundary == BoundaryCondition::Periodic
                        ? circular_shift(values, size, -1)
                        : boundary_shift(values, size, -1, values[0]);
  const auto right = boundary == BoundaryCondition::Periodic
                         ? circular_shift(values, size, 1)
                         : boundary_shift(values, size, 1, values[size - 1]);

  std::vector<double> difference(size);
  for (int index = 0; index < size; ++index)
    difference[index] = 0.5 * (right[index] - left[index]);
  return difference;
}

std::vector<double>
Numerics::limited_flux(const std::vector<double> &values, int size, double dx,
                       double dt, int velocity_sign,
                       const std::vector<double> &boundary_values,
                       BoundaryCondition boundary) const {
  if (size < 0 || static_cast<std::size_t>(size) != values.size())
    throw std::invalid_argument(
        "Numerics::limited_flux size does not match values");
  if (size == 0)
    return {};
  if (dx == 0.0)
    throw std::invalid_argument("Numerics::limited_flux requires dx != 0");
  if (velocity_sign != 1 && velocity_sign != -1)
    throw std::invalid_argument(
        "Numerics::limited_flux velocity sign must be +/-1");
  if (boundary != BoundaryCondition::Periodic &&
      boundary != BoundaryCondition::Reflective)
    throw std::invalid_argument("Numerics::limited_flux boundary is invalid");
  if (boundary == BoundaryCondition::Reflective &&
      boundary_values.size() != values.size())
    throw std::invalid_argument(
        "Numerics::limited_flux boundary values must match the grid size");

  const auto left =
      boundary == BoundaryCondition::Periodic
          ? circular_shift(values, size, -1)
          : boundary_shift(values, size, -1, boundary_values.front());
  const auto right =
      boundary == BoundaryCondition::Periodic
          ? circular_shift(values, size, 1)
          : boundary_shift(values, size, 1, boundary_values.back());

  std::vector<double> difference(size);
  for (int index = 0; index < size; ++index) {
    const double one_sided_difference = velocity_sign > 0
                                            ? values[index] - left[index]
                                            : right[index] - values[index];
    difference[index] = dt / dx * one_sided_difference;
  }
  return difference;
}

void Numerics::advance_kinetic_euler(const std::vector<double> &density,
                                     const std::vector<double> &velocity,
                                     const std::vector<double> &temperature,
                                     int size, double dx, double dt,
                                     BoundaryCondition boundary,
                                     std::vector<double> &density_change,
                                     std::vector<double> &momentum_change,
                                     std::vector<double> &energy_change) const {
  if (size < 0 || static_cast<std::size_t>(size) != density.size() ||
      velocity.size() != density.size() || temperature.size() != density.size())
    throw std::invalid_argument(
        "Numerics::advance_kinetic_euler input size mismatch");
  if (density_change.size() != density.size() ||
      momentum_change.size() != density.size() ||
      energy_change.size() != density.size())
    throw std::invalid_argument(
        "Numerics::advance_kinetic_euler output size mismatch");
  if (size == 0)
    return;
  if (dx == 0.0)
    throw std::invalid_argument(
        "Numerics::advance_kinetic_euler requires dx != 0");

  std::vector<double> g1_positive(size), g1_negative(size);
  std::vector<double> g2_positive(size), g2_negative(size);
  std::vector<double> g3_positive(size), g3_negative(size);
  constexpr double euler_pi = 3.141592653589793238462643383279502884;
  constexpr double lambda_state = 1.0; // corresponding to Dim = 3

  for (int index = 0; index < size; ++index) {
    if (!(temperature[index] > 0.0))
      throw std::invalid_argument(
          "Numerics::advance_kinetic_euler requires positive temperature");
    const double positive_fraction =
        0.5 * (1.0 + error_function(velocity[index] /
                                    std::sqrt(2.0 * temperature[index])));
    const double negative_fraction = 1.0 - positive_fraction;
    const double gaussian_tail =
        std::sqrt(temperature[index] / (2.0 * euler_pi)) *
        std::exp(-0.5 * velocity[index] * velocity[index] / temperature[index]);

    g1_positive[index] =
        density[index] * (velocity[index] * positive_fraction + gaussian_tail);
    g1_negative[index] =
        density[index] * (velocity[index] * negative_fraction - gaussian_tail);
    g2_positive[index] =
        density[index] *
        ((temperature[index] + velocity[index] * velocity[index]) *
             positive_fraction +
         velocity[index] * gaussian_tail);
    g2_negative[index] =
        density[index] *
        ((temperature[index] + velocity[index] * velocity[index]) *
             negative_fraction -
         velocity[index] * gaussian_tail);
    g3_positive[index] =
        density[index] *
        ((1.5 * temperature[index] + 0.5 * velocity[index] * velocity[index]) *
             velocity[index] * positive_fraction +
         (temperature[index] + 0.5 * velocity[index] * velocity[index]) *
             gaussian_tail);
    g3_negative[index] =
        density[index] *
        ((1.5 * temperature[index] + 0.5 * velocity[index] * velocity[index]) *
             velocity[index] * negative_fraction -
         (temperature[index] + 0.5 * velocity[index] * velocity[index]) *
             gaussian_tail);
    g3_positive[index] +=
        lambda_state * temperature[index] * g1_positive[index];
    g3_negative[index] +=
        lambda_state * temperature[index] * g1_negative[index];
  }

  const auto d1_positive =
      limited_flux(g1_positive, size, dx, dt, 1, g1_positive, boundary);
  const auto d1_negative =
      limited_flux(g1_negative, size, dx, dt, -1, g1_negative, boundary);
  const auto d2_positive =
      limited_flux(g2_positive, size, dx, dt, 1, g2_positive, boundary);
  const auto d2_negative =
      limited_flux(g2_negative, size, dx, dt, -1, g2_negative, boundary);
  const auto d3_positive =
      limited_flux(g3_positive, size, dx, dt, 1, g3_positive, boundary);
  const auto d3_negative =
      limited_flux(g3_negative, size, dx, dt, -1, g3_negative, boundary);
  for (int index = 0; index < size; ++index) {
    density_change[index] = d1_positive[index] + d1_negative[index];
    momentum_change[index] = d2_positive[index] + d2_negative[index];
    energy_change[index] = d3_positive[index] + d3_negative[index];
  }
}

} // namespace coulomb
