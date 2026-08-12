#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace coulomb::resampling {

class ResamplingNumerics {
public:
  std::vector<double> frequencies(std::size_t count);
  std::vector<std::size_t> augmented_locations(std::size_t count,
                                               std::size_t augmentation_factor);
  std::vector<std::complex<double>> imaginary_frequencies(std::size_t count);
  double evaluate_quadratic_taylor(double delta_x, double delta_y,
                                   double delta_z,
                                   const std::vector<double> &derivatives);
};

} // namespace coulomb::resampling
