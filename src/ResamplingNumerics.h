#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace coulomb::resampling {

class ResamplingNumerics {
public:
  static std::vector<double> frequencies(std::size_t count);
  static std::vector<std::size_t>
  augmented_locations(std::size_t count, std::size_t augmentation_factor);
  static std::vector<std::complex<double>>
  imaginary_frequencies(std::size_t count);
  static double
  evaluate_quadratic_taylor(double delta_x, double delta_y, double delta_z,
                            const std::vector<double> &derivatives);
};

} // namespace coulomb::resampling
