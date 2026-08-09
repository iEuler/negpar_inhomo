#include "ResamplingNumerics.h"

#include <limits>
#include <stdexcept>

namespace coulomb::resampling {
namespace {

int checked_count(std::size_t count) {
  if (count == 0 ||
      count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument(
        "resampling frequency count must fit in a positive int");
  return static_cast<int>(count);
}

int frequency(std::size_t index, std::size_t count) {
  const int count_int = checked_count(count);
  if (index >= count)
    throw std::out_of_range("resampling frequency index is outside the grid");
  const int index_int = static_cast<int>(index);
  return index_int >= count_int / 2 + 1 ? index_int - count_int : index_int;
}

}  // namespace

std::vector<double> frequencies(std::size_t count) {
  checked_count(count);
  std::vector<double> result(count);
  for (std::size_t index = 0; index < count; ++index)
    result[index] = static_cast<double>(frequency(index, count));
  return result;
}

std::vector<std::size_t> augmented_locations(
    std::size_t count, std::size_t augmentation_factor) {
  checked_count(count);
  if (augmentation_factor == 0 ||
      count > static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                  augmentation_factor)
    throw std::invalid_argument(
        "resampling augmentation must produce a positive int-sized grid");

  const auto augmented_count = count * augmentation_factor;
  std::vector<std::size_t> result(count);
  for (std::size_t index = 0; index < count; ++index) {
    const int mode = frequency(index, count);
    result[index] = mode < 0
                        ? static_cast<std::size_t>(
                              mode + static_cast<int>(augmented_count))
                        : static_cast<std::size_t>(mode);
  }
  return result;
}

std::vector<std::complex<double>> imaginary_frequencies(std::size_t count) {
  checked_count(count);
  std::vector<std::complex<double>> result(count);
  for (std::size_t index = 0; index < count; ++index)
    result[index] = {0.0, static_cast<double>(frequency(index, count))};
  return result;
}

double evaluate_quadratic_taylor(
    double delta_x, double delta_y, double delta_z,
    const std::vector<double>& derivatives) {
  if (derivatives.size() != 10)
    throw std::invalid_argument(
        "quadratic 3-D Taylor evaluation requires 10 derivatives");

  const double f = derivatives[0];
  const double fx = derivatives[1];
  const double fy = derivatives[2];
  const double fz = derivatives[3];
  const double fxx = derivatives[4];
  const double fyy = derivatives[5];
  const double fzz = derivatives[6];
  const double fxy = derivatives[7];
  const double fxz = derivatives[8];
  const double fyz = derivatives[9];
  return f + fx * delta_x + fy * delta_y + fz * delta_z +
         0.5 * fxx * delta_x * delta_x +
         0.5 * fyy * delta_y * delta_y +
         0.5 * fzz * delta_z * delta_z + fxy * delta_x * delta_y +
         fxz * delta_x * delta_z + fyz * delta_y * delta_z;
}

}  // namespace coulomb::resampling
