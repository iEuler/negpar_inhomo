#include "ResamplingNumerics.h"

#include <limits>
#include <stdexcept>

namespace coulomb::resampling {
namespace {

int checkedCount(std::size_t count) {
	if (count == 0 ||
		count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
		throw std::invalid_argument(
			"resampling frequency count must fit in a positive int");
	return static_cast<int>(count);
}

int frequency(std::size_t index, std::size_t count) {
	const int countInt = checkedCount(count);
	if (index >= count)
		throw std::out_of_range(
			"resampling frequency index is outside the grid");
	const int indexInt = static_cast<int>(index);
	return indexInt >= countInt / 2 + 1 ? indexInt - countInt : indexInt;
}

} // namespace

std::vector<double> ResamplingNumerics::frequencies(std::size_t count) {
	checkedCount(count);
	std::vector<double> result(count);
	for (std::size_t index = 0; index < count; ++index)
		result[index] = static_cast<double>(frequency(index, count));
	return result;
}

std::vector<std::size_t>
ResamplingNumerics::augmentedLocations(std::size_t count,
									   std::size_t augmentationFactor) {
	checkedCount(count);
	if (augmentationFactor == 0 ||
		count > static_cast<std::size_t>(std::numeric_limits<int>::max()) /
					augmentationFactor)
		throw std::invalid_argument(
			"resampling augmentation must produce a positive int-sized grid");

	const auto augmentedCount = count * augmentationFactor;
	std::vector<std::size_t> result(count);
	for (std::size_t index = 0; index < count; ++index) {
		const int mode = frequency(index, count);
		result[index] = mode < 0 ? static_cast<std::size_t>(
									   mode + static_cast<int>(augmentedCount))
								 : static_cast<std::size_t>(mode);
	}
	return result;
}

std::vector<std::complex<double>>
ResamplingNumerics::imaginaryFrequencies(std::size_t count) {
	checkedCount(count);
	std::vector<std::complex<double>> result(count);
	for (std::size_t index = 0; index < count; ++index)
		result[index] = {0.0, static_cast<double>(frequency(index, count))};
	return result;
}

double ResamplingNumerics::evaluateQuadraticTaylor(
	double deltaX, double deltaY, double deltaZ,
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
	return f + fx * deltaX + fy * deltaY + fz * deltaZ +
		   0.5 * fxx * deltaX * deltaX + 0.5 * fyy * deltaY * deltaY +
		   0.5 * fzz * deltaZ * deltaZ + fxy * deltaX * deltaY +
		   fxz * deltaX * deltaZ + fyz * deltaY * deltaZ;
}

} // namespace coulomb::resampling
