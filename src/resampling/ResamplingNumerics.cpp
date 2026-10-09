#include "ResamplingNumerics.h"

#include <limits>
#include <cmath>
#include <stdexcept>
#include "Constants.h"
#include "RandomSampling.h"

namespace coulomb::resampling {
int StratifiedProposalAllocator::allocate(double expected, RandomContext& random) {
	if (!std::isfinite(expected) || expected < 0.0 || expected >= std::numeric_limits<int>::max())
		throw std::invalid_argument("Invalid stratified proposal expectation");
	if (expected == 0.0) return 0;
	const double end = cumulative + expected;
	// Keep unit strata representable and detect lost increments.
	if (!std::isfinite(end) || end >= 0x1p52 || end == cumulative)
		throw std::overflow_error("Stratified proposal accumulation overflow");
	if (!initialized) {
		boundaryUniform = RandomSampling(random).uniform();
		initialized = true;
	}
	const double oldFloor = std::floor(cumulative), newFloor = std::floor(end);
	if (newFloor != oldFloor) boundaryUniform = RandomSampling(random).uniform();
	const bool included = boundaryUniform < end - newFloor;
	const double count = newFloor - oldFloor - (boundaryIncluded ? 1.0 : 0.0) + (included ? 1.0 : 0.0);
	if (count < 0.0 || count > std::numeric_limits<int>::max())
		throw std::overflow_error("Stratified proposal count overflow");
	cumulative = end;
	boundaryIncluded = included;
	return static_cast<int>(count);
}


std::vector<double> ResamplingNumerics::wrapPeriodicSample(const std::vector<double>& sample) {
	if (sample.size() != 3) throw std::invalid_argument("Periodic sample requires three coordinates");
	auto result = sample;
	for (double& v : result) {
		if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite periodic sample");
		v = std::fmod(v, 2.0*pi);
		if (v < 0.0) v += 2.0*pi;
		if (v >= 2.0*pi) v = 0.0;
	}
	return result;
}

double ResamplingNumerics::quadraticEnvelope(double h, const std::vector<double>& d) {
	if (d.size() != 10 || !std::isfinite(h) || h < 0.0)
		throw std::invalid_argument("Invalid quadratic envelope input");
	for (double v : d)
		if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite quadratic derivative");
	const double bound = std::abs(d[0]) + h * (std::abs(d[1]) + std::abs(d[2]) + std::abs(d[3]))
		+ h*h * (0.5*(std::abs(d[4])+std::abs(d[5])+std::abs(d[6]))
			+ std::abs(d[7])+std::abs(d[8])+std::abs(d[9]));
	if (!std::isfinite(bound)) throw std::invalid_argument("Quadratic envelope overflow");
	const double padded = bound * (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
	if (!std::isfinite(padded)) throw std::invalid_argument("Quadratic envelope padding overflow");
	return padded;
}

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
