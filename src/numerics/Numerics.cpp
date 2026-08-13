#include "Numerics.h"

#include "SimulationTypes.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace coulomb {

double Numerics::errorFunction(double x) const {
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
	const double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) *
							   t * std::exp(-x * x);
	return sign * y;
}

std::vector<double> Numerics::circularShift(const std::vector<double>& values,
											int size, int shiftDistance) const {
	if (size < 0 || static_cast<std::size_t>(size) != values.size())
		throw std::invalid_argument(
			"Numerics::circularShift size does not match values");
	if (size == 0)
		return {};
	std::vector<double> shifted(size);
	shiftDistance %= size;
	for (int index = 0; index < size; ++index) {
		int oldIndex = index + shiftDistance;
		if (oldIndex >= size)
			oldIndex -= size;
		if (oldIndex < 0)
			oldIndex += size;
		shifted[index] = values[oldIndex];
	}
	return shifted;
}

std::vector<double> Numerics::boundaryShift(const std::vector<double>& values,
											int size, int shiftDistance,
											double boundary) const {
	if (size < 0 || static_cast<std::size_t>(size) != values.size())
		throw std::invalid_argument(
			"Numerics::boundaryShift size does not match values");
	if (size == 0)
		return {};
	std::vector<double> shifted(size);
	for (int index = 0; index < size; ++index) {
		const int oldIndex = index + shiftDistance;
		shifted[index] =
			(oldIndex >= size || oldIndex < 0) ? boundary : values[oldIndex];
	}
	return shifted;
}

std::vector<double>
Numerics::centralDifference(const std::vector<double>& values, int size,
							BoundaryCondition boundary) const {
	if (size < 0 || static_cast<std::size_t>(size) != values.size())
		throw std::invalid_argument(
			"Numerics::centralDifference size does not match values");
	if (boundary != BoundaryCondition::Periodic &&
		boundary != BoundaryCondition::Reflective)
		throw std::invalid_argument(
			"Numerics::centralDifference boundary is invalid");
	if (size == 0)
		return {};
	const auto left = boundary == BoundaryCondition::Periodic
						  ? circularShift(values, size, -1)
						  : boundaryShift(values, size, -1, values[0]);
	const auto right = boundary == BoundaryCondition::Periodic
						   ? circularShift(values, size, 1)
						   : boundaryShift(values, size, 1, values[size - 1]);

	std::vector<double> difference(size);
	for (int index = 0; index < size; ++index)
		difference[index] = 0.5 * (right[index] - left[index]);
	return difference;
}

std::vector<double>
Numerics::limitedFlux(const std::vector<double>& values, int size, double dx,
					  double dt, int velocitySign,
					  const std::vector<double>& boundaryValues,
					  BoundaryCondition boundary) const {
	if (size < 0 || static_cast<std::size_t>(size) != values.size())
		throw std::invalid_argument(
			"Numerics::limitedFlux size does not match values");
	if (size == 0)
		return {};
	if (dx == 0.0)
		throw std::invalid_argument("Numerics::limitedFlux requires dx != 0");
	if (velocitySign != 1 && velocitySign != -1)
		throw std::invalid_argument(
			"Numerics::limitedFlux velocity sign must be +/-1");
	if (boundary != BoundaryCondition::Periodic &&
		boundary != BoundaryCondition::Reflective)
		throw std::invalid_argument(
			"Numerics::limitedFlux boundary is invalid");
	if (boundary == BoundaryCondition::Reflective &&
		boundaryValues.size() != values.size())
		throw std::invalid_argument(
			"Numerics::limitedFlux boundary values must match the grid size");

	const auto left =
		boundary == BoundaryCondition::Periodic
			? circularShift(values, size, -1)
			: boundaryShift(values, size, -1, boundaryValues.front());
	const auto right =
		boundary == BoundaryCondition::Periodic
			? circularShift(values, size, 1)
			: boundaryShift(values, size, 1, boundaryValues.back());

	std::vector<double> difference(size);
	for (int index = 0; index < size; ++index) {
		const double oneSidedDifference = velocitySign > 0
											  ? values[index] - left[index]
											  : right[index] - values[index];
		difference[index] = dt / dx * oneSidedDifference;
	}
	return difference;
}

void Numerics::advanceKineticEuler(const std::vector<double>& density,
								   const std::vector<double>& velocity,
								   const std::vector<double>& temperature,
								   int size, double dx, double dt,
								   BoundaryCondition boundary,
								   std::vector<double>& densityChange,
								   std::vector<double>& momentumChange,
								   std::vector<double>& energyChange) const {
	if (size < 0 || static_cast<std::size_t>(size) != density.size() ||
		velocity.size() != density.size() ||
		temperature.size() != density.size())
		throw std::invalid_argument(
			"Numerics::advanceKineticEuler input size mismatch");
	if (densityChange.size() != density.size() ||
		momentumChange.size() != density.size() ||
		energyChange.size() != density.size())
		throw std::invalid_argument(
			"Numerics::advanceKineticEuler output size mismatch");
	if (size == 0)
		return;
	if (dx == 0.0)
		throw std::invalid_argument(
			"Numerics::advanceKineticEuler requires dx != 0");

	std::vector<double> g1Positive(size), g1Negative(size);
	std::vector<double> g2Positive(size), g2Negative(size);
	std::vector<double> g3Positive(size), g3Negative(size);
	constexpr double eulerPi = 3.141592653589793238462643383279502884;
	constexpr double lambdaState = 1.0; // corresponding to Dim = 3

	for (int index = 0; index < size; ++index) {
		if (!(temperature[index] > 0.0))
			throw std::invalid_argument("Numerics::advanceKineticEuler "
										"requires positive temperature");
		const double positiveFraction =
			0.5 * (1.0 + errorFunction(velocity[index] /
									   std::sqrt(2.0 * temperature[index])));
		const double negativeFraction = 1.0 - positiveFraction;
		const double gaussianTail =
			std::sqrt(temperature[index] / (2.0 * eulerPi)) *
			std::exp(-0.5 * velocity[index] * velocity[index] /
					 temperature[index]);

		g1Positive[index] = density[index] *
							(velocity[index] * positiveFraction + gaussianTail);
		g1Negative[index] = density[index] *
							(velocity[index] * negativeFraction - gaussianTail);
		g2Positive[index] =
			density[index] *
			((temperature[index] + velocity[index] * velocity[index]) *
				 positiveFraction +
			 velocity[index] * gaussianTail);
		g2Negative[index] =
			density[index] *
			((temperature[index] + velocity[index] * velocity[index]) *
				 negativeFraction -
			 velocity[index] * gaussianTail);
		g3Positive[index] =
			density[index] *
			((1.5 * temperature[index] +
			  0.5 * velocity[index] * velocity[index]) *
				 velocity[index] * positiveFraction +
			 (temperature[index] + 0.5 * velocity[index] * velocity[index]) *
				 gaussianTail);
		g3Negative[index] =
			density[index] *
			((1.5 * temperature[index] +
			  0.5 * velocity[index] * velocity[index]) *
				 velocity[index] * negativeFraction -
			 (temperature[index] + 0.5 * velocity[index] * velocity[index]) *
				 gaussianTail);
		g3Positive[index] +=
			lambdaState * temperature[index] * g1Positive[index];
		g3Negative[index] +=
			lambdaState * temperature[index] * g1Negative[index];
	}

	const auto d1Positive =
		limitedFlux(g1Positive, size, dx, dt, 1, g1Positive, boundary);
	const auto d1Negative =
		limitedFlux(g1Negative, size, dx, dt, -1, g1Negative, boundary);
	const auto d2Positive =
		limitedFlux(g2Positive, size, dx, dt, 1, g2Positive, boundary);
	const auto d2Negative =
		limitedFlux(g2Negative, size, dx, dt, -1, g2Negative, boundary);
	const auto d3Positive =
		limitedFlux(g3Positive, size, dx, dt, 1, g3Positive, boundary);
	const auto d3Negative =
		limitedFlux(g3Negative, size, dx, dt, -1, g3Negative, boundary);
	for (int index = 0; index < size; ++index) {
		densityChange[index] = d1Positive[index] + d1Negative[index];
		momentumChange[index] = d2Positive[index] + d2Negative[index];
		energyChange[index] = d3Positive[index] + d3Negative[index];
	}
}

} // namespace coulomb
