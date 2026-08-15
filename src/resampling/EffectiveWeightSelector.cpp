#include "EffectiveWeightSelector.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace coulomb {
namespace {

double clampFinite(double value, double lower, double upper, double fallback) {
	if (!std::isfinite(value))
		return fallback;
	return std::clamp(value, lower, upper);
}

double positiveOr(double value, double fallback) {
	return std::isfinite(value) && value > 0.0 ? value : fallback;
}

} // namespace

double EffectiveWeightSelector::interpolateQuadratic(const double* x,
													 const double* y,
													 double xValue) {
	double value = 0.0;
	for (int index = 0; index < 3; ++index) {
		double basis = 1.0;
		for (int other = 0; other < 3; ++other) {
			if (other != index)
				basis *= (xValue - x[other]) / (x[index] - x[other]);
		}
		value += y[index] * basis;
	}
	return value;
}

EffectiveWeightSelection EffectiveWeightSelector::select(
	double currentSignedWeight, double currentFullWeight, double signedWeightMin,
	double signedWeightMax, double fullWeightMin, double fullWeightMax,
	double cpuCostConstant, double cpuCostCollisionCoefficient, double dt,
	double collisionCoefficient, const std::vector<EffectiveWeightCell>& cells) const {
	const auto boundedSigned = [&](double value) {
		return clampFinite(value, signedWeightMin, signedWeightMax,
						   currentSignedWeight);
	};
	const auto boundedFull = [&](double value) {
		return clampFinite(value, fullWeightMin, fullWeightMax, currentFullWeight);
	};
	if (!(currentSignedWeight > 0.0) || !(currentFullWeight > 0.0) ||
		!(signedWeightMin > 0.0) || signedWeightMin > signedWeightMax ||
		!(fullWeightMin > 0.0) || fullWeightMin > fullWeightMax ||
		cells.empty()) {
		return {boundedSigned(currentSignedWeight),
				boundedFull(currentFullWeight)};
	}

	double minimumOmega = 1.0;
	double maximumOmega = 0.0;
	double fullDensity = 0.0;
	double signedDensity = 0.0;
	double maximumCountRatio = 1.0;
	bool hasFiniteCell = false;
	for (const auto& cell : cells) {
		if (!std::isfinite(cell.omega))
			continue;
		const double omega = std::clamp(cell.omega, 0.0, 1.0);
		minimumOmega = std::min(minimumOmega, omega);
		maximumOmega = std::max(maximumOmega, omega);
		if (cell.fullParticleCount > 0 && cell.signedParticleCount > 0)
			maximumCountRatio = std::max(
				maximumCountRatio,
				static_cast<double>(cell.signedParticleCount) /
					static_cast<double>(cell.fullParticleCount));
		signedDensity +=
			currentSignedWeight *
			static_cast<double>(std::max(0, cell.signedParticleCount));
		fullDensity +=
			currentFullWeight *
			static_cast<double>(std::max(0, cell.fullParticleCount));
		hasFiniteCell = true;
	}
	if (!hasFiniteCell || !(signedDensity > 0.0) || !(fullDensity >= 0.0))
		return {boundedSigned(currentSignedWeight),
				boundedFull(currentFullWeight)};

	const double safeOneMinusMin = std::max(1e-12, 1.0 - minimumOmega);
	const double safeOneMinusMax = std::max(1e-12, 1.0 - maximumOmega);
	const double minimumSlope =
		currentFullWeight / currentSignedWeight * minimumOmega /
		safeOneMinusMin;
	const double maximumSlope =
		currentFullWeight / currentSignedWeight * maximumOmega /
		safeOneMinusMax;
	const double alphaCpu =
		1.0 + std::max(0.0, cpuCostConstant) +
		std::max(0.0, cpuCostCollisionCoefficient) * std::max(0.0, dt) *
			std::max(0.0, collisionCoefficient);
	const double objectiveSlope =
		0.5 / positiveOr(alphaCpu, 1.0) * fullDensity /
		positiveOr(signedDensity, 1.0);

	double x[3]{};
	double y[3]{};
	x[0] = 1.0 / fullWeightMax;
	const double minimumFullInverse =
		1.0 / currentFullWeight -
		(1.0 / fullWeightMax - 1.0 / currentSignedWeight) /
			positiveOr(minimumSlope, 1e-12);
	y[0] = positiveOr(minimumFullInverse, 1.0 / fullWeightMin);
	x[1] = 1.0 / currentSignedWeight;
	y[1] = 1.0 / currentFullWeight;
	const double maximumFullDensity =
		std::max(maximumCountRatio, 1e-12);
	const double maximumFullInverse =
		(1.0 / currentFullWeight +
		 1.0 / currentSignedWeight / positiveOr(maximumSlope, 1e-12)) /
		(1.0 + 1.0 / maximumFullDensity /
				   positiveOr(maximumSlope, 1e-12));
	y[2] = positiveOr(maximumFullInverse, 1.0 / fullWeightMax);
	x[2] = 1.0 / positiveOr(1.0 / y[2] * maximumFullDensity,
							 currentSignedWeight);

	double targetX = x[1];
	double targetY = y[1];
	const double coefficient =
		2.0 * y[0] / ((x[0] - x[1]) * (x[0] - x[2])) +
		2.0 * y[1] / ((x[1] - x[0]) * (x[1] - x[2])) +
		2.0 * y[2] / ((x[2] - x[0]) * (x[2] - x[1]));
	const double constant =
		-y[0] * (x[1] + x[2]) / ((x[0] - x[1]) * (x[0] - x[2])) -
		y[1] * (x[0] + x[2]) / ((x[1] - x[0]) * (x[1] - x[2])) -
		y[2] * (x[0] + x[1]) / ((x[2] - x[0]) * (x[2] - x[1]));
	if (std::isfinite(coefficient) && std::abs(coefficient) > 1e-12 &&
		std::isfinite(objectiveSlope) && objectiveSlope > 1e-12) {
		targetX = (-(1.0 / objectiveSlope) - constant) / coefficient;
		targetX = std::clamp(targetX,
							 std::min({x[0], x[1], x[2]}),
							 std::max({x[0], x[1], x[2]}));
		targetY = interpolateQuadratic(x, y, targetX);
	}

	double proposedSigned = 1.0 / positiveOr(targetX, 1.0 / currentSignedWeight);
	double proposedFull = 1.0 / positiveOr(targetY, 1.0 / currentFullWeight);
	proposedSigned = boundedSigned(proposedSigned);
	proposedFull = boundedFull(proposedFull);

	// A cell must retain enough full particles for the proposed signed weight.
	for (const auto& cell : cells) {
		if (cell.signedParticleCount <= 0 || !(cell.macroMass > 0.0) ||
			!std::isfinite(cell.macroMass))
			continue;
		const double expectedSignedCount =
			cell.signedParticleCount * currentSignedWeight / proposedSigned;
		const double stability =
			cell.macroMass / std::max(1.0, 1.1 * expectedSignedCount);
		if (std::isfinite(stability) && stability > currentFullWeight)
			proposedFull = std::min(proposedFull, stability);
	}
	proposedFull = std::max(proposedFull, currentFullWeight);
	proposedFull = boundedFull(proposedFull);
	return {proposedSigned, proposedFull};
}

} // namespace coulomb
