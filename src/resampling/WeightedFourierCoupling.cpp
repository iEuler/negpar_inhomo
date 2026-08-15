#include "WeightedFourierCoupling.h"

#include <algorithm>
#include <cmath>

#include "Constants.h"

namespace coulomb::resampling {

double WeightedFourierCoupling::particleVariance(
	const std::complex<double>& coefficient, double effectiveWeight,
	int particleCount) {
	if (!(effectiveWeight > 0.0) || particleCount < 0 ||
		!std::isfinite(effectiveWeight))
		return 0.0;
	const double count = static_cast<double>(particleCount) + 1.0;
	const double cubic2Pi = 8.0 * pi * pi * pi;
	const double rho = effectiveWeight * count;
	const double variance =
		rho * rho / count -
		(cubic2Pi * cubic2Pi / count) * std::norm(coefficient);
	return std::max(0.0, std::isfinite(variance) ? variance : 0.0);
}

double WeightedFourierCoupling::optimalWeight(
	const std::complex<double>& fullCoefficient,
	const std::complex<double>& positiveCoefficient,
	const std::complex<double>& negativeCoefficient,
	double fullEffectiveWeight, double signedEffectiveWeight,
	int fullParticleCount, int positiveParticleCount, int negativeParticleCount) {
	const double fullVariance =
		particleVariance(fullCoefficient, fullEffectiveWeight, fullParticleCount);
	const double signedVariance =
		particleVariance(positiveCoefficient, signedEffectiveWeight,
						 positiveParticleCount) +
		particleVariance(negativeCoefficient, signedEffectiveWeight,
						 negativeParticleCount);
	const double denominator = fullVariance + signedVariance;
	if (!(denominator > 1e-12) || !std::isfinite(denominator))
		return 0.0;
	return std::clamp(signedVariance / denominator, 0.0, 1.0);
}

std::complex<double> WeightedFourierCoupling::blend(
	const std::complex<double>& fullResidual,
	const std::complex<double>& signedCoefficient, double weight) {
	const double boundedWeight =
		std::clamp(std::isfinite(weight) ? weight : 0.0, 0.0, 1.0);
	return boundedWeight * fullResidual +
		   (1.0 - boundedWeight) * signedCoefficient;
}

} // namespace coulomb::resampling
