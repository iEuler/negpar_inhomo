#pragma once

#include <complex>

namespace coulomb::resampling {

class WeightedFourierCoupling {
  public:
	static double particleVariance(const std::complex<double>& coefficient,
								   double effectiveWeight,
								   int particleCount);
	static double optimalWeight(const std::complex<double>& fullCoefficient,
								const std::complex<double>& positiveCoefficient,
								const std::complex<double>& negativeCoefficient,
								double fullEffectiveWeight,
								double signedEffectiveWeight,
								int fullParticleCount,
								int positiveParticleCount,
								int negativeParticleCount);
	static std::complex<double>
	blend(const std::complex<double>& fullResidual,
		  const std::complex<double>& signedCoefficient, double weight);
};

} // namespace coulomb::resampling
