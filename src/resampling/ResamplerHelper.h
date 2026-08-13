#pragma once

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

class ResamplerHelper {
  public:
	explicit ResamplerHelper(RandomContext& random) : randomContext(random) {}

	VectorBool3D filterFourierCoefficients(VectorComplex3D& coefficients);
	Vector3D upperBound(const Vector3D& values);
	std::vector<double> valuesAt(const std::vector<Vector3D>& values, int x,
								 int y, int z);
	double valueFromFft(const std::vector<double>& sample,
						const VectorComplex3D& coefficients,
						const std::vector<std::complex<double>>& frequencyX,
						const std::vector<std::complex<double>>& frequencyY,
						const std::vector<std::complex<double>>& frequencyZ,
						const VectorBool3D& coefficientMask);
	void acceptSample(const std::vector<double>& sample, NeParticleGroup& group,
					  double value, double& maximum);

  private:
	RandomContext& randomContext;
};

} // namespace coulomb::resampling
