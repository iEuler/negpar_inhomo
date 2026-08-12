#pragma once

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

class ResamplerHelper {
  public:
	explicit ResamplerHelper(RandomContext& random) : random_(random) {}

	VectorBool3D filter_fourier_coefficients(VectorComplex3D& coefficients);
	Vector3D upper_bound(const Vector3D& values);
	std::vector<double> values_at(const std::vector<Vector3D>& values, int x,
								  int y, int z);
	double value_from_fft(const std::vector<double>& sample,
						  const VectorComplex3D& coefficients,
						  const std::vector<std::complex<double>>& frequency_x,
						  const std::vector<std::complex<double>>& frequency_y,
						  const std::vector<std::complex<double>>& frequency_z,
						  const VectorBool3D& coefficient_mask);
	void accept_sample(const std::vector<double>& sample,
					   NeParticleGroup& group, double value, double& maximum);

  private:
	RandomContext& random_;
};

} // namespace coulomb::resampling
