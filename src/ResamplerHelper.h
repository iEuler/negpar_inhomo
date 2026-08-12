#pragma once

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

class ResamplerHelper {
public:
  static VectorBool3D
  filter_fourier_coefficients(VectorComplex3D &coefficients);
  static Vector3D upper_bound(const Vector3D &values);
  static std::vector<double> values_at(const std::vector<Vector3D> &values,
                                       int x, int y, int z);
  static double
  value_from_fft(const std::vector<double> &sample,
                 const VectorComplex3D &coefficients,
                 const std::vector<std::complex<double>> &frequency_x,
                 const std::vector<std::complex<double>> &frequency_y,
                 const std::vector<std::complex<double>> &frequency_z,
                 const VectorBool3D &coefficient_mask);
  static void accept_sample(const std::vector<double> &sample,
                            NeParticleGroup &group, double value,
                            double &maximum, RandomContext &random);
};

} // namespace coulomb::resampling
