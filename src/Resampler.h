#pragma once
#include <complex>
#include <vector>

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

struct FourierResamplerConfig {
  double effective_particle_weight{1.0};
  size_t frequency_count{30};
  bool use_approximation{true};
  size_t max_sampling_attempts{1'000'000};
};

class FourierResampler {
 public:
  explicit FourierResampler(const NeParticleGroup& particles,
                            FourierResamplerConfig config = {});

  void reinit(const NeParticleGroup& particles) { particles_ = particles; }

  NeParticleGroup resample(RandomContext& random) const;

 private:
  NeParticleGroup particles_;
  double Neff_;
  size_t Nfreq_;
  bool useApproximation_;
  size_t augFactor_ = 2;
  size_t maxSamplingAttempts_;

  VectorComplex3D fft3d(NeParticleGroup& S_x) const;
  VectorComplex3D fft3dApprox(NeParticleGroup& S_x) const;
  std::vector<Vector3D> derivativesFromFFT(
      const VectorComplex3D& Fouriercoeff) const;

  Vector3D funcOnAugGrid(const VectorComplex3D& Fouriercoeff) const;
  Vector3D derivativesFromFFTOneTerm(const VectorComplex3D& Fouriercoeff,
                                     int orderx, int ordery, int orderz) const;

  VectorComplex3D fft3dApproxOneterm(const Vector3D& f, int orderx, int ordery,
                                     int orderz) const;
};
}  // namespace coulomb::resampling
