#pragma once
#include <complex>
#include <memory>
#include <vector>

#include "Classes.h"

namespace coulomb {

class Resampler {
 public:
  void reinit(const NeParticleGroup& negParGroup) {
    negParGroup_ = std::make_shared<NeParticleGroup>(negParGroup);
  };

  std::shared_ptr<NeParticleGroup> resample();

  void resetNeff(double Neff, double Neff_F) {
    grid_.Neff = Neff;
    grid_.Neff_F = Neff_F;
  };

 private:
  std::shared_ptr<NeParticleGroup> negParGroup_;
  NumericGridClass grid_;
  ParaClass para_;
  bool useApproximation_;
  double Neff_;
  size_t Nfreq_;
  size_t augFactor_ = 2;

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
}  // namespace coulomb