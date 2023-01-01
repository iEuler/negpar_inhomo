#pragma once
#include <complex>
#include <memory>
#include <vector>

#include "Classes.h"
#include "Resampler.h"

namespace coulomb {

class InhomoResampler {
 public:
  InhomoResampler(double Neff, size_t Nfreq, bool useApproximation)
      : Neff_(Neff), Nfreq_(Nfreq), useApproximation_(useApproximation){};
  void reinit(const NeParticleGroup& negParGroup) {
    negParGroup_ = std::make_shared<NeParticleGroup>(negParGroup);
  };

  void resample(std::vector<NeParticleGroup>& S_x);

 private:
  NumericGridClass grid_;
  ParaClass para_;
  std::shared_ptr<NeParticleGroup> negParGroup_;
  double Neff_;
  size_t Nfreq_;
  bool useApproximation_;
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