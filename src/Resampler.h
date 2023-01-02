#pragma once
#include <complex>
#include <memory>
#include <vector>

#include "Classes.h"

namespace coulomb {

class Resampler {
 public:
  Resampler(double Neff, double NeffF, size_t Nfreq, bool useApproximation,
            double dxSpace)
      : Neff_(Neff),
        NeffF_(NeffF),
        Nfreq_(Nfreq),
        useApproximation_(useApproximation),
        dxSpace_(dxSpace){};

  Resampler(const NeParticleGroup& negParGroup, double Neff = 1.0,
            double NeffF = 1.0, size_t Nfreq = 30, bool useApproximation = true,
            double dxSpace = 1.0)
      : Resampler(Neff, NeffF, Nfreq, useApproximation, dxSpace) {
    negParGroup_ = std::make_shared<NeParticleGroup>(negParGroup);
  };

  void reinit(const NeParticleGroup& negParGroup) {
    negParGroup_ = std::make_shared<NeParticleGroup>(negParGroup);
  };

  NeParticleGroup resample(bool sampleFromFullDistribution = false) const;

 private:
  std::shared_ptr<NeParticleGroup> negParGroup_;
  double Neff_, NeffF_;  // effective number for deviational particles and
                         // F particles (i.e. coarse particles)
  size_t Nfreq_;
  bool useApproximation_;
  double dxSpace_ = 1.0;  // to calculate mass from densitiy rho
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