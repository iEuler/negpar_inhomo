#pragma once

#include <complex>
#include <vector>

#include "fftw3.h"

namespace coulomb {

class FFT1D {
 public:
  FFT1D(size_t n);
  ~FFT1D();

  std::vector<std::complex<double>> fft(const std::vector<double>& func);
  std::vector<double> ifft(const std::vector<std::complex<double>>& funcFFT);

 private:
  size_t n_;
  fftw_complex *func_, *funcFFT_;  // func for FFT, funcFFT for inverse FFT
  fftw_plan fftPlan_, ifftPlan_;
};

class FFT3D {
 public:
  using Vector3D = std::vector<std::vector<std::vector<double>>>;
  using VectorComplex3D =
      std::vector<std::vector<std::vector<std::complex<double>>>>;
  FFT3D(size_t n1, size_t n2, size_t n3);
  ~FFT3D();

  VectorComplex3D fft(const Vector3D& func);
  Vector3D ifft(const VectorComplex3D& funcFFT);

 private:
  size_t n1_, n2_, n3_;
  fftw_complex *func_, *funcFFT_;  // func for FFT, funcFFT for inverse FFT
  fftw_plan fftPlan_, ifftPlan_;
};

}  // namespace coulomb