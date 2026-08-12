#pragma once

#include <complex>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "fftw3.h"

namespace coulomb {

struct FFTWBufferDeleter {
  void operator()(fftw_complex *buffer) const noexcept;
};

struct FFTWPlanDeleter {
  void operator()(std::remove_pointer_t<fftw_plan> *plan) const noexcept;
};

using FFTWBuffer = std::unique_ptr<fftw_complex[], FFTWBufferDeleter>;
using FFTWPlan =
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, FFTWPlanDeleter>;

using Vector3D = std::vector<std::vector<std::vector<double>>>;
using VectorComplex3D =
    std::vector<std::vector<std::vector<std::complex<double>>>>;

class TensorReshape {
public:
  template <typename T>
  std::vector<std::vector<std::vector<T>>>
  to_3d(const std::vector<T> &vec1d, size_t n1, size_t n2, size_t n3) {
    if (vec1d.size() != n1 * n2 * n3)
      throw std::invalid_argument(
          "reshape1dTo3d - output vector sizes do not match input size");

    auto vec3d = std::vector(n1, std::vector(n2, std::vector<T>(n3)));
    for (size_t kk1 = 0; kk1 < n1; kk1++) {
      for (size_t kk2 = 0; kk2 < n2; kk2++) {
        for (size_t kk3 = 0; kk3 < n3; kk3++) {
          const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
          vec3d[kk1][kk2][kk3] = vec1d[kk];
        }
      }
    }
    return vec3d;
  }

  template <typename T>
  std::vector<T> to_1d(const std::vector<std::vector<std::vector<T>>> &vec3d) {
    if (vec3d.empty() || vec3d.front().empty() || vec3d.front().front().empty())
      throw std::invalid_argument("reshape3dTo1d - input vector is empty");
    const auto n1 = vec3d.size();
    const auto n2 = vec3d.front().size();
    const auto n3 = vec3d.front().front().size();

    for (const auto &plane : vec3d) {
      if (plane.size() != n2)
        throw std::invalid_argument("reshape3dTo1d - input is not rectangular");
      for (const auto &row : plane) {
        if (row.size() != n3)
          throw std::invalid_argument(
              "reshape3dTo1d - input is not rectangular");
      }
    }

    auto vec1d = std::vector<T>(n1 * n2 * n3);
    for (size_t kk1 = 0; kk1 < n1; kk1++) {
      for (size_t kk2 = 0; kk2 < n2; kk2++) {
        for (size_t kk3 = 0; kk3 < n3; kk3++) {
          const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
          vec1d[kk] = vec3d[kk1][kk2][kk3];
        }
      }
    }
    return vec1d;
  }
};

class FFT1D {
public:
  FFT1D(size_t n);
  ~FFT1D();

  FFT1D(const FFT1D &) = delete;
  FFT1D &operator=(const FFT1D &) = delete;
  FFT1D(FFT1D &&) = delete;
  FFT1D &operator=(FFT1D &&) = delete;

  std::vector<std::complex<double>> fft(const std::vector<double> &func);
  std::vector<double> ifft(const std::vector<std::complex<double>> &funcFFT);

private:
  size_t n_;
  FFTWBuffer func_, funcFFT_; // func for FFT, funcFFT for inverse FFT
  FFTWPlan fftPlan_, ifftPlan_;
};

class FFT3D {
public:
  FFT3D(size_t n1, size_t n2, size_t n3);
  ~FFT3D();

  FFT3D(const FFT3D &) = delete;
  FFT3D &operator=(const FFT3D &) = delete;
  FFT3D(FFT3D &&) = delete;
  FFT3D &operator=(FFT3D &&) = delete;

  VectorComplex3D fft(const Vector3D &func);
  Vector3D ifft(const VectorComplex3D &funcFFT);

private:
  size_t n1_, n2_, n3_;
  FFTWBuffer func_, funcFFT_; // func for FFT, funcFFT for inverse FFT
  FFTWPlan fftPlan_, ifftPlan_;
};

} // namespace coulomb
