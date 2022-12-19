#include "FFT.h"

namespace coulomb {

FFT1D::FFT1D(size_t n) : n_(n) {
  func_ = new fftw_complex[n_];
  funcFFT_ = new fftw_complex[n_];

  fftPlan_ = fftw_plan_dft_1d(n_, func_, funcFFT_, FFTW_FORWARD, FFTW_ESTIMATE);
  ifftPlan_ =
      fftw_plan_dft_1d(n_, funcFFT_, func_, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFT1D::~FFT1D() {
  fftw_destroy_plan(fftPlan_);
  fftw_destroy_plan(ifftPlan_);
  delete[] func_;
  delete[] funcFFT_;
}

std::vector<std::complex<double>> FFT1D::fft(const std::vector<double>& func) {
  for (int kx = 0; kx < n_; kx++) {
    func_[kx][0] = func[kx];
    func_[kx][1] = 0;
  }

  fftw_execute(fftPlan_);

  std::vector<std::complex<double>> funcFFT(n_);
  for (int kx = 0; kx < n_; kx++) {
    funcFFT[kx] = {funcFFT_[kx][0], funcFFT_[kx][1]};
  }
  return funcFFT;
}

std::vector<double> FFT1D::ifft(
    const std::vector<std::complex<double>>& funcFFT) {
  for (int kx = 0; kx < n_; kx++) {
    funcFFT_[kx][0] = funcFFT[kx].real();
    funcFFT_[kx][1] = funcFFT[kx].imag();
  }

  fftw_execute(ifftPlan_);

  std::vector<double> func(n_);
  for (int kx = 0; kx < n_; kx++) {
    func[kx] = func_[kx][0];
  }
  return func;
}

FFT3D::FFT3D(size_t n1, size_t n2, size_t n3) : n1_(n1), n2_(n2), n3_(n3) {
  fftw_complex temp;
  func_ = (fftw_complex*)fftw_malloc(n1_ * n2_ * n3_ * sizeof(temp));
  funcFFT_ = (fftw_complex*)fftw_malloc(n1_ * n2_ * n3_ * sizeof(temp));
  // TODO -  check whether we can use the following instead
  // func_ = new fftw_complex[n_*n_*n_];

  fftPlan_ = fftw_plan_dft_3d(n1_, n2_, n3_, func_, funcFFT_, FFTW_FORWARD,
                              FFTW_ESTIMATE);
  ifftPlan_ = fftw_plan_dft_3d(n1_, n2_, n3_, funcFFT_, func_, FFTW_BACKWARD,
                               FFTW_ESTIMATE);
}

FFT3D::~FFT3D() {
  fftw_destroy_plan(fftPlan_);
  fftw_destroy_plan(ifftPlan_);
  fftw_free(func_);
  fftw_free(funcFFT_);
}

FFT3D::VectorComplex3D FFT3D::fft(const Vector3D& func) {
  // the (i,j,k)-th element of the array with size (Nx,Ny,Nz), you would use the
  // expression an_array[k + Nz * (j + Ny * i)].

  for (size_t kk1 = 0; kk1 < n1_; kk1++) {
    for (size_t kk2 = 0; kk2 < n2_; kk2++) {
      for (size_t kk3 = 0; kk3 < n3_; kk3++) {
        const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
        func_[kk][0] = func[kk1][kk2][kk3];
        func_[kk][1] = 0;
      }
    }
  }

  fftw_execute(fftPlan_);
  VectorComplex3D funcFFT = std::vector(
      n1_, std::vector(n2_, std::vector<std::complex<double>>(n3_)));

  for (size_t kk1 = 0; kk1 < n1_; kk1++) {
    for (size_t kk2 = 0; kk2 < n2_; kk2++) {
      for (size_t kk3 = 0; kk3 < n3_; kk3++) {
        const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
        funcFFT[kk1][kk2][kk3] = {funcFFT_[kk][0], funcFFT_[kk][1]};
      }
    }
  }

  return funcFFT;
}

FFT3D::Vector3D FFT3D::ifft(const VectorComplex3D& funcFFT) {
  for (size_t kk1 = 0; kk1 < n1_; kk1++) {
    for (size_t kk2 = 0; kk2 < n2_; kk2++) {
      for (size_t kk3 = 0; kk3 < n3_; kk3++) {
        const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
        funcFFT_[kk][0] = funcFFT[kk1][kk2][kk3].real();
        funcFFT_[kk][1] = funcFFT[kk1][kk2][kk3].imag();
      }
    }
  }

  fftw_execute(ifftPlan_);

  Vector3D func = std::vector(n1_, std::vector(n2_, std::vector<double>(n3_)));
  for (size_t kk1 = 0; kk1 < n1_; kk1++) {
    for (size_t kk2 = 0; kk2 < n2_; kk2++) {
      for (size_t kk3 = 0; kk3 < n3_; kk3++) {
        const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
        func[kk1][kk2][kk3] = func_[kk][0];
      }
    }
  }

  return func;
}

}  // namespace coulomb