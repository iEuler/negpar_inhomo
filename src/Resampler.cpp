#include "Resampler.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "FFT.h"
#include "ParticleGroupOperations.h"
#include "ResamplerHelper.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"
#include "Constants.h"
#include "RandomSampling.h"

namespace coulomb::resampling {

using std::abs;
using std::complex;
using std::exp;
using std::max;
using std::to_string;

FourierResampler::FourierResampler(const NeParticleGroup& particles,
                                   FourierResamplerConfig config)
    : particles_(particles),
      Neff_(config.effective_particle_weight),
      Nfreq_(config.frequency_count),
      useApproximation_(config.use_approximation),
      maxSamplingAttempts_(config.max_sampling_attempts) {
  if (!(Neff_ > 0.0))
    throw std::invalid_argument(
        "Resampler effective particle weight must be positive");
  if (Nfreq_ < 2 || Nfreq_ % 2 != 0)
    throw std::invalid_argument(
        "Resampler frequency count must be an even value of at least 2");
  if (Nfreq_ > static_cast<size_t>(std::numeric_limits<int>::max()) /
                   augFactor_)
    throw std::invalid_argument("Resampler frequency grid is too large");
  if (maxSamplingAttempts_ == 0)
    throw std::invalid_argument(
        "Resampler sampling-attempt budget must be positive");
}

Vector3D FourierResampler::funcOnAugGrid(
    const VectorComplex3D &Fouriercoeff) const {
  auto fft3d =
      FFT3D(Nfreq_ * augFactor_, Nfreq_ * augFactor_, Nfreq_ * augFactor_);
  return fft3d.ifft(Fouriercoeff);
}

Vector3D FourierResampler::derivativesFromFFTOneTerm(
    const VectorComplex3D &Fouriercoeff, int orderx, int ordery,
    int orderz) const {
  const auto augmented_size = augFactor_ * Nfreq_;
  auto FSaug = VectorComplex3D(
      augmented_size,
      std::vector(augmented_size,
                  std::vector<std::complex<double>>(augmented_size)));

  // 1i *freq
  const auto freq1 = resampling::frequencies(Nfreq_);
  const auto freq2 = resampling::frequencies(Nfreq_);
  const auto freq3 = resampling::frequencies(Nfreq_);

  const auto loc1 = resampling::augmented_locations(Nfreq_, augFactor_);
  const auto loc2 = resampling::augmented_locations(Nfreq_, augFactor_);
  const auto loc3 = resampling::augmented_locations(Nfreq_, augFactor_);

  for (int kk1 = 0; kk1 < Nfreq_; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq_; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq_; kk3++) {
        size_t kk1aug = loc1[kk1];
        size_t kk2aug = loc2[kk2];
        size_t kk3aug = loc3[kk3];

        double freq = 1.0;
        for (int kx = 0; kx < orderx; kx++) freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++) freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++) freq *= freq3[kk3];

        FSaug[kk1aug][kk2aug][kk3aug] = freq * Fouriercoeff[kk1][kk2][kk3];
      }
    }
  }

  return funcOnAugGrid(FSaug);
}

std::vector<Vector3D> FourierResampler::derivativesFromFFT(
    const VectorComplex3D &Fouriercoeff) const {
  const auto orders = std::vector<std::vector<int>>{
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
      {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

  std::vector<Vector3D> fDerivatives;
  for (const auto &order : orders) {
    fDerivatives.push_back(
        derivativesFromFFTOneTerm(Fouriercoeff, order[0], order[1], order[2]));
  }
  return fDerivatives;
}

VectorComplex3D FourierResampler::fft3dApproxOneterm(
    const Vector3D &f, int orderx, int ordery, int orderz) const {
  const auto n = augFactor_ * Nfreq_;
  auto Fouriercoeff =
      std::vector(n, std::vector(n, std::vector<std::complex<double>>(n)));

  auto fft3d = FFT3D(n, n, n);

  const auto FSaug = fft3d.fft(f);

  // 1i *freq
  const auto freq1 = resampling::frequencies(Nfreq_);
  const auto freq2 = resampling::frequencies(Nfreq_);
  const auto freq3 = resampling::frequencies(Nfreq_);

  const auto loc1 = resampling::augmented_locations(Nfreq_, augFactor_);
  const auto loc2 = resampling::augmented_locations(Nfreq_, augFactor_);
  const auto loc3 = resampling::augmented_locations(Nfreq_, augFactor_);

  for (size_t kk1 = 0; kk1 < Nfreq_; kk1++) {
    for (size_t kk2 = 0; kk2 < Nfreq_; kk2++) {
      for (size_t kk3 = 0; kk3 < Nfreq_; kk3++) {
        size_t kk1aug = loc1[kk1];
        size_t kk2aug = loc2[kk2];
        size_t kk3aug = loc3[kk3];

        double freq = Neff_;
        for (int kx = 0; kx < orderx; kx++) freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++) freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++) freq *= freq3[kk3];

        if ((orderx + ordery + orderz) == 2) {
          if ((orderx == 2) || (ordery == 2) || (orderz == 2))
            freq *= -.5;
          else
            freq *= -1;
        }

        const auto augOne = FSaug[kk1aug][kk2aug][kk3aug];
        if ((orderx + ordery + orderz) == 1) {
          Fouriercoeff[kk1][kk2][kk3] =
              freq * complex<double>(augOne.imag(), -augOne.real());
        } else {
          Fouriercoeff[kk1][kk2][kk3] = freq * augOne;
        }
      }
    }
  }

  return Fouriercoeff;
}

NeParticleGroup FourierResampler::resample(RandomContext& random) const {
  NeParticleGroup S_x_new;
  auto S_x = particles_;

  const auto parTypes = std::vector<ParticleKind>{ParticleKind::Positive,
                                                  ParticleKind::Negative};

  /* Normalize particle velocity to [0 2*pi] */
  S_x.set_xyzrange({ParticleKind::Positive, ParticleKind::Negative});
  auto S_x_renormalized = resampling::normalize_signed_velocities(S_x);

  /* Prepare the grids in physical space and frequence space */
  // double dx = 2.0*pi/Nfreq;

  const auto ifreq = resampling::imaginary_frequencies(Nfreq_);
  std::vector<double> interp_x(Nfreq_);
  for (int kx = 0; kx < Nfreq_; kx++) interp_x[kx] = kx * 2 * pi / Nfreq_;

  /* Compute the Fourier coefficient */
  VectorComplex3D Fouriercoeff = useApproximation_
                                     ? fft3dApprox(S_x_renormalized)
                                     : fft3d(S_x_renormalized);

  // Apply the filter on Fourier coefficients
  auto flag_Fouriercoeff = filterFourierCoeff(Fouriercoeff);

  // cout << " F coeff computed " << endl;

  /* Compute a coarse interpolation in physical space */
  //  const auto fcoarse = interp3d_fcoarse(Fouriercoeff, Nfreq, Nfreq, Nfreq);

  auto fDerivatives = derivativesFromFFT(Fouriercoeff);

  /* evaluate the upperbound of f */
  const auto f_up = upperBoundFunc(fDerivatives[0]);

  /* refined x grid */
  double dxaug = 2.0 * pi / Nfreq_ / augFactor_;
  std::vector<double> interp_xaug(Nfreq_ * augFactor_);
  for (int kx = 0; kx < Nfreq_ * augFactor_; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq_ / augFactor_;

  const auto &f = fDerivatives[0];
  /* create a NeParticleGroup to host the P and N particles in current cell */

  /* Start sampling */

  size_t sampling_attempts = 0;

  for (int kx = 0; kx < augFactor_ * Nfreq_; kx++) {
    for (int ky = 0; ky < augFactor_ * Nfreq_; ky++) {
      for (int kz = 0; kz < augFactor_ * Nfreq_; kz++) {
        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        double fcc = f_up[kx][ky][kz];

        if (fcc < std::abs(f[kx][ky][kz])) throw std::exception("small bound!");

        double maxf = 1.5 * fcc;
        int N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff_, random);

        int k_virtual = 0;
        NeParticleGroup S_x_incell;
        const auto fDeriv = getValuesByLoc(fDerivatives, kx, ky, kz);

        while (k_virtual < N_incell) {
          if (++sampling_attempts > maxSamplingAttempts_)
            throw std::runtime_error(
                "Resampler exceeded its sampling-attempt budget at cell (" +
                std::to_string(kx) + ", " + std::to_string(ky) + ", " +
                std::to_string(kz) + ") with target " +
                std::to_string(N_incell) + " and envelope " +
                std::to_string(maxf));
          // create a particle in the cell
          // Sample offsets from the explicit RandomContext below.
          double deltax = myrand(random) * dxaug - 0.5 * dxaug;
          double deltay = myrand(random) * dxaug - 0.5 * dxaug;
          double deltaz = myrand(random) * dxaug - 0.5 * dxaug;
          std::vector<double> Sf{xc + deltax, yc + deltay, zc + deltaz};

          // compute f at this point
          double fval =
              useApproximation_
                  ? resampling::evaluate_quadratic_taylor(
                        deltax, deltay, deltaz, fDeriv)
                  : fvalueFromFFT(Sf, Fouriercoeff, ifreq, ifreq, ifreq,
                                  flag_Fouriercoeff);

          // reset current cell if fval>maxf, otherwise continue sampling
          // in current cell
          acceptSampled(Sf, S_x_incell, fval, maxf, random);

          // reset N_incell if maxf is changed
          N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff_, random);
          k_virtual++;
        }

        ::coulomb::mergeNeParticleGroup(S_x_new, S_x_incell, parTypes);
      }
    }
  }

  // cout << "Resampled." << endl;

  // rescale to the original coordinates

  const auto &xyz_minmax = S_x.xyz_minmax;
  for (const auto parType : parTypes) {
    auto &Sp_sampled = S_x_new.list(parType);
    resampling::restore_velocities(Sp_sampled, xyz_minmax);
  }

  // cout << "Rescaled." << endl;

  return S_x_new;
  // return std::make_shared<NeParticleGroup>(S_x_new);
}

VectorComplex3D FourierResampler::fft3d(NeParticleGroup &S_x) const {
  auto Fouriercoeff = std::vector(
      Nfreq_, std::vector(Nfreq_, std::vector<std::complex<double>>(Nfreq_)));

  int Np = S_x.size(ParticleKind::Positive);
  int Nn = S_x.size(ParticleKind::Negative);

  auto &Sp = S_x.list(ParticleKind::Positive);
  auto &Sn = S_x.list(ParticleKind::Negative);

  double Neff = Neff_;

  // double Neff_temp = 1./Np;

  const auto ifreq1 = resampling::imaginary_frequencies(Nfreq_);
  const auto ifreq2 = resampling::imaginary_frequencies(Nfreq_);
  const auto ifreq3 = resampling::imaginary_frequencies(Nfreq_);

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi * Nfreq_ * Nfreq_ * Nfreq_;
  double maxFS = 0.0;

  // the (i,j,k)-th element of the array with size (Nx,Ny,Nz), you would use the
  // expression an_array[k + Nz * (j + Ny * i)].

  for (int kk1 = 0; kk1 < Nfreq_; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq_; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq_; kk3++) {
        auto coeff = complex<double>(0., 0.);
        for (int kp = 0; kp < Np; kp++) {
          auto &vp = Sp[kp].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          coeff += exp(-expterm);
        }
        for (int kn = 0; kn < Nn; kn++) {
          auto &vp = Sn[kn].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          coeff -= exp(-expterm);
        }
        // Fouriercoeff[kk] *= Neff * coeff_fft;
        coeff *= Neff * coeff_fft / (Nfreq_ * Nfreq_ * Nfreq_);
        maxFS = max(maxFS, abs(coeff));
        Fouriercoeff[kk1][kk2][kk3] = coeff;
      }
    }
  }

  return Fouriercoeff;
}

VectorComplex3D FourierResampler::fft3dApprox(NeParticleGroup &S_x) const {
  auto fouriercoeff =
      std::vector(Nfreq_, std::vector(Nfreq_, std::vector<std::complex<double>>(
                                                  Nfreq_, {0.0, 0.0})));

  size_t augFactor = augFactor_;

  int Np = S_x.size(ParticleKind::Positive);
  int Nn = S_x.size(ParticleKind::Negative);

  auto &Sp = S_x.list(ParticleKind::Positive);
  auto &Sn = S_x.list(ParticleKind::Negative);

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi;

  // for the (i,j,k)-th element of the array with size (Nx,Ny,Nz), use the
  // expression an_array[k + Nz * (j + Ny * i)].

  // create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
  double dx = 2.0 * pi / augFactor / Nfreq_;
  double dy = 2.0 * pi / augFactor / Nfreq_;
  double dz = 2.0 * pi / augFactor / Nfreq_;

  const auto n = augFactor_ * Nfreq_;
  Vector3D f = std::vector(n, std::vector(n, std::vector<double>(n, 0.0)));
  auto fx = f;
  auto fy = f;
  auto fz = f;
  auto fxx = f;
  auto fyy = f;
  auto fzz = f;
  auto fxy = f;
  auto fxz = f;
  auto fyz = f;

  auto sizeF = n * n * n;

  for (int kp = 0; kp < Np; kp++) {
    double x0 = Sp[kp].velocity(0);
    double y0 = Sp[kp].velocity(1);
    double z0 = Sp[kp].velocity(2);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= n) xloc--;
    if (yloc >= n) yloc--;
    if (zloc >= n) zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    size_t loc = zloc + n * (yloc + n * xloc);
    if ((loc >= sizeF) || (loc < 0)) {
      std::string errMsg =
          "error: in approximation. Particle moved out of range. x =  (" +
          to_string(x0) + ", " + to_string(y0) + ", " + to_string(z0) +
          "), dx = " + to_string(dx);
      throw std::exception(errMsg.c_str());
    }

    f[xloc][yloc][zloc]++;
    fx[xloc][yloc][zloc] += xdelta;
    fy[xloc][yloc][zloc] += ydelta;
    fz[xloc][yloc][zloc] += zdelta;
    fxx[xloc][yloc][zloc] += xdelta * xdelta;
    fyy[xloc][yloc][zloc] += ydelta * ydelta;
    fzz[xloc][yloc][zloc] += zdelta * zdelta;
    fxy[xloc][yloc][zloc] += xdelta * ydelta;
    fyz[xloc][yloc][zloc] += ydelta * zdelta;
    fxz[xloc][yloc][zloc] += zdelta * xdelta;
  }

  // cout << "Approx 2" << endl;

  for (int kp = 0; kp < Nn; kp++) {
    double x0 = Sn[kp].velocity(0);
    double y0 = Sn[kp].velocity(1);
    double z0 = Sn[kp].velocity(2);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= n) xloc--;
    if (yloc >= n) yloc--;
    if (zloc >= n) zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    size_t loc = zloc + n * (yloc + n * xloc);
    if ((loc >= sizeF) || (loc < 0)) {
      std::string errMsg =
          "error: in approximation. Particle moved out of range. x =  (" +
          to_string(x0) + ", " + to_string(y0) + ", " + to_string(z0) +
          "), dx = " + to_string(dx);
      throw std::exception(errMsg.c_str());
    }

    f[xloc][yloc][zloc]--;
    fx[xloc][yloc][zloc] -= xdelta;
    fy[xloc][yloc][zloc] -= ydelta;
    fz[xloc][yloc][zloc] -= zdelta;
    fxx[xloc][yloc][zloc] -= xdelta * xdelta;
    fyy[xloc][yloc][zloc] -= ydelta * ydelta;
    fzz[xloc][yloc][zloc] -= zdelta * zdelta;
    fxy[xloc][yloc][zloc] -= xdelta * ydelta;
    fyz[xloc][yloc][zloc] -= ydelta * zdelta;
    fxz[xloc][yloc][zloc] -= zdelta * xdelta;
  }

  const auto orders = std::vector<std::vector<int>>{
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
      {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

  const auto fDerivs = std::vector{
      std::move(f),   std::move(fx),  std::move(fy),  std::move(fz),
      std::move(fxx), std::move(fyy), std::move(fzz), std::move(fxy),
      std::move(fxz), std::move(fyz),
  };
  for (size_t korder = 0; korder < orders.size(); ++korder) {
    const auto order = orders[korder];
    const auto fourierCoeffOneTerm =
        fft3dApproxOneterm(fDerivs[korder], order[0], order[1], order[2]);
    for (size_t kx = 0; kx < Nfreq_; kx++) {
      for (size_t ky = 0; ky < Nfreq_; ky++) {
        for (size_t kz = 0; kz < Nfreq_; kz++) {
          fouriercoeff[kx][ky][kz] += fourierCoeffOneTerm[kx][ky][kz];
        }
      }
    }
  }

  for (size_t kx = 0; kx < Nfreq_; kx++) {
    for (size_t ky = 0; ky < Nfreq_; ky++) {
      for (size_t kz = 0; kz < Nfreq_; kz++) {
        fouriercoeff[kx][ky][kz] *= coeff_fft;
      }
    }
  }

  return fouriercoeff;
}

}  // namespace coulomb::resampling
