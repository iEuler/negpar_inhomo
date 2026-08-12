#include "FullParticleFourier.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>

#include "Constants.h"
#include "FFT.h"
#include "ParticleGroup.h"
#include "ResamplingNumerics.h"

namespace coulomb {
using std::abs;
using std::complex;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::min;
using std::pow;
using std::sqrt;
using std::vector;

namespace {

/* ======================================================== *\
        Use Fourier transform for 3D interpolation
\* ======================================================== */

// xyz_minmax = [xmin, xmax, ymin, ymax, zmin, zmax];

void interp3d_fft_approx_terms(std::vector<std::complex<double>> &Fouriercoeff,
                               vector<double> &f, int Nfreq1, int Nfreq2,
                               int Nfreq3, int augFactor, int orderx,
                               int ordery, int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;

  // 1i *freq
  const auto freq1 = resampling::ResamplingNumerics::frequencies(Nfreq1);
  const auto freq2 = resampling::ResamplingNumerics::frequencies(Nfreq2);
  const auto freq3 = resampling::ResamplingNumerics::frequencies(Nfreq3);

  const auto loc1 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq1, augFactor);
  const auto loc2 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq2, augFactor);
  const auto loc3 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq3, augFactor);

  FFTWBuffer fin(static_cast<fftw_complex *>(
      fftw_malloc(static_cast<std::size_t>(sizeF) * sizeof(fftw_complex))));
  FFTWBuffer FSaug(static_cast<fftw_complex *>(
      fftw_malloc(static_cast<std::size_t>(sizeF) * sizeof(fftw_complex))));
  if (!fin || !FSaug)
    throw std::runtime_error("Unable to allocate interpolated FFTW buffers");

  FFTWPlan plan3d_ft(fftw_plan_dft_3d(
      augFactor * Nfreq1, augFactor * Nfreq2, augFactor * Nfreq3, fin.get(),
      FSaug.get(), FFTW_FORWARD, FFTW_ESTIMATE));
  if (!plan3d_ft)
    throw std::runtime_error("Unable to create interpolated FFTW plan");

  for (int kfc = 0; kfc < sizeF; kfc++) {
    fin[kfc][0] = f[kfc];
    fin[kfc][1] = 0.;
  }

  fftw_execute(plan3d_ft.get());

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);

        const std::size_t kk1aug = loc1[kk1];
        const std::size_t kk2aug = loc2[kk2];
        const std::size_t kk3aug = loc3[kk3];
        const std::size_t kkaug =
            kk3aug + static_cast<std::size_t>(augFactor) * Nfreq3 *
                         (kk2aug + static_cast<std::size_t>(augFactor) *
                                       Nfreq2 * kk1aug);

        double freq = 1.0;
        for (int kx = 0; kx < orderx; kx++)
          freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++)
          freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++)
          freq *= freq3[kk3];

        if ((orderx + ordery + orderz) == 2) {
          if ((orderx == 2) || (ordery == 2) || (orderz == 2))
            freq *= -.5;
          else
            freq *= -1;
        }

        if ((orderx + ordery + orderz) == 1) {
          Fouriercoeff[kk] +=
              freq * complex<double>(FSaug[kkaug][1], -FSaug[kkaug][0]);
        } else {
          Fouriercoeff[kk] +=
              freq * complex<double>(FSaug[kkaug][0], FSaug[kkaug][1]);
        }
      }
    }
  }
}

} // namespace

std::vector<std::complex<double>>
FullParticleFourier::approximate_transform(NeParticleGroup &S_x, int Nfreq1,
                                           int Nfreq2, int Nfreq3) {
  std::vector<std::complex<double>> Fouriercoeff(Nfreq1 * Nfreq2 * Nfreq3,
                                                 {0., 0.});
  int augFactor = 2;

  int Np = S_x.size(ParticleKind::Positive);
  int Nn = S_x.size(ParticleKind::Negative);

  auto &Sp = S_x.list(ParticleKind::Positive);
  auto &Sn = S_x.list(ParticleKind::Negative);

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi;

  // for the (i,j,k)-th element of the array with size (Nx,Ny,Nz), use the
  // expression an_array[k + Nz * (j + Ny * i)].

  // create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  double dx = 2.0 * pi / augFactor / Nfreq1;
  double dy = 2.0 * pi / augFactor / Nfreq2;
  double dz = 2.0 * pi / augFactor / Nfreq3;

  vector<double> f(sizeF);
  vector<double> fx(sizeF);
  vector<double> fy(sizeF);
  vector<double> fz(sizeF);
  vector<double> fxx(sizeF);
  vector<double> fyy(sizeF);
  vector<double> fzz(sizeF);
  vector<double> fxy(sizeF);
  vector<double> fxz(sizeF);
  vector<double> fyz(sizeF);

  // cout << "Approx 1" << endl;

  for (int kk = 0; kk < sizeF; kk++) {
    f[kk] = 0.;
    fx[kk] = 0.;
    fy[kk] = 0.;
    fz[kk] = 0.;
    fxx[kk] = 0.;
    fyy[kk] = 0.;
    fzz[kk] = 0.;
    fxy[kk] = 0.;
    fxz[kk] = 0.;
    fyz[kk] = 0.;
  }

  for (int kp = 0; kp < Np; kp++) {
    double x0 = Sp[kp].velocity(0);
    double y0 = Sp[kp].velocity(1);
    double z0 = Sp[kp].velocity(2);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= augFactor * Nfreq1)
      xloc--;
    if (yloc >= augFactor * Nfreq2)
      yloc--;
    if (zloc >= augFactor * Nfreq3)
      zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

    if ((loc >= sizeF) || (loc < 0))
      cout << x0 << ' ' << y0 << ' ' << z0 << ' ' << dx << ' ' << loc << endl;

    f[loc]++;
    fx[loc] += xdelta;
    fy[loc] += ydelta;
    fz[loc] += zdelta;
    fxx[loc] += xdelta * xdelta;
    fyy[loc] += ydelta * ydelta;
    fzz[loc] += zdelta * zdelta;
    fxy[loc] += xdelta * ydelta;
    fyz[loc] += ydelta * zdelta;
    fxz[loc] += zdelta * xdelta;
  }

  // cout << "Approx 2" << endl;

  for (int kp = 0; kp < Nn; kp++) {
    double x0 = Sn[kp].velocity(0);
    double y0 = Sn[kp].velocity(1);
    double z0 = Sn[kp].velocity(2);
    // int xloc = floor(x0/dx);
    // int yloc = floor(y0/dy);
    // int zloc = floor(z0/dz);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= augFactor * Nfreq1)
      xloc--;
    if (yloc >= augFactor * Nfreq2)
      yloc--;
    if (zloc >= augFactor * Nfreq3)
      zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

    if ((loc >= sizeF) || (loc < 0)) {
      cout << "error: in approximation. Particle moved out of range. kx = "
           << xloc << ' ' << yloc << ' ' << zloc << ' ' << loc << endl;
      exit(0);
    }

    f[loc]--;
    fx[loc] -= xdelta;
    fy[loc] -= ydelta;
    fz[loc] -= zdelta;
    fxx[loc] -= xdelta * xdelta;
    fyy[loc] -= ydelta * ydelta;
    fzz[loc] -= zdelta * zdelta;
    fxy[loc] -= xdelta * ydelta;
    fyz[loc] -= ydelta * zdelta;
    fxz[loc] -= zdelta * xdelta;
  }
  // cout << "Approx 3" << endl;

  interp3d_fft_approx_terms(Fouriercoeff, f, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 0, 0);

  interp3d_fft_approx_terms(Fouriercoeff, fx, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            1, 0, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fy, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 1, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fz, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 0, 1);
  interp3d_fft_approx_terms(Fouriercoeff, fxx, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 2, 0, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fyy, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 2, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fzz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 0, 2);
  interp3d_fft_approx_terms(Fouriercoeff, fxy, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 1, 1, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fxz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 1, 0, 1);
  interp3d_fft_approx_terms(Fouriercoeff, fyz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 1, 1);

  for (int kk = 0; kk < Nfreq1 * Nfreq2 * Nfreq3; kk++)
    Fouriercoeff[kk] *= coeff_fft;

  return Fouriercoeff;
  // cout << "Approx finished." << endl;
}

void FullParticleFourier::filter(
    std::vector<std::complex<double>> & /*Fouriercoeff*/,
    vector<int> &flag_Fouriercoeff, int size_FC) {
  // double thres = 10.0;
  for (int k = 0; k < size_FC; k++) {
    flag_Fouriercoeff[k] = 1;
    /*
    double abs_FC = abs(Fouriercoeff[k]);
    if (abs_FC < thres) {
      Fouriercoeff[k] *= 0.;
      flag_Fouriercoeff[k] = 0;
    }
    */
  }
}

/******************************************************************/
/* ---------- Use Fourier transform for 3D interpolation -------- */
/******************************************************************/

/*
  Find the coarse approximation with the given Fourier coefficients
  Need to include 'fftw3.f'
*/
vector<double> FullParticleFourier::interpolate_coarse(
    const std::vector<std::complex<double>> &Fouriercoeff, int Nfreq1,
    int Nfreq2, int Nfreq3) {
  // double Lcubic = (double) (Nfreq1*Nfreq2*Nfreq3);

  vector<double> fcoarse(Nfreq1 * Nfreq2 * Nfreq3);

  const auto element_count = static_cast<std::size_t>(Nfreq1) * Nfreq2 * Nfreq3;
  FFTWBuffer FC(static_cast<fftw_complex *>(
      fftw_malloc(element_count * sizeof(fftw_complex))));
  FFTWBuffer fcoarse_c(static_cast<fftw_complex *>(
      fftw_malloc(element_count * sizeof(fftw_complex))));
  if (!FC || !fcoarse_c)
    throw std::runtime_error("Unable to allocate coarse FFTW buffers");

  for (int kfc = 0; kfc < Nfreq1 * Nfreq2 * Nfreq3; kfc++) {
    FC[kfc][0] = Fouriercoeff[kfc].real();
    FC[kfc][1] = Fouriercoeff[kfc].imag();
  }

  // use fftw to obtain an estimation of f_p - f_n
  // // cout << " resample 1.3" << endl;
  FFTWPlan plan3d_ift(fftw_plan_dft_3d(Nfreq1, Nfreq2, Nfreq3, FC.get(),
                                       fcoarse_c.get(), FFTW_BACKWARD,
                                       FFTW_ESTIMATE));
  if (!plan3d_ift)
    throw std::runtime_error("Unable to create coarse FFTW plan");

  fftw_execute(plan3d_ift.get());

  for (int kfc = 0; kfc < Nfreq1 * Nfreq2 * Nfreq3; kfc++) {
    // fcoarse[kfc] = fcoarse_c[kfc][0] / Lcubic;
    fcoarse[kfc] = fcoarse_c[kfc][0];
  }
  return fcoarse;
}

/** The derivatives of f on grids
 */
// CONTINUE HERE

namespace {

vector<double>
interp3d_fxyz_terms(const std::vector<std::complex<double>> &Fouriercoeff,
                    int Nfreq1, int Nfreq2, int Nfreq3, int augFactor,
                    int orderx, int ordery, int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  std::vector<complex<double>> FSaug(sizeF, {0., 0.});

  // 1i *freq
  const auto freq1 = resampling::ResamplingNumerics::frequencies(Nfreq1);
  const auto freq2 = resampling::ResamplingNumerics::frequencies(Nfreq2);
  const auto freq3 = resampling::ResamplingNumerics::frequencies(Nfreq3);

  const auto loc1 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq1, augFactor);
  const auto loc2 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq2, augFactor);
  const auto loc3 =
      resampling::ResamplingNumerics::augmented_locations(Nfreq3, augFactor);

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);

        const std::size_t kk1aug = loc1[kk1];
        const std::size_t kk2aug = loc2[kk2];
        const std::size_t kk3aug = loc3[kk3];
        const std::size_t kkaug =
            kk3aug + static_cast<std::size_t>(augFactor) * Nfreq3 *
                         (kk2aug + static_cast<std::size_t>(augFactor) *
                                       Nfreq2 * kk1aug);

        double freq = 1.0;
        for (int kx = 0; kx < orderx; kx++)
          freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++)
          freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++)
          freq *= freq3[kk3];

        FSaug[kkaug] = freq * Fouriercoeff[kk];
      }
    }
  }

  return FullParticleFourier::interpolate_coarse(
      FSaug, augFactor * Nfreq1, augFactor * Nfreq2, augFactor * Nfreq3);
}

} // namespace

std::vector<std::vector<double>> FullParticleFourier::interpolate_derivatives(
    const std::vector<std::complex<double>> &Fouriercoeff, int Nfreq1,
    int Nfreq2, int Nfreq3, int augFactor) {
  const auto f = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                     augFactor, 0, 0, 0);
  const auto fx = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                      augFactor, 1, 0, 0);
  const auto fy = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                      augFactor, 0, 1, 0);
  const auto fz = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                      augFactor, 0, 0, 1);
  const auto fxx = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 2, 0, 0);
  const auto fyy = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 0, 2, 0);
  const auto fzz = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 0, 0, 2);
  const auto fxy = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 1, 1, 0);
  const auto fxz = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 1, 0, 1);
  const auto fyz = interp3d_fxyz_terms(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3,
                                       augFactor, 0, 1, 1);
  return {f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz};
}

// void interp3d_fft_eachlevel(NeParticleGroup * S_x, MultlLevelGroup * MLsol,
// int Nlevel); void interp3d_fft_ml(complex<double> *Fouriercoeff, int *
// flag_Fouriercoeff, MultlLevelGroup * MLsol, int Nlevel);
/*
  Sampled particles are stored in particles_Sp_sampled and particles_Sn_sampled
  with size particles%Np_sampled and particles%Nn_sampled
*/

std::vector<double>
FullParticleFourier::values_at(const std::vector<std::vector<double>> &fvecs,
                               int k) {
  std::vector<double> result;
  result.reserve(fvecs.size());

  for (const auto &fvec : fvecs) {
    result.push_back(fvec[k]);
  }
  return result;
}

vector<double> FullParticleFourier::upper_bound(int N,
                                                const vector<double> &fc) {
  vector<double> f_up(N * N * N);
  // Find f_up>abs(fc)
  double f_all[8];
  int kk[8];

  for (int kx = 0; kx < N; kx++) {
    int xr = kx + 1;
    if (kx == N - 1)
      xr = 0;

    for (int ky = 0; ky < N; ky++) {
      int yr = ky + 1;
      if (ky == N - 1)
        yr = 0;

      for (int kz = 0; kz < N; kz++) {
        int zr = kz + 1;
        if (kz == N - 1)
          zr = 0;

        kk[0] = kz + N * (ky + N * kx);
        kk[1] = zr + N * (ky + N * kx);
        kk[2] = kz + N * (yr + N * kx);
        kk[3] = zr + N * (yr + N * kx);
        kk[4] = kz + N * (ky + N * xr);
        kk[5] = zr + N * (ky + N * xr);
        kk[6] = kz + N * (yr + N * xr);
        kk[7] = zr + N * (yr + N * xr);

        double max_f_all = 0.;
        for (int k = 0; k < 8; k++) {
          f_all[k] = abs(fc[kk[k]]);
          max_f_all = max(max_f_all, f_all[k]);
        }

        f_up[kk[0]] = max_f_all;
      }
    }
  }
  return f_up;
}

namespace {

void addMaxwellian_terms(double rhoM, vector<double> uM, vector<double> TM,
                         double Neff, vector<double> &f, int Nfreq,
                         int augFactor, int orderx, int ordery, int orderz) {
  double Mcc_coe = rhoM / sqrt(8.0 * pi * pi * pi * TM[0] * TM[1] * TM[2]);

  vector<double> interp_xaug(Nfreq * augFactor);
  vector<double> exp_x(Nfreq * augFactor);
  vector<double> exp_y(Nfreq * augFactor);
  vector<double> exp_z(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++) {
    double xk = kx * 2 * pi / Nfreq / augFactor;
    interp_xaug[kx] = xk;
    exp_x[kx] = exp(-(xk - uM[0]) * (xk - uM[0]) / 2 / TM[0]);
    exp_y[kx] = exp(-(xk - uM[1]) * (xk - uM[1]) / 2 / TM[1]);
    exp_z[kx] = exp(-(xk - uM[2]) * (xk - uM[2]) / 2 / TM[2]);
  }

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double Mcc = Mcc_coe * exp_x[kx] * exp_y[ky] * exp_z[kz];
        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        if (orderx == 1)
          Mcc *= -(xc - uM[0]) / TM[0];
        if (orderx == 2)
          Mcc *= ((xc - uM[0]) * (xc - uM[0]) - TM[0]) / TM[0] / TM[0];
        if (ordery == 1)
          Mcc *= -(yc - uM[1]) / TM[1];
        if (ordery == 2)
          Mcc *= ((yc - uM[1]) * (yc - uM[1]) - TM[1]) / TM[1] / TM[1];
        if (orderz == 1)
          Mcc *= -(zc - uM[2]) / TM[2];
        if (orderz == 2)
          Mcc *= ((zc - uM[2]) * (zc - uM[2]) - TM[2]) / TM[2] / TM[2];

        f[kk] = Neff * f[kk] + Mcc;
      }
    }
  }
}

} // namespace

void FullParticleFourier::add_maxwellian(
    double rhoM, vector<double> uM, vector<double> TM, double Neff,
    std::vector<vector<double>> &fDerivatives, int Nfreq, int augFactor) {
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[0], Nfreq, augFactor, 0,
                      0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[1], Nfreq, augFactor, 1,
                      0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[2], Nfreq, augFactor, 0,
                      1, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[3], Nfreq, augFactor, 0,
                      0, 1);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[4], Nfreq, augFactor, 2,
                      0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[5], Nfreq, augFactor, 0,
                      2, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[6], Nfreq, augFactor, 0,
                      0, 2);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[7], Nfreq, augFactor, 1,
                      1, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[8], Nfreq, augFactor, 1,
                      0, 1);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fDerivatives[9], Nfreq, augFactor, 0,
                      1, 1);
}

} // namespace coulomb
