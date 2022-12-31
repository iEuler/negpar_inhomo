#include "Resampler.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "FFT.h"
#include "_global_variables.h"
#include "utils.h"

namespace coulomb {

void merge_NeParticleGroup(NeParticleGroup &S_x,
                           const NeParticleGroup &S_x_new);

NeParticleGroup interp3dRenormalize(NeParticleGroup &S_x) {
  NeParticleGroup S_x_new;
  int Np = S_x.size('p');
  int Nn = S_x.size('n');

  Particle1d3d S_one;

  auto &Sp = S_x.list('p');
  auto &Sn = S_x.list('n');

  const auto &xyz_minmax = S_x.xyz_minmax;

  // interp3d_xyzminmax(S_x, xyz_minmax);
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }
  // for (int k2 = 0; k2 < 6; k2 ++)
  // // cout << xyz_minmax[k2] << ' ';
  // // cout << endl;

  // renormalizaed value
  double v1[3];
  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v1[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v1);
    S_x_new.push_back(S_one, 'p');
  }

  for (int kn = 0; kn < Nn; kn++) {
    auto &v0 = Sn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v1[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v1);
    S_x_new.push_back(S_one, 'n');
  }

  // renormalize Maxwellian
  S_x_new.u1M = (S_x.u1M - xyz_minmax[0]) * 2.0 * pi / Lxyz[0];
  S_x_new.u2M = (S_x.u2M - xyz_minmax[2]) * 2.0 * pi / Lxyz[1];
  S_x_new.u3M = (S_x.u3M - xyz_minmax[4]) * 2.0 * pi / Lxyz[2];
  S_x_new.T1M = S_x.TprtM * (4.0 * pi * pi / Lxyz[0] / Lxyz[0]);
  S_x_new.T2M = S_x.TprtM * (4.0 * pi * pi / Lxyz[1] / Lxyz[1]);
  S_x_new.T3M = S_x.TprtM * (4.0 * pi * pi / Lxyz[2] / Lxyz[2]);
}

void interp3dInvertRenormalize(std::vector<Particle1d3d> &Sp, int Np,
                               const std::vector<double> &xyz_minmax) {
  // rescale to the original coordinates stored in xyz_minmax
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }

  for (int kp = 0; kp < Np; kp++) {
    auto v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v0[k2] = xyz_minmax[2 * k2] + v0[k2] * Lxyz[k2] / (2.0 * pi);
    }
    Sp[kp].set_velocity(v0);
  }
}

/**
 the frequency grids
*/

int frequency(int kth, int Nfreq) {
  int kfreq = kth;
  if (kth >= Nfreq / 2 + 1) kfreq = kth - Nfreq;
  return kfreq;
}

int frequencyInverse(int kfreq, int Nfreq) {
  int kth = kfreq;
  if (kfreq < 0) kth = kfreq + Nfreq;
  return kth;
}

vector<double> interpFrequencies(int Nfreq) {
  vector<double> freq(Nfreq);
  for (int j = 0; j < Nfreq; j++) freq[j] = (double)(frequency(j, Nfreq));
  return freq;
}

vector<int> interpFrequenciesAugmented(int Nfreq, int augFactor) {
  vector<int> loc(Nfreq);
  for (int j = 0; j < Nfreq; j++) {
    int kfreq = frequency(j, Nfreq);
    loc[j] = frequencyInverse(kfreq, augFactor * Nfreq);
  }
  return loc;
}

std::vector<std::complex<double>> interpFrequenciesComplex(int Nfreq) {
  std::vector<std::complex<double>> ifreq(Nfreq);
  for (int j = 0; j < Nfreq / 2 + 1; j++) {
    ifreq[j] = complex<double>(0., (double)j);
  }
  for (int j = Nfreq / 2 + 1; j < Nfreq; j++) {
    ifreq[j] = complex<double>(0., (double)(j - Nfreq));
  }
  return ifreq;
}

VectorBool3D filterFourierCoeff(VectorComplex3D &Fouriercoeff) {
  // double thres = 10.0;
  const auto n1 = Fouriercoeff.size();
  const auto n2 = Fouriercoeff.front().size();
  const auto n3 = Fouriercoeff.front().front().size();
  auto flag_Fouriercoeff =
      std::vector(n1, std::vector(n2, std::vector<bool>(n3, true)));

  // for (size_t kk1 = 0; kk1 < n1; kk1++) {
  //   for (size_t kk2 = 0; kk2 < n2; kk2++) {
  //     for (size_t kk3 = 0; kk3 < n3; kk3++) {
  //       const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
  //       flag_Fouriercoeff[kk1][kk2][kk3] = true;
  //     }
  //   }
  // }
  // for (int k = 0; k < size_FC; k++) {
  //   /*
  //   double abs_FC = abs(Fouriercoeff[k]);
  //   if (abs_FC < thres) {
  //     Fouriercoeff[k] *= 0.;
  //     flag_Fouriercoeff[k] = 0;
  //   }
  //   */
  // }
  return flag_Fouriercoeff;
}

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
Vector3D upperBoundFunc(const Vector3D &fc) {
  const auto n = fc.size();
  auto fUp = std::vector(n, std::vector(n, std::vector<double>(n, 0.0)));

  for (int kx = 0; kx < n; kx++) {
    int xr = kx + 1;
    if (kx == n - 1) xr = 0;

    for (int ky = 0; ky < n; ky++) {
      int yr = ky + 1;
      if (ky == n - 1) yr = 0;

      for (int kz = 0; kz < n; kz++) {
        int zr = kz + 1;
        if (kz == n - 1) zr = 0;

        double max_f_all = std::abs(fc[kx][ky][kz]);
        max_f_all = std::max(max_f_all, std::abs(fc[kx][ky][zr]));
        max_f_all = std::max(max_f_all, std::abs(fc[kx][yr][kz]));
        max_f_all = std::max(max_f_all, std::abs(fc[xr][ky][kz]));
        max_f_all = std::max(max_f_all, std::abs(fc[kx][yr][zr]));
        max_f_all = std::max(max_f_all, std::abs(fc[xr][yr][kz]));
        max_f_all = std::max(max_f_all, std::abs(fc[xr][ky][zr]));
        max_f_all = std::max(max_f_all, std::abs(fc[xr][yr][zr]));

        fUp[kx][ky][zr] = max_f_all;
      }
    }
  }
  return fUp;
}

std::vector<double> getValuesByLoc(const std::vector<Vector3D> &fvecs, int kx,
                                   int ky, int kz) {
  std::vector<double> result;
  result.reserve(fvecs.size());

  for (const auto &fvec : fvecs) {
    result.push_back(fvec[kx][ky][kz]);
  }
  return result;
}

double fvalueApproxFromDeriv(double deltax, double deltay, double deltaz,
                             const std::vector<double> &fDeriv) {
  double f0;
  double f = fDeriv[0];
  double fx = fDeriv[1], fy = fDeriv[2], fz = fDeriv[3];
  double fxx = fDeriv[4], fyy = fDeriv[5], fzz = fDeriv[6];
  double fxy = fDeriv[7], fxz = fDeriv[8], fyz = fDeriv[9];

  f0 = f + fx * deltax + fy * deltay + fx * deltay +
       .5 * fxx * deltax * deltax + .5 * fyy * deltay * deltay +
       .5 * fzz * deltaz * deltaz + fxy * deltax * deltay +
       fxz * deltax * deltaz + fyz * deltay * deltaz;
  return f0;
}

double fvalueFromFFT(const std::vector<double> &Sf,
                     const VectorComplex3D &Fouriercoeff,
                     const std::vector<complex<double>> &ifreq1,
                     const std::vector<complex<double>> &ifreq2,
                     const std::vector<complex<double>> &ifreq3,
                     const VectorBool3D &flag_Fouriercoeff, int Nfreq1,
                     int Nfreq2, int Nfreq3) {
  complex<double> fval_c(0., 0.);

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        if (flag_Fouriercoeff[kk1][kk2][kk3]) {
          fval_c += Fouriercoeff[kk1][kk2][kk3] *
                    exp(ifreq1[kk1] * Sf[0] + ifreq2[kk2] * Sf[1] +
                        ifreq3[kk3] * Sf[2]);
        }
      }
    }
  }

  // return real(fval_c)/(Nfreq1*Nfreq2*Nfreq3);
  return fval_c.real();
}

void acceptSampled(const std::vector<double> &Sf,
                   NeParticleGroup &ptr_S_x_incell, double fval, double &maxf) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);

    int Np_remove = myfloor((1 - keeprate) * ptr_S_x_incell.size('p'));
    int Nn_remove = myfloor((1 - keeprate) * ptr_S_x_incell.size('n'));

    for (int kp = 0; kp < Np_remove; kp++) {
      int k_remove = (int)(myrand() * ptr_S_x_incell.size('p'));
      ptr_S_x_incell.erase(k_remove, 'p');
    }

    for (int kn = 0; kn < Nn_remove; kn++) {
      int k_remove = (int)(myrand() * ptr_S_x_incell.size('n'));
      ptr_S_x_incell.erase(k_remove, 'n');
    }
  }

  // accept this particle with rate abs(fval/maxf)
  if (myrand() < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      if (fval > 0) {
        ptr_S_x_incell.push_back(S_one, 'p');
      } else {
        ptr_S_x_incell.push_back(S_one, 'n');
      }
    }
  }
}

Vector3D Resampler::funcOnAugGrid(const VectorComplex3D &Fouriercoeff) const {
  auto fft3d =
      FFT3D(Nfreq_ * augFactor_, Nfreq_ * augFactor_, Nfreq_ * augFactor_);
  return fft3d.ifft(Fouriercoeff);
}

Vector3D Resampler::derivativesFromFFTHelper(
    const VectorComplex3D &Fouriercoeff, int orderx, int ordery,
    int orderz) const {
  auto FSaug = Fouriercoeff;

  // 1i *freq
  const auto freq1 = interpFrequencies(Nfreq_);
  const auto freq2 = interpFrequencies(Nfreq_);
  const auto freq3 = interpFrequencies(Nfreq_);

  const auto loc1 = interpFrequenciesAugmented(Nfreq_, augFactor_);
  const auto loc2 = interpFrequenciesAugmented(Nfreq_, augFactor_);
  const auto loc3 = interpFrequenciesAugmented(Nfreq_, augFactor_);

  for (int kk1 = 0; kk1 < Nfreq_; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq_; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq_; kk3++) {
        int kk1aug = loc1[kk1];
        int kk2aug = loc2[kk2];
        int kk3aug = loc3[kk3];

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

std::vector<Vector3D> Resampler::derivativesFromFFT(
    const VectorComplex3D &Fouriercoeff) const {
  const auto orders = std::vector<std::vector<int>>{
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
      {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

  std::vector<Vector3D> fDerivatives;
  for (const auto &order : orders) {
    fDerivatives.push_back(
        derivativesFromFFTHelper(Fouriercoeff, order[0], order[1], order[2]));
  }
  return fDerivatives;
}

std::shared_ptr<NeParticleGroup> Resampler::resample() {
  NeParticleGroup S_x_new;
  auto &S_x = *negParGroup_;
  auto Nfreq = Nfreq_;
  double Neff = Neff_;

  /* Normalize particle velocity to [0 2*pi] */
  S_x.set_xyzrange();

  auto S_x_renormalized = interp3dRenormalize(S_x);

  /* Prepare the grids in physical space and frequence space */
  // double dx = 2.0*pi/Nfreq;

  const auto ifreq = interpFrequenciesComplex(Nfreq);  // 1i *freq
  std::vector<double> interp_x(Nfreq);
  for (int kx = 0; kx < Nfreq; kx++) interp_x[kx] = kx * 2 * pi / Nfreq;

  /* Compute the Fourier coefficient */
  VectorComplex3D Fouriercoeff;

  if (useApproximation_)
    Fouriercoeff = fft3dApprox(S_x_renormalized, Nfreq, Nfreq, Nfreq);
  else
    Fouriercoeff = fft3d(S_x_renormalized, Nfreq, Nfreq, Nfreq);

  // Apply the filter on Fourier coefficients
  auto flag_Fouriercoeff = filterFourierCoeff(Fouriercoeff);

  // cout << " F coeff computed " << endl;

  /* Compute a coarse interpolation in physical space */
  //  const auto fcoarse = interp3d_fcoarse(Fouriercoeff, Nfreq, Nfreq, Nfreq);

  auto augFactor = augFactor_;
  auto sizeF = augFactor * augFactor * augFactor * Nfreq * Nfreq * Nfreq;

  const auto fDerivatives = derivativesFromFFT(Fouriercoeff);

  /* evaluate the upperbound of f */
  const auto f_up = upperBoundFunc(fDerivatives[0]);

  /* refined x grid */
  double dxaug = 2.0 * pi / Nfreq / augFactor;
  std::vector<double> interp_xaug(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq / augFactor;

  const auto &f = fDerivatives[0];
  /* create a NeParticleGroup to host the P and N particles in current cell */

  /* Start sampling */

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        double fcc = f_up[kx][ky][kz];

        if (fcc < std::abs(f[kx][ky][kz])) throw std::exception("small bound!");

        double maxf = 1.5 * fcc;
        int N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff);

        int k_virtual = 0;
        NeParticleGroup S_x_incell;

        while (k_virtual < N_incell) {
          // create a particle in the cell
          // double Sf[3] = {xc+myrand()*dx, yc+myrand()*dx, zc+myrand()*dx};
          double deltax = myrand() * dxaug - 0.5 * dxaug;
          double deltay = myrand() * dxaug - 0.5 * dxaug;
          double deltaz = myrand() * dxaug - 0.5 * dxaug;
          std::vector<double> Sf{xc + deltax, yc + deltay, zc + deltaz};

          // compute f at this point
          double fval = 0;
          if (useApproximation_) {
            const auto fDeriv = getValuesByLoc(fDerivatives, kx, ky, kz);
            fval = fvalueApproxFromDeriv(deltax, deltay, deltaz, fDeriv);
          } else
            fval = fvalueFromFFT(Sf, Fouriercoeff, ifreq, ifreq, ifreq,
                                 flag_Fouriercoeff, Nfreq, Nfreq, Nfreq);

          // reset current cell if fval>maxf, otherwise continue sampling in
          // current cell
          acceptSampled(Sf, S_x_incell, fval, maxf);

          // reset N_incell if maxf is changed
          N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff);
          k_virtual++;
        }

        merge_NeParticleGroup(S_x_new, S_x_incell);
      }
    }
  }

  // cout << "Resampled." << endl;

  // rescale to the original coordinates
  auto &Sp_sampled = S_x_new.list('p');
  auto &Sn_sampled = S_x_new.list('n');
  const auto &xyz_minmax = S_x.xyz_minmax;
  interp3dInvertRenormalize(Sp_sampled, S_x_new.size('p'), xyz_minmax);
  interp3dInvertRenormalize(Sn_sampled, S_x_new.size('n'), xyz_minmax);

  // cout << "Rescaled." << endl;

  return std::make_shared<NeParticleGroup>(S_x_new);
}

VectorComplex3D Resampler::fft3d(NeParticleGroup &S_x, int Nfreq1, int Nfreq2,
                                 int Nfreq3) const {
  std::vector<std::complex<double>> Fouriercoeff(Nfreq1 * Nfreq2 * Nfreq3);

  int Np = S_x.size('p');
  int Nn = S_x.size('n');

  auto &Sp = S_x.list('p');
  auto &Sn = S_x.list('n');

  double Neff = 1.0;

  // double Neff_temp = 1./Np;

  const auto ifreq1 = interpFrequencies(Nfreq1);
  const auto ifreq2 = interpFrequencies(Nfreq2);
  const auto ifreq3 = interpFrequencies(Nfreq3);

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi * Nfreq1 * Nfreq2 * Nfreq3;
  double maxFS = 0.0;

  // the (i,j,k)-th element of the array with size (Nx,Ny,Nz), you would use the
  // expression an_array[k + Nz * (j + Ny * i)].

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);
        Fouriercoeff[kk] = complex<double>(0., 0.);
        for (int kp = 0; kp < Np; kp++) {
          auto &vp = Sp[kp].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          Fouriercoeff[kk] += exp(-expterm);
        }
        for (int kn = 0; kn < Nn; kn++) {
          auto &vp = Sn[kn].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          Fouriercoeff[kk] -= exp(-expterm);
        }
        // Fouriercoeff[kk] *= Neff * coeff_fft;
        Fouriercoeff[kk] *= Neff * coeff_fft / (Nfreq1 * Nfreq2 * Nfreq3);
        maxFS = max(maxFS, abs(Fouriercoeff[kk]));
      }
    }
  }

  return reshape1dTo3d(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3);
}

VectorComplex3D Resampler::fft3dApprox(NeParticleGroup &S_x, int Nfreq1,
                                       int Nfreq2, int Nfreq3) const {
  std::vector<std::complex<double>> Fouriercoeff(Nfreq1 * Nfreq2 * Nfreq3,
                                                 {0., 0.});
  int augFactor = 2;

  int Np = S_x.size('p');
  int Nn = S_x.size('n');

  auto &Sp = S_x.list('p');
  auto &Sn = S_x.list('n');

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
    if (xloc >= augFactor * Nfreq1) xloc--;
    if (yloc >= augFactor * Nfreq2) yloc--;
    if (zloc >= augFactor * Nfreq3) zloc--;
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
    if (xloc >= augFactor * Nfreq1) xloc--;
    if (yloc >= augFactor * Nfreq2) yloc--;
    if (zloc >= augFactor * Nfreq3) zloc--;
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

  return reshape1dTo3d(Fouriercoeff, Nfreq1, Nfreq2, Nfreq3);
  // cout << "Approx finished." << endl;
}
}  // namespace coulomb