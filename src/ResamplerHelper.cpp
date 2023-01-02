#include "ResamplerHelper.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "_global_variables.h"
#include "utils.h"

namespace coulomb {

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
  S_x_new.rhoM = S_x.rhoM;
  S_x_new.u1M = (S_x.u1M - xyz_minmax[0]) * 2.0 * pi / Lxyz[0];
  S_x_new.u2M = (S_x.u2M - xyz_minmax[2]) * 2.0 * pi / Lxyz[1];
  S_x_new.u3M = (S_x.u3M - xyz_minmax[4]) * 2.0 * pi / Lxyz[2];
  S_x_new.T1M = S_x.TprtM * (4.0 * pi * pi / Lxyz[0] / Lxyz[0]);
  S_x_new.T2M = S_x.TprtM * (4.0 * pi * pi / Lxyz[1] / Lxyz[1]);
  S_x_new.T3M = S_x.TprtM * (4.0 * pi * pi / Lxyz[2] / Lxyz[2]);
}

void interp3dInvertRenormalize(std::vector<Particle1d3d> &Sp,
                               const std::vector<double> &xyz_minmax) {
  // rescale to the original coordinates stored in xyz_minmax
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }

  for (int kp = 0; kp < Sp.size(); kp++) {
    auto v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v0[k2] = xyz_minmax[2 * k2] + v0[k2] * Lxyz[k2] / (2.0 * pi);
    }
    Sp[kp].set_velocity(v0);
  }
}

int frequency(int kth, size_t Nfreq) {
  int kfreq = kth;
  if (kth >= Nfreq / 2 + 1) kfreq = kth - Nfreq;
  return kfreq;
}

int frequencyInverse(int kfreq, size_t Nfreq) {
  int kth = kfreq;
  if (kfreq < 0) kth = kfreq + Nfreq;
  return kth;
}

std::vector<double> interpFrequencies(size_t Nfreq) {
  std::vector<double> freq(Nfreq);
  for (size_t j = 0; j < Nfreq; j++) freq[j] = (double)(frequency(j, Nfreq));
  return freq;
}

std::vector<size_t> augmentedLocation(size_t Nfreq, size_t augFactor) {
  std::vector<size_t> loc(Nfreq);
  for (size_t j = 0; j < Nfreq; j++) {
    int kfreq = frequency(j, Nfreq);
    loc[j] = frequencyInverse(kfreq, augFactor * Nfreq);
  }
  return loc;
}

std::vector<std::complex<double>> interpFrequenciesComplex(size_t Nfreq) {
  std::vector<std::complex<double>> ifreq(Nfreq);
  for (size_t j = 0; j < Nfreq / 2 + 1; j++) {
    ifreq[j] = complex<double>(0., (double)j);
  }
  for (size_t j = Nfreq / 2 + 1; j < Nfreq; j++) {
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
                     const std::vector<std::complex<double>> &ifreq1,
                     const std::vector<std::complex<double>> &ifreq2,
                     const std::vector<std::complex<double>> &ifreq3,
                     const VectorBool3D &flag_Fouriercoeff) {
  std::complex<double> fval_c(0., 0.);

  for (int kk1 = 0; kk1 < ifreq1.size(); kk1++) {
    for (int kk2 = 0; kk2 < ifreq2.size(); kk2++) {
      for (int kk3 = 0; kk3 < ifreq3.size(); kk3++) {
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

void acceptSampled(const std::vector<double> &Sf, NeParticleGroup &S_x_incell,
                   double fval, double &maxf, bool sampleFromFullDistribution) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    const auto parTypes = sampleFromFullDistribution ? "f"s : "pn"s;
    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);

    for (const auto parType : parTypes) {
      int Np_remove = myfloor((1 - keeprate) * S_x_incell.size(parType));
      for (int kp = 0; kp < Np_remove; kp++) {
        int k_remove = (int)(myrand() * S_x_incell.size(parType));
        S_x_incell.erase(k_remove, parType);
      }
    }
  }

  // accept this particle with rate abs(fval/maxf)
  if (myrand() < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      const auto parType =
          sampleFromFullDistribution ? 'f' : (fval > 0 ? 'p' : 'n');
      S_x_incell.push_back(S_one, parType);
    }
  }
}

void addMaxwellian_terms(double rhoM, vector<double> uM, vector<double> TM,
                         double Neff, Vector3D &f, int Nfreq, int augFactor,
                         int orderx, int ordery, int orderz);

void addMaxwellian(const NeParticleGroup &S_x, double Neff,
                   std::vector<Vector3D> &fDerivatives, int Nfreq,
                   int augFactor, double dxSpace) {
  double rhoM = S_x.rhoM * dxSpace;
  std::vector<double> uM{S_x.u1M, S_x.u2M, S_x.u3M};
  std::vector<double> TM{S_x.T1M, S_x.T2M, S_x.T3M};

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

}  // namespace coulomb