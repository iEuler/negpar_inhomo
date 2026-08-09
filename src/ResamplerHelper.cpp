#include "ResamplerHelper.h"

#include <cmath>
#include <complex>
#include <vector>

#include "Constants.h"
#include "utils.h"

namespace coulomb::resampling {

using std::abs;
using std::complex;
using std::exp;
using std::sqrt;
using std::vector;

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

        fUp[kx][ky][kz] = max_f_all;
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
                   double fval, double &maxf, RandomContext& random) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);

    const auto remove_particles = [&](ParticleKind kind) {
      const int count =
          myfloor((1 - keeprate) * S_x_incell.size(kind), random);
      for (int particle = 0; particle < count; ++particle) {
        const int index =
            static_cast<int>(myrand(random) * S_x_incell.size(kind));
        S_x_incell.erase(index, kind);
      }
    };
    remove_particles(ParticleKind::Positive);
    remove_particles(ParticleKind::Negative);
  }

  // accept this particle with rate abs(fval/maxf)
  if (myrand(random) < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      const auto kind = fval > 0 ? ParticleKind::Positive
                                 : ParticleKind::Negative;
      S_x_incell.push_back(S_one, kind);
    }
  }
}

}  // namespace coulomb::resampling
