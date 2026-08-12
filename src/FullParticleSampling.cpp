#include "FullParticleSampling.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "FullParticleFourier.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"

namespace coulomb {
using std::abs;
using std::sqrt;
using std::vector;

namespace {

void resampleF_acceptsampled(const std::vector<double> &Sf,
                             NeParticleGroup &ptr_S_x_incell, double fval,
                             double &maxf, RandomContext &random) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new
    double keeprate = maxf / (1.5 * abs(fval));
    maxf = 1.5 * abs(fval);

    int Np_remove = RandomSampling(random).stochastic_floor(
        (1 - keeprate) * ptr_S_x_incell.size(ParticleKind::Full));

    for (int kp = 0; kp < Np_remove; kp++) {
      int k_remove = (int)(RandomSampling(random).uniform() *
                           ptr_S_x_incell.size(ParticleKind::Full));
      ptr_S_x_incell.erase(k_remove, ParticleKind::Full);
    }
  }

  // accept this particle with rate abs(fval/maxf)
  if (RandomSampling(random).uniform() < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      ptr_S_x_incell.push_back(S_one, ParticleKind::Full);
    }
  }
}

} // namespace

NeParticleGroup FullParticleSampling::resample(NeParticleGroup &S_x, int Nfreq,
                                               double Neff, double Neff_F,
                                               double dx_space,
                                               RandomContext &random) {
  NeParticleGroup S_x_new;
  /* Normalize particle velocity to [0 2*pi] */
  S_x.set_xyzrange();

  auto S_x_renormalized =
      resampling::ResamplingVelocity{}.normalize_signed(S_x);

  const auto ifreq =
      resampling::ResamplingNumerics{}.imaginary_frequencies(Nfreq);
  vector<double> interp_x(Nfreq);
  for (int kx = 0; kx < Nfreq; kx++)
    interp_x[kx] = kx * 2 * pi / Nfreq;

  vector<int> flag_Fouriercoeff(Nfreq * Nfreq * Nfreq);

  auto Fouriercoeff = FullParticleFourier{}.approximate_transform(
      S_x_renormalized, Nfreq, Nfreq, Nfreq);
  FullParticleFourier{}.filter(Fouriercoeff, flag_Fouriercoeff,
                               Nfreq * Nfreq * Nfreq);

  const auto fcoarse = FullParticleFourier{}.interpolate_coarse(
      Fouriercoeff, Nfreq, Nfreq, Nfreq);

  int augFactor = 2;
  auto fDerivatives = FullParticleFourier{}.interpolate_derivatives(
      Fouriercoeff, Nfreq, Nfreq, Nfreq, augFactor);

  vector<double> uM(3);
  vector<double> TM(3);
  double rhoM = S_x.rhoM * dx_space;
  uM[0] = S_x_renormalized.u1M;
  uM[1] = S_x_renormalized.u2M;
  uM[2] = S_x_renormalized.u3M;
  TM[0] = S_x_renormalized.T1M;
  TM[1] = S_x_renormalized.T2M;
  TM[2] = S_x_renormalized.T3M;

  FullParticleFourier{}.add_maxwellian(rhoM, uM, TM, Neff, fDerivatives, Nfreq,
                                       augFactor);
  const auto f = fDerivatives[0];

  const auto f_up = FullParticleFourier{}.upper_bound(augFactor * Nfreq, f);

  double dxaug = 2.0 * pi / Nfreq / augFactor;
  vector<double> interp_xaug(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq / augFactor;

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];
        double fcc = f_up[kk];

        double maxf = 1.5 * abs(fcc);
        int N_incell = RandomSampling(random).stochastic_floor(
            maxf * dxaug * dxaug * dxaug / Neff_F);

        int k_virtual = 0;
        NeParticleGroup S_x_incell;

        while (k_virtual < N_incell) {
          double deltax =
              RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
          double deltay =
              RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
          double deltaz =
              RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
          std::vector<double> Sf{xc + deltax, yc + deltay, zc + deltaz};

          const auto fDeriv = FullParticleFourier{}.values_at(fDerivatives, kk);
          double fval =
              resampling::ResamplingNumerics{}.evaluate_quadratic_taylor(
                  deltax, deltay, deltaz, fDeriv);

          resampleF_acceptsampled(Sf, S_x_incell, fval, maxf, random);

          N_incell = RandomSampling(random).stochastic_floor(
              maxf / (Neff_F / (dxaug * dxaug * dxaug)));
          k_virtual++;
        }

        ParticleGroupOperations{}.merge_full(S_x_new, S_x_incell);
      }
    }
  }

  auto &Sp_sampled = S_x_new.list(ParticleKind::Full);
  const auto &xyz_minmax = S_x.xyz_minmax;
  resampling::ResamplingVelocity{}.restore(Sp_sampled, xyz_minmax);

  std::cout << "# resampled F = " << S_x_new.size(ParticleKind::Full)
            << std::endl;

  return S_x_new;
}

} // namespace coulomb
