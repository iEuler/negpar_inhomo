#include "InhomoResampler.h"

#include <iostream>
namespace coulomb {

void merge_NeParticleGroup(NeParticleGroup& S_x,
                           const NeParticleGroup& S_x_new);

bool resampleHomo(NeParticleGroup& S_x, double Neff, size_t Nfreq,
                  bool useApproximation) {
  Resampler resampler(S_x, Neff, Nfreq, useApproximation);

  // NUM_RESAMPLE++;

  // resample particles
  auto S_x_new = resampler.resample();
  S_x.isResampled = true;

  int Np_old = S_x.size('p'), Nn_old = S_x.size('n');
  int Np_new = S_x_new.size('p'), Nn_new = S_x_new.size('n');

  if ((Np_new < Np_old) && (Nn_new < Nn_old)) {
    // Replace the original particles by new sampled particles
    S_x_new.setPositionRangeAndRandomizeValues(S_x.get_xmin(), S_x.get_xmax());
    S_x.clear('p');
    S_x.clear('n');
    merge_NeParticleGroup(S_x, S_x_new);
    return true;
  } else {
    std::cout << "New sampled particles rejected." << std::endl;
    return false;
  }
}
void InhomoResampler::resample(std::vector<NeParticleGroup>& S_x) {
  bool needGlobalResample = false;

  bool flag_resample_success = true;

  for (const auto& S_kx : S_x) {
    if ((S_kx.size('p') + S_kx.size('n')) >= S_kx.size('f'))
      needGlobalResample = true;
  }
  if (needGlobalResample) {
    double resample_spatial_ratio = 0;
    for (const auto& S_kx : S_x) {
      resample_spatial_ratio +=
          (S_kx.size('p') + S_kx.size('n')) / (S_kx.size('f'));
    }
    resample_spatial_ratio /= grid_.Nx;

    resample_spatial_ratio = para_.resample_spatial_ratio;

    // for (int kx = 0; kx < grid.Nx; kx ++) {
    int kx = 0;
    while ((flag_resample_success) && (kx < grid_.Nx)) {
      if ((S_x[kx].size('p') + S_x[kx].size('n')) >=
          resample_spatial_ratio * S_x[kx].size('f')) {
        std::cout << "Particles resampling: ( " << S_x[kx].size('p') << ", "
                  << S_x[kx].size('n') << ", " << S_x[kx].size('f') << ") "
                  << std::endl;

        flag_resample_success = resampleHomo(S_x[kx], grid_.Neff, para_.Nfreq,
                                             para_.useApproximation);

        std::cout << "After resampling: ( " << S_x[kx].size('p') << ", "
                  << S_x[kx].size('n') << ", " << S_x[kx].size('f') << ") "
                  << std::endl;
      }
      kx++;
    }
  }

  if (!flag_resample_success) {
    resampleF_inhomo(S_x, grid_.Neff_F / 2, grid_, para_.Nfreq);

    int Nx = grid_.Nx;
    std::vector<double> rho(Nx), rho_F(Nx);

    for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

    for (int kx = 0; kx < Nx; kx++)
      rho[kx] =
          S_x[kx].rhoM + (S_x[kx].m0P - S_x[kx].m0N) * grid_.Neff / grid_.dx;

    for (int kx = 0; kx < Nx; kx++)
      rho_F[kx] = S_x[kx].m0F * grid_.Neff_F / grid_.dx;
  }
}
}  // namespace coulomb