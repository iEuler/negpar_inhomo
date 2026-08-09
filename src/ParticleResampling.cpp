#include "ParticleResampling.h"

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationConfig.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include "Diagnostics.h"
#include "FullParticleSampling.h"
#include "ParticleGroupOperations.h"
#include "Output.h"
#include "Resampler.h"
#include "utils.h"

namespace coulomb {
using std::cout;
using std::endl;
using std::vector;

bool particleresample_homo(NeParticleGroup &S_x, const ParaClass &para,
                           SimulationState& state) {
  // // cout << " resample 0" << endl;

  state.resampleCount++;

  int Np_old = S_x.size(ParticleKind::Positive);
  int Nn_old = S_x.size(ParticleKind::Negative);

  // int Nmax = 2*max(S_x . size('p'), S_x . size('n'));

  // cout << " resample 1" << endl;

  // resample particles
  resampling::FourierResamplerConfig config;
  config.frequency_count = static_cast<std::size_t>(para.Nfreq);
  config.use_approximation = true;
  const resampling::FourierResampler resampler(S_x, config);
  auto S_x_new = resampler.resample(state.random);

  int Np_new = S_x_new.size(ParticleKind::Positive);
  int Nn_new = S_x_new.size(ParticleKind::Negative);

  // cout << " Resample finished." << endl;
  // cout << "After resampling N = (" << ptr_S_x_new .size('p') << ", " <<
  // ptr_S_x_new .size('n') << ");" << endl;
  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax(), state.random);

  S_x.isResampled = true;

  // replace old particles by new sampled particles

  // save_particles(S_x, ptr_S_x_new);

  // Replace the original particles by new sampled particles
  if ((Np_new < Np_old) && (Nn_new < Nn_old)) {
    // cout << "Replace by new sampled particles" << endl;
    S_x.clear(ParticleKind::Positive);
    S_x.clear(ParticleKind::Negative);
    merge_NeParticleGroup(S_x, S_x_new);
    return true;
  } else {
    cout << "New sampled particles rejected." << endl;
    return false;
  }
}

// void particleresample_inhomo(NeParticleGroup * S_x, NumericGridClass & grid,
// ParaClass & para, MultlLevelGroup * MLsol) {
void particleresample_inhomo(std::vector<NeParticleGroup> &S_x,
                             NumericGridClass &grid, ParaClass &para,
                             SimulationState& state) {
  bool needGlobalResample = false;

  bool flag_resample_success = true;

  for (int kx = 0; kx < grid.Nx; kx++) {
    if ((S_x[kx].size(ParticleKind::Positive) +
         S_x[kx].size(ParticleKind::Negative)) >=
        S_x[kx].size(ParticleKind::Full))
      needGlobalResample = true;
  }
  if (needGlobalResample) {
    double resample_spatial_ratio = 0;
    for (int kx = 0; kx < grid.Nx; kx++) {
      resample_spatial_ratio +=
          (S_x[kx].size(ParticleKind::Positive) +
           S_x[kx].size(ParticleKind::Negative)) /
          S_x[kx].size(ParticleKind::Full);
    }
    resample_spatial_ratio /= grid.Nx;

    resample_spatial_ratio = para.resample_spatial_ratio;

    // for (int kx = 0; kx < grid.Nx; kx ++) {
    int kx = 0;
    while ((flag_resample_success) && (kx < grid.Nx)) {
      if ((S_x[kx].size(ParticleKind::Positive) +
           S_x[kx].size(ParticleKind::Negative)) >=
          resample_spatial_ratio * S_x[kx].size(ParticleKind::Full)) {
        cout << "Particles resampling: ( "
             << S_x[kx].size(ParticleKind::Positive) << ", "
             << S_x[kx].size(ParticleKind::Negative) << ", "
             << S_x[kx].size(ParticleKind::Full) << ") " << endl;

        flag_resample_success = particleresample_homo(S_x[kx], para, state);

        cout << "After resampling: ( "
             << S_x[kx].size(ParticleKind::Positive) << ", "
             << S_x[kx].size(ParticleKind::Negative) << ", "
             << S_x[kx].size(ParticleKind::Full) << ") " << endl;
      }
      kx++;
    }
  }

  if (!flag_resample_success) {
    resampleF_inhomo(S_x, grid.Neff_F / 2, grid, para.Nfreq, state);

    int Nx = grid.Nx;
    vector<double> rho(Nx), rho_F(Nx);

    for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

    for (int kx = 0; kx < Nx; kx++)
      rho[kx] =
          S_x[kx].rhoM + (S_x[kx].positive_moments.m0 - S_x[kx].negative_moments.m0) * grid.Neff / grid.dx;

    for (int kx = 0; kx < Nx; kx++)
      rho_F[kx] = S_x[kx].full_moments.m0 * grid.Neff_F / grid.dx;

    save_macro<double>(rho, "rho_test", state);
    save_macro<double>(rho_F, "rhoF_test", state);
  }
}

void resampleF_homo(NeParticleGroup &S_x, double Neff_F_new, double Neff,
                    int Nfreq, double dx_space, RandomContext& random) {
  // resample particles
  auto S_x_new = resample_F_from_MPN(S_x, Nfreq, Neff, Neff_F_new, dx_space,
                                     random);

  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax(), random);

  // replace old particles by new sampled particles

  S_x.clear(ParticleKind::Full);
  mergeF_NeParticleGroup(S_x, S_x_new);
}

void resampleF_inhomo(std::vector<NeParticleGroup> &S_x, double Neff_F_new,
                      NumericGridClass &grid, int Nfreq,
                      SimulationState& state) {
  for (int kx = 0; kx < grid.Nx; kx++) {
    resampleF_homo(S_x[kx], Neff_F_new, grid.Neff, Nfreq, grid.dx,
                   state.random);
  }

  grid.Neff_F = Neff_F_new;
  state.syncTime = 0;

  cout << "F particle resampled." << endl;
}

void resampleF_keeptotalmass(std::vector<NeParticleGroup> &S_x,
                             NumericGridClass &grid, int Nf_old,
                             RandomContext& random) {
  int Nf_new = count_particle_number(S_x, grid.Nx, ParticleKind::Full);
  if (Nf_new > Nf_old) {
    double Neff_F_new = grid.Neff_F;
    double totalmass = 0;

    for (int kx = 0; kx < grid.Nx; kx++) {
      double mass = S_x[kx].rhoM * grid.dx;
      totalmass += mass;
      int Nk = (int)(mass / Neff_F_new);

      int Nk_remove = S_x[kx].size(ParticleKind::Full) - Nk;

      for (int kp = 0; kp < Nk_remove; kp++) {
        int k_remove = (int)(myrand(random) * S_x[kx].size(ParticleKind::Full));
        S_x[kx].erase(k_remove, ParticleKind::Full);
      }
    }

    grid.Neff_F =
        totalmass / count_particle_number(S_x, grid.Nx, ParticleKind::Full);
  }
}

void sync_coarse(std::vector<NeParticleGroup> &S_x, NumericGridClass &grid,
                 ParaClass &para, SimulationState& state) {
  if (para.collisionType == CollisionType::Coulomb) {
    if (state.syncTime > para.sync_time_interval) {
      cout << "Start resample F" << endl;

      // cout << "First resample P and N" << endl;
      state.syncTime = 0;
      particleresample_inhomo(S_x, grid, para, state);

      cout << "P and N resampled" << endl;

      /*
      double totalmass = 0;
      for (int kx = 0; kx < grid.Nx; kx ++)  totalmass += (S_x+kx) . rhoM;
      totalmass *= grid.dx;

      int Nd = 0;
      for (int kx = 0; kx < grid.Nx; kx ++)  Nd += ( (S_x+kx) . size('p') +
      (S_x+kx) . size('n') ); int Nf = (int)(1.2*Nd);

      double Neff_F_new = totalmass/Nf;
      if (Neff_F_new < grid.Neff_F) Neff_F_new = grid.Neff_F;
      */

      double Neff_F_new = 100;
      for (int kx = 0; kx < grid.Nx; kx++) {
        int N_one = (S_x[kx].size(ParticleKind::Positive) +
                     S_x[kx].size(ParticleKind::Negative));
        double Neff_F_one = (S_x[kx].rhoM) * grid.dx / N_one / 1.1;
        if (Neff_F_new > Neff_F_one) Neff_F_new = Neff_F_one;
      }

      if (Neff_F_new < grid.Neff_F) Neff_F_new = grid.Neff_F;

      // cout << "s resample F" << endl;

      int Nf_old = count_particle_number(S_x, grid.Nx, ParticleKind::Full);

      resampleF_inhomo(S_x, Neff_F_new, grid, para.Nfreq, state);
      cout << "F resampled" << endl;
      resampleF_keeptotalmass(S_x, grid, Nf_old, state.random);
    }
  }
}

}  // namespace coulomb
