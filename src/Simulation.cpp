
#include <fftw3.h>
#include <omp.h>

#include <cmath>
#include <complex>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "_global_variables.h"
#include "Advection.h"
#include "Diagnostics.h"
#include "ElectricField.h"
#include "Initialization.h"
#include "LegacyResampling.h"
#include "Moments.h"
#include "NegativeParticle.h"
#include "Numerics.h"
#include "Output.h"
#include "RunOptions.h"
#include "Simulation.h"

//#define I complex<double>(0.,1.)

// g++ *.cpp -o out -lfftw3 -std=c++17 -fopenmp

// ===============================================================
// Main program

using namespace coulomb;
using std::cout;
using std::endl;

void coulomb::save_macro_evolution(std::vector<NeParticleGroup>& groups,
                                   const NumericGridClass& grid,
                                   const SimulationState& state) {
  for (int cell = 0; cell < grid.Nx; ++cell) groups[cell].computemoments();
  update_macro(groups, grid);
  save_rhouT(groups, grid, state);
  save_rhouT_F(groups, grid, state);
  save_E(groups, grid, state);
  save_dist(groups, grid, state);
  save_rhouT_maxwellian(groups, grid, state);

  std::ofstream file(output_path("numOfDist.txt", state));
  if (!file) throw std::runtime_error("Unable to open distribution count output");
  file << state.saveIndex;
}

int run_simulation_impl(const coulomb::RunOptions& options,
                        coulomb::SimulationState& state) {
  // cout << "Programe start time: " << asctime(localtime(&t_start)) << endl;

  // clock_t t0;
  // t0 = clock();

  // save_homo_rdist();

  // return 0;

  ParaClass para;
  NumericGridClass grid(100, para.method);
  if (options.steps) grid.Nt = *options.steps;

  // int kt_convtest = 2;

  // grid.Nt = 100*kt_convtest;
  // grid.dt = grid.dt / kt_convtest;

  // grid.Nt = 10*kt_convtest;
  // grid.dt = 0.2 / kt_convtest;
  // grid.Nt = 0;

  // grid.Nt = 2*grid.Nt;
  // grid.dt = .5*grid.dt;

  para.dt = grid.dt;
  grid.lambda_Poisson = para.lambda_Poisson;

  save_grids(grid, state);
  saveparameter(para, grid, state);

  std::vector<NeParticleGroup> S_x(grid.Nx);
  NeParticleGroup* ptr_S_x = &S_x[0];

  initialize_distri_Negpar(grid, S_x, state);
  para.lambda_Poisson = grid.lambda_Poisson;

  cout << "method = " << method_name(para.method) << endl;

  state.filenameWithNumber = false;

  int kt100 = 0, kt10 = 0;

  std::vector<double> elec_energy;
  std::vector<double> elec_energy_F;
  std::vector<double> total_energy;
  std::vector<double> total_energy_F;
  std::vector<double> Neff_F_rec;

  std::vector<double> time_dist;

  std::vector<double> cputime_adve, cputime_coll, cputime_all, cputime_resamp;

  std::vector<int> Np_rec;
  std::vector<int> Nn_rec;
  std::vector<int> Nf_rec;
  std::vector<int> num_resample_rec;

  update_macro(S_x, grid);
  updateelecfiled(S_x, grid);

  state.syncTime = 0;

  for (int kt = 0; kt < grid.Nt; kt++) {
    std::cout << "step " << kt << '/' << grid.Nt << endl;

    if (kt >= kt10) {
      state.filenameWithNumber = true;
      save_macro_evolution(S_x, grid, state);
      state.saveIndex++;
      kt10 += 40;
      state.filenameWithNumber = false;
      time_dist.push_back(kt * grid.dt);
    }

    // std::cout << "a" << std::endl;
    elec_energy.push_back(compute_elec_energy(S_x, grid));
    elec_energy_F.push_back(compute_elec_energy_F(S_x, grid));
    total_energy.push_back(compute_total_energy(S_x, grid));
    total_energy_F.push_back(compute_total_energy_F(S_x, grid));
    Np_rec.push_back(count_particle_number(S_x, grid.Nx,
                                            ParticleKind::Positive));
    Nn_rec.push_back(count_particle_number(S_x, grid.Nx,
                                            ParticleKind::Negative));
    Nf_rec.push_back(count_particle_number(S_x, grid.Nx, ParticleKind::Full));
    Neff_F_rec.push_back(grid.Neff_F);

    num_resample_rec.push_back(state.resampleCount);
    state.resampleCount = 0;
    // cout << "Energy = (" << elec_energy[kt] << ", " << elec_energy_F[kt] <<
    // ", "
    //      << total_energy[kt] << ", " << total_energy_F[kt] << ")" << endl;

    state.t0All = clock();

    if (para.method == SimulationMethod::HDP)
      Negpar_inhomo_onestep(S_x, grid, para, state);
    else
      Negpar_inhomo_onestep_PIC(S_x, grid, para, state);

    state.syncTime += grid.dt;

    // Negpar_inhomo_onestep_PIC(ptr_S_x, grid, para);

    state.t1All = clock();

    cputime_all.push_back(((float)(state.t1All - state.t0All)) / CLOCKS_PER_SEC);
    cputime_adve.push_back(((float)(state.t1Advection - state.t0Advection)) / CLOCKS_PER_SEC);
    cputime_coll.push_back(((float)(state.t1Collision - state.t0Collision)) / CLOCKS_PER_SEC);
    cputime_resamp.push_back(((float)(state.t1Resampling - state.t0Resampling)) / CLOCKS_PER_SEC);

    if (kt >= kt100) {
      state.filenameWithNumber = false;
      save_macro<double>(elec_energy, "elec_energy", state);
      save_macro<double>(elec_energy_F, "elec_energy_F", state);
      save_macro<double>(total_energy, "total_energy", state);
      save_macro<double>(total_energy_F, "total_energy_F", state);
      save_macro<double>(cputime_all, "cputime_all", state);
      save_macro<double>(cputime_adve, "cputime_adve", state);
      save_macro<double>(cputime_coll, "cputime_coll", state);
      save_macro<double>(cputime_resamp, "cputime_resamp", state);
      save_macro<int>(Np_rec, "Np_rec", state);
      save_macro<int>(Nn_rec, "Nn_rec", state);
      save_macro<int>(Nf_rec, "Nf_rec", state);
      save_macro<double>(Neff_F_rec, "Neff_F_rec", state);
      save_macro<double>(time_dist, "time_dist", state);
      save_macro<int>(num_resample_rec, "num_resample", state);

      state.filenameWithNumber = true;

      kt100 += 50;
    }

    // if (kt == 4000) return 0;
  }

  state.filenameWithNumber = false;
  save_macro<double>(elec_energy, "elec_energy", state);
  save_macro<double>(elec_energy_F, "elec_energy_F", state);
  save_macro<int>(Np_rec, "Np_rec", state);
  save_macro<int>(Nn_rec, "Nn_rec", state);
  save_macro<double>(Neff_F_rec, "Neff_F_rec", state);
  save_macro<int>(num_resample_rec, "num_resample", state);

  save_macro_evolution(S_x, grid, state);

  state.filenameWithNumber = false;
  save_macro<double>(elec_energy, "elec_energy", state);
  save_macro<double>(elec_energy_F, "elec_energy_F", state);
  save_macro<double>(total_energy, "total_energy", state);
  save_macro<double>(total_energy_F, "total_energy_F", state);
  save_macro<double>(cputime_all, "cputime_all", state);
  save_macro<double>(cputime_adve, "cputime_adve", state);
  save_macro<double>(cputime_coll, "cputime_coll", state);
  save_macro<double>(cputime_resamp, "cputime_resamp", state);
  save_macro<double>(time_dist, "time_dist", state);
  save_macro<int>(Np_rec, "Np_rec", state);
  save_macro<int>(Nn_rec, "Nn_rec", state);
  save_macro<int>(Nf_rec, "Nf_rec", state);
  save_macro<int>(num_resample_rec, "num_resample", state);

  cout << "Finished" << endl;
  // state.t1All = clock();
  // cout << "Elapsed time = " << ((float)(state.t1All - state.t0All))/CLOCKS_PER_SEC <<
  // endl;

  // getchar();

  return 0;
}

coulomb::Simulation::Simulation(coulomb::RunOptions options,
                                coulomb::SimulationState& state)
    : options_(options), state_(state) {}

int coulomb::Simulation::run() { return run_simulation_impl(options_, state_); }

int coulomb::run_simulation(const RunOptions& options,
                            SimulationState& state) {
  apply_run_options(options, state);
  return Simulation(options, state).run();
}
