
#include <fftw3.h>
#include <omp.h>

#include <cmath>
#include <complex>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "coulomb_class.h"
#include "coulomb_function.h"
#include "coulomb_initialize.h"
#include "coulomb_negpar.h"
#include "coulomb_resample.h"
#include "coulomb_save.h"

//#define I complex<double>(0.,1.)

// g++ -Wall -fopenmp inhomo_neg_coulomb.cpp -o prog -std=c++11 -lfftw3

// ===============================================================
// Main program

int main() {
  int num_threads = omp_get_max_threads();
  omp_set_num_threads(num_threads - 2);

  time_t t_start = time(NULL);
  // cout << "Programe start time: " << asctime(localtime(&t_start)) << endl;

  // clock_t t0;
  // t0 = clock();

  // save_homo_rdist();

  // return 0;

  ParaClass para;
  NumericGridClass grid(400, para.method);

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

  save_grids(grid);
  saveparameter(para, grid);

  vector<NeParticleGroup> S_x(grid.Nx);
  NeParticleGroup* ptr_S_x = &S_x[0];

  initialize_distri_Negpar(grid, ptr_S_x);
  para.lambda_Poisson = grid.lambda_Poisson;

  cout << "method = " << para.method << endl;

  FLAG_FILENAME_WITH_NUM = false;

  int kt100 = 0, kt10 = 0;

  vector<double> elec_energy;
  vector<double> elec_energy_F;
  vector<double> total_energy;
  vector<double> total_energy_F;
  vector<double> Neff_F_rec;

  vector<double> time_dist;

  vector<double> cputime_adve, cputime_coll, cputime_all, cputime_resamp;

  vector<int> Np_rec;
  vector<int> Nn_rec;
  vector<int> Nf_rec;
  vector<int> num_resample_rec;

  update_macro(ptr_S_x, grid);
  updateelecfiled(ptr_S_x, grid);

  SYNC_TIME = 0;

  for (int kt = 0; kt < grid.Nt; kt++) {
    cout << kt << ' ' << grid.Nt << endl;

    if (kt >= kt10) {
      FLAG_FILENAME_WITH_NUM = true;
      save_macro_evolution(ptr_S_x, grid);
      K_SAVE_TIME++;
      kt10 += 40;
      FLAG_FILENAME_WITH_NUM = false;
      time_dist.push_back(kt * grid.dt);
    }

    std::cout << "a" << std::endl;
    elec_energy.push_back(compute_elec_energy(ptr_S_x, grid));
    elec_energy_F.push_back(compute_elec_energy_F(ptr_S_x, grid));
    total_energy.push_back(compute_total_energy(ptr_S_x, grid));
    total_energy_F.push_back(compute_total_energy_F(ptr_S_x, grid));
    Np_rec.push_back(count_particle_number(ptr_S_x, grid.Nx, 'p'));
    Nn_rec.push_back(count_particle_number(ptr_S_x, grid.Nx, 'n'));
    Nf_rec.push_back(count_particle_number(ptr_S_x, grid.Nx, 'f'));
    Neff_F_rec.push_back(grid.Neff_F);

    num_resample_rec.push_back(NUM_RESAMPLE);
    NUM_RESAMPLE = 0;
    cout << "Energy = (" << elec_energy[kt] << ", " << elec_energy_F[kt] << ", "
         << total_energy[kt] << ", " << total_energy_F[kt] << ")" << endl;

    t0_all = clock();

    if (para.method == "HDP")
      Negpar_inhomo_onestep(ptr_S_x, grid, para);
    else
      Negpar_inhomo_onestep_PIC(ptr_S_x, grid, para);

    SYNC_TIME += grid.dt;

    // Negpar_inhomo_onestep_PIC(ptr_S_x, grid, para);

    t1_all = clock();

    cputime_all.push_back(((float)(t1_all - t0_all)) / CLOCKS_PER_SEC);
    cputime_adve.push_back(((float)(t1_adve - t0_adve)) / CLOCKS_PER_SEC);
    cputime_coll.push_back(((float)(t1_coll - t0_coll)) / CLOCKS_PER_SEC);
    cputime_resamp.push_back(((float)(t1_resamp - t0_resamp)) / CLOCKS_PER_SEC);

    if (kt >= kt100) {
      FLAG_FILENAME_WITH_NUM = false;
      save_macro<double>(elec_energy, "elec_energy");
      save_macro<double>(elec_energy_F, "elec_energy_F");
      save_macro<double>(total_energy, "total_energy");
      save_macro<double>(total_energy_F, "total_energy_F");
      save_macro<double>(cputime_all, "cputime_all");
      save_macro<double>(cputime_adve, "cputime_adve");
      save_macro<double>(cputime_coll, "cputime_coll");
      save_macro<double>(cputime_resamp, "cputime_resamp");
      save_macro<int>(Np_rec, "Np_rec");
      save_macro<int>(Nn_rec, "Nn_rec");
      save_macro<int>(Nf_rec, "Nf_rec");
      save_macro<double>(Neff_F_rec, "Neff_F_rec");
      save_macro<double>(time_dist, "time_dist");
      save_macro<int>(num_resample_rec, "num_resample");

      FLAG_FILENAME_WITH_NUM = true;

      kt100 += 50;
    }

    // if (kt == 4000) return 0;
  }

  FLAG_FILENAME_WITH_NUM = false;
  save_macro<double>(elec_energy, "elec_energy");
  save_macro<double>(elec_energy_F, "elec_energy_F");
  save_macro<int>(Np_rec, "Np_rec");
  save_macro<int>(Nn_rec, "Nn_rec");
  save_macro<double>(Neff_F_rec, "Neff_F_rec");
  save_macro<int>(num_resample_rec, "num_resample");

  save_macro_evolution(ptr_S_x, grid);

  FLAG_FILENAME_WITH_NUM = false;
  save_macro<double>(elec_energy, "elec_energy");
  save_macro<double>(elec_energy_F, "elec_energy_F");
  save_macro<double>(total_energy, "total_energy");
  save_macro<double>(total_energy_F, "total_energy_F");
  save_macro<double>(cputime_all, "cputime_all");
  save_macro<double>(cputime_adve, "cputime_adve");
  save_macro<double>(cputime_coll, "cputime_coll");
  save_macro<double>(cputime_resamp, "cputime_resamp");
  save_macro<double>(time_dist, "time_dist");
  save_macro<int>(Np_rec, "Np_rec");
  save_macro<int>(Nn_rec, "Nn_rec");
  save_macro<int>(Nf_rec, "Nf_rec");
  save_macro<int>(num_resample_rec, "num_resample");

  cout << "Finished" << endl;
  // t1_all = clock();
  // cout << "Elapsed time = " << ((float)(t1_all - t0_all))/CLOCKS_PER_SEC <<
  // endl;

  // getchar();

  return 0;
}
