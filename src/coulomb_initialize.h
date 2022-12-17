#include <fstream>

#include "Classes.h"

namespace coulomb {
void save_initial(IniValClass &inidata);

// specify the initial macro
void initialize_Negpar(NeParticleGroup &S_x, const IniValClass &inidata,
                       double Neff, double Neff_F, double dx);
void initialize_TwoStreamInstab(IniValClass &inidata);
void update_macro(std::vector<NeParticleGroup> &, const NumericGridClass &grid);

/** specify the initial macroscopic quantities for NegPar method
    @param grid contains numerical paramters on grids
    @param S_x contains the P, N, F particles information
*/

void initialize_distri_Negpar(NumericGridClass &grid,
                              std::vector<NeParticleGroup> &S_x) {
  int Nx = grid.Nx;
  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> u2(Nx);
  vector<double> u3(Nx);
  vector<double> Tprt(Nx);
  vector<double> x = grid.x;
  double xk;
  double xcenter = (x[0] + x[Nx - 1]) / 2.0;

  double vol_x = grid.xmax - grid.xmin;

  // ofstream file0;
  // file0.open("IniParticle.txt");

  IniValClass inidata;

  // int Np, Nn, Nf;
  double dx = x[1] - x[0];

  inidata.probname = "LandauDamping";
  // inidata.probname = "TwoStreamInstab";
  // inidata.probname = "test1";
  // inidata.probname = "BumpOnTail";
  // inidata.probname = "Analytic";
  // inidata.probname = "Efficiency";	// for efficiency test

  if (inidata.probname == "Efficiency") {
    grid.lambda_Poisson = 0;
  }

  if (inidata.probname == "TwoStreamInstab") {
    initialize_TwoStreamInstab(inidata);
  }

  cout << "Problem name = " << inidata.probname << endl;

  inidata.totalmass = 0;
  // grid.bdry_x = 'n';
  for (int kx = 0; kx < Nx; kx++) {
    xk = x[kx];

    if ((inidata.probname == "LandauDamping") ||
        (inidata.probname == "Efficiency")) {  // Landau dampling
      double LD_alpha = inidata.LD_alpha;
      rho[kx] = 1 + LD_alpha * sin(xk);
      // rho[kx] = 1 + 0.4*sin(xk);
      u1[kx] = 0.;
      u2[kx] = 0.;
      u3[kx] = 0.;
      Tprt[kx] = 1.0;

      inidata.rho = rho[kx];
      inidata.velocity[0] = u1[kx];
      inidata.velocity[1] = u2[kx];
      inidata.velocity[2] = u3[kx];
      inidata.Tprt = Tprt[kx];

      inidata.totalmass += rho[kx] * dx;

    } else if (inidata.probname == "TwoStreamInstab") {
      double alpha_TSI = 0.5;
      double k_TSI = 0.5;

      inidata.TSI_coe =
          1. / 12 / pi * (1.0 + alpha_TSI * cos(k_TSI * xk)) / 40.;
    } else if (inidata.probname == "BumpOnTail") {
      double BOT_beta = inidata.BOT_beta;
      double BOT_rho0 = inidata.BOT_rho0;
      double BOT_Tprt = inidata.BOT_Tprt;
      double BOT_dTprt = inidata.BOT_dTprt;
      double BOT_Tx = inidata.BOT_Tx;
      double BOT_ub = inidata.BOT_ub;

      double rho0 = BOT_rho0 / vol_x * BOT_beta;
      // double rho1 = BOT_rho0 *
      // (1-BOT_beta)/sqrt(2*pi*BOT_Tx)*exp(-(xk-xcenter)*(xk-xcenter)/2/BOT_Tx);
      double rho1 = BOT_rho0 * (1 - BOT_beta) / vol_x *
                    exp(-(xk - xcenter) * (xk - xcenter) / 2 / BOT_Tx);
      double rho_tot = rho0 + rho1;
      double u0 = 0.;
      double u1_tot = BOT_ub;
      double m1 = rho1 * BOT_ub;
      double u = m1 / rho_tot;
      double Tprt0 = BOT_Tprt;
      double Tprt1 = BOT_dTprt;
      double energy =
          1.5 * rho0 * Tprt0 + rho1 * (.5 * u1_tot * u1_tot + 1.5 * Tprt1);
      double Tprt_tot = (energy - .5 * rho_tot * u * u) / 1.5 / rho_tot;

      rho[kx] = rho_tot;
      u1[kx] = u;
      u2[kx] = 0.;
      u3[kx] = 0.;
      Tprt[kx] = Tprt_tot;

      inidata.rho = rho[kx];
      inidata.velocity[0] = u1[kx];
      inidata.velocity[1] = u2[kx];
      inidata.velocity[2] = u3[kx];
      inidata.Tprt = Tprt[kx];

      inidata.rho1 = rho0;
      inidata.velocity1[0] = u0;
      inidata.velocity1[1] = 0.;
      inidata.velocity1[2] = 0.;
      inidata.Tprt1 = Tprt0;

      inidata.rho2 = rho1;
      inidata.velocity2[0] = u1_tot;
      inidata.velocity2[1] = 0.;
      inidata.velocity2[2] = 0.;
      inidata.Tprt2 = Tprt1;

      inidata.totalmass += rho[kx] * dx;
    }

    S_x[kx].set_xrange(xk - dx / 2, xk + dx / 2);
    S_x[kx].reset_flag_resampled();
    initialize_Negpar(S_x[kx], inidata, grid.Neff, grid.Neff_F, grid.dx);
  }

  save_initial(inidata);

  update_macro(S_x, grid);
}

void initialize_distri_Negpar_test(NumericGridClass &grid,
                                   std::vector<NeParticleGroup> &S_x) {
  int Nx = grid.Nx;
  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> u2(Nx);
  vector<double> u3(Nx);
  vector<double> Tprt(Nx);
  vector<double> x = grid.x;
  double xk;

  // ofstream file0;
  // file0.open("IniParticle.txt");

  IniValClass inidata;

  // int Np, Nn, Nf;
  double dx = x[1] - x[0];

  inidata.probname = "Delta";

  // grid.bdry_x = 'n';
  for (int kx = 0; kx < Nx; kx++) {
    xk = x[kx];

    // rho[kx] = 1 + 0.4*sin(xk);
    rho[kx] = 1.;
    if (xk > pi) rho[kx] = 0.5;

    u1[kx] = 1.;
    u2[kx] = 0.;
    u3[kx] = 0.;
    Tprt[kx] = 1.0;

    inidata.rho = rho[kx];
    inidata.velocity[0] = u1[kx];
    inidata.velocity[1] = u2[kx];
    inidata.velocity[2] = u3[kx];
    inidata.Tprt = Tprt[kx];

    S_x[kx].set_xrange(xk - dx / 2, xk + dx / 2);
    S_x[kx].reset_flag_resampled();
    initialize_Negpar(S_x[kx], inidata, grid.Neff, grid.Neff_F, grid.dx);
  }

  update_macro(S_x, grid);
}

// ========================================================================

/** Precompute in the Two Stream Instablity problem:
   rhof = ||f||_1 is the mass
   rhop = 1/2 * ||f-M||_1 is the deviational mass
   Tprt is the Temperature
   max_f_minus_M = max_v f/M
   Here f = exp(-|v|^2 / 2) * (1+5 v_x^2)
*/

void enforce_conservation(double m0, double m11, double m12, double m13,
                          double m21, double m22, double m23,
                          NeParticleGroup *S_new, double Neff,
                          bool flag_conserve_energyvector);

void initialize_TwoStreamInstab(IniValClass &inidata) {
  double vmax = 6.;
  int Nv = 200;
  double dv = 2.0 * vmax / Nv;
  vector<double> v(Nv);
  for (int kv = 0; kv < Nv; kv++) v[kv] = (kv + 0.5) * dv - vmax;

  double v1, v2, v3, energyf = 0.;
  vector<double> M0(Nv * Nv * Nv);
  vector<double> f0(Nv * Nv * Nv);
  double rhof = 0.;
  for (int kv1 = 0; kv1 < Nv; kv1++) {
    v1 = v[kv1];
    for (int kv2 = 0; kv2 < Nv; kv2++) {
      v2 = v[kv2];
      for (int kv3 = 0; kv3 < Nv; kv3++) {
        v3 = v[kv3];
        int kk = kv3 + Nv * (kv2 + Nv * kv1);
        double vsq = v1 * v1 + v2 * v2 + v3 * v3;
        f0[kk] = exp(-vsq / 2.) * (1.0 + 5.0 * v1 * v1);
        rhof += f0[kk];
        energyf += .5 * vsq * f0[kk];
      }
    }
  }
  rhof *= dv * dv * dv;
  energyf *= dv * dv * dv;
  double Tprt = energyf / (1.5 * rhof);

  double rhop = 0.;
  // max_f_minus_M = 0.;
  double max_f_over_M = 0.;
  double m21 = 0., m22 = 0., m23 = 0.;
  for (int kv1 = 0; kv1 < Nv; kv1++) {
    v1 = v[kv1];
    for (int kv2 = 0; kv2 < Nv; kv2++) {
      v2 = v[kv2];
      for (int kv3 = 0; kv3 < Nv; kv3++) {
        v3 = v[kv3];
        int kk = kv3 + Nv * (kv2 + Nv * kv1);
        double vsq = v1 * v1 + v2 * v2 + v3 * v3;
        M0[kk] = rhof / pow(sqrt(2. * pi * Tprt), 3) * exp(-vsq / 2. / Tprt);
        // max_f_minus_M = max(max_f_minus_M, abs(f0[kk] - M0[kk]) );
        max_f_over_M = max(max_f_over_M, abs(f0[kk] - M0[kk]) / M0[kk]);
        m21 += v1 * v1 * f0[kk];
        m22 += v2 * v2 * f0[kk];
        m23 += v3 * v3 * f0[kk];
        if (f0[kk] > M0[kk]) rhop += f0[kk] - M0[kk];
      }
    }
  }

  rhop *= dv * dv * dv;
  m21 *= dv * dv * dv;
  m22 *= dv * dv * dv;
  m23 *= dv * dv * dv;

  inidata.TSI_rhof = rhof;
  inidata.TSI_rhop = rhop;
  inidata.TSI_Tprt = Tprt;
  inidata.TSI_max_f_over_M = max_f_over_M;

  inidata.TSI_m21 = m21;
  inidata.TSI_m22 = m22;
  inidata.TSI_m23 = m23;

  cout << "Initially " << m21 + m22 + m23 << " vs " << 3 * rhof * Tprt << endl;
}

// Generate a P, N, F particle list with designated distribution

void assign_positions(NeParticleGroup &S_new, double xmin, double xmax);

void initialize_Negpar(NeParticleGroup &S_x, const IniValClass &inidata,
                       double Neff, double Neff_F, double dx) {
  // initialize_Negpar_size(int &Np, int &Nn, int &Nf);

  string probname = inidata.probname;

  double rhof, rhop, Tprt, max_f_over_M;
  rhof = inidata.TSI_rhof;
  rhop = inidata.TSI_rhop;
  Tprt = inidata.TSI_Tprt;
  max_f_over_M = inidata.TSI_max_f_over_M;

  /*
  int Np = S_x->size('p');
  int Nn = S_x->size('n');
  int Nf = S_x->size('f');
  */

  double x1, x2;
  x1 = S_x.get_xmin();
  x2 = S_x.get_xmax();

  if ((probname == "LandauDamping") || (probname == "Efficiency")) {
    // decide the size
    // int Np = 0, Nn = 0;
    int Nf = myfloor(inidata.rho * dx / Neff_F);

    // Particle1d3d * Sp = S_x->list('p');
    // Particle1d3d * Sn = S_x->list('n');
    // Particle1d3d * Sf = S_x->list('f');

    // update maxwellian part

    S_x.rhoM = inidata.rho;
    S_x.u1M = inidata.velocity[0];
    S_x.u2M = inidata.velocity[1];
    S_x.u3M = inidata.velocity[2];
    S_x.TprtM = inidata.Tprt;

    // create F particles

    std::vector<double> vf(3);

    double sqrtT = sqrt(inidata.Tprt);
    for (int kf = 0; kf < Nf; kf++) {
      for (int k = 0; k < 3; k++)
        vf[k] = inidata.velocity[k] + sqrtT * myrandn();

      Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, 'f');
    }

    // create P and N particles

  } else if (probname == "Delta") {
    // decide the size
    // int Np = 0, Nn = 0;
    int Nf = myfloor(inidata.rho * dx / Neff_F);

    // Particle1d3d * Sp = S_x->list('p');
    // Particle1d3d * Sn = S_x->list('n');
    // Particle1d3d * Sf = S_x->list('f');

    // update maxwellian part

    S_x.rhoM = inidata.rho;
    S_x.u1M = inidata.velocity[0];
    S_x.u2M = inidata.velocity[1];
    S_x.u3M = inidata.velocity[2];
    S_x.TprtM = inidata.Tprt;

    // create F particles

    std::vector<double> vf(3);

    for (int kf = 0; kf < Nf; kf++) {
      for (int k = 0; k < 3; k++) vf[k] = inidata.velocity[k];
      // for (int k=0;k<3;k++) vf[k] = inidata.velocity[k] +
      // sqrt(inidata.Tprt)*myrandn();

      Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, 'f');
    }

  } else if (probname == "TwoStreamInstab") {
    // decide the size
    int Np, Nn, Nf;
    Np = myfloor(rhop * inidata.TSI_coe * dx / Neff);
    Nn = Np;
    Nf = myfloor(rhof * inidata.TSI_coe * dx / Neff_F);

    cout << Np << ' ' << Nn << ' ' << Nf << endl;

    // Particle1d3d * Sp = S_x->list('p');
    // Particle1d3d * Sn = S_x->list('n');
    // Particle1d3d * Sf = S_x->list('f');

    // update maxwellian part

    S_x.rhoM = rhof * inidata.TSI_coe;
    S_x.u1M = 0.;
    S_x.u2M = 0.;
    S_x.u3M = 0.;
    S_x.TprtM = Tprt;

    // create F particles

    double v_sq = 1.8;
    double maxf0 = exp(-v_sq / 2) * (1 + 5 * v_sq);
    std::vector<double> vp(3);
    double vmax = 6.;

    // int kf = 0;
    while (S_x.size('f') < Nf) {
      double v1 = (myrand() - .5) * 2 * vmax;
      if (myrand() < (exp(-v1 * v1 / 2) * (1 + 5 * v1 * v1) / maxf0)) {
        vp[0] = v1;
        for (int k = 1; k < 3; k++) vp[k] = myrandn();

        Particle1d3d S_one(myrand() * (x2 - x1) + x1, vp);
        S_x.push_back(S_one, 'f');
      }
    }

    // create P and N particles
    double coe_m0 = rhof / pow(sqrt(2. * pi * Tprt), 3);
    double sqrtT = sqrt(Tprt);
    while ((S_x.size('p') < Np) || (S_x.size('n') < Nn)) {
      // double vp[3] = { (myrand()-.5)*2*vmax, (myrand()-.5)*2*vmax,
      // (myrand()-.5)*2*vmax};
      std::vector<double> vp{myrandn() * sqrtT, myrandn() * sqrtT,
                             myrandn() * sqrtT};
      double vsq = vp[0] * vp[0] + vp[1] * vp[1] + vp[2] * vp[2];
      double f0 = exp(-vsq / 2) * (1 + 5 * vp[0] * vp[0]);
      double m0 = coe_m0 * exp(-vsq / 2 / Tprt);
      // cout << (f0/m0) /max_f_over_M << endl;
      if (myrand() < abs((f0 - m0) / m0 / max_f_over_M)) {
        if (f0 > m0) {
          if (S_x.size('p') < Np) {
            Particle1d3d S_one(myrand() * (x2 - x1) + x1, vp);
            S_x.push_back(S_one, 'p');
          }
        } else {
          if (S_x.size('n') < Nn) {
            Particle1d3d S_one(myrand() * (x2 - x1) + x1, vp);
            S_x.push_back(S_one, 'n');
          }
        }
      }
    }

    assign_positions(S_x, x1, x2);

    double m21 = inidata.TSI_coe * (inidata.TSI_m21 - rhof * Tprt);
    double m22 = inidata.TSI_coe * (inidata.TSI_m22 - rhof * Tprt);
    double m23 = inidata.TSI_coe * (inidata.TSI_m23 - rhof * Tprt);
    // cout << "Now " <<  m21 << ' ' << m22 << ' ' << m23 << endl;
    enforce_conservation(0., 0., 0., 0., m21, m22, m23, &S_x, Neff, true);

  } else if (probname == "BumpOnTail") {
    // decide the size
    int Np = myfloor(inidata.rho * dx / Neff);
    // int Nn = Np;
    int Nf = myfloor(inidata.rho * dx / Neff_F);

    // Particle1d3d * Sp = S_x->list('p');
    // Particle1d3d * Sn = S_x->list('n');
    // Particle1d3d * Sf = S_x->list('f');
    // Particle1d3d * Sp = S_x->list('p');
    // Particle1d3d * Sn = S_x->list('n');

    // update maxwellian part

    S_x.rhoM = inidata.rho;
    S_x.u1M = inidata.velocity[0];
    S_x.u2M = inidata.velocity[1];
    S_x.u3M = inidata.velocity[2];
    S_x.TprtM = inidata.Tprt;

    double rho = inidata.rho;
    double rho1 = inidata.rho1;
    double rho2 = inidata.rho2;
    double u = inidata.velocity[0];
    double u1 = inidata.velocity1[0];
    double u2 = inidata.velocity2[0];
    double Tprt = inidata.Tprt;
    double Tprt1 = inidata.Tprt1;
    double Tprt2 = inidata.Tprt2;

    // create F particles

    std::vector<double> vf(3);
    double center[3] = {0, 0, 0}, sigma;
    for (int kf = 0; kf < Nf; kf++) {
      if (myrand() < rho1 / rho) {
        center[0] = u1;
        sigma = sqrt(Tprt1);
      } else {
        center[0] = u2;
        sigma = sqrt(Tprt2);
      }
      for (int k = 0; k < 3; k++) vf[k] = center[k] + sigma * myrandn();
      Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, 'f');
    }

    // create P and N particles
    // int kp = 0, kn = 0;
    double coeT0T = sqrt(Tprt / Tprt1);
    coeT0T = coeT0T * coeT0T * coeT0T;
    for (int kf = 0; kf < Np; kf++) {
      if (myrand() < rho1 / rho) {
        center[0] = u;
        sigma = sqrt(Tprt);
        for (int k = 0; k < 3; k++) vf[k] = center[k] + sigma * myrandn();
        double rho_temp =
            rho1 * coeT0T *
                exp(-(vf[0] * vf[0] + vf[1] * vf[1] + vf[2] * vf[2]) / 2 /
                        Tprt1 +
                    ((vf[0] - u) * (vf[0] - u) + vf[1] * vf[1] +
                     vf[2] * vf[2]) /
                        2 / Tprt) -
            rho1 - rho2;
        if (myrand() < (abs(rho_temp) / rho1)) {
          if (rho_temp > 0) {
            Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
            S_x.push_back(S_one, 'p');
          } else {
            Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
            S_x.push_back(S_one, 'n');
          }
        }

      } else {
        // generate one particle in the high bump
        center[0] = u2;
        sigma = sqrt(Tprt2);
        for (int k = 0; k < 3; k++) vf[k] = center[k] + sigma * myrandn();
        Particle1d3d S_one(myrand() * (x2 - x1) + x1, vf);
        S_x.push_back(S_one, 'p');
      }
    }

    // cout << rho1 << " " << rho2 << " " << u << " " << Tprt1 << " " << Tprt2
    // << " " << Tprt << endl; cout << "[ "<< kp << ", " << kn << " ]" << endl;
  }

  // cout << "bounds " <<  S_x->get_xmin() << ' ' << S_x->get_xmax() << endl;

  S_x.computemoments();
}
}  // namespace coulomb