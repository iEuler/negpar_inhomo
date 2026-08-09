#include "ParticleInitialization.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Constants.h"
#include "ParticleConservation.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "SimulationConfig.h"
#include "RandomSampling.h"

namespace coulomb {
using std::abs;
using std::cout;
using std::endl;
using std::exp;
using std::pow;
using std::sqrt;
using std::string;
using std::vector;
// Generate a P, N, F particle list with designated distribution

void initialize_Negpar(NeParticleGroup &S_x, const IniValClass &inidata,
                       double Neff, double Neff_F, double dx,
                       RandomContext& random) {
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
    int Nf = myfloor(inidata.rho * dx / Neff_F, random);

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
        vf[k] = inidata.velocity[k] + sqrtT * myrandn(random);

      Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, ParticleKind::Full);
    }

    // create P and N particles

  } else if (probname == "Delta") {
    // decide the size
    // int Np = 0, Nn = 0;
    int Nf = myfloor(inidata.rho * dx / Neff_F, random);

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
      // A deterministic thermal draw would use myrandn(random).

      Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, ParticleKind::Full);
    }

  } else if (probname == "TwoStreamInstab") {
    // decide the size
    int Np, Nn, Nf;
    Np = myfloor(rhop * inidata.TSI_coe * dx / Neff, random);
    Nn = Np;
    Nf = myfloor(rhof * inidata.TSI_coe * dx / Neff_F, random);

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
    while (S_x.size(ParticleKind::Full) < Nf) {
      double v1 = (myrand(random) - .5) * 2 * vmax;
      if (myrand(random) <
          (exp(-v1 * v1 / 2) * (1 + 5 * v1 * v1) / maxf0)) {
        vp[0] = v1;
        for (int k = 1; k < 3; k++) vp[k] = myrandn(random);

        Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vp);
        S_x.push_back(S_one, ParticleKind::Full);
      }
    }

    // create P and N particles
    double coe_m0 = rhof / pow(sqrt(2. * pi * Tprt), 3);
    double sqrtT = sqrt(Tprt);
    while ((S_x.size(ParticleKind::Positive) < Np) ||
           (S_x.size(ParticleKind::Negative) < Nn)) {
      // Velocity draws use the explicit RandomContext parameter.
      std::vector<double> vp_sample{myrandn(random) * sqrtT,
                                    myrandn(random) * sqrtT,
                                    myrandn(random) * sqrtT};
      double vsq = vp_sample[0] * vp_sample[0] +
                   vp_sample[1] * vp_sample[1] +
                   vp_sample[2] * vp_sample[2];
      double f0 = exp(-vsq / 2) * (1 + 5 * vp_sample[0] * vp_sample[0]);
      double m0 = coe_m0 * exp(-vsq / 2 / Tprt);
      // cout << (f0/m0) /max_f_over_M << endl;
      if (myrand(random) < abs((f0 - m0) / m0 / max_f_over_M)) {
        if (f0 > m0) {
          if (S_x.size(ParticleKind::Positive) < Np) {
            Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vp_sample);
            S_x.push_back(S_one, ParticleKind::Positive);
          }
        } else {
          if (S_x.size(ParticleKind::Negative) < Nn) {
            Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vp_sample);
            S_x.push_back(S_one, ParticleKind::Negative);
          }
        }
      }
    }

    assign_positions(S_x, x1, x2, random);

    double m21 = inidata.TSI_coe * (inidata.TSI_m21 - rhof * Tprt);
    double m22 = inidata.TSI_coe * (inidata.TSI_m22 - rhof * Tprt);
    double m23 = inidata.TSI_coe * (inidata.TSI_m23 - rhof * Tprt);
    // cout << "Now " <<  m21 << ' ' << m22 << ' ' << m23 << endl;
    enforce_conservation(0., 0., 0., 0., m21, m22, m23, S_x, Neff, true,
                         random);

  } else if (probname == "BumpOnTail") {
    // decide the size
    int Np = myfloor(inidata.rho * dx / Neff, random);
    // int Nn = Np;
    int Nf = myfloor(inidata.rho * dx / Neff_F, random);

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
    double Tprt_bump = inidata.Tprt;
    double Tprt1 = inidata.Tprt1;
    double Tprt2 = inidata.Tprt2;

    // create F particles

    std::vector<double> vf(3);
    double center[3] = {0, 0, 0}, sigma;
    for (int kf = 0; kf < Nf; kf++) {
      if (myrand(random) < rho1 / rho) {
        center[0] = u1;
        sigma = sqrt(Tprt1);
      } else {
        center[0] = u2;
        sigma = sqrt(Tprt2);
      }
      for (int k = 0; k < 3; k++)
        vf[k] = center[k] + sigma * myrandn(random);
      Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
      S_x.push_back(S_one, ParticleKind::Full);
    }

    // create P and N particles
    // int kp = 0, kn = 0;
    double coeT0T = sqrt(Tprt_bump / Tprt1);
    coeT0T = coeT0T * coeT0T * coeT0T;
    for (int kf = 0; kf < Np; kf++) {
      if (myrand(random) < rho1 / rho) {
        center[0] = u;
        sigma = sqrt(Tprt_bump);
        for (int k = 0; k < 3; k++)
          vf[k] = center[k] + sigma * myrandn(random);
        double rho_temp =
            rho1 * coeT0T *
                exp(-(vf[0] * vf[0] + vf[1] * vf[1] + vf[2] * vf[2]) / 2 /
                        Tprt1 +
                    ((vf[0] - u) * (vf[0] - u) + vf[1] * vf[1] +
                     vf[2] * vf[2]) /
                        2 / Tprt_bump) -
            rho1 - rho2;
        if (myrand(random) < (abs(rho_temp) / rho1)) {
          if (rho_temp > 0) {
            Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
            S_x.push_back(S_one, ParticleKind::Positive);
          } else {
            Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
            S_x.push_back(S_one, ParticleKind::Negative);
          }
        }

      } else {
        // generate one particle in the high bump
        center[0] = u2;
        sigma = sqrt(Tprt2);
        for (int k = 0; k < 3; k++)
          vf[k] = center[k] + sigma * myrandn(random);
        Particle1d3d S_one(myrand(random) * (x2 - x1) + x1, vf);
        S_x.push_back(S_one, ParticleKind::Positive);
      }
    }

    // cout << rho1 << " " << rho2 << " " << u << " " << Tprt1 << " " << Tprt2
    // << " " << Tprt << endl; cout << "[ "<< kp << ", " << kn << " ]" << endl;
  }

  // cout << "bounds " <<  S_x->get_xmin() << ' ' << S_x->get_xmax() << endl;

  S_x.computemoments();
}
}  // namespace coulomb
