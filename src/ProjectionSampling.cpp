#include "ProjectionSampling.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "ParticleGroupOperations.h"
#include "Constants.h"
#include "Grid.h"
#include "ParticleGroup.h"
#include "utils.h"

namespace coulomb {
using std::abs;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::pow;
using std::sqrt;
using std::vector;

void sample_from_P3M_coeff_ver3(NeParticleGroup &S_x, double dt, double dx,
                                double &a0, double &a11, double &a2,
                                double &a21, double &a31) {
  double rhoM = S_x.rhoM;
  double u1M = S_x.u1M;
  double TprtM = S_x.TprtM;

  // double dx_rhoM = S_x . dx_rhoM;
  double dx_u1M = S_x.dx_u1M;
  double dx_TprtM = S_x.dx_TprtM;

  const double dimen = 3.;

  // coefficients from -\Delta t (I-\Pi_M) (v\cdot\nabla_x  + E \cdot\nabla_v) M

  double sqrtT = sqrt(TprtM);

  a0 = 0.;
  a11 = -dt * rhoM * (-(dimen + 2) / 2. / sqrtT * dx_TprtM);
  a21 = -dt * rhoM * dx_u1M;
  a2 = -dt * rhoM * (-dx_u1M / dimen);
  a31 = -dt * rhoM * (dx_TprtM / sqrtT / 2.);

  // coefficients from \Delta t \Pi_M (v\cdot\nabla_x  + E \cdot\nabla_v) (f_p -
  // f_n)
  S_x.computemoments();

  // inner product with 1
  double drho = S_x.drho_g;
  double coe_0 = drho;
  // inner product with (v_1 - u_1)
  double dm1 = S_x.dm1_g;  // inner product with v_1
  double coe_1 = dm1 - u1M * drho;
  // inner product with ( |v - u|^2/T - d )
  double denergy = S_x.denergy_g;  // inner product with |v|^2
  double coe_2 = 2. / TprtM * denergy - 2. * u1M / TprtM * dm1 +
                 (u1M * u1M / TprtM - dimen) * drho;

  // cout <<"a11: " << a11 << ' ' << coe_1 / sqrt(TprtM) << endl;
  // cout <<"a2: " << a2 << ' ' << coe_2 / 2. / dimen << endl;

  a0 += coe_0 - .5 * coe_2;
  a11 += coe_1 / sqrt(TprtM);
  a2 += coe_2 / 2. / dimen;

  // coefficients need to multiply the grid size dx

  a0 *= dx;
  a11 *= dx;
  a21 *= dx;
  a2 *= dx;
  a31 *= dx;

  // rhoM *= u1M; // non sense
}

//  Step 1, Determine the number of particles to be sampled

int sample_from_P3M_getsize(double a0, double a11, double a2, double a21,
                            double a31, double Neff, RandomContext& random) {
  double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
                    (abs(a2) + abs(a21)) * 4 * exp(-1.) +
                    abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
  maxratio = maxratio * pow(sqrt(2), 3);
  return myfloor(maxratio / Neff, random);
}

//  Step2, sample.

NeParticleGroup sample_from_P3M_sample(double a0, double a11, double a2,
                                       double a21, double a31, int Ntotal,
                                       RandomContext& random) {
  NeParticleGroup S_new;
  double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
                    (abs(a2) + abs(a21)) * 4 * exp(-1.) +
                    abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
  // double maxratio = abs(a0) + abs(b) * sqrt(2.)*exp(-0.5) +
  //                  abs(c) * 4*exp(-1.) + abs(d) * (6*sqrt(6.) +
  //                  4*sqrt(2.))*exp(-1.5);
  maxratio = maxratio * pow(sqrt(2), 3);

  std::array<double, 3> v{};
  double coe_M0 = 1.0 / pow(sqrt(2. * pi), 3);
  double coe_M1 = 1.0 / pow(sqrt(4. * pi), 3);

  Particle1d3d S_one;

  double sqrt2 = sqrt(2.0);

  for (int k = 0; k < Ntotal; k++) {
    double vsq = 0.;
    for (int kv = 0; kv < 3; kv++) {
      v[kv] = sqrt2 * myrandn(random);
      vsq += v[kv] * v[kv];
    }

    // double M0 = coe_M0 * (a + b * v[0] + c * vsq + d * v[0]*vsq) *
    // exp(-vsq/2.);
    double M0 =
        coe_M0 *
        (a0 + a11 * v[0] + a2 * vsq + a21 * v[0] * v[0] + a31 * v[0] * vsq) *
        exp(-vsq / 2.);
    double M1 = coe_M1 * exp(-vsq / 4.);

    if (myrand(random) < (abs(M0) / M1 / maxratio)) {
      if (M0 > 0) {
        S_one.set_velocity(v);
        S_new.push_back(S_one, ParticleKind::Positive);
      } else {
        S_one.set_velocity(v);
        S_new.push_back(S_one, ParticleKind::Negative);
      }
    }
  }

  return S_new;

  // cout << (kp+kn)/( (double) Ntotal) << endl;
}

//  Step3, enforce conservation
NeParticleGroup sample_from_P3M_rescale(const NeParticleGroup &S_new, double u1,
                                        double Tprt) {
  NeParticleGroup S_rescaled;
  const auto &Sp = S_new.list(ParticleKind::Positive);
  const auto &Sn = S_new.list(ParticleKind::Negative);

  std::vector<double> v_rescale(3);

  std::vector<double> u_center{u1, 0., 0.};
  double sqrtT = sqrt(Tprt);

  for (const auto& Sone : Sp) {
    const auto &v_normalized = Sone.velocity();
    for (int kv = 0; kv < 3; kv++)
      v_rescale[kv] = u_center[kv] + sqrtT * v_normalized[kv];
    S_rescaled.push_back(Particle1d3d(Sone.position(), v_rescale),
                         ParticleKind::Positive);
  }

  for (const auto& Sone : Sn) {
    const auto &v_normalized = Sone.velocity();
    for (int kv = 0; kv < 3; kv++)
      v_rescale[kv] = u_center[kv] + sqrtT * v_normalized[kv];
    S_rescaled.push_back(Particle1d3d(Sone.position(), v_rescale),
                         ParticleKind::Negative);
  }

  return S_rescaled;
}

/**
  The whole algorithm of sampling from the micro-macro projection of advection
  part Sample P and N particles from - (I-\Pi) T M + \Pi T g
*/

// in one grid
void sample_from_MMprojection_homo(NeParticleGroup &S_x,
                                   const NumericGridClass &grid,
                                   RandomContext& random) {
  double a0, a11, a2, a21, a31;
  int Ntotal;

  sample_from_P3M_coeff_ver3(S_x, grid.dt, grid.dx, a0, a11, a2, a21,
                             a31);
  // sample_from_P3M_coeff_nog(S_x, grid.dt, grid.Neff, a0, a11, a2, a21,
  // a31);
  Ntotal = sample_from_P3M_getsize(a0, a11, a2, a21, a31, grid.Neff,
                                   random);

  if (S_x.TprtM < 0) {
    cout << " (" << S_x.rhoM << ' ' << S_x.u1M << ' ' << S_x.TprtM << ") ";
    cout << a0 << ' ' << a11 << ' ' << a2 << ' ' << a21 << ' ' << a31 << ' '
         << Ntotal << endl;
  }

  auto S_x_new =
      sample_from_P3M_sample(a0, a11, a2, a21, a31, Ntotal, random);

  S_x_new = sample_from_P3M_rescale(S_x_new, S_x.u1M, S_x.TprtM);

  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax(), random);
  // cout << "( " << ptr_S_x_new.size('p') << ", " << ptr_S_x_new.size('n')
  // <<
  // ") ";

  merge_NeParticleGroup(S_x, S_x_new);

  if ((S_x.size(ParticleKind::Positive) +
       S_x.size(ParticleKind::Negative)) > 200) {
    // enforce_conservation_zero(S_x, grid.Neff);
  }

  /*
  cout << "Np, Nn = " << S_x.size(ParticleKind::Positive) << ' '
       << S_x.size(ParticleKind::Negative) << endl;
  cout << ", after cons 2d = " <<  S_x . positive_moments.m0 - S_x . negative_moments.m0 << endl;
  S_x . computemoments();
  cout << ", after cons 2e = " <<  S_x . positive_moments.m0 - S_x . negative_moments.m0 << endl;
  if (abs(S_x . positive_moments.m0 - S_x. negative_moments.m0)>.5) {
    cout << "conserv m0, out, " <<  S_x . positive_moments.m0 << ' ' << S_x . negative_moments.m0 << endl;
    exit(0);
  }
  */
}

// over all grids
void sample_from_MMprojection(std::vector<NeParticleGroup> &S_x,
                              const NumericGridClass &grid,
                              RandomContext& random) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) {
    sample_from_MMprojection_homo(S_x[kx], grid, random);
  }
}


}  // namespace coulomb
