#include "InitialConditions.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "Grid.h"
#include "SimulationConfig.h"

namespace coulomb {
using std::abs;
using std::cos;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::pow;
using std::sin;
using std::sqrt;
using std::vector;

IniValClass InitialConditions::create(NumericGridClass &grid) {
  IniValClass initial_data;

  initial_data.probname = "LandauDamping";
  // initial_data.probname = "TwoStreamInstab";
  // initial_data.probname = "test1";
  // initial_data.probname = "BumpOnTail";
  // initial_data.probname = "Analytic";
  // initial_data.probname = "Efficiency";

  if (initial_data.probname == "Efficiency")
    grid.lambda_Poisson = 0;
  if (initial_data.probname == "TwoStreamInstab")
    InitialConditions::configure_two_stream(initial_data);

  cout << "Problem name = " << initial_data.probname << endl;
  initial_data.totalmass = 0;
  return initial_data;
}

void InitialConditions::configure(IniValClass &initial_data,
                                  const NumericGridClass &grid, int cell) {
  const double x = grid.x[cell];
  const double dx = grid.x[1] - grid.x[0];
  const double spatial_volume = grid.xmax - grid.xmin;
  const double spatial_center = (grid.x[0] + grid.x[grid.Nx - 1]) / 2.0;

  if ((initial_data.probname == "LandauDamping") ||
      (initial_data.probname == "Efficiency")) {
    initial_data.rho = 1 + initial_data.LD_alpha * sin(x);
    initial_data.velocity[0] = 0.;
    initial_data.velocity[1] = 0.;
    initial_data.velocity[2] = 0.;
    initial_data.Tprt = 1.0;
    initial_data.totalmass += initial_data.rho * dx;
  } else if (initial_data.probname == "TwoStreamInstab") {
    constexpr double alpha = 0.5;
    constexpr double wave_number = 0.5;
    initial_data.TSI_coe =
        1. / 12 / pi * (1.0 + alpha * cos(wave_number * x)) / 40.;
  } else if (initial_data.probname == "BumpOnTail") {
    const double rho0 =
        initial_data.BOT_rho0 / spatial_volume * initial_data.BOT_beta;
    const double rho1 = initial_data.BOT_rho0 * (1 - initial_data.BOT_beta) /
                        spatial_volume *
                        exp(-(x - spatial_center) * (x - spatial_center) /
                            (2 * initial_data.BOT_Tx));
    const double rho_total = rho0 + rho1;
    const double background_velocity = 0.;
    const double bump_velocity = initial_data.BOT_ub;
    const double momentum = rho1 * bump_velocity;
    const double velocity = momentum / rho_total;
    const double background_temperature = initial_data.BOT_Tprt;
    const double bump_temperature = initial_data.BOT_dTprt;
    const double energy =
        1.5 * rho0 * background_temperature +
        rho1 * (.5 * bump_velocity * bump_velocity + 1.5 * bump_temperature);
    const double temperature =
        (energy - .5 * rho_total * velocity * velocity) / 1.5 / rho_total;

    initial_data.rho = rho_total;
    initial_data.velocity[0] = velocity;
    initial_data.velocity[1] = 0.;
    initial_data.velocity[2] = 0.;
    initial_data.Tprt = temperature;

    initial_data.rho1 = rho0;
    initial_data.velocity1[0] = background_velocity;
    initial_data.velocity1[1] = 0.;
    initial_data.velocity1[2] = 0.;
    initial_data.Tprt1 = background_temperature;

    initial_data.rho2 = rho1;
    initial_data.velocity2[0] = bump_velocity;
    initial_data.velocity2[1] = 0.;
    initial_data.velocity2[2] = 0.;
    initial_data.Tprt2 = bump_temperature;
    initial_data.totalmass += initial_data.rho * dx;
  }
}

void InitialConditions::configure_two_stream(IniValClass &initial_data) {
  double vmax = 6.;
  int Nv = 200;
  double dv = 2.0 * vmax / Nv;
  vector<double> v(Nv);
  for (int kv = 0; kv < Nv; kv++)
    v[kv] = (kv + 0.5) * dv - vmax;

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
        max_f_over_M = max(max_f_over_M, abs(f0[kk] - M0[kk]) / M0[kk]);
        m21 += v1 * v1 * f0[kk];
        m22 += v2 * v2 * f0[kk];
        m23 += v3 * v3 * f0[kk];
        if (f0[kk] > M0[kk])
          rhop += f0[kk] - M0[kk];
      }
    }
  }

  rhop *= dv * dv * dv;
  m21 *= dv * dv * dv;
  m22 *= dv * dv * dv;
  m23 *= dv * dv * dv;

  initial_data.TSI_rhof = rhof;
  initial_data.TSI_rhop = rhop;
  initial_data.TSI_Tprt = Tprt;
  initial_data.TSI_max_f_over_M = max_f_over_M;
  initial_data.TSI_m21 = m21;
  initial_data.TSI_m22 = m22;
  initial_data.TSI_m23 = m23;

  cout << "Initially " << m21 + m22 + m23 << " vs " << 3 * rhof * Tprt << endl;
}

} // namespace coulomb
