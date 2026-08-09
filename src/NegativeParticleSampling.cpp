#include "NegativeParticleSampling.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include "Constants.h"
#include "Grid.h"
#include "Numerics.h"
#include "ParticleGroup.h"
#include "SimulationConfig.h"
#include "utils.h"

namespace coulomb {
using std::abs;
using std::cos;
using std::exp;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::sin;
using std::sqrt;
using std::vector;

/* ======================================================================== *\
        Sample from the source term: the change in Maxwellian due to collisions
        with P and N particles
\* ======================================================================== */

/**
  Evaluate M(v)
  evaluate M(v0) where the moments of Maxwellian is given in moment
*/

double evaluateM(const std::vector<double> &v0, const NeParticleGroup &S_x) {
  double rho = S_x.rhoM;
  double u1 = S_x.u1M;
  double u2 = S_x.u2M;
  double u3 = S_x.u3M;
  double Tprt = S_x.TprtM;

  double usq = (v0[0] - u1) * (v0[0] - u1) + (v0[1] - u2) * (v0[1] - u2) +
               (v0[2] - u3) * (v0[2] - u3);
  return rho / pow(sqrt(2.0 * pi * Tprt), 3) * exp(-usq / 2.0 / Tprt);
}

// ========================================================================

/**
  Evaluate h(v;v1)
  evaluate h(v0;v1) = h(v0), with a source particle at v1
  mode = 0, calculate
  mode = 1, calculate c0
  mode = 2, calculate c1
  mode = 3, calculate c2
*/

double evaluateH(const std::vector<double> &v0, const std::vector<double> &v1,
                 const NeParticleGroup &S_x, const ParaClass &para, int mode) {
  // double rho = S_x.rhoM;
  double u1 = S_x.u1M;
  double u2 = S_x.u2M;
  double u3 = S_x.u3M;
  double Tprt = S_x.TprtM;

  double h = 0;

  if (para.method_binarycoll == BinaryCollisionMethod::TA) {
    double u[3], sqrt_u = 0;
    for (int k = 0; k < 3; k++) {
      u[k] = v0[k] - v1[k];
      sqrt_u += u[k] * u[k];
    }

    sqrt_u = sqrt(sqrt_u);
    double sigma2_delta = para.coeff_binarycoll * para.dt / pow(sqrt_u, 3);

    double r = 0.5 / sigma2_delta;
    double logr = log(r);
    double maxdelta;
    if (r > 1.0) {
      maxdelta = exp(-0.0026 * logr * logr - 0.4150 * logr + 0.9154);
    } else {
      maxdelta = exp(-0.0006 * logr * logr - 0.2242 * logr + 0.4711);
    }

    double v0_perp[3], v0_perp_sq = 0, v0_perp_dot_u = 0;
    v0_perp[0] = v0[0] - u1;
    v0_perp[1] = v0[1] - u2;
    v0_perp[2] = v0[2] - u3;
    for (int k = 0; k < 3; k++) v0_perp_dot_u += v0_perp[k] * u[k];
    for (int k = 0; k < 3; k++)
      v0_perp[k] -= v0_perp_dot_u * u[k] / (sqrt_u * sqrt_u);
    for (int k = 0; k < 3; k++) v0_perp_sq += v0_perp[k] * v0_perp[k];

    double vperp2 = v0_perp_sq / (2.0 * Tprt);

    int Ndelta = 16;
    vector<double> delta_all(Ndelta);
    for (int j = 0; j < Ndelta; j++)
      delta_all[j] = ((double)j) / Ndelta * maxdelta;

    vector<double> coeff_sum(Ndelta);
    for (int j = 0; j < Ndelta; j++) {
      double zeta2p1_1p5 = pow(sqrt(1.0 + delta_all[j] * delta_all[j]), 3);
      coeff_sum[j] = sqrt(r / pi) * pow(sqrt(zeta2p1_1p5), 3) *
                     exp(-r * delta_all[j] * delta_all[j] * zeta2p1_1p5);
    }

    double dx_delta = delta_all[1] - delta_all[0];

    coeff_sum[0] = coeff_sum[0] / 2.0;
    coeff_sum[Ndelta - 1] = coeff_sum[Ndelta - 1] / 2.0;

    vector<double> eps2(Ndelta);
    for (int j = 0; j < Ndelta; j++)
      eps2[j] = sqrt_u * sqrt_u * delta_all[j] * delta_all[j] / 2.0 / Tprt;

    double M = evaluateM(v0, S_x);
    double hM = 0;

    for (int j = 0; j < Ndelta; j++)
      hM += exp(-eps2[j]) *
            (1.0 + (eps2[j] * vperp2) +
             .25 * (eps2[j] * vperp2 * eps2[j] * vperp2)) *
            coeff_sum[j];

    hM = hM * dx_delta * 2.0;

    h = M * hM;

    double hh = 0;
    if (mode == 1) {  // this h gives c0
      for (int j = 0; j < Ndelta; j++) hh += exp(-eps2[j]) * coeff_sum[j];
      h = sqrt_u * sqrt_u * (hh * dx_delta * 2.0 - 1.0);
    } else if (mode == 2) {  // this h gives c1
      for (int j = 0; j < Ndelta; j++)
        hh += exp(-eps2[j]) * eps2[j] * coeff_sum[j];
      h = sqrt_u * sqrt_u * hh * dx_delta * 2.0;
    } else if (mode == 3) {  // this h gives c2
      for (int j = 0; j < Ndelta; j++)
        hh += exp(-eps2[j]) * eps2[j] * eps2[j] * coeff_sum[j];
      h = sqrt_u * sqrt_u * hh * dx_delta * 2.0;
    }

  } else {
    h = 1.0;
  }

  return h;
}

/**
  Find the lower/upper bound of delta m (v;v1)
  for the delta source at v1, find lower/upper bound of delta m(v;v1) in the
  following form \delta m_n(v;v1) > - alpha_neg * M(v) |v-v1|^2 delta m_p(v;v1)
  < alpha_pos * rho_m \detla m = \delta m_p - \delta m_n, with \delta m_p(m) = 0
  if |\delta m_p(m)|<(>) alpha_neg * M(v)
  |\delta m_n(v;v1)| < alpha_neg * M(v)
  |v-v1|^2 delta m_p(v;v1) < alpha_pos * rho_m
*/

void finddeltambound(NeParticleGroup &S_x, const ParaClass &para) {
  double Tprt = S_x.TprtM;
  std::vector<double> v1{sqrt(Tprt), 0, 0};

  int Neps_in = 20;
  vector<double> eps_all(Neps_in + 1);
  for (int k = 0; k <= Neps_in; k++) eps_all[k] = 0.1 / (1 << k);

  double vrange = 3.0 * sqrt(Tprt);

  int Neps_out = 40;
  int length_v_all = 2 * Neps_in + 2 * Neps_out;

  vector<double> v_all_1(length_v_all);

  double dv1 = (v1[0] - eps_all[0] + vrange) / Neps_out;
  double dv2 = (vrange - v1[0] - eps_all[0]) / Neps_out;

  for (int k = 0; k < Neps_out; k++) v_all_1[k] = -vrange + (k + 1) * dv1;
  for (int k = Neps_out; k < (Neps_out + Neps_in); k++)
    v_all_1[k] = v1[0] - eps_all[k - Neps_out + 1];
  for (int k = Neps_out + Neps_in; k < (Neps_out + 2 * Neps_in); k++)
    v_all_1[k] = v1[0] + eps_all[k - Neps_out - Neps_in + 1];
  for (int k = Neps_out + 2 * Neps_in; k < length_v_all; k++)
    v_all_1[k] = vrange - (k - Neps_out - 2 * Neps_in + 1) * dv2;

  vector<double> hh(length_v_all);
  vector<double> MM(length_v_all);

  // Look for lower bound

  std::vector<double> v0{0, 0, 0};

  for (int kv = 0; kv < length_v_all; kv++) {
    v0[0] = v_all_1[kv];
    hh[kv] = evaluateH(v0, v1, S_x, para);
    MM[kv] = evaluateM(v0, S_x);
  }

  // save_macro<double>(hh, "hh");

  // delta m = h - m
  vector<double> hhMM(length_v_all);
  for (int kv = 0; kv < length_v_all; kv++) hhMM[kv] = hh[kv] / MM[kv] - 1;
  double alpha_neg = -minval(hhMM);

  // Look for upper bound

  vector<double> hh0(length_v_all);
  vector<double> hh1(length_v_all);
  vector<double> hh2(length_v_all);
  for (int kv = 0; kv < length_v_all; kv++) {
    v0[0] = v_all_1[kv];

    hh0[kv] = evaluateH(v0, v1, S_x, para, 1);
    hh1[kv] = evaluateH(v0, v1, S_x, para, 2);
    hh2[kv] = evaluateH(v0, v1, S_x, para, 3);
  }

  double beta = 3.0;

  double alpha_pos =
      maxval(hh0) + maxval(hh1) * beta * beta + maxval(hh2) * pow(beta, 4);

  S_x.alpha_neg = alpha_neg;
  S_x.alpha_pos = alpha_pos;
  S_x.rmax = 6.0 * sqrt(2 * Tprt);

  // cout << "alpha = ( " << alpha_neg << ", " << alpha_pos << ", " << S_x .
  // rmax << " )" << endl;
}

void finddeltambound_inhomo(std::vector<NeParticleGroup> &S_x,
                            const NumericGridClass &grid,
                            const ParaClass &para) {
  double minTprt = S_x.front().TprtM;
  int kx_minTprt = 0;
  for (int kx = 1; kx < grid.Nx; kx++) {
    if (minTprt > S_x[kx].TprtM) {
      minTprt = S_x[kx].TprtM;
      kx_minTprt = kx;
    }
  }

  finddeltambound(S_x[kx_minTprt], para);

  double alpha_neg = S_x[kx_minTprt].alpha_neg;
  double alpha_pos = S_x[kx_minTprt].alpha_pos;

  for (int kx = 0; kx < grid.Nx; kx++) {
    S_x[kx].alpha_neg = alpha_neg;
    S_x[kx].alpha_pos = alpha_pos;
    S_x[kx].rmax = 6.0 * sqrt(2 * (S_x[kx].TprtM));
  }
}

// ========================================================================

struct NegativeParticleSample {
  std::vector<double> velocity;
  ParticleKind kind;
};

/** Sample one accepted particle from the negative part of Delta M. */
std::optional<NegativeParticleSample> samplefromh_neg(
    NeParticleGroup &S_x, const ParaClass &para, double Neff,
    RandomContext& random) {
  std::vector<double> v0{0.0, 0.0, 0.0};

  double alpha_neg = S_x.alpha_neg;

  double rhof = S_x.rho;
  double rhop = S_x.positive_moments.m0 * Neff;
  double rhon = S_x.negative_moments.m0 * Neff;

  int Np, Nn;
  Np = S_x.size(ParticleKind::Positive);
  Nn = S_x.size(ParticleKind::Negative);

  int Npickup = para.Npickup_neg;

  int NNp = min(Npickup, Np);
  int NNn = min(Npickup, Nn);

  const auto idp = myrandperm(Np, NNp, random);
  const auto idn = myrandperm(Nn, NNn, random);

  v0[0] = myrandn(random) * sqrt(S_x.TprtM) + S_x.u1M;
  v0[1] = myrandn(random) * sqrt(S_x.TprtM) + S_x.u2M;
  v0[2] = myrandn(random) * sqrt(S_x.TprtM) + S_x.u3M;

  double M0 = evaluateM(v0, S_x);

  double hp = 0, hn = 0;

  auto &Sp = S_x.list(ParticleKind::Positive);
  auto &Sn = S_x.list(ParticleKind::Negative);

  for (int kp = 0; kp < NNp; kp++) {
    auto &v1 = Sp[idp[kp] - 1].velocity();
    double h0 = evaluateH(v0, v1, S_x, para) - M0;
    if (h0 < (alpha_neg * M0)) hp += h0;
  }
  for (int kn = 0; kn < NNn; kn++) {
    auto &v1 = Sn[idn[kn] - 1].velocity();
    double h0 = evaluateH(v0, v1, S_x, para) - M0;
    if (h0 < (alpha_neg * M0)) hn += h0;
  }
  double h = hp * Np / (NNp + 1.0e-15) - hn * Nn / (NNn + 1.0e-15);
  h = h * Neff / rhof;
  double hbar = max(rhop, rhon) / rhof * M0 * alpha_neg;
  double r0 = myrand(random);
  if (r0 < (abs(h) / hbar)) {
    return NegativeParticleSample{
        v0, h > 0 ? ParticleKind::Positive : ParticleKind::Negative};
  }
  return std::nullopt;
}

// ========================================================================

/**
  the number of virtual particles sampled from h_+
  determine the number of virtual particles sampled from \Delta m_+
  Max_p(j) is the upper bound for h due to source at Sp(idp(j))
*/

int samplefromDeltamp_Npv(const NeParticleGroup &S_x, double Neff,
                           RandomContext& random) {
  int Np = S_x.size(ParticleKind::Positive);
  int Nn = S_x.size(ParticleKind::Negative);

  double rhom = S_x.rhoM;
  double Tprtm = S_x.TprtM;

  double rho = rhom + Neff * (Np - Nn);
  return myfloor(4.0 * pi * S_x.rmax * S_x.alpha_pos * rhom /
                 pow(sqrt(2.0 * pi * Tprtm), 3) / rho * (Np + Nn), random);
}

// ========================================================================

/**
  Sample particles from \Delta M, in homogeneous case
*/

void samplefromDeltam(NeParticleGroup &S_x, NeParticleGroup &S_x_new,
                      const ParaClass &para, double Neff,
                      RandomContext& random) {
  // Sample particles from Delta m
  double alpha_neg = S_x.alpha_neg;
  double alpha_pos = S_x.alpha_pos;
  double rmax = S_x.rmax;

  int Np = S_x.size(ParticleKind::Positive);
  int Nn = S_x.size(ParticleKind::Negative);

  auto &Sp = S_x.list(ParticleKind::Positive);
  auto &Sn = S_x.list(ParticleKind::Negative);

  double rhof = S_x.rho;
  double rhom = S_x.rhoM;
  double Tprtm = S_x.TprtM;
  double maxm = rhom / pow(sqrt(2.0 * pi * Tprtm), 3);

  // Sample from negative part

  double Nneg_f = max(Np, Nn) * alpha_neg * rhom / rhof;
  int Nneg = myfloor(Nneg_f, random);  // Number of virtual particles

  // Particle1d3d * Sp_new = S_x_new . list('p');
  // Particle1d3d * Sn_new = S_x_new . list('n');

  Particle1d3d S_one;

  for (int kneg = 0; kneg < Nneg; kneg++) {
    const auto sample = samplefromh_neg(S_x, para, Neff, random);
    if (sample) {
      S_one.set_velocity(sample->velocity);
      S_x_new.push_back(S_one, sample->kind);
    }
  }

  // cout << "negpart " << COUNT_MYRAND << endl;

  int Npos = samplefromDeltamp_Npv(S_x, Neff, random);

  // cout << "Npos = " << Npos << endl;

  double rate_P = ((double)Np) / (Np + Nn);

  // int kk_test = 418;
  // cout << Npos << ' ' << rate_P << ' ' << COUNT_MYRAND << endl;
  for (int kpos = 0; kpos < Npos; kpos++) {
    double rrr = myrand(random);

    /*
    if (FLAG_CHECK == 1) {
      if (kpos == 52) {
        cout << COUNT_MYRAND << endl;
        cout << rrr << ' ' << rate_P << endl;
        // std::exit(0);
      }
    }
                */

    if (rrr < rate_P) {
      // Sample positve particles
      // choose the source particle
      int kp = (int)(Np * myrand(random));
      auto &v1 = Sp[kp].velocity();
      // cout << v1[0] << ' '<< v1[1] << ' '<< v1[2] << endl;
      // sample a particle from the nearby
      double r1 = myrand(random) * rmax;
      double costheta = 2.0 * myrand(random) - 1.0;
      double sintheta = sqrt(1.0 - costheta * costheta);
      double phi = myrand(random) * pi * 2.0;
      std::vector<double> v0{v1[0] + r1 * sintheta * cos(phi),
                             v1[1] + r1 * sintheta * sin(phi),
                             v1[2] + r1 * costheta};

      double M0 = evaluateM(v0, S_x);

      // if (kpos == kk_test) cout << "test " << COUNT_MYRAND <<' '<< v0[0] << '
      // '<< v0[1] << ' '<< v0[2] << ' ' << M0 << endl;

      if (myrand(random) < (M0 / maxm)) {
        double H0 = evaluateH(v0, v1, S_x, para);
        double Hbar0 = H0 - M0 - alpha_neg * M0;
        if (Hbar0 > 0) {
          // check v0 is in the pos zone
          double r2h0 = r1 * r1 * Hbar0;
          double rr = myrand(random);
          if (rr < (r2h0 / (alpha_pos * M0))) {
            // accept the virtual particle v0 with suitable rate
            S_one.set_velocity(v0);
        S_x_new.push_back(S_one, ParticleKind::Positive);
            // cout << "pos " <<  COUNT_MYRAND << ' ' << kpos << endl;
            // cout << v0[0] << ' '<< v0[1] << ' '<< v0[2] << endl;
          }
        }
      }
    } else {
      // Sample negative particles
      // choose the source particle
      int kn = (int)(Nn * myrand(random));
      auto &v1 = Sn[kn].velocity();
      // sample a particle from the nearby
      double r1 = myrand(random) * rmax;
      double costheta = 2.0 * myrand(random) - 1.0;
      double sintheta = sqrt(1.0 - costheta * costheta);
      double phi = myrand(random) * pi * 2.0;
      std::vector<double> v0{v1[0] + r1 * sintheta * cos(phi),
                             v1[1] + r1 * sintheta * sin(phi),
                             v1[2] + r1 * costheta};

      double M0 = evaluateM(v0, S_x);

      rrr = myrand(random);
      if (rrr < (M0 / maxm)) {
        double H0 = evaluateH(v0, v1, S_x, para);
        double Hbar0 = H0 - M0 - alpha_neg * M0;
        if (Hbar0 > 0) {
          // check v0 is in the pos zone
          double r2h0 = r1 * r1 * Hbar0;
          double rr = myrand(random);
          if (rr < (r2h0 / (alpha_pos * M0))) {
            // accept the virtual particle v0 with suitable rate
            S_one.set_velocity(v0);
        S_x_new.push_back(S_one, ParticleKind::Negative);
          }
        }
      }
    }
  }
}


}  // namespace coulomb
