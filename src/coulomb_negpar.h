#include "Classes.h"

namespace coulomb {
void particleresample_inhomo(NeParticleGroup *S_x, NumericGridClass &grid,
                             ParaClass &para);
void sync_coarse(std::vector<NeParticleGroup> &S_x, NumericGridClass &grid,
                 ParaClass &para);
// ========================================================================

/**
  Perform coulomb collisions in homogeneous case between P/N and F particles
*/

void coulomb_collision_homo_PFNF(NeParticleGroup *S_x, const ParaClass &para) {
  int Nf = S_x->size('f');
  int Np = S_x->size('p');
  int Nn = S_x->size('n');
  if (Nf < (Np + Nn)) {
    cout << "Too few F particles." << endl;
    cout << "(" << Np << ", " << Nn << ", " << Nf << ") " << endl;
    // particleresample_homo(S_x, para);
  }

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');
  auto &Sf = S_x->list('f');

  const auto p = myrandperm(Nf, Np + Nn);
  int kf;

  for (int kp = 0; kp < Np; kp++) {
    kf = p[kp] - 1;
    const auto v1 = Sp[kp].velocity();
    const auto v2 = Sf[kf].velocity();

    const auto vp = coulombBinary3d(v1, v2, para);

    Sp[kp].set_velocity(vp.first);
  }
  for (int kn = 0; kn < Nn; kn++) {
    kf = p[kn + Np] - 1;
    const auto v1 = Sn[kn].velocity();
    const auto v2 = Sf[kf].velocity();

    const auto vp = coulombBinary3d(v1, v2, para);

    Sn[kn].set_velocity(vp.first);
  }
}

/* ======================================================================== *\
        Sample from the source term: the change in Maxwellian due to collisions
        with P and N particles
\* ======================================================================== */

/**
  Evaluate M(v)
  evaluate M(v0) where the moments of Maxwellian is given in moment
*/

double evaluateM(double *v0, NeParticleGroup *S_x) {
  double rho = S_x->rhoM;
  double u1 = S_x->u1M;
  double u2 = S_x->u2M;
  double u3 = S_x->u3M;
  double Tprt = S_x->TprtM;

  double usq = (*v0 - u1) * (*v0 - u1) + (*(v0 + 1) - u2) * (*(v0 + 1) - u2) +
               (*(v0 + 2) - u3) * (*(v0 + 2) - u3);
  return rho / pow(sqrt(2.0 * pi * Tprt), 3) * exp(-usq / 2.0 / Tprt);
}

// ========================================================================

/**
  Evaluate h(v;v1)
  evaluate h(v0;v1) = h(v0), with a source particle at v1
*/

double evaluateH(double *v0, double *v1, NeParticleGroup *S_x,
                 const ParaClass &para) {
  // double rho = S_x->rhoM;
  double u1 = S_x->u1M;
  double u2 = S_x->u2M;
  double u3 = S_x->u3M;
  double Tprt = S_x->TprtM;

  double h = 0;

  if (para.method_binarycoll.compare("TA") == 0) {
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
    if (FLAG_PRECOMPUTE_ALPHA_U == 100) {  // this h gives c0
      for (int j = 0; j < Ndelta; j++) hh += exp(-eps2[j]) * coeff_sum[j];
      h = sqrt_u * sqrt_u * (hh * dx_delta * 2.0 - 1.0);
    } else if (FLAG_PRECOMPUTE_ALPHA_U == 101) {  // this h gives c1
      for (int j = 0; j < Ndelta; j++)
        hh += exp(-eps2[j]) * eps2[j] * coeff_sum[j];
      h = sqrt_u * sqrt_u * hh * dx_delta * 2.0;
    } else if (FLAG_PRECOMPUTE_ALPHA_U == 102) {  // this h gives c2
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

void finddeltambound(NeParticleGroup *S_x, const ParaClass &para) {
  double Tprt = S_x->TprtM;
  double v1[3] = {sqrt(Tprt), 0, 0};

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

  double v0[3] = {0, 0, 0};

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

    FLAG_PRECOMPUTE_ALPHA_U = 100;
    hh0[kv] = evaluateH(v0, v1, S_x, para);
    FLAG_PRECOMPUTE_ALPHA_U = 101;
    hh1[kv] = evaluateH(v0, v1, S_x, para);
    FLAG_PRECOMPUTE_ALPHA_U = 102;
    hh2[kv] = evaluateH(v0, v1, S_x, para);
  }

  double beta = 3.0;

  double alpha_pos =
      maxval(hh0) + maxval(hh1) * beta * beta + maxval(hh2) * pow(beta, 4);

  S_x->alpha_neg = alpha_neg;
  S_x->alpha_pos = alpha_pos;
  S_x->rmax = 6.0 * sqrt(2 * Tprt);

  FLAG_PRECOMPUTE_ALPHA_U = 0;

  // cout << "alpha = ( " << alpha_neg << ", " << alpha_pos << ", " << S_x ->
  // rmax << " )" << endl;
}

void finddeltambound_inhomo(NeParticleGroup *S_x, const NumericGridClass &grid,
                            const ParaClass &para) {
  double minTprt = S_x->TprtM;
  int kx_minTprt = 0;
  for (int kx = 1; kx < grid.Nx; kx++) {
    if (minTprt > (S_x + kx)->TprtM) {
      minTprt = (S_x + kx)->TprtM;
      kx_minTprt = kx;
    }
  }

  finddeltambound(S_x + kx_minTprt, para);

  double alpha_neg = (S_x + kx_minTprt)->alpha_neg;
  double alpha_pos = (S_x + kx_minTprt)->alpha_pos;

  for (int kx = 0; kx < grid.Nx; kx++) {
    (S_x + kx)->alpha_neg = alpha_neg;
    (S_x + kx)->alpha_pos = alpha_pos;
    (S_x + kx)->rmax = 6.0 * sqrt(2 * ((S_x + kx)->TprtM));
  }
}

// ========================================================================

/**
  Sample one particle from negative part of Delta M
*/

void samplefromh_neg(double *v0, int &signv, bool &flag_accept,
                     NeParticleGroup *S_x, const ParaClass &para, double Neff) {
  double alpha_neg = S_x->alpha_neg;

  double rhof = S_x->rho;
  double rhop = S_x->m0P * Neff;
  double rhon = S_x->m0N * Neff;

  int Np, Nn;
  Np = S_x->size('p');
  Nn = S_x->size('n');

  flag_accept = false;
  signv = 1;

  for (int k = 0; k < 3; k++) v0[k] = 0.;

  int Npickup = para.Npickup_neg;

  int NNp = min(Npickup, Np);
  int NNn = min(Npickup, Nn);

  const auto idp = myrandperm(Np, NNp);
  const auto idn = myrandperm(Nn, NNn);

  v0[0] = myrandn() * sqrt(S_x->TprtM) + S_x->u1M;
  v0[1] = myrandn() * sqrt(S_x->TprtM) + S_x->u2M;
  v0[2] = myrandn() * sqrt(S_x->TprtM) + S_x->u3M;

  double M0 = evaluateM(v0, S_x);

  double hp = 0, hn = 0;

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  for (int kp = 0; kp < NNp; kp++) {
    auto &v1 = Sp[idp[kp] - 1].velocity();
    double h0 = evaluateH(v0, &v1[0], S_x, para) - M0;
    if (h0 < (alpha_neg * M0)) hp += h0;
  }
  for (int kn = 0; kn < NNn; kn++) {
    auto &v1 = Sn[idn[kn] - 1].velocity();
    double h0 = evaluateH(v0, &v1[0], S_x, para) - M0;
    if (h0 < (alpha_neg * M0)) hn += h0;
  }
  double h = hp * Np / (NNp + 1.0e-15) - hn * Nn / (NNn + 1.0e-15);
  h = h * Neff / rhof;
  double hbar = max(rhop, rhon) / rhof * M0 * alpha_neg;
  double r0 = myrand();
  if (r0 < (abs(h) / hbar)) {
    flag_accept = true;
    if (h > 0) {
      signv = 1;
    } else {
      signv = -1;
    }
  }
}

// ========================================================================

/**
  the number of virtual particles sampled from h_+
  determine the number of virtual particles sampled from \Delta m_+
  Max_p(j) is the upper bound for h due to source at Sp(idp(j))
*/

int samplefromDeltamp_Npv(NeParticleGroup *S_x, double Neff) {
  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  double rhom = S_x->rhoM;
  double Tprtm = S_x->TprtM;

  double rho = rhom + Neff * (Np - Nn);
  return myfloor(4.0 * pi * S_x->rmax * S_x->alpha_pos * rhom /
                 pow(sqrt(2.0 * pi * Tprtm), 3) / rho * (Np + Nn));
}

// ========================================================================

/**
  Sample particles from \Delta M, in homogeneous case
*/

void samplefromDeltam(NeParticleGroup *S_x, NeParticleGroup *S_x_new,
                      const ParaClass &para, double Neff) {
  // Sample particles from Delta m
  double alpha_neg = S_x->alpha_neg;
  double alpha_pos = S_x->alpha_pos;
  double rmax = S_x->rmax;

  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  double rhof = S_x->rho;
  double rhom = S_x->rhoM;
  double Tprtm = S_x->TprtM;
  double maxm = rhom / pow(sqrt(2.0 * pi * Tprtm), 3);

  // Sample from negative part

  double Nneg_f = max(Np, Nn) * alpha_neg * rhom / rhof;
  int Nneg = myfloor(Nneg_f);  // Number of virtual particles

  // Particle1d3d * Sp_new = S_x_new -> list('p');
  // Particle1d3d * Sn_new = S_x_new -> list('n');

  Particle1d3d S_one;

  int signv;
  bool flag_accept;
  double v0[3];

  for (int kneg = 0; kneg < Nneg; kneg++) {
    samplefromh_neg(v0, signv, flag_accept, S_x, para, Neff);
    if (flag_accept) {
      if (signv > 0) {
        S_one.set_velocity(v0);
        S_x_new->push_back(S_one, 'p');
      } else {
        S_one.set_velocity(v0);
        S_x_new->push_back(S_one, 'n');
      }
    }
  }

  // cout << "negpart " << COUNT_MYRAND << endl;

  int Npos = samplefromDeltamp_Npv(S_x, Neff);

  // cout << "Npos = " << Npos << endl;

  double rate_P = ((double)Np) / (Np + Nn);

  // int kk_test = 418;
  // cout << Npos << ' ' << rate_P << ' ' << COUNT_MYRAND << endl;
  for (int kpos = 0; kpos < Npos; kpos++) {
    double rrr = myrand();

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
      int kp = (int)(Np * myrand());
      auto &v1 = Sp[kp].velocity();
      // cout << v1[0] << ' '<< v1[1] << ' '<< v1[2] << endl;
      // sample a particle from the nearby
      double r1 = myrand() * rmax;
      double costheta = 2.0 * myrand() - 1.0;
      double sintheta = sqrt(1.0 - costheta * costheta);
      double phi = myrand() * pi * 2.0;
      double v0[3] = {v1[0] + r1 * sintheta * cos(phi),
                      v1[1] + r1 * sintheta * sin(phi), v1[2] + r1 * costheta};

      double M0 = evaluateM(v0, S_x);

      // if (kpos == kk_test) cout << "test " << COUNT_MYRAND <<' '<< v0[0] << '
      // '<< v0[1] << ' '<< v0[2] << ' ' << M0 << endl;

      if (myrand() < (M0 / maxm)) {
        double H0 = evaluateH(v0, &v1[0], S_x, para);
        double Hbar0 = H0 - M0 - alpha_neg * M0;
        if (Hbar0 > 0) {
          // check v0 is in the pos zone
          double r2h0 = r1 * r1 * Hbar0;
          double rr = myrand();
          if (rr < (r2h0 / (alpha_pos * M0))) {
            // accept the virtual particle v0 with suitable rate
            S_one.set_velocity(v0);
            S_x_new->push_back(S_one, 'p');
            // cout << "pos " <<  COUNT_MYRAND << ' ' << kpos << endl;
            // cout << v0[0] << ' '<< v0[1] << ' '<< v0[2] << endl;
          }
        }
      }
    } else {
      // Sample negative particles
      // choose the source particle
      int kn = (int)(Nn * myrand());
      auto &v1 = Sn[kn].velocity();
      // sample a particle from the nearby
      double r1 = myrand() * rmax;
      double costheta = 2.0 * myrand() - 1.0;
      double sintheta = sqrt(1.0 - costheta * costheta);
      double phi = myrand() * pi * 2.0;
      double v0[3] = {v1[0] + r1 * sintheta * cos(phi),
                      v1[1] + r1 * sintheta * sin(phi), v1[2] + r1 * costheta};

      double M0 = evaluateM(v0, S_x);

      rrr = myrand();
      if (rrr < (M0 / maxm)) {
        double H0 = evaluateH(v0, &v1[0], S_x, para);
        double Hbar0 = H0 - M0 - alpha_neg * M0;
        if (Hbar0 > 0) {
          // check v0 is in the pos zone
          double r2h0 = r1 * r1 * Hbar0;
          double rr = myrand();
          if (rr < (r2h0 / (alpha_pos * M0))) {
            // accept the virtual particle v0 with suitable rate
            S_one.set_velocity(v0);
            S_x_new->push_back(S_one, 'n');
          }
        }
      }
    }
  }
}

// ========================================================================

/**
  Merge the new P and N particles to the exising ones
*/

void merge_NeParticleGroup(NeParticleGroup *S_x, NeParticleGroup *S_x_new) {
  int Np_new = S_x_new->size('p');
  int Nn_new = S_x_new->size('n');
  // int Nf_new = S_x_new->size('f');
  auto &Sp = S_x_new->list('p');
  auto &Sn = S_x_new->list('n');
  // Particle1d3d * Sf = S_x_new -> list('f');
  for (int kp = 0; kp < Np_new; kp++) S_x->push_back(&Sp[kp], 'p');
  for (int kn = 0; kn < Nn_new; kn++) S_x->push_back(&Sn[kn], 'n');
  // for (int kf=0; kf<Nf_new; kf++)	S_x->push_back(Sf+kf, 'f');
  // S_x -> computemoments();
}

void mergeF_NeParticleGroup(NeParticleGroup *S_x, NeParticleGroup *S_x_new) {
  int Nf_new = S_x_new->size('f');
  auto &Sf = S_x_new->list('f');
  for (int kf = 0; kf < Nf_new; kf++) S_x->push_back(&Sf[kf], 'f');
  // S_x -> computemoments();
}

// ========================================================================

/**
  Assign positions to the new particles
*/
void assign_positions(NeParticleGroup *S_new, double xmin, double xmax) {
  double x1 = xmin, x2 = xmax;
  auto &Sp = S_new->list('p');
  auto &Sn = S_new->list('n');
  auto &Sf = S_new->list('f');
  for (int kp = 0; kp < S_new->size('p'); kp++)
    Sp[kp].set_position(myrand() * (x2 - x1) + x1);
  for (int kp = 0; kp < S_new->size('n'); kp++)
    Sn[kp].set_position(myrand() * (x2 - x1) + x1);
  for (int kp = 0; kp < S_new->size('f'); kp++)
    Sf[kp].set_position(myrand() * (x2 - x1) + x1);
}

// ========================================================================

/**
  evolve one step, in homo case
*/

void NegPar_collision_homo(NeParticleGroup *S_x, const ParaClass &para,
                           double Neff) {
  // NeParticleGroup S_x_new(S_x->size('p'), S_x->size('n'), 0);
  NeParticleGroup S_x_new;
  NeParticleGroup *ptr_S_x_new = &S_x_new;

  // sample from the change of maxwellian due to M-P and M-N collisions
  // finddeltambound(S_x, para);

  // cout << "before sample " << COUNT_MYRAND << endl;
  samplefromDeltam(S_x, ptr_S_x_new, para, Neff);
  // cout << "after sample " << COUNT_MYRAND << endl;

  assign_positions(ptr_S_x_new, S_x->get_xmin(), S_x->get_xmax());

  // perform P-F and N-F collisions
  coulomb_collision_homo_PFNF(S_x, para);

  // merge the new sampled particles to the post-collisional particles
  merge_NeParticleGroup(S_x, ptr_S_x_new);

  // perform F-F collisions
  auto &Sf = S_x->list('f');
  coulomb_collision_homo(&Sf[0], S_x->size('f'), para);
}

/**
  evolve one step in inhomo case
*/

void NegPar_collision(NeParticleGroup *S_x, const NumericGridClass &grid,
                      const ParaClass &para) {
  finddeltambound_inhomo(S_x, grid, para);

  for (int kx = 0; kx < grid.Nx; kx++) {
    // cout << "kx = " << kx << " N = " << (S_x+kx)->size('p')<< ' ' <<
    // (S_x+kx)->size('n') << endl;
    NegPar_collision_homo(S_x + kx, para, grid.Neff);
  }
}

void NegPar_collision_openmp(NeParticleGroup *S_x, const NumericGridClass &grid,
                             const ParaClass &para) {
  finddeltambound_inhomo(S_x, grid, para);
#pragma omp parallel if (para.FLAG_USE_OPENMP)
  {
#pragma omp for
    for (int kx = 0; kx < grid.Nx; kx++) {
      // cout << " " << kx << " ";
      NegPar_collision_homo(S_x + kx, para, grid.Neff);
      /*
          NeParticleGroup S_x_one = *(S_x+kx);
            NegPar_collision_homo(&S_x_one, para, grid.Neff);
          #pragma omp critical
            *(S_x+kx) = S_x_one;
        */
    }
  }
}

/* ======================================================== *\
        Enforce conservation
\* ======================================================== */

/**
  Enforce conservation
  Input: the desired moments: m0, m11, m12, m13, m21, m22, m23
         Neg particle group: S_new
         effecitve number: Neff
  Output: update particles in S_new
*/

void enforce_conservation(double m0, double m11, double m12, double m13,
                          double m21, double m22, double m23,
                          NeParticleGroup *S_new, double Neff,
                          bool flag_conserve_energyvector) {
  // enforce m0
  double m0_need = m0;
  S_new->computemoments();
  double m0_actual = Neff * (S_new->m0P - S_new->m0N);
  // cout << "before cons = " <<  S_new -> m0P - S_new -> m0N;
  int N_remove;
  if (m0_actual < m0_need) {
    N_remove = myfloor((m0_need - m0_actual) / Neff);
    for (int kp = 0; kp < N_remove; kp++) {
      int k_remove = (int)(myrand() * S_new->size('n'));
      S_new->erase(k_remove, 'n');
    }
  } else {
    N_remove = myfloor((m0_actual - m0_need) / Neff);
    for (int kp = 0; kp < N_remove; kp++) {
      int k_remove = (int)(myrand() * S_new->size('p'));
      S_new->erase(k_remove, 'p');
    }
  }

  int Np = S_new->size('p');
  int Nn = S_new->size('n');

  // enforce m11, m12, m13

  S_new->computemoments();
  double m1_actual[3], m1_need[3] = {m11, m12, m13};
  m1_actual[0] = Neff * (S_new->m11P - S_new->m11N);
  m1_actual[1] = Neff * (S_new->m12P - S_new->m12N);
  m1_actual[2] = Neff * (S_new->m13P - S_new->m13N);

  double v0[3], m1_mod[3];
  for (int kv = 0; kv < 3; kv++) m1_mod[kv] = -m1_actual[kv] + m1_need[kv];

  if (Np > Nn) {
    auto &Sp = S_new->list('p');
    for (int kp = 0; kp < Np; kp++) {
      auto &vkp = Sp[kp].velocity();
      for (int kv = 0; kv < 3; kv++) v0[kv] = vkp[kv] + m1_mod[kv] / Neff / Np;
      Sp[kp].set_velocity(v0);
    }
  } else {
    auto &Sn = S_new->list('n');
    for (int kn = 0; kn < Nn; kn++) {
      auto &vkn = Sn[kn].velocity();
      for (int kv = 0; kv < 3; kv++) v0[kv] = vkn[kv] - m1_mod[kv] / Neff / Nn;
      Sn[kn].set_velocity(v0);
    }
  }

  /*
  {
  Particle1d3d *Snn = S_new -> list('n');
  for (int kn = 0; kn < Nn; kn ++) cout << Snn -> velocity(0) << ' ';
  cout << endl;
  }
  */

  // enforce m21, m22, m23
  // mu2p*Tp - mu2n*Tn + Np*cp^2 - Nn*cn^2 = Ep_need - En_need

  S_new->computemoments();

  double cp[3], cn[3], Tp[3], Tn[3], RHS[3];
  double m2p_actual[3], m2n_actual[3];
  double m2_need[3] = {m21 / Neff, m22 / Neff, m23 / Neff};

  cp[0] = S_new->m11P;
  cp[1] = S_new->m12P;
  cp[2] = S_new->m13P;
  cn[0] = S_new->m11N;
  cn[1] = S_new->m12N;
  cn[2] = S_new->m13N;
  for (int kv = 0; kv < 3; kv++) {
    cp[kv] /= Np;
    cn[kv] /= Nn;
  }

  m2p_actual[0] = S_new->m21P;
  m2p_actual[1] = S_new->m22P;
  m2p_actual[2] = S_new->m23P;
  m2n_actual[0] = S_new->m21N;
  m2n_actual[1] = S_new->m22N;
  m2n_actual[2] = S_new->m23N;

  /*
  cout << m2p_actual[0] - m2n_actual[0] << " vs " << m2_need[0] << ", "
      << m2p_actual[1] - m2n_actual[1] << " vs " << m2_need[1] << ", "
      << m2p_actual[2] - m2n_actual[2] << " vs " << m2_need[2] << endl;
  */

  double mu2p[3], mu2n[3];
  for (int kv = 0; kv < 3; kv++) {
    Tp[kv] = m2p_actual[kv] - Np * cp[kv] * cp[kv];
    Tn[kv] = m2n_actual[kv] - Nn * cn[kv] * cn[kv];
    RHS[kv] = m2_need[kv] - Np * cp[kv] * cp[kv] + Nn * cn[kv] * cn[kv];
  }

  if (flag_conserve_energyvector) {
    for (int kv = 0; kv < 3; kv++) {
      if (Np > Nn) {
        mu2n[kv] = 1.0;
        mu2p[kv] = 1.0 / Tp[kv] * (mu2n[kv] * Tn[kv] + RHS[kv]);
        if (mu2p[kv] < 0) {
          mu2p[kv] = 1.0;
          cout << "ERROR NOT CONSERVATIVE" << endl;
        }
      } else {
        mu2p[kv] = 1.0;
        mu2n[kv] = 1.0 / Tn[kv] * (mu2p[kv] * Tp[kv] - RHS[kv]);
        if (mu2n[kv] < 0) {
          mu2n[kv] = 1.0;
          cout << "ERROR NOT CONSERVATIVE" << endl;
        }
      }
    }
  } else {
    double sum_RHS = 0., sum_Tp = 0., sum_Tn = 0.;
    double mu2n_all, mu2p_all;
    for (int kv = 0; kv < 3; kv++) {
      sum_RHS += RHS[kv];
      sum_Tp += Tp[kv];
      sum_Tn += Tn[kv];
    }
    if (Np > Nn) {
      mu2n_all = 1.0;
      mu2p_all = 1.0 / sum_Tp * (mu2n_all * sum_Tn + sum_RHS);
      if (mu2p_all < 0) mu2p_all = 1.0;
      // if (mu2p_all<0) { mu2p_all = 1.0; cout << "ERROR NOT CONSERVATIVE TOO"
      // << endl;}
    } else {
      mu2p_all = 1.0;
      mu2n_all = 1.0 / sum_Tn * (mu2p_all * sum_Tp - sum_RHS);
      if (mu2n_all < 0) mu2n_all = 1.0;
      // if (mu2n_all<0) { mu2n_all = 1.0; cout << "ERROR NOT CONSERVATIVE TOO"
      // << endl;}
    }
    for (int kv = 0; kv < 3; kv++) {
      mu2n[kv] = mu2n_all;
      mu2p[kv] = mu2p_all;
    }
  }

  if (Np > Nn) {
    auto &Sp = S_new->list('p');
    for (int kp = 0; kp < Np; kp++) {
      auto &vkp = Sp[kp].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = sqrt(mu2p[kv]) * (vkp[kv] - cp[kv]) + cp[kv];
      Sp[kp].set_velocity(v0);
    }
  } else {
    auto &Sn = S_new->list('n');
    for (int kn = 0; kn < Nn; kn++) {
      auto &vkn = Sn[kn].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = sqrt(mu2n[kv]) * (vkn[kv] - cn[kv]) + cn[kv];
      Sn[kn].set_velocity(v0);
    }
  }
  S_new->computemoments();
  /*
  {
  Particle1d3d *Snn = S_new -> list('n');
  for (int kn = 0; kn < Nn; kn ++) cout << Snn -> velocity(0) << ' ';
  cout << endl;
  }
  */
}

void enforce_conservation_zero(NeParticleGroup *S_new, double Neff) {
  enforce_conservation(0., 0., 0., 0., 0., 0., 0., S_new, Neff, false);
}

/* ======================================================== *\
        Sample from the micro-macro projection in advection step
\* ======================================================== */

/**
  Sample P and N particles from (a0 + a11 v1 + a21 v1^2 + a2 |v|^2 + a31 v
  |v|^2) 1/(2*pi)^(3/2) exp(-|v|^2/2) For this function b = (b, 0, 0) Denote M0
  = 1/(2*pi)^(3/2) * exp(-|v|^2/2), M1 = 1/(2*pi*2)^(3/2) * exp(-|v|^2/4)
  Precompute  max_v M0/M1 = 2^(3/2)
              max_v v M0/M1 = 2^(3/2) sqrt(2) exp(-1/2)
              max_v v_1^2 M0/M1 = 2^(3/2) 4 exp(-1)
              max_v |v|^2 M0/M1 = 2^(3/2) 4 exp(-1)
              max_v v|v|^2 M0/M1 = 2^(3/2) (6 sqrt(6) + 4 sqrt(2)) exp(-3/2)
*/

//  Step 0, Determine the coefficients. Return a0, a11, a2, a21, a31

void sample_from_P3M_coeff(NeParticleGroup *S_x, double dt, double &a0,
                           double &a11, double &a2, double &a21, double &a31) {
  double rhoM = S_x->rhoM;
  double u1M = S_x->u1M;
  double TprtM = S_x->TprtM;

  double dx_rhoM = S_x->dx_rhoM;
  double dx_u1M = S_x->dx_u1M;
  double dx_TprtM = S_x->dx_TprtM;

  double elecfield = S_x->elecfield;

  const double dimen = 3.;
  double sqrtT = sqrt(TprtM);

  // coefficients from (v\cdot\nabla_x  + E \cdot\nabla_v) M

  a0 = -dt * rhoM * u1M * (dx_rhoM / rhoM - dimen * dx_TprtM / TprtM / 2.);
  a11 = -dt * rhoM *
        (u1M * dx_u1M / sqrtT - elecfield / sqrtT +
         sqrtT * (dx_rhoM / rhoM - dimen * dx_TprtM / TprtM / 2.));
  a21 = -dt * rhoM * dx_u1M;
  a2 = -dt * rhoM * u1M * dx_TprtM / TprtM / 2.;
  a31 = -dt * rhoM * dx_TprtM / sqrtT / 2.;

  // coefficients from \Pi_M (v\cdot\nabla_x  + E \cdot\nabla_v) f

  // inner product with 1
  double coe_0 = S_x->drho;
  // inner product with (v_1 - u_1)
  double coe_1 = S_x->dm1 - u1M * S_x->drho;
  // inner product with ( |v - u|^2/T - d )
  double coe_2 = 2. / TprtM * S_x->denergy - 2. * u1M / TprtM * S_x->dm1 +
                 (u1M * u1M / TprtM - dimen) * S_x->drho;

  // cout << a0 << ' ' << coe_0 << endl;
  a0 += coe_0 - .5 * coe_2;
  a11 += coe_1 / sqrtT;
  a2 += coe_2 / 2. / dimen;
}
/*
void sample_from_P3M_coeff_ver2(NeParticleGroup * S_x, double dt, double Neff,
double &a0, double &a11, double &a2, double &a21, double &a31) { double rhoM =
S_x -> rhoM; double u1M = S_x -> u1M; double TprtM = S_x -> TprtM;

  // double dx_rhoM = S_x -> dx_rhoM;
  double dx_u1M = S_x -> dx_u1M;
  double dx_TprtM = S_x -> dx_TprtM;

  const double dimen = 3.;

  // coefficients from -\Delta t (I-\Pi_M) (v\cdot\nabla_x  + E \cdot\nabla_v) M

  double sqrtT = sqrt(TprtM);

  a0  = 0.;
  a11 = -dt * rhoM * ( - (dimen + 2) / 2. / sqrtT * dx_TprtM );
  a21 = -dt * rhoM * dx_u1M;
  a2  = -dt * rhoM * ( - dx_u1M / dimen );
  a31 = -dt * rhoM * ( dx_TprtM / sqrtT / 2.);



  // coefficients from \Delta t \Pi_M (v\cdot\nabla_x  + E \cdot\nabla_v) (f_p -
f_n)
  // The particles have been advected
  S_x -> computemoments();

  // inner product with 1
  double drho = Neff/grid.dx * (S_x->m0P - S_x->m0N);
  double coe_0 = drho;
  // inner product with (v_1 - u_1)
  double dm1 = Neff/grid.dx * ( S_x -> m11P - S_x -> m11N ); // inner product
with v_1 double coe_1 = dm1 - u1M * drho;
  // inner product with ( |v - u|^2/T - d )
  double denergy = 0.5 * Neff * ( S_x -> m2P - S_x -> m2N ); // inner product
with |v|^2 double coe_2 = 2./TprtM * denergy - 2.*u1M/TprtM * dm1 +
(u1M*u1M/TprtM - dimen) * drho;

  cout <<"a11: " << a11 << ' ' << coe_1 / sqrt(TprtM) << endl;
  cout <<"a2: " << a2 << ' ' << coe_2 / 2. / dimen << endl;

  a0  -= coe_0 - .5*coe_2;
  a11 -= coe_1 / sqrt(TprtM);
  a2  -= coe_2 / 2. / dimen;



  // rhoM *= u1M; // non sense
}

*/

void sample_from_P3M_coeff_ver3(NeParticleGroup *S_x, double dt, double Neff,
                                double dx, double &a0, double &a11, double &a2,
                                double &a21, double &a31) {
  double rhoM = S_x->rhoM;
  double u1M = S_x->u1M;
  double TprtM = S_x->TprtM;

  // double dx_rhoM = S_x -> dx_rhoM;
  double dx_u1M = S_x->dx_u1M;
  double dx_TprtM = S_x->dx_TprtM;

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
  S_x->computemoments();

  // inner product with 1
  double drho = S_x->drho_g;
  double coe_0 = drho;
  // inner product with (v_1 - u_1)
  double dm1 = S_x->dm1_g;  // inner product with v_1
  double coe_1 = dm1 - u1M * drho;
  // inner product with ( |v - u|^2/T - d )
  double denergy = S_x->denergy_g;  // inner product with |v|^2
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
                            double a31, double Neff) {
  double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
                    (abs(a2) + abs(a21)) * 4 * exp(-1.) +
                    abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
  maxratio = maxratio * pow(sqrt(2), 3);
  return myfloor(maxratio / Neff);
}

//  Step2, sample.

void sample_from_P3M_sample(double a0, double a11, double a2, double a21,
                            double a31, int Ntotal, NeParticleGroup *S_new) {
  double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
                    (abs(a2) + abs(a21)) * 4 * exp(-1.) +
                    abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
  // double maxratio = abs(a0) + abs(b) * sqrt(2.)*exp(-0.5) +
  //                  abs(c) * 4*exp(-1.) + abs(d) * (6*sqrt(6.) +
  //                  4*sqrt(2.))*exp(-1.5);
  maxratio = maxratio * pow(sqrt(2), 3);

  double v[3];
  double coe_M0 = 1.0 / pow(sqrt(2. * pi), 3);
  double coe_M1 = 1.0 / pow(sqrt(4. * pi), 3);

  Particle1d3d S_one;

  double sqrt2 = sqrt(2.0);

  for (int k = 0; k < Ntotal; k++) {
    double vsq = 0.;
    for (int kv = 0; kv < 3; kv++) {
      v[kv] = sqrt2 * myrandn();
      vsq += v[kv] * v[kv];
    }

    // double M0 = coe_M0 * (a + b * v[0] + c * vsq + d * v[0]*vsq) *
    // exp(-vsq/2.);
    double M0 =
        coe_M0 *
        (a0 + a11 * v[0] + a2 * vsq + a21 * v[0] * v[0] + a31 * v[0] * vsq) *
        exp(-vsq / 2.);
    double M1 = coe_M1 * exp(-vsq / 4.);

    if (myrand() < (abs(M0) / M1 / maxratio)) {
      if (M0 > 0) {
        S_one.set_velocity(v);
        S_new->push_back(S_one, 'p');
      } else {
        S_one.set_velocity(v);
        S_new->push_back(S_one, 'n');
      }
    }
  }

  // cout << (kp+kn)/( (double) Ntotal) << endl;
}

//  Step3, enforce conservation
void sample_from_P3M_conserve(double a0, double a11, double a2, double a21,
                              double a31, NeParticleGroup *S_new, double Neff) {
  double m0_need = a0 + 3. * a2 + a21;
  double m11_need = a11 + 5. * a31;
  double m1k_need = 0.;  // k = 2, 3
  double m21_need = a0 + 5. * a2 + 3. * a21;
  double m2k_need = a0 + 5. * a2 + a21;  // k = 2, 3
  /*
  S_new -> computemoments();
  cout << S_new->size('p') << ' ' << S_new->size('n') << endl;
  cout  << Neff * ( S_new->m0P - S_new->m0N ) << ' '
        << Neff * ( S_new->m11P - S_new->m11N ) << ' '
        << Neff * ( S_new->m21P - S_new->m21N ) << ' ' << endl;
  cout << m0_need << ' ' << m11_need << ' ' << m21_need << endl;
  */
  enforce_conservation(m0_need, m11_need, m1k_need, m1k_need, m21_need,
                       m2k_need, m2k_need, S_new, Neff, true);
  /*
  S_new -> computemoments();
  cout  << Neff * ( S_new->m0P - S_new->m0N ) << ' '
        << Neff * ( S_new->m11P - S_new->m11N ) << ' '
        << Neff * ( S_new->m21P - S_new->m21N ) << ' ' << endl << endl;
  */
}
/*
void sample_from_P3M_conserve_aftermerge(NeParticleGroup * S_x, double Neff) {
  enforce_conservation_zero(S_x, Neff);
}
*/

//  Step4, rescale.

void sample_from_P3M_rescale(NeParticleGroup *S_new, double u1, double Tprt) {
  int Np = S_new->size('p');
  int Nn = S_new->size('n');
  auto &Sp = S_new->list('p');
  auto &Sn = S_new->list('n');

  std::vector<double> v_rescale(3);

  std::vector<double> u_center{u1, 0., 0.};
  double sqrtT = sqrt(Tprt);

  for (int kp = 0; kp < Np; kp++) {
    auto &v_normalized = Sp[kp].velocity();
    for (int kv = 0; kv < 3; kv++)
      v_rescale[kv] = u_center[kv] + sqrtT * v_normalized[kv];
    Sp[kp].set_velocity(v_rescale);
  }

  for (int kn = 0; kn < Nn; kn++) {
    auto &v_normalized = Sn[kn].velocity();
    for (int kv = 0; kv < 3; kv++)
      v_rescale[kv] = u_center[kv] + sqrtT * v_normalized[kv];
    Sn[kn].set_velocity(v_rescale);
  }
}

/**
  The whole algorithm of sampling from the micro-macro projection of advection
  part Sample P and N particles from - (I-\Pi) T M + \Pi T g
*/

// in one grid
void sample_from_MMprojection_homo(NeParticleGroup *S_x,
                                   const NumericGridClass &grid) {
  double a0, a11, a2, a21, a31;
  int Ntotal;

  sample_from_P3M_coeff_ver3(S_x, grid.dt, grid.Neff, grid.dx, a0, a11, a2, a21,
                             a31);
  // sample_from_P3M_coeff_nog(S_x, grid.dt, grid.Neff, a0, a11, a2, a21, a31);
  Ntotal = sample_from_P3M_getsize(a0, a11, a2, a21, a31, grid.Neff);

  if (S_x->TprtM < 0) {
    cout << " (" << S_x->rhoM << ' ' << S_x->u1M << ' ' << S_x->TprtM << ") ";
    cout << a0 << ' ' << a11 << ' ' << a2 << ' ' << a21 << ' ' << a31 << ' '
         << Ntotal << endl;
  }

  NeParticleGroup S_x_new;
  NeParticleGroup *ptr_S_x_new = &S_x_new;

  sample_from_P3M_sample(a0, a11, a2, a21, a31, Ntotal, ptr_S_x_new);

  // sample_from_P3M_conserve(a0, a11, a2, a21, a31, ptr_S_x_new, grid.Neff);

  sample_from_P3M_rescale(ptr_S_x_new, S_x->u1M, S_x->TprtM);

  assign_positions(ptr_S_x_new, S_x->get_xmin(), S_x->get_xmax());
  // cout << "( " << ptr_S_x_new->size('p') << ", " << ptr_S_x_new->size('n') <<
  // ") ";

  merge_NeParticleGroup(S_x, ptr_S_x_new);

  if ((S_x->size('p') + S_x->size('n')) > 200) {
    // enforce_conservation_zero(S_x, grid.Neff);
  }

  /*
  cout << "Np, Nn = " << S_x->size('p') << ' ' << S_x->size('n') << endl;
  cout << ", after cons 2d = " <<  S_x -> m0P - S_x -> m0N << endl;
  S_x -> computemoments();
  cout << ", after cons 2e = " <<  S_x -> m0P - S_x -> m0N << endl;
  if (abs(S_x -> m0P - S_x-> m0N)>.5) {
    cout << "conserv m0, out, " <<  S_x -> m0P << ' ' << S_x -> m0N << endl;
    exit(0);
  }
  */
}

// over all grids
void sample_from_MMprojection(NeParticleGroup *S_x,
                              const NumericGridClass &grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) {
    sample_from_MMprojection_homo(S_x + kx, grid);
  }
}

// =================================================================================

/** Update all macro quantities in S_x including
  (1) dx_rhoM, dx_u1M, dx_TprtM; // derivative in x direction, used in sampling
  from source part (2) rho, u1, Tprt; // the moments of all distributions (3)
  m0P, m11P, m12P, m13P, m2P, m21P, m22P, m23P, m31P, m32P, m33P; // the moments
  of P particles, with Neff = 1 m0N, m11N, m12N, m13N, m2N, m21N, m22N, m23N,
  m31N, m32N, m33N; // the moments of N particles, with Neff = 1 m0F, m11F,
  m12F, m13F, m2F, m21F, m22F, m23F, m31F, m32F, m33F; // the moments of F
  particles, with Neff = 1
*/

void update_macro(std::vector<NeParticleGroup> &S_x,
                  const NumericGridClass &grid) {
  int Nx = grid.Nx;

  // update group (3)
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // update group (2)
  update_rhouT(&S_x[0], grid);

  // update group (1)
  update_dx_rhouT_M(&S_x[0], grid);
}

// =================================================================================

void NegPar_BGK_collision_homo(NeParticleGroup *S_x, ParaClass &para) {
  int Np = S_x->size('p');
  int Nn = S_x->size('n');
  int Nf = S_x->size('f');

  int Np_remove = myfloor(Np * (para.dt * para.coeff_binarycoll));
  int Nn_remove = myfloor(Nn * (para.dt * para.coeff_binarycoll));

  for (int kp = 0; kp < Np_remove; kp++) {
    int k_remove = (int)(myrand() * S_x->size('p'));
    S_x->erase(k_remove, 'p');
  }

  for (int kp = 0; kp < Nn_remove; kp++) {
    int k_remove = (int)(myrand() * S_x->size('n'));
    S_x->erase(k_remove, 'n');
  }

  double rate_change = para.dt * para.coeff_binarycoll;

  double vf[3];

  double sqrtT = sqrt(S_x->TprtM);
  auto &Sf = S_x->list('f');
  for (int kf = 0; kf < Nf; kf++) {
    if (myrand() < rate_change) {
      vf[0] = S_x->u1M + sqrtT * myrandn();
      vf[1] = S_x->u2M + sqrtT * myrandn();
      vf[2] = S_x->u3M + sqrtT * myrandn();
      Sf[kf].set_velocity(vf);
    }
  }
}

void NegPar_BGK_collision(NeParticleGroup *S_x, NumericGridClass &grid,
                          ParaClass &para) {
  for (int kx = 0; kx < grid.Nx; kx++) {
    // cout << '(' << (S_x+kx) -> size('p') << ", " << (S_x+kx) -> size('n') <<
    // ") -> (";
    NegPar_BGK_collision_homo(S_x + kx, para);
    // cout << (S_x+kx) -> size('p') << ", " << (S_x+kx) -> size('n') << ")\n";
  }
}

/* ======================================================== *\
        Forward one step
\* ======================================================== */

int count_particle_number(const std::vector<NeParticleGroup> &S_x, int Nx,
                          char partype);

/** One step
  Forward one step in time, with time splitting
*/

// void Negpar_inhomo_onestep(NeParticleGroup * S_x, NumericGridClass & grid,
// ParaClass & para, MultlLevelGroup * MLsol) {
void Negpar_inhomo_onestep(std::vector<NeParticleGroup> &S_x,
                           NumericGridClass &grid, ParaClass &para) {
  int Nplast = count_particle_number(S_x, grid.Nx, 'p');
  int Nnlast = count_particle_number(S_x, grid.Nx, 'n');

  cout << "step start" << endl;

  // Step 1, collision.

  // Step 1.0 update all macro quantities
  update_macro(S_x, grid);

  // Step 1.0 perform negative collisions

  t0_coll = clock();

  if (para.flag_collision == 1)
    // NegPar_collision(S_x, grid, para);
    NegPar_collision_openmp(&S_x[0], grid, para);
  else if (para.flag_collision == 2)
    NegPar_BGK_collision(&S_x[0], grid, para);

  cout << "step 1" << endl;

  int Npcoll = count_particle_number(S_x, grid.Nx, 'p');
  int Nncoll = count_particle_number(S_x, grid.Nx, 'n');

  t1_coll = clock();

  // step 2, advection

  t0_adve = t1_coll;

  // Step 2.0 update all macro quantities and electric field
  update_macro(S_x, grid);
  updateelecfiled(S_x, grid);

  for (int kx = 0; kx < grid.Nx; kx++) S_x[kx].copymoments();

  cout << "step 2.0" << endl;

  // Switch 2.1 and 2.2

  // Step 2.1, compute moment change: S_x->drho, dm1, denergy
  compute_change_in_macro(&S_x[0], grid);
  cout << "step 2.1" << endl;

  // Step 2.2, advect P N F particles.
  particleadvection(S_x, grid);
  cout << "step 2.2" << endl;

  // Step 2.3, Sample P and N particles from micro-macro projection
  sample_from_MMprojection(&S_x[0], grid);
  cout << "step 2.3" << endl;

  // Step 2.4, update maxwellian part:S_x->rhoM, u1M, TprtM
  update_maxwellian(&S_x[0], grid);
  cout << "step 2.4" << endl;

  int Npadve = count_particle_number(S_x, grid.Nx, 'p');
  int Nnadve = count_particle_number(S_x, grid.Nx, 'n');

  cout << "d(Np, Nn) = (" << Npcoll - Nplast << ", " << Nncoll - Nnlast
       << "), (" << Npadve - Npcoll << ", " << Nnadve - Nncoll << ")" << endl;

  t1_adve = clock();

  // Step 3, resampling particles when needed
  // particleresample_inhomo(S_x, grid, para, MLsol);

  t0_resamp = t1_adve;
  if (para.flag_collision == 1) {
    particleresample_inhomo(&S_x[0], grid, para);
  }
  t1_resamp = clock();

  sync_coarse(S_x, grid, para);

  cout << "Np = " << count_particle_number(S_x, grid.Nx, 'p')
       << "; Nn = " << count_particle_number(S_x, grid.Nx, 'n')
       << "; Nf = " << count_particle_number(S_x, grid.Nx, 'f') << endl;
}

void Negpar_inhomo_onestep_ver2(std::vector<NeParticleGroup> &S_x,
                                NumericGridClass &grid, ParaClass &para) {
  int Nplast = count_particle_number(S_x, grid.Nx, 'p');
  int Nnlast = count_particle_number(S_x, grid.Nx, 'n');

  cout << "step start" << endl;

  // Step 1, collision.

  // Step 1.0 update all macro quantities
  update_macro(S_x, grid);

  // Step 1.0 perform negative collisions

  if (para.flag_collision == 1)
    NegPar_collision(&S_x[0], grid, para);
  else if (para.flag_collision == 2)
    NegPar_BGK_collision(&S_x[0], grid, para);

  cout << "step 1" << endl;

  int Npcoll = count_particle_number(S_x, grid.Nx, 'p');
  int Nncoll = count_particle_number(S_x, grid.Nx, 'n');

  // step 2, advection

  // Step 2.0 update all macro quantities and electric field
  update_macro(S_x, grid);
  updateelecfiled(S_x, grid);
  cout << "step 2.0" << endl;

  // Switch 2.1 and 2.2

  // Step 2.1, advect P N F particles.
  particleadvection(S_x, grid);

  // Step 2.1, compute moment change: S_x->drho, dm1, denergy
  compute_change_in_macro(&S_x[0], grid);

  // Step 2.3, Sample P and N particles from micro-macro projection
  sample_from_MMprojection(&S_x[0], grid);

  // Step 2.4, update maxwellian part:S_x->rhoM, u1M, TprtM
  update_maxwellian(&S_x[0], grid);
  cout << "step 2.4" << endl;

  int Npadve = count_particle_number(S_x, grid.Nx, 'p');
  int Nnadve = count_particle_number(S_x, grid.Nx, 'n');

  cout << "d(Np, Nn) = (" << Npcoll - Nplast << ", " << Nncoll - Nnlast
       << "), (" << Npadve - Npcoll << ", " << Nnadve - Nncoll << ")" << endl;

  // Step 3, resampling particles when needed
  // particleresample_inhomo(S_x, grid, para, MLsol);
  particleresample_inhomo(&S_x[0], grid, para);

  cout << "Np = " << count_particle_number(S_x, grid.Nx, 'p')
       << "; Nn = " << count_particle_number(S_x, grid.Nx, 'n')
       << "; Nf = " << count_particle_number(S_x, grid.Nx, 'f') << endl;
}

void Negpar_inhomo_onestep_PIC(std::vector<NeParticleGroup> &S_x,
                               NumericGridClass &grid, ParaClass &para) {
  t0_coll = clock();

  for (int kx = 0; kx < grid.Nx; kx++) {
    auto &Sf = S_x[kx].list('f');
    coulomb_collision_homo(&Sf[0], S_x[kx].size('f'), para);
  }

  t1_coll = clock();

  updateelecfiled_PIC(S_x, grid);

  t0_adve = clock();
  particleadvection(S_x, 'f', grid);
  t1_adve = clock();

  cout << "Np = " << count_particle_number(S_x, grid.Nx, 'p')
       << "; Nn = " << count_particle_number(S_x, grid.Nx, 'n')
       << "; Nf = " << count_particle_number(S_x, grid.Nx, 'f') << endl;
}

void Negpar_inhomo_onestep_stop(std::vector<NeParticleGroup> &S_x,
                                NumericGridClass &grid, ParaClass &para,
                                int flag_stop) {
  int Nplast = count_particle_number(S_x, grid.Nx, 'p');
  int Nnlast = count_particle_number(S_x, grid.Nx, 'n');

  cout << "step start" << endl;

  // Step 1, collision.

  // Step 1.0 update all macro quantities
  update_macro(S_x, grid);

  // Step 1.0 perform negative collisions

  t0_coll = clock();

  if (para.flag_collision == 1)
    NegPar_collision(&S_x[0], grid, para);
  else if (para.flag_collision == 2)
    NegPar_BGK_collision(&S_x[0], grid, para);

  cout << "step 1" << endl;

  int Npcoll = count_particle_number(S_x, grid.Nx, 'p');
  int Nncoll = count_particle_number(S_x, grid.Nx, 'n');

  t1_coll = clock();

  // step 2, advection

  t0_adve = t1_coll;

  if (flag_stop == 0) {
    // Step 2.0 update all macro quantities and electric field
    update_macro(S_x, grid);
    updateelecfiled(S_x, grid);

    for (int kx = 0; kx < grid.Nx; kx++) S_x[kx].copymoments();

    cout << "step 2.0" << endl;

    // Switch 2.1 and 2.2

    // Step 2.1, compute moment change: S_x->drho, dm1, denergy
    compute_change_in_macro(&S_x[0], grid);
    cout << "step 2.1" << endl;

    // Step 2.2, advect P N F particles.
    particleadvection(S_x, grid);
    cout << "step 2.2" << endl;

    // Step 2.3, Sample P and N particles from micro-macro projection
    sample_from_MMprojection(&S_x[0], grid);
    cout << "step 2.3" << endl;

    // Step 2.4, update maxwellian part:S_x->rhoM, u1M, TprtM
    update_maxwellian(&S_x[0], grid);
    cout << "step 2.4" << endl;
  }

  int Npadve = count_particle_number(S_x, grid.Nx, 'p');
  int Nnadve = count_particle_number(S_x, grid.Nx, 'n');

  cout << "d(Np, Nn) = (" << Npcoll - Nplast << ", " << Nncoll - Nnlast
       << "), (" << Npadve - Npcoll << ", " << Nnadve - Nncoll << ")" << endl;

  t1_adve = clock();

  // Step 3, resampling particles when needed
  // particleresample_inhomo(S_x, grid, para, MLsol);

  t0_resamp = t1_adve;
  particleresample_inhomo(&S_x[0], grid, para);
  t1_resamp = clock();

  cout << "Np = " << count_particle_number(S_x, grid.Nx, 'p')
       << "; Nn = " << count_particle_number(S_x, grid.Nx, 'n')
       << "; Nf = " << count_particle_number(S_x, grid.Nx, 'f') << endl;
}
}  // namespace coulomb