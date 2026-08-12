#include "ParticleConservation.h"

#include <array>
#include <cmath>
#include <iostream>

#include "ParticleGroup.h"
#include "RandomSampling.h"

namespace coulomb {

void ParticleConservation::enforce(double m0, double m11, double m12,
                                   double m13, double m21, double m22,
                                   double m23, NeParticleGroup &S_new,
                                   double Neff, bool flag_conserve_energyvector,
                                   RandomContext &random) {
  // enforce m0
  double m0_need = m0;
  S_new.computemoments();
  double m0_actual =
      Neff * (S_new.positive_moments.m0 - S_new.negative_moments.m0);
  // cout << "before cons = " <<  S_new . positive_moments.m0 - S_new .
  // negative_moments.m0;
  int N_remove;
  if (m0_actual < m0_need) {
    N_remove =
        RandomSampling(random).stochastic_floor((m0_need - m0_actual) / Neff);
    for (int kp = 0; kp < N_remove; kp++) {
      int k_remove = (int)(RandomSampling(random).uniform() *
                           S_new.size(ParticleKind::Negative));
      S_new.erase(k_remove, ParticleKind::Negative);
    }
  } else {
    N_remove =
        RandomSampling(random).stochastic_floor((m0_actual - m0_need) / Neff);
    for (int kp = 0; kp < N_remove; kp++) {
      int k_remove = (int)(RandomSampling(random).uniform() *
                           S_new.size(ParticleKind::Positive));
      S_new.erase(k_remove, ParticleKind::Positive);
    }
  }

  int Np = S_new.size(ParticleKind::Positive);
  int Nn = S_new.size(ParticleKind::Negative);

  // enforce m11, m12, m13

  S_new.computemoments();
  double m1_actual[3], m1_need[3] = {m11, m12, m13};
  m1_actual[0] =
      Neff * (S_new.positive_moments.m11 - S_new.negative_moments.m11);
  m1_actual[1] =
      Neff * (S_new.positive_moments.m12 - S_new.negative_moments.m12);
  m1_actual[2] =
      Neff * (S_new.positive_moments.m13 - S_new.negative_moments.m13);

  std::array<double, 3> v0{};
  double m1_mod[3];
  for (int kv = 0; kv < 3; kv++)
    m1_mod[kv] = -m1_actual[kv] + m1_need[kv];

  if (Np > Nn) {
    auto &Sp = S_new.list(ParticleKind::Positive);
    for (int kp = 0; kp < Np; kp++) {
      auto &vkp = Sp[kp].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = vkp[kv] + m1_mod[kv] / Neff / Np;
      Sp[kp].set_velocity(v0);
    }
  } else {
    auto &Sn = S_new.list(ParticleKind::Negative);
    for (int kn = 0; kn < Nn; kn++) {
      auto &vkn = Sn[kn].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = vkn[kv] - m1_mod[kv] / Neff / Nn;
      Sn[kn].set_velocity(v0);
    }
  }

  // enforce m21, m22, m23
  // mu2p*Tp - mu2n*Tn + Np*cp^2 - Nn*cn^2 = Ep_need - En_need

  S_new.computemoments();

  double cp[3], cn[3], Tp[3], Tn[3], RHS[3];
  double m2p_actual[3], m2n_actual[3];
  double m2_need[3] = {m21 / Neff, m22 / Neff, m23 / Neff};

  cp[0] = S_new.positive_moments.m11;
  cp[1] = S_new.positive_moments.m12;
  cp[2] = S_new.positive_moments.m13;
  cn[0] = S_new.negative_moments.m11;
  cn[1] = S_new.negative_moments.m12;
  cn[2] = S_new.negative_moments.m13;
  for (int kv = 0; kv < 3; kv++) {
    cp[kv] /= Np;
    cn[kv] /= Nn;
  }

  m2p_actual[0] = S_new.positive_moments.m21;
  m2p_actual[1] = S_new.positive_moments.m22;
  m2p_actual[2] = S_new.positive_moments.m23;
  m2n_actual[0] = S_new.negative_moments.m21;
  m2n_actual[1] = S_new.negative_moments.m22;
  m2n_actual[2] = S_new.negative_moments.m23;

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
          std::cout << "ERROR NOT CONSERVATIVE" << std::endl;
        }
      } else {
        mu2p[kv] = 1.0;
        mu2n[kv] = 1.0 / Tn[kv] * (mu2p[kv] * Tp[kv] - RHS[kv]);
        if (mu2n[kv] < 0) {
          mu2n[kv] = 1.0;
          std::cout << "ERROR NOT CONSERVATIVE" << std::endl;
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
      if (mu2p_all < 0)
        mu2p_all = 1.0;
    } else {
      mu2p_all = 1.0;
      mu2n_all = 1.0 / sum_Tn * (mu2p_all * sum_Tp - sum_RHS);
      if (mu2n_all < 0)
        mu2n_all = 1.0;
    }
    for (int kv = 0; kv < 3; kv++) {
      mu2n[kv] = mu2n_all;
      mu2p[kv] = mu2p_all;
    }
  }

  if (Np > Nn) {
    auto &Sp = S_new.list(ParticleKind::Positive);
    for (int kp = 0; kp < Np; kp++) {
      auto &vkp = Sp[kp].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = std::sqrt(mu2p[kv]) * (vkp[kv] - cp[kv]) + cp[kv];
      Sp[kp].set_velocity(v0);
    }
  } else {
    auto &Sn = S_new.list(ParticleKind::Negative);
    for (int kn = 0; kn < Nn; kn++) {
      auto &vkn = Sn[kn].velocity();
      for (int kv = 0; kv < 3; kv++)
        v0[kv] = std::sqrt(mu2n[kv]) * (vkn[kv] - cn[kv]) + cn[kv];
      Sn[kn].set_velocity(v0);
    }
  }
  S_new.computemoments();
}

void ParticleConservation::enforce_zero(NeParticleGroup &S_new, double Neff,
                                        RandomContext &random) {
  ParticleConservation::enforce(0., 0., 0., 0., 0., 0., 0., S_new, Neff, false,
                                random);
}

} // namespace coulomb
