#include "Classes.h"

#include "_global_variables.h"

namespace coulomb {

ParaClass::ParaClass() {
  method = "HDP";
  // method = "PIC";
  FLAG_USE_OPENMP = true;
  dt = 0.01;
  coeff_binarycoll = 10.0;
  method_binarycoll = "TA";
  resample_ratio = 1.2;
  Npickup_neg = 100;
  Nfreq = 30;
  // Nlevel = 4;
  //  Num_grids = (3 * Nlevel * Nlevel - 3 * Nlevel + 2) / 2;
  //  Num_gridpoints = (Nlevel * (Nlevel + 1) * (1 << (Nlevel + 2)) +
  //                    Nlevel * (Nlevel - 1) * (1 << (Nlevel + 1)) +
  //                    (Nlevel - 1) * (Nlevel - 2) * (1 << Nlevel)) /
  // 2;
  collisionType = COULOMB_COLLISION;
  lambda_Poisson = 10.0;
  resample_spatial_ratio = 0.9;
  sync_time_interval = 0.5;
  resample_sync_ratio = 1.1;
}

IniValClass::IniValClass() {
  probname = "TwoPeaks";
  probname_ext = "BumpOnTail";
  totalmass = 0;
  rho1 = .9;
  rho2 = .1;
  velocity1[0] = 0.0;
  velocity1[1] = 0.0;
  velocity1[2] = 0.0;
  velocity2[0] = 5.0;
  velocity2[1] = 0.0;
  velocity2[2] = 0.0;

  Tprt1 = 1.0;
  Tprt2 = .01;

  LD_alpha = 0.4;
  // LD_alpha = 0.4;

  BOT_beta = 0.9;
  BOT_rho0 = 1.0;
  BOT_Tprt = 1.0;
  BOT_dTprt = 0.01;
  BOT_ub = 5.0;
  BOT_Tx = 0.25;
}

NumericGridClass::NumericGridClass(int n_x, std::string method) {
  lambda_Poisson = 10.0;

  xmax = 2 * pi / 0.5;
  xmin = 0;
  vmax = 6.0;
  vmin = -vmax;
  tmax = 10;
  // tmax = 8.0;
  // tmax = 0.2;

  Nx = n_x;
  dx = (xmax - xmin) / Nx;
  dt = dx / 2 / vmax;
  if (Nx == 1) dt = 0.001;
  Nt = (int)(tmax / dt);

  if (method == "HDP") {
    Neff = 2.5e-7;
    Neff_F = 1e-5;
    // Neff_F = Neff*20;
  } else {
    Neff = 1e-4;
    Neff_F = 5e-7;
  }

  x.resize(Nx);
  for (int kx = 0; kx < Nx; kx++) x[kx] = xmin + (kx + 0.5) * dx;

  Nv = 200;
  dv = (vmax - vmin) / Nv;
  vx.resize(Nv);
  for (int kv = 0; kv < Nv; kv++) vx[kv] = vmin + (kv + 0.5) * dv;

  bdry_x = 'p';
  bdry_v = 'p';
}

void ParticleGroup::set_xrange(double x1, double x2) {
  xmin = x1;
  xmax = x2;
}

void ParticleGroup::push_back(const Particle1d3d &Spnew) {
  vS.push_back(Spnew);
}

void ParticleGroup::push_back(Particle1d3d *Spnew) { vS.push_back(*Spnew); }

void ParticleGroup::erase(int kp) {
  vS[kp] = vS[vS.size() - 1];
  vS.pop_back();
}

void ParticleGroup::computemoments() {
  int Np = static_cast<int>(vS.size());
  m0 = (double)Np;
  m11 = 0.0;
  m12 = 0.0;
  m13 = 0.0;
  m21 = 0.0;
  m22 = 0.0;
  m23 = 0.0;

  double v1, v2, v3;
  for (int kp = 0; kp < Np; kp++) {
    v1 = vS[kp].velocity(0);
    v2 = vS[kp].velocity(1);
    v3 = vS[kp].velocity(2);

    m11 += v1;
    m12 += v2;
    m13 += v3;
    m21 += v1 * v1;
    m22 += v2 * v2;
    m23 += v3 * v3;
  }

  // m11 = m11/Np;	m12 = m12/Np;	m13 = m13/Np;
  // m21 = m21/Np;	m22 = m22/Np;	m23 = m23/Np;

  m2 = m21 + m22 + m23;
}

void NeParticleGroup::set_xrange(double x1, double x2) {
  xmin = x1;
  xmax = x2;
}

int NeParticleGroup::size(char partype) const {
  int n0 = 0;
  if (partype == 'p') {
    n0 = static_cast<int>(vSp.size());
  } else if (partype == 'n') {
    n0 = static_cast<int>(vSn.size());
  } else if (partype == 'f') {
    n0 = static_cast<int>(vSf.size());
  }

  return n0;
}

void NeParticleGroup::push_back(const Particle1d3d &Snew, char partype) {
  if (partype == 'p') {
    vSp.push_back(Snew);
  } else if (partype == 'n') {
    vSn.push_back(Snew);
  } else if (partype == 'f') {
    vSf.push_back(Snew);
  }
}

void NeParticleGroup::push_back(Particle1d3d *Snew, char partype) {
  if (partype == 'p') {
    vSp.push_back(*Snew);
  } else if (partype == 'n') {
    vSn.push_back(*Snew);
  } else if (partype == 'f') {
    vSf.push_back(*Snew);
  }
}

void NeParticleGroup::erase(int k, char partype) {
  if (partype == 'p') {
    vSp[k] = vSp[vSp.size() - 1];
    vSp.pop_back();
  } else if (partype == 'n') {
    vSn[k] = vSn[vSn.size() - 1];
    vSn.pop_back();
  } else if (partype == 'f') {
    vSf[k] = vSf[vSf.size() - 1];
    vSf.pop_back();
  }
}

void NeParticleGroup::clear(char partype) {
  if (partype == 'p') {
    vSp.clear();
  } else if (partype == 'n') {
    vSn.clear();
  } else if (partype == 'f') {
    vSf.clear();
  }
}

void NeParticleGroup::computemoments() {
  // the moments of P particles,
  int Np = static_cast<int>(vSp.size());
  m0P = 1.0 * Np;
  m11P = 0.0;
  m12P = 0.0;
  m13P = 0.0;
  m21P = 0.0;
  m22P = 0.0;
  m23P = 0.0;
  m31P = 0.0;
  m32P = 0.0;
  m33P = 0.0;

  double v1, v2, v3, vsq;
  for (int kp = 0; kp < Np; kp++) {
    v1 = vSp[kp].velocity(0);
    v2 = vSp[kp].velocity(1);
    v3 = vSp[kp].velocity(2);

    vsq = v1 * v1 + v2 * v2 + v3 * v3;

    m11P += v1;
    m12P += v2;
    m13P += v3;
    m21P += v1 * v1;
    m22P += v2 * v2;
    m23P += v3 * v3;
    m31P += v1 * vsq;
    m32P += v2 * vsq;
    m33P += v3 * vsq;
  }

  m2P = m21P + m22P + m23P;

  // the moments of N particles,
  int Nn = static_cast<int>(vSn.size());
  m0N = 1.0 * Nn;
  m11N = 0.0;
  m12N = 0.0;
  m13N = 0.0;
  m21N = 0.0;
  m22N = 0.0;
  m23N = 0.0;
  m31N = 0.0;
  m32N = 0.0;
  m33N = 0.0;

  for (int kn = 0; kn < Nn; kn++) {
    v1 = vSn[kn].velocity(0);
    v2 = vSn[kn].velocity(1);
    v3 = vSn[kn].velocity(2);

    vsq = (v1 * v1 + v2 * v2 + v3 * v3);

    m11N += v1;
    m12N += v2;
    m13N += v3;
    m21N += v1 * v1;
    m22N += v2 * v2;
    m23N += v3 * v3;
    m31N += v1 * vsq;
    m32N += v2 * vsq;
    m33N += v3 * vsq;
  }

  m2N = m21N + m22N + m23N;

  // the moments of F particles,
  int Nf = static_cast<int>(vSf.size());
  m0F = 1.0 * Nf;
  m11F = 0.0;
  m12F = 0.0;
  m13F = 0.0;
  m21F = 0.0;
  m22F = 0.0;
  m23F = 0.0;
  m31F = 0.0;
  m32F = 0.0;
  m33F = 0.0;

  for (int kf = 0; kf < Nf; kf++) {
    v1 = vSf[kf].velocity(0);
    v2 = vSf[kf].velocity(1);
    v3 = vSf[kf].velocity(2);

    vsq = (v1 * v1 + v2 * v2 + v3 * v3);

    m11F += v1;
    m12F += v2;
    m13F += v3;
    m21F += v1 * v1;
    m22F += v2 * v2;
    m23F += v3 * v3;
    m31F += v1 * vsq;
    m32F += v2 * vsq;
    m33F += v3 * vsq;
  }

  m2F = m21F + m22F + m23F;
}

void NeParticleGroup::copymoments() {
  m0P_o = m0P;
  m11P_o = m11P;
  m12P_o = m12P;
  m13P_o = m13P;
  m2P_o = m2P;
  m21P_o = m21P;
  m22P_o = m22P;
  m23P_o = m23P;
  m31P_o = m31P;
  m32P_o = m32P;
  m33P_o = m33P;

  m0N_o = m0N;
  m11N_o = m11N;
  m12N_o = m12N;
  m13N_o = m13N;
  m2N_o = m2N;
  m21N_o = m21N;
  m22N_o = m22N;
  m23N_o = m23N;
  m31N_o = m31N;
  m32N_o = m32N;
  m33N_o = m33N;

  rho_o = rho;
  u1_o = u1;
  Tprt_o = Tprt;
}

void NeParticleGroup::set_xyzrange() {
  int Np = static_cast<int>(vSp.size());
  int Nn = static_cast<int>(vSn.size());
  for (int k = 0; k < 6; k++) {
    xyz_minmax[k] = 0.;
  }

  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = vSp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = std::min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = std::max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }
  for (int kn = 0; kn < Nn; kn++) {
    auto &v0 = vSn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = std::min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = std::max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }

  for (int k2 = 0; k2 < 3; k2++) {
    xyz_minmax[2 * k2] -= 1e-6;
    xyz_minmax[2 * k2 + 1] += 1e-6;
  }

  // xyz_minmax[0] = 0.; xyz_minmax[2] = 0.; xyz_minmax[4] = 0.;
  // xyz_minmax[1] = 2*pi; xyz_minmax[3] = 2*pi; xyz_minmax[5] = 2*pi;
}

std::vector<Particle1d3d> &NeParticleGroup::list(char partype) {
  if (partype == 'p') {
    return vSp;
  }

  if (partype == 'n') {
    return vSn;
  }

  return vSf;
}

const std::vector<Particle1d3d> &NeParticleGroup::list(char partype) const {
  if (partype == 'p') {
    return vSp;
  }

  if (partype == 'n') {
    return vSn;
  }

  return vSf;
}

Particle1d3d &NeParticleGroup::list(int k, char partype) {
  if (partype == 'p') {
    return vSp[k];
  }

  if (partype == 'n') {
    return vSn[k];
  }
  // else  if (partype == 'f')

  return vSf[k];
}

const Particle1d3d &NeParticleGroup::list(int k, char partype) const {
  if (partype == 'p') {
    return vSp[k];
  }

  if (partype == 'n') {
    return vSn[k];
  }
  // else  if (partype == 'f')

  return vSf[k];
}

}  // namespace coulomb