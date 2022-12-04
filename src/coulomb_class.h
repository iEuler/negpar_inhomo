
using namespace std;

// Global variables
#define pi 3.1415926535897932
int FLAG_PRECOMPUTE_ALPHA_U = 1;
int K_SAVE_TIME = 0;                  // record the file numbers
bool FLAG_FILENAME_WITH_NUM = false;  // append number to the filename

bool FLAG_SAVEFLUX = false;
int NUM_MOVED = 0;
int NUM_RESAMPLE = 0;

double resample_spatial_ratio = 1.0;

double SYNC_TIME = 0.0;

// global variables on time
clock_t t0_all, t1_all, t0_coll, t1_coll, t0_adve, t1_adve, t0_resamp,
    t1_resamp;

template <class T>
void save_macro(const vector<T> &macro, string filename);

// ========================================================================
// group parameters in a class
class ParaClass {
 public:
  int FLAG_USE_OPENMP;
  double dt, coeff_binarycoll, resample_ratio, resample_spatial_ratio,
      sync_time_interval, resample_sync_ratio;
  string method_binarycoll, method;  // method = HDP or PIC
  int Npickup_neg, Nfreq, Nlevel, Num_grids, Num_gridpoints;
  // Npickup_neg: number of particles used in the evaluation of neg part of
  // Delta M Nfreq: parameters used in resampling. Number of frequeces in
  // Fourier interpolation Nlevel: parameters used in resampling. Number of
  // levels in Multi Level version of Fourier interpolation Num_grids =
  // (3*Nlevel*Nlevel - 3*Nlevel + 2)/2: parameters used in resampling.
  int flag_collision;     // 0 = no collision, 1 = Coulomb collison, 2 = BGK
                          // collision
  double lambda_Poisson;  // the lambda used in the Poission equation: 1/\lambda
                          // * \nabla^2 \phi = \rho

  int flag_resample;  // 1 = resample current cell when Nf>=Np+Nn is violated; 2
                      // = resample all cells when it's violated in one cell
  int flag_source;    // 1 = include a source term in Bump On Tail problem

  ParaClass();
};

ParaClass::ParaClass() {
  method = "HDP";
  // method = "PIC";
  FLAG_USE_OPENMP = 1;
  dt = 0.01;
  coeff_binarycoll = 10.0;
  method_binarycoll = "TA";
  resample_ratio = 1.2;
  Npickup_neg = 100;
  Nfreq = 30;
  Nlevel = 4;
  Num_grids = (3 * Nlevel * Nlevel - 3 * Nlevel + 2) / 2;
  Num_gridpoints = (Nlevel * (Nlevel + 1) * (1 << (Nlevel + 2)) +
                    Nlevel * (Nlevel - 1) * (1 << (Nlevel + 1)) +
                    (Nlevel - 1) * (Nlevel - 2) * (1 << Nlevel)) /
                   2;
  flag_collision =
      1;  // 0 = no collision, 1 = Coulomb collison, 2 = BGK collision
  lambda_Poisson = 10.0;
  flag_resample = 2;  // NOT USED
  flag_source = 1;
  resample_spatial_ratio = 0.9;
  sync_time_interval = 0.5;
  resample_sync_ratio = 1.1;
}

// ========================================================================
// specify the parameters in initial distribution in velocity
class IniValClass {
 public:
  string probname, probname_ext;
  double totalmass;

  // for probname=='TwoPeaks'
  double rho1, rho2, Tprt1, Tprt2;
  double velocity1[3], velocity2[3];

  // for probname=='TwoTprtMaxwellian' and 'Maxwellian'
  double rho, Tprt, dTprt;
  double velocity[3];

  // for probname=='TwoStream'
  double TSI_alpha, TSI_coe, TSI_rhof, TSI_rhop, TSI_Tprt, TSI_max_f_over_M,
      TSI_m21, TSI_m22, TSI_m23;

  // for probname=='Rosenbluth'
  // f(v,0) = a*exp( -(|v|-b)^2/c)
  // typical choice: f(v,0) = 0.01*exp( -10* ((|v|-0.3)/0.3)^2)
  double a, b, c;

  // for probname=='LandauDamping'
  // f(v,0) = M, with
  // rho = 1 + alpha*sin(xk);
  // u = 0.;
  // T = 1.;
  double LD_alpha;

  // for probname=='BumpOnTail'
  // f(v,0) = beta / |X| * rho0 / (2*pi*T)^3/2 * exp(-|v|^2/2/T) +
  //          (1-\beta) * rho0 / |X| * exp(-|x-x_c|^2/2/Tx) * 1/(2*pi*dT)^3/2 *
  //          exp(-|v-u_b|^2/2/dT)
  double BOT_beta, BOT_rho0, BOT_Tprt, BOT_dTprt, BOT_Tx, BOT_ub;

  IniValClass();
};

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

// ========================================================================
// specify  parameters for numerical grid
class NumericGridClass {
 public:
  int Nx, Nt, Nv;
  double xmax, xmin, vmax, vmin, tmax, dx, dv, dt, Neff,
      Neff_F;  // Neff is the effective number
  // double *x, *vx;
  vector<double> x, vx;
  char bdry_x, bdry_v;
  double lambda_Poisson;

  // NumericGridClass();
  NumericGridClass();         // default n_x = 100
  NumericGridClass(int n_x);  // n_x is the number of grid in x direction
  NumericGridClass(int n_x,
                   string method);  // n_x is the number of grid in x direction
};

NumericGridClass::NumericGridClass(int n_x, string method) {
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

NumericGridClass::NumericGridClass(int n_x) { NumericGridClass(n_x, "HDP"); }

NumericGridClass::NumericGridClass() { NumericGridClass(100); }

// ========================================================================
// define particle class in 1D x and 3D v

class Particle1d3d {
 public:
  Particle1d3d();
  Particle1d3d(double *);          // with v
  Particle1d3d(double, double *);  // with x and v
  void set_position(double);
  void set_velocity(double *);

  double position() { return x; }
  double *velocity() { return v; }
  double velocity(int k) { return *(v + k); }

  bool flag_moved;

 private:
  double x, v[3];
};

Particle1d3d::Particle1d3d() {
  x = 0;
  *v = 0;
  *(v + 1) = 0;
  *(v + 2) = 0;
  flag_moved = false;
}

Particle1d3d::Particle1d3d(double *v0) {
  x = 0;
  *v = *v0;
  *(v + 1) = *(v0 + 1);
  *(v + 2) = *(v0 + 2);
  flag_moved = false;
}

Particle1d3d::Particle1d3d(double x0, double *v0) {
  x = x0;
  *v = *v0;
  *(v + 1) = *(v0 + 1);
  *(v + 2) = *(v0 + 2);
  flag_moved = false;
}

void Particle1d3d::set_position(double x0) { x = x0; }

void Particle1d3d::set_velocity(double *v0) {
  *v = *v0;
  *(v + 1) = *(v0 + 1);
  *(v + 2) = *(v0 + 2);
}

// ========================================================================
// define particle class in 1D x and 3D v
// ========================================================================
// define particle class in 1D x and 3D v

class ParticleGroup {
 public:
  double m0, m11, m12, m13, m2, m21, m22, m23;  // the moments, with Neff = 1
  double elecfield;

  ParticleGroup() : xmin(0.), xmax(0.){};

  void set_xrange(double, double);

  double get_xmin() { return xmin; }
  double get_xmax() { return xmax; }

  int size() { return static_cast<int>(vS.size()); }
  void push_back(Particle1d3d *);        // add one more particle
  void push_back(const Particle1d3d &);  // add one more particle
  void erase(int);                       // remove the int particle

  auto &list() { return vS; }
  Particle1d3d &list(int nn) { return vS[nn]; }
  void computemoments();  // compute the moments

 private:
  double xmin, xmax;  // the range of particle positions in x space
  vector<Particle1d3d> vS;
};

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

// ========================================================================
// define P, N, F particle class in 1D x and 3D v

class NeParticleGroup {
 public:
  // public variables
  double rhoM, u1M, u2M, u3M, TprtM;  // the moments of maxwellian part
  double T1M, T2M, T3M;               // Used in resampling coarse particles
  double rhoF, u1F, u2F, u3F, TprtF;  // the moments of maxwellian part

  double dx_rhoM, dx_u1M,
      dx_TprtM;  // derivative in x direction, used in sampling from source part

  double rho, u1, Tprt;  // the moments of all distributions
  // double dx_rho, dx_u1, dx_Tprt; // derivative in x direction, used in
  // sampling from source part

  double rho_o, u1_o,
      Tprt_o;  // the moments of all distributions before advection

  double drho, dm1, denergy;  // change in moments: rho^{n+1} = rho^n - drho
  double drho_g, dm1_g, denergy_g;  // change in moments due to g = f_p - f_n,
                                    // used in sampling from source part

  double m0P, m11P, m12P, m13P, m2P, m21P, m22P, m23P, m31P, m32P,
      m33P;  // the moments of P particles, with Neff = 1
  double m0N, m11N, m12N, m13N, m2N, m21N, m22N, m23N, m31N, m32N,
      m33N;  // the moments of N particles, with Neff = 1
  double m0F, m11F, m12F, m13F, m2F, m21F, m22F, m23F, m31F, m32F,
      m33F;  // the moments of F particles, with Neff = 1

  double m0P_o, m11P_o, m12P_o, m13P_o, m2P_o, m21P_o, m22P_o, m23P_o, m31P_o,
      m32P_o,
      m33P_o;  // the moments of P particles before advection, with Neff = 1
  double m0N_o, m11N_o, m12N_o, m13N_o, m2N_o, m21N_o, m22N_o, m23N_o, m31N_o,
      m32N_o,
      m33N_o;  // the moments of N particles before advection, with Neff = 1

  double elecfield, elecfield_F;

  double alpha_neg, alpha_pos, rmax;

  double xyz_minmax[6];  // record the information on [xmin, xmax, ymin, ymax,
                         // zmin, zmax]

  int flag_resampled;

  // public functions

  NeParticleGroup() : xmin(0.), xmax(0.){};

  void set_xrange(double, double);
  double get_xmin() { return xmin; }
  double get_xmax() { return xmax; }

  int size(char);
  void push_back(Particle1d3d *Snew, char partype);
  void push_back(const Particle1d3d &,
                 char);   // add one more particle in char type
  void erase(int, char);  // remove the int particle in char type
  void clear(char);       // remove every particle in char type

  vector<Particle1d3d> &list(char);
  Particle1d3d &list(int, char);
  void computemoments();   // compute the moments
  void computexyzrange();  // update xyz_minmax = [xmin, xmax, ymin, ymax, zmin,
                           // zmax]
  void copymoments();      // copy m0P ... to m0P_o

  void reset_flag_resampled() { flag_resampled = 0; }

 private:
  double xmin, xmax;  // the range of particle positions in x space
  vector<Particle1d3d> vSp, vSn, vSf;
};

void NeParticleGroup::set_xrange(double x1, double x2) {
  xmin = x1;
  xmax = x2;
}

int NeParticleGroup::size(char partype) {
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

void NeParticleGroup::computexyzrange() {
  int Np = static_cast<int>(vSp.size());
  int Nn = static_cast<int>(vSn.size());
  for (int k = 0; k < 6; k++) {
    xyz_minmax[k] = 0.;
  }

  for (int kp = 0; kp < Np; kp++) {
    double *v0 = vSp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }
  for (int kn = 0; kn < Nn; kn++) {
    double *v0 = vSn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }

  for (int k2 = 0; k2 < 3; k2++) {
    xyz_minmax[2 * k2] -= 1e-6;
    xyz_minmax[2 * k2 + 1] += 1e-6;
  }

  // xyz_minmax[0] = 0.; xyz_minmax[2] = 0.; xyz_minmax[4] = 0.;
  // xyz_minmax[1] = 2*pi; xyz_minmax[3] = 2*pi; xyz_minmax[5] = 2*pi;
}

vector<Particle1d3d> &NeParticleGroup::list(char partype) {
  if (partype == 'p') {
    return vSp;
  }

  if (partype == 'n') {
    return vSn;
  }

  return vSf;
}

Particle1d3d &NeParticleGroup::list(int k, char partype) {
  Particle1d3d S0;
  if (partype == 'p') {
    S0 = vSp[k];
  } else if (partype == 'n') {
    S0 = vSn[k];
  } else if (partype == 'f') {
    S0 = vSf[k];
  }

  return S0;
}