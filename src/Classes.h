#pragma once

#include <string>
#include <vector>

namespace coulomb {
// ========================================================================
// group parameters in a class
class ParaClass {
 public:
  int FLAG_USE_OPENMP;
  double dt, coeff_binarycoll, resample_ratio, resample_spatial_ratio,
      sync_time_interval, resample_sync_ratio;
  std::string method_binarycoll, method;  // method = HDP or PIC
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

// ========================================================================
// specify the parameters in initial distribution in velocity
class IniValClass {
 public:
  std::string probname, probname_ext;
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

// ========================================================================
// specify  parameters for numerical grid
class NumericGridClass {
 public:
  int Nx, Nt, Nv;
  double xmax, xmin, vmax, vmin, tmax, dx, dv, dt, Neff,
      Neff_F;  // Neff is the effective number
  // double *x, *vx;
  std::vector<double> x, vx;
  char bdry_x, bdry_v;
  double lambda_Poisson;

  NumericGridClass(
      int n_x,
      std::string method);  // n_x is the number of grid in x direction

  NumericGridClass(int n_x) {
    NumericGridClass(n_x, "HDP");
  };  // n_x is the number of grid in x direction

  NumericGridClass() { NumericGridClass(100); };  // default n_x = 100
};

// ========================================================================
// define particle class in 1D x and 3D v

class Particle1d3d {
 public:
  Particle1d3d() : flag_moved(false){};
  Particle1d3d(const std::vector<double> &v0)
      : x(0.0), v(v0), flag_moved(false){};  // with v
  Particle1d3d(double x0, const std::vector<double> &v0)
      : x(x0), v(v0), flag_moved(false){};  // with x and v
  void set_position(double x0) { x = x0; };
  void set_velocity(const std::vector<double> &v0) { v = v0; };
  void set_velocity(double *v0) {
    v[0] = *v0;
    v[1] = *(v0 + 1);
    v[2] = *(v0 + 2);
  };

  double position() { return x; }
  auto &velocity() { return v; }
  double velocity(int k) { return v[k]; }

  bool flag_moved;

 private:
  double x{0.0};
  std::vector<double> v{0.0, 0.0, 0.0};
};

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
  std::vector<Particle1d3d> vS;
};

// ========================================================================
// define P, N, F particle class in 1D x and 3D v

class NeParticleGroup {
 public:
  // public variables
  double rhoM, u1M, u2M, u3M, TprtM;  // the moments of maxwellian part
  double T1M, T2M, T3M;               // Used in resampling coarse particles
  double rhoF, u1F, u2F, u3F, TprtF;  // the moments of maxwellian part

  double dx_rhoM, dx_u1M,
      dx_TprtM;  // derivative in x direction, used in sampling from source
                 // part

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

  std::vector<Particle1d3d> &list(char);
  Particle1d3d &list(int, char);
  void computemoments();   // compute the moments
  void computexyzrange();  // update xyz_minmax = [xmin, xmax, ymin, ymax,
                           // zmin, zmax]
  void copymoments();      // copy m0P ... to m0P_o

  void reset_flag_resampled() { flag_resampled = 0; }

 private:
  double xmin, xmax;  // the range of particle positions in x space
  std::vector<Particle1d3d> vSp, vSn, vSf;
};

}  // namespace coulomb