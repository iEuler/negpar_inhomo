#include "Classes.h"
#include "coulomb_fwd_declare.h"
#include "utils.h"
// ========================================================================

namespace coulomb {
// perform Coulomb collisions between two particle velocities

std::pair<std::vector<double>, std::vector<double>> coulombBinary3d(
    const std::vector<double>& v1, const std::vector<double>& v2,
    const ParaClass& para) {
  // need para.dt, para.method_binarycoll and para.coeff_binarycoll
  double dt = para.dt;
  double coeff = para.coeff_binarycoll;

  std::vector<double> v1p(3), v2p(3);

  if (para.method_binarycoll.compare("TA") == 0) {
    std::vector<double> u(3), du(3);
    for (int k = 0; k < 3; k++) u[k] = v1[k] - v2[k];

    double g = 0;
    for (int k = 0; k < 3; k++) g += u[k] * u[k];
    g = sqrt(g);

    double sigma2_delta = coeff * dt / (g * g * g);

    double deltak = sqrt(sigma2_delta) * myrandn();
    double phi = 2.0 * pi * myrand();
    double sintheta = 2.0 * deltak / (1.0 + deltak * deltak);
    double costheta = 1.0 - 2.0 * deltak * deltak / (1.0 + deltak * deltak);

    double uperp = sqrt(u[0] * u[0] + u[1] * u[1]) + 1e-10;
    du[0] = (u[0] / uperp) * u[2] * sintheta * cos(phi) -
            (u[1] / uperp) * g * sintheta * sin(phi) - u[0] * (1.0 - costheta);
    du[1] = (u[1] / uperp) * u[2] * sintheta * cos(phi) +
            (u[0] / uperp) * g * sintheta * sin(phi) - u[1] * (1.0 - costheta);
    du[2] = -uperp * sintheta * cos(phi) - u[2] * (1.0 - costheta);

    for (int k = 0; k < 3; k++) {
      v1p[k] = v1[k] + 0.5 * du[k];
      v2p[k] = v2[k] - 0.5 * du[k];
    }
  }

  return {std::move(v1p), std::move(v2p)};
}

// ========================================================================

// Perform coulomb collisions in homogeneous case

void coulomb_collision_homo(std::vector<Particle1d3d>& Sp, int Np,
                            const ParaClass& para) {
  const auto p = myrandperm(Np, Np);
  int kp1, kp2;

  for (int kp = 0; kp < Np / 2; kp++) {
    kp1 = p[2 * kp] - 1;
    kp2 = p[2 * kp + 1] - 1;
    auto& v1 = Sp[kp1].velocity();
    auto& v2 = Sp[kp2].velocity();

    const auto& vp = coulombBinary3d(v1, v2, para);

    Sp[kp1].set_velocity(vp.first);
    Sp[kp2].set_velocity(vp.second);
  }
}

// ========================================================================

// given rho, (u, T) -> ( momentum, energy)

void uT2mE(const double& rho, double* u, const double& Tprt, double* momentum,
           double& energy) {
  for (int k = 0; k < 3; k++) momentum[k] = rho * u[k];
  double e0 = 0.0;
  for (int k = 0; k < 3; k++) e0 += u[k] * u[k];
  energy = 0.5 * rho * (e0 + 3.0 * Tprt);
}

// given rho, ( momentum, energy) -> (u, T)

void mE2uT(const double& rho, double* momentum, const double& energy, double* u,
           double& Tprt) {
  for (int k = 0; k < 3; k++) u[k] = momentum[k] / rho;
  double e0 = 0.0;
  for (int k = 0; k < 3; k++) e0 += u[k] * u[k];
  Tprt = (energy * 2.0 / rho - e0) / 3.0;
}

// ========================================================================

/**
  given rho, (u, T) -> ( momentum, energy) in 1d x, 3d v
*/

void uT2mE_x1v3(int Nx, const vector<double>& rho, const vector<double>& u1,
                const vector<double>& Tprt, vector<double>& m1,
                vector<double>& energy) {
  double u_temp[3] = {0., 0., 0.};
  double momentum_temp[3] = {0., 0., 0.};
  double energy_temp;
  for (int kx = 0; kx < Nx; kx++) {
    u_temp[0] = u1[kx];
    uT2mE(rho[kx], u_temp, Tprt[kx], momentum_temp, energy_temp);
    m1[kx] = momentum_temp[0];
    energy[kx] = energy_temp;
  }
}
// ========================================================================

/**
  given rho, ( momentum, energy) -> (u, T) in 1d x, 3d v
*/

void mE2uT_x1v3(int Nx, const vector<double>& rho, const vector<double>& m1,
                const vector<double>& energy, vector<double>& u1,
                vector<double>& Tprt) {
  double u_temp[3] = {0., 0., 0.};
  double momentum_temp[3] = {0., 0., 0.};
  double Tprt_temp;
  for (int kx = 0; kx < Nx; kx++) {
    momentum_temp[0] = m1[kx];
    mE2uT(rho[kx], momentum_temp, energy[kx], u_temp, Tprt_temp);
    u1[kx] = u_temp[0];
    Tprt[kx] = Tprt_temp;
  }
}

// ========================================================================

// compute the density, momentum and energy from the particle group
void particle2rhomE(ParticleGroup* Sp_x, double Neff, double& rho,
                    double* momentum, double& energy) {
  rho = (Sp_x->m0) * Neff;
  *momentum = (Sp_x->m11) * Neff;
  *(momentum + 1) = (Sp_x->m12) * Neff;
  *(momentum + 2) = (Sp_x->m13) * Neff;
  energy = 0.5 * (Sp_x->m2) * Neff;
}

// ========================================================================

// compute the density, momentum and energy from the particle group in the space
// inhomo grid
void compute_rhouT(int Nx, ParticleGroup* Sp_x, double Neff,
                   vector<double>& rho, vector<double>& u1, vector<double>& u2,
                   vector<double>& u3, vector<double>& Tprt) {
  double rhok, Tprtk, energyk;

  double momentumk[3];
  double uk[3];

  for (int kx = 0; kx < Nx; kx++) {
    particle2rhomE(Sp_x + kx, Neff, rhok, momentumk, energyk);
    mE2uT(rhok, momentumk, energyk, uk, Tprtk);
    rho[kx] = rhok;
    u1[kx] = *uk;
    u2[kx] = *(uk + 1);
    u3[kx] = *(uk + 2);
    Tprt[kx] = Tprtk;
  }
}

// ========================================================================
// compute electric filed

/**
        A 1-D Poisson solver
        Input: \rho, Nx, xdomainsize
        Output: E
        Solves: - \nabla^2 \Phi = \rho
                        E =  lambda \nabla \Phi

*/
void PoissonSolver(const vector<double>& rho, vector<double>& elecfield, int Nx,
                   double xdomainsize, double lambda) {
  fftw_complex* rho_c = new fftw_complex[Nx];
  fftw_complex* Frho = new fftw_complex[Nx];
  fftw_complex* Felec = new fftw_complex[Nx];
  fftw_complex* IFelec = new fftw_complex[Nx];
  fftw_plan p_fft, p_ifft;

  p_fft = fftw_plan_dft_1d(Nx, rho_c, Frho, FFTW_FORWARD, FFTW_ESTIMATE);
  p_ifft = fftw_plan_dft_1d(Nx, Felec, IFelec, FFTW_BACKWARD, FFTW_ESTIMATE);

  // forward FFT
  for (int kx = 0; kx < Nx; kx++) {
    rho_c[kx][0] = rho[kx];
    rho_c[kx][1] = 0;
  }

  fftw_execute(p_fft);

  // compute Felec
  vector<double> lx(Nx);
  for (int kx = 0; kx < Nx / 2 + 1; kx++) lx[kx] = kx * 2.0 * pi / xdomainsize;
  for (int kx = Nx / 2 + 1; kx < Nx; kx++)
    lx[kx] = (kx - Nx) * 2.0 * pi / xdomainsize;

  Felec[0][0] = 0;
  Felec[0][1] = 0;
  for (int kx = 1; kx < Nx; kx++) {
    Felec[kx][0] = Frho[kx][1] / lx[kx] / Nx;
    Felec[kx][1] = -Frho[kx][0] / lx[kx] / Nx;
  }

  // backward FFT
  fftw_execute(p_ifft);

  for (int kx = 0; kx < Nx; kx++) elecfield[kx] = lambda * IFelec[kx][0];

  fftw_destroy_plan(p_fft);
  fftw_destroy_plan(p_ifft);
  // fftw_free(rho_c); fftw_free(Frho);	fftw_free(Felec); fftw_free(IFelec);
  delete[] rho_c;
  delete[] Frho;
  delete[] Felec;
  delete[] IFelec;
}

// Update electric filed, without negative particle

void updateelecfiled(std::vector<ParticleGroup>& Sp_x,
                     const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) Sp_x[kx].computemoments();

  // compute rho
  vector<double> rho(Nx);
  vector<double> elecfield(Nx);

  for (int kx = 0; kx < Nx; kx++) rho[kx] = Sp_x[kx].m0 * grid.Neff / grid.dx;

  double lambda = 10.0;
  PoissonSolver(rho, elecfield, grid.Nx, grid.xmax - grid.xmin, lambda);

  for (int kx = 0; kx < Nx; kx++) Sp_x[kx].elecfield = elecfield[kx];
}

// Update electric filed, with negative particle

void updateelecfiled(std::vector<NeParticleGroup>& S_x,
                     const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // compute rho
  vector<double> rho(Nx);
  vector<double> elecfield(Nx);
  for (int kx = 0; kx < Nx; kx++)
    rho[kx] = S_x[kx].rhoM + (S_x[kx].m0P - S_x[kx].m0N) * grid.Neff / grid.dx;

  // double lambda = 10.0;
  double lambda = grid.lambda_Poisson;
  PoissonSolver(rho, elecfield, grid.Nx, grid.xmax - grid.xmin, lambda);

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield = elecfield[kx];

  // for the F particles
  // for (int kx = 0; kx<Nx; kx++)	(S_x+kx)->elecfield_F = *(elecfield+kx);
  for (int kx = 0; kx < Nx; kx++) rho[kx] = S_x[kx].m0F * grid.Neff_F / grid.dx;
  PoissonSolver(rho, elecfield, grid.Nx, grid.xmax - grid.xmin, lambda);
  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield_F = elecfield[kx];
}

void updateelecfiled_PIC(std::vector<NeParticleGroup>& S_x,
                         const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // compute rho
  vector<double> rho(Nx);
  vector<double> elecfield(Nx);

  double lambda = grid.lambda_Poisson;

  // for the F particles
  // for (int kx = 0; kx<Nx; kx++)	(S_x+kx)->elecfield_F = *(elecfield+kx);
  for (int kx = 0; kx < Nx; kx++) rho[kx] = S_x[kx].m0F * grid.Neff_F / grid.dx;
  PoissonSolver(rho, elecfield, grid.Nx, grid.xmax - grid.xmin, lambda);
  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield_F = elecfield[kx];
}

void updateelecfiled_fromcoarse(std::vector<NeParticleGroup>& S_x,
                                const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // compute rho
  vector<double> rho(Nx);
  vector<double> elecfield(Nx);

  for (int kx = 0; kx < Nx; kx++) rho[kx] = S_x[kx].m0F * grid.Neff_F / grid.dx;

  // double lambda = 10.0;
  double lambda = grid.lambda_Poisson;
  PoissonSolver(rho, elecfield, grid.Nx, grid.xmax - grid.xmin, lambda);

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield = elecfield[kx];
}
void updateelecfiled_zero(std::vector<NeParticleGroup>& S_x,
                          const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // compute rho

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield = 0.;

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield_F = 0.;
}

void updateelecfiled_rho(std::vector<NeParticleGroup>& S_x,
                         const NumericGridClass& grid) {
  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

  // compute rho

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield = S_x[kx].rho;

  for (int kx = 0; kx < Nx; kx++) S_x[kx].elecfield_F = S_x[kx].rhoF;
}

// move particles according to velocity

Particle1d3d moveparticle(const Particle1d3d& Sp, double elecfield,
                          const NumericGridClass& grid) {
  Particle1d3d SpMoved;
  double xnew = Sp.position() + grid.dt * Sp.velocity(0);

  auto vnew = Sp.velocity();
  double vxnew = vnew[0] + grid.dt * elecfield;
  // cout << elecfield << endl;

  if (grid.bdry_x == 'p') {
    double xperiod = grid.xmax - grid.xmin;
    while (xnew >= grid.xmax) xnew -= xperiod;
    while (xnew < grid.xmin) xnew += xperiod;
  } else if (grid.bdry_x == 'n') {
    while ((xnew >= grid.xmax) || (xnew < grid.xmin)) {
      if (xnew >= grid.xmax) {
        xnew = 2 * grid.xmax - xnew;
        vxnew = -vxnew;
      }
      if (xnew < grid.xmin) {
        xnew = 2 * grid.xmin - xnew;
        vxnew = -vxnew;
      }
    }
  }

  // cout << *vnew << ' ' << vxnew << endl;

  vnew[0] = vxnew;
  NUM_MOVED++;

  return Particle1d3d(xnew, vnew, true);
}

// find the group according to particle position

int findparticlegroup(Particle1d3d* Sp, const NumericGridClass& grid) {
  int kx_after = (int)((Sp->position() - grid.xmin) / grid.dx);
  if (kx_after >= grid.Nx) {
    cout << "error: Particle moved out of range. kx = " << kx_after << endl;
    kx_after = grid.Nx;
  }
  if (kx_after < 0) {
    cout << "error: Particle moved out of range. kx = " << kx_after << endl;
    kx_after = 0;
    exit(0);
  }
  return kx_after;
}

// place particle in the correct group
// move the kp-th particle in group kx_before to group kx_after

void relocateparticle(ParticleGroup* Sp_x, int kx_before, int kp,
                      int kx_after) {
  if (kx_before != kx_after) {
    (Sp_x + kx_after)->push_back((Sp_x + kx_before)->list().at(kp));
    (Sp_x + kx_before)->erase(kp);
  }
}

void relocateparticle(NeParticleGroup* S_x, char partype, int kx_before, int kp,
                      int kx_after) {
  if (kx_before != kx_after) {
    (S_x + kx_after)
        ->push_back((S_x + kx_before)->list(partype).at(kp), partype);
    (S_x + kx_before)->erase(kp, partype);
  }
}

// reset the flag_moved

void reset_flag_moved(ParticleGroup* Sp_x, int Nx) {
  for (int kx = 0; kx < Nx; kx++) {
    auto& Sp = (Sp_x + kx)->list();
    for (int kp = 0; kp < (Sp_x + kx)->size(); kp++) Sp[kp].flag_moved = false;
  }
}

void reset_flag_moved(NeParticleGroup* S_x, char partype, int Nx) {
  for (int kx = 0; kx < Nx; kx++) {
    auto& Sp = (S_x + kx)->list(partype);
    for (int kp = 0; kp < (S_x + kx)->size(partype); kp++) {
      if (!(Sp[kp].flag_moved)) cout << "NOT MOVED" << endl;
      Sp[kp].flag_moved = false;
    }
  }
}

// main program on particles advection
void particleadvection(std::vector<ParticleGroup>& Sp_x,
                       const NumericGridClass& grid) {
  updateelecfiled(Sp_x, grid);

  for (int kx = 0; kx < grid.Nx; kx++) {
    // cout << "move kx = " << kx << " [" <<(Sp_x+kx)->get_xmin() << ','
    // <<(Sp_x+kx)->get_xmax() << "]" << endl;
    auto& Sp = Sp_x[kx].list();
    double elecfield = Sp_x[kx].elecfield;
    int kp = 0;
    while (kp < Sp_x[kx].size()) {
      // cout << kp+1 << " / " << (Sp_x+kx)->size() << ' ' <<(Sp +
      // kp)->flag_moved << ' ' << Sp[kp].position();
      if (!(Sp[kp].flag_moved)) {
        // cout << " old x = " << Sp[kp].position() << " vx = " << (Sp +
        // kp)->velocity(0);
        Sp[kp] = moveparticle(Sp[kp], elecfield, grid);
        // cout << " new x = " << Sp[kp].position() << endl;
        // cout << " new x = " << Sp[kp].position();
        int kx_after = findparticlegroup(&Sp[kp], grid);
        relocateparticle(&Sp_x[0], kx, kp, kx_after);
        // cout << " new kx = " << kx_after << endl;
      } else {
        kp++;
      }
      // cout << endl;
    }
  }
  // cout << "a4" << endl;
  reset_flag_moved(&Sp_x[0], grid.Nx);
}

// main program on particles advection, with negative particles

void particleadvection(std::vector<NeParticleGroup>& S_x, char partype,
                       const NumericGridClass& grid) {
  // cout << "a1" << endl;

  for (int kx = 0; kx < grid.Nx; kx++) {
    // cout << "move kx = " << kx << " [" <<(S_x+kx)->get_xmin() << ','
    // <<(S_x+kx)->get_xmax() << "]" << endl;
    auto& Sp = S_x[kx].list(partype);
    double elecfield = S_x[kx].elecfield;
    if (partype == 'f') elecfield = S_x[kx].elecfield_F;
    int kp = 0;
    while (kp < S_x[kx].size(partype)) {
      if (!(Sp[kp].flag_moved)) {
        Sp[kp] = moveparticle(Sp[kp], elecfield, grid);
        int kx_after = findparticlegroup(&Sp[kp], grid);
        relocateparticle(&S_x[0], partype, kx, kp, kx_after);
      } else {
        kp++;
      }
      // cout << endl;
    }
  }
  // cout << "a7" << endl;
  // cout << "a4" << endl;
  reset_flag_moved(&S_x[0], partype, grid.Nx);
  cout << "Number of particle moved = " << NUM_MOVED << endl;
  NUM_MOVED = 0;
}

void particleadvection(std::vector<NeParticleGroup>& S_x,
                       const NumericGridClass& grid) {
  // cout << "advect p" << endl;
  particleadvection(S_x, 'p', grid);
  // cout << "advect n" << endl;
  particleadvection(S_x, 'n', grid);
  // cout << "advect f" << endl;
  particleadvection(S_x, 'f', grid);
}

// ========================================================================

// move particles according to velocity
void limiter_x1_o2(const vector<double>& f, int Nx, double dx, double dt,
                   int sgn_v, const vector<double>& f_bd, char bdry,
                   vector<double>& df);

double myerf(double x);

void Euler_kinetic_x1(const vector<double>& rho, const vector<double>& u1,
                      const vector<double>& Tprt, int Nx, double dx, double dt,
                      char bdry, vector<double>& drho, vector<double>& dm1,
                      vector<double>& denergy) {
  // kinetic scheme for 1D compressible Euler system U_t + (F(U))_x = 0

  vector<double> Ap(Nx);
  vector<double> Am(Nx);
  vector<double> B(Nx);
  vector<double> G1p(Nx);
  vector<double> G1m(Nx);
  vector<double> G2p(Nx);
  vector<double> G2m(Nx);
  vector<double> G3p(Nx);
  vector<double> G3m(Nx);
  vector<double> dG1p(Nx);
  vector<double> dG1m(Nx);
  vector<double> dG2p(Nx);
  vector<double> dG2m(Nx);
  vector<double> dG3p(Nx);
  vector<double> dG3m(Nx);

  vector<double> bd_G1p(Nx);
  vector<double> bd_G1m(Nx);
  vector<double> bd_G3p(Nx);
  vector<double> bd_G3m(Nx);

  const double lambda_state = 1.;  // corresponding to Dim = 3
  // const double lambda_state = 0; // corresponding to Dim = 1

  for (int kx = 0; kx < Nx; kx++) {
    Ap[kx] = .5 * (1.0 + myerf(u1[kx] / sqrt(2.0 * Tprt[kx])));
    Am[kx] = 1.0 - Ap[kx];
    B[kx] = sqrt(Tprt[kx] / 2.0 / pi) * exp(-.5 * u1[kx] * u1[kx] / Tprt[kx]);
    G1p[kx] = rho[kx] * (u1[kx] * Ap[kx] + B[kx]);
    G1m[kx] = rho[kx] * (u1[kx] * Am[kx] - B[kx]);
    G2p[kx] =
        rho[kx] * ((Tprt[kx] + u1[kx] * u1[kx]) * Ap[kx] + u1[kx] * B[kx]);
    G2m[kx] =
        rho[kx] * ((Tprt[kx] + u1[kx] * u1[kx]) * Am[kx] - u1[kx] * B[kx]);
    G3p[kx] =
        rho[kx] * ((1.5 * Tprt[kx] + .5 * u1[kx] * u1[kx]) * u1[kx] * Ap[kx] +
                   (Tprt[kx] + .5 * u1[kx] * u1[kx]) * B[kx]);
    G3m[kx] =
        rho[kx] * ((1.5 * Tprt[kx] + .5 * u1[kx] * u1[kx]) * u1[kx] * Am[kx] -
                   (Tprt[kx] + .5 * u1[kx] * u1[kx]) * B[kx]);

    G3p[kx] += lambda_state * Tprt[kx] * G1p[kx];
    G3m[kx] += lambda_state * Tprt[kx] * G1m[kx];

    bd_G1p[kx] = -G1m[kx];
    bd_G1m[kx] = -G1p[kx];
    bd_G3p[kx] = -G3m[kx];
    bd_G3m[kx] = -G3p[kx];
  }
  /*
  limiter_x1_o2(G1p, Nx, dx, dt,  1, bd_G1p, bdry, dG1p);
  limiter_x1_o2(G1m, Nx, dx, dt, -1, bd_G1m, bdry, dG1m);
  limiter_x1_o2(G2p, Nx, dx, dt,  1, G2m, bdry, dG2p);
  limiter_x1_o2(G2m, Nx, dx, dt, -1, G2p, bdry, dG2m);
  limiter_x1_o2(G3p, Nx, dx, dt,  1, bd_G3p, bdry, dG3p);
  limiter_x1_o2(G3m, Nx, dx, dt, -1, bd_G3m, bdry, dG3m);
  */
  limiter_x1_o2(G1p, Nx, dx, dt, 1, G1p, bdry, dG1p);
  limiter_x1_o2(G1m, Nx, dx, dt, -1, G1m, bdry, dG1m);
  limiter_x1_o2(G2p, Nx, dx, dt, 1, G2p, bdry, dG2p);
  limiter_x1_o2(G2m, Nx, dx, dt, -1, G2m, bdry, dG2m);
  limiter_x1_o2(G3p, Nx, dx, dt, 1, G3p, bdry, dG3p);
  limiter_x1_o2(G3m, Nx, dx, dt, -1, G3m, bdry, dG3m);

  for (int kx = 0; kx < Nx; kx++) {
    drho[kx] = dG1p[kx] + dG1m[kx];
    dm1[kx] = dG2p[kx] + dG2m[kx];
    denergy[kx] = dG3p[kx] + dG3m[kx];
  }
}

// ========================================================================

/**
  slope limiter
*/

void cshift_1d(const vector<double>& f, int size, int shiftdistance,
               vector<double>& f_shift);
void eoshift_1d(const vector<double>& f, int size, int shiftdistance,
                vector<double>& f_shift, double boundary);

void limiter_x1_o2(const vector<double>& f, int Nx, double dx, double dt,
                   int sgn_v, const vector<double>& f_bd, char bdry,
                   vector<double>& df) {
  // sgn_v = 1 or -1

  vector<double> f_L1(Nx);
  vector<double> f_R1(Nx);

  if (bdry == 'p') {
    cshift_1d(f, Nx, -1, f_L1);
    cshift_1d(f, Nx, 1, f_R1);
  } else {
    eoshift_1d(f, Nx, -1, f_L1, f_bd[0]);
    eoshift_1d(f, Nx, 1, f_R1, f_bd[Nx - 1]);
  }

  vector<double> df_upwind(Nx);
  vector<double> df_downwind(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    df_upwind[kx] = f[kx] - f_L1[kx];
    df_downwind[kx] = f_R1[kx] - f[kx];
  }

  if (sgn_v > 0) {
    for (int kx = 0; kx < Nx; kx++) df[kx] = dt / dx * df_upwind[kx];
    // df = dt/dx * (f-f_L1 + 0.5d0 * (f_sigma - f_sigman));
  } else {
    for (int kx = 0; kx < Nx; kx++) df[kx] = dt / dx * df_downwind[kx];
    // df = dt/dx * (f_R1-f - 0.5d0 * (f_sigmap - f_sigma));
  }
}

// ========================================================================
/**
  1d cshift function as in Fortran
*/

void cshift_1d(const vector<double>& f, int size, int shiftdistance,
               vector<double>& f_shift) {
  // shiftdistance>0, shift to the left
  // shiftdistance<0, shift to the write
  for (int kx = 0; kx < size; kx++) {
    int kx_old = kx + shiftdistance;
    if (kx_old >= size) kx_old -= size;
    if (kx_old < 0) kx_old += size;
    f_shift[kx] = f[kx_old];
  }
}

// ========================================================================
/**
  1d eoshift function as in Fortran
*/

void eoshift_1d(const vector<double>& f, int size, int shiftdistance,
                vector<double>& f_shift, double boundary) {
  // shiftdistance>0, shift to the left
  // shiftdistance<0, shift to the write
  for (int kx = 0; kx < size; kx++) {
    int kx_old = kx + shiftdistance;
    double fnew;
    if ((kx_old >= size) || (kx_old < 0))
      fnew = boundary;
    else
      fnew = f[kx_old];

    f_shift[kx] = fnew;
  }
}

// ========================================================================

/**
  compute the central difference of rho
*/

void diff_1d_central(const vector<double>& rho, int Nx, char bdry,
                     vector<double>& drho) {
  vector<double> rho_L1(Nx);
  vector<double> rho_R1(Nx);

  if (bdry == 'p') {
    cshift_1d(rho, Nx, -1, rho_L1);
    cshift_1d(rho, Nx, 1, rho_R1);
  } else {
    eoshift_1d(rho, Nx, -1, rho_L1, rho[0]);
    eoshift_1d(rho, Nx, 1, rho_R1, rho[Nx - 1]);
  }

  for (int kx = 0; kx < Nx; kx++) drho[kx] = .5 * (rho_R1[kx] - rho_L1[kx]);
}

// ========================================================================

/**
  Update rho, u, T from M, P and N
  input: S_x, grid
  output: none
*/

void update_rhouT(NeParticleGroup* S_x, const NumericGridClass& grid) {
  int Nx = grid.Nx;
  double Neff = grid.Neff;
  double dx = grid.dx;

  vector<double> rho_m(Nx);
  vector<double> u1_m(Nx);
  vector<double> Tprt_m(Nx);

  vector<double> rho_pn(Nx);
  vector<double> m1_pn(Nx);
  vector<double> energy_pn(Nx);

  vector<double> rho(Nx);
  vector<double> m1(Nx);
  vector<double> energy(Nx);

  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho_m[kx] = (S_x + kx)->rhoM;
    u1_m[kx] = (S_x + kx)->u1M;
    Tprt_m[kx] = (S_x + kx)->TprtM;

    rho_pn[kx] = Neff / dx * ((S_x + kx)->m0P - (S_x + kx)->m0N);
    m1_pn[kx] = Neff / dx * ((S_x + kx)->m11P - (S_x + kx)->m11N);
    energy_pn[kx] = .5 * Neff / dx * ((S_x + kx)->m2P - (S_x + kx)->m2N);

    rho[kx] = rho_m[kx];
  }

  // convert rho_m, u1_m, Tprt_m to momentum and energy
  uT2mE_x1v3(Nx, rho_m, u1_m, Tprt_m, m1, energy);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] += rho_pn[kx];
    m1[kx] += m1_pn[kx];
    energy[kx] += energy_pn[kx];
  }

  mE2uT_x1v3(Nx, rho, m1, energy, u1, Tprt);

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->rho = rho[kx];
    (S_x + kx)->u1 = u1[kx];
    (S_x + kx)->Tprt = Tprt[kx];
  }
}

void update_rhouT_F(NeParticleGroup* S_x, const NumericGridClass& grid) {
  int Nx = grid.Nx;
  double Neff_F = grid.Neff_F;

  vector<double> rho(Nx);
  vector<double> m1(Nx);
  vector<double> energy(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = Neff_F / grid.dx * (S_x + kx)->m0F;
    m1[kx] = Neff_F / grid.dx * (S_x + kx)->m11F;
    energy[kx] = .5 * Neff_F / grid.dx * (S_x + kx)->m2F;
  }

  mE2uT_x1v3(Nx, rho, m1, energy, u1, Tprt);

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->rhoF = rho[kx];
    (S_x + kx)->u1F = u1[kx];
    (S_x + kx)->TprtF = Tprt[kx];
  }
}

// ========================================================================

/**
  Compute the derivative of rho, u, T
  input: S_x, grid
  output: none
*/

void update_dx_rhouT_M(NeParticleGroup* S_x, const NumericGridClass& grid) {
  int Nx = grid.Nx;
  double dx = grid.dx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);
  vector<double> dx_rho(Nx);
  vector<double> dx_u1(Nx);
  vector<double> dx_Tprt(Nx);

  // the derivative of maxwellian moments
  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->rhoM;
    u1[kx] = (S_x + kx)->u1M;
    Tprt[kx] = (S_x + kx)->TprtM;
  }

  diff_1d_central(rho, Nx, grid.bdry_x, dx_rho);
  diff_1d_central(u1, Nx, grid.bdry_x, dx_u1);
  diff_1d_central(Tprt, Nx, grid.bdry_x, dx_Tprt);

  for (int kx = 0; kx < Nx; kx++) {
    dx_rho[kx] /= dx;
    dx_u1[kx] /= dx;
    dx_Tprt[kx] /= dx;
  }

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->dx_rhoM = dx_rho[kx];
    (S_x + kx)->dx_u1M = dx_u1[kx];
    (S_x + kx)->dx_TprtM = dx_Tprt[kx];
  }
}

// ========================================================================
/**
  the moments of \Delta t (v\cdot\nabla_x g)
  It's also the moments of \Delta t (v\cdot\nabla_x g + E\cdot\nabla_v g)
  input: S_x, grid
  output: drho, dm1, denergy
*/
void momentchange_g(NeParticleGroup* S_x, const NumericGridClass& grid,
                    vector<double>& drho, vector<double>& dm1,
                    vector<double>& denergy) {
  int Nx = grid.Nx;
  double Neff = grid.Neff;
  double dx = grid.dx;

  vector<double> rho_flux(Nx);
  vector<double> m1_flux(Nx);
  vector<double> energy_flux(Nx);

  vector<double> rho(Nx);
  vector<double> m1(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho_flux[kx] = Neff / dx * ((S_x + kx)->m11P - (S_x + kx)->m11N);
    m1_flux[kx] = Neff / dx * ((S_x + kx)->m21P - (S_x + kx)->m21N);
    energy_flux[kx] = .5 * Neff / dx * ((S_x + kx)->m31P - (S_x + kx)->m31N);

    rho[kx] = Neff / dx * ((S_x + kx)->m0P - (S_x + kx)->m0N);
    m1[kx] = Neff / dx * ((S_x + kx)->m11P - (S_x + kx)->m11N);
  }

  diff_1d_central(rho_flux, Nx, grid.bdry_x, drho);
  diff_1d_central(m1_flux, Nx, grid.bdry_x, dm1);
  diff_1d_central(energy_flux, Nx, grid.bdry_x, denergy);

  double dtdx = grid.dt / grid.dx;
  for (int kx = 0; kx < Nx; kx++) {
    drho[kx] *= dtdx;
    dm1[kx] *= dtdx;
    denergy[kx] *= dtdx;
  }

  for (int kx = 0; kx < Nx; kx++) {
    double elec = (S_x + kx)->elecfield;
    dm1[kx] -= grid.dt * rho[kx] * elec;
    denergy[kx] -= grid.dt * m1[kx] * elec;
  }
}

void momentchange_g_ver2(NeParticleGroup* S_x, const NumericGridClass& grid,
                         vector<double>& drho, vector<double>& dm1,
                         vector<double>& denergy) {
  int Nx = grid.Nx;
  double Neff = grid.Neff;
  double dx = grid.dx;

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->computemoments();

    double rho_old =
        Neff / dx *
        ((S_x + kx)->m0P_o - (S_x + kx)->m0N_o);  // inner product with 1
    double m1_old =
        Neff / dx *
        ((S_x + kx)->m11P_o - (S_x + kx)->m11N_o);  // inner product with v_1
    double energy_old =
        0.5 * Neff / dx *
        ((S_x + kx)->m2P_o - (S_x + kx)->m2N_o);  // inner product with .5*|v|^2

    double rho_new =
        Neff / dx *
        ((S_x + kx)->m0P - (S_x + kx)->m0N);  // inner product with 1
    double m1_new =
        Neff / dx *
        ((S_x + kx)->m11P - (S_x + kx)->m11N);  // inner product with v_1
    double energy_new =
        0.5 * Neff / dx *
        ((S_x + kx)->m2P - (S_x + kx)->m2N);  // inner product with .5*|v|^2

    drho[kx] = -(rho_new - rho_old);
    dm1[kx] = -(m1_new - m1_old);
    denergy[kx] = -(energy_new - energy_old);
  }
}

// ===============================================================
// compute the change in the macro part
// update S_x -> drho, dm1, denergy
// update S_x -> drho_g, dm1_g, denergy_g
void compute_change_in_macro(std::vector<NeParticleGroup>& S_x,
                             const NumericGridClass& grid) {
  int Nx = grid.Nx;

  vector<double> rho_m(Nx);
  vector<double> u1_m(Nx);
  vector<double> Tprt_m(Nx);

  vector<double> drho(Nx);
  vector<double> dm1(Nx);
  vector<double> denergy(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho_m[kx] = S_x[kx].rhoM;
    u1_m[kx] = S_x[kx].u1M;
    Tprt_m[kx] = S_x[kx].TprtM;
  }

  // change of macro part due to v\cdot\nabla_x M

  Euler_kinetic_x1(rho_m, u1_m, Tprt_m, Nx, grid.dx, grid.dt, grid.bdry_x, drho,
                   dm1, denergy);

  if (FLAG_SAVEFLUX) {
    save_macro<double>(drho, "drho_euler");
    save_macro<double>(dm1, "dm1_euler");
    save_macro<double>(denergy, "denergy_euler");
  }

  // for (int kx = 0; kx<Nx; kx ++) cout << "dm1 = " <<  dm1[kx] << endl;

  for (int kx = 0; kx < Nx; kx++) {
    S_x[kx].drho = drho[kx];
    S_x[kx].dm1 = dm1[kx];
    S_x[kx].denergy = denergy[kx];
  }

  // cout << "denergy =  " << denergy[1];
  // save_macro<double>( drho, "drho_a");
  // save_macro<double>( denergy, "denergy_a");

  // change of macro part due to v\cdot\nabla_x g + E\cdot\nabla_v g

  momentchange_g(&S_x[0], grid, drho, dm1, denergy);
  // momentchange_g_ver2(S_x, grid, drho, dm1, denergy);

  if (FLAG_SAVEFLUX) {
    save_macro<double>(drho, "drho_g");
    save_macro<double>(dm1, "dm1_g");
    save_macro<double>(denergy, "denergy_g");
  }

  for (int kx = 0; kx < Nx; kx++) {
    S_x[kx].drho_g = drho[kx];
    S_x[kx].dm1_g = dm1[kx];
    S_x[kx].denergy_g = denergy[kx];

    S_x[kx].drho += drho[kx];
    S_x[kx].dm1 += dm1[kx];
    S_x[kx].denergy += denergy[kx];
  }

  // cout << ", " << denergy[1];
  // save_macro<double>( drho, "drho_b");
  // save_macro<double>( denergy, "denergy_b");

  // change of macro part due to E\cdot\nabla_v (M)

  for (int kx = 0; kx < Nx; kx++) {
    S_x[kx].dm1 -= grid.dt * S_x[kx].rho_o * S_x[kx].elecfield;
    S_x[kx].denergy -=
        grid.dt * S_x[kx].rho_o * S_x[kx].u1_o * S_x[kx].elecfield;
    // (S_x + kx) -> dm1 -= grid.dt * (S_x+kx)-> rhoM * (S_x+kx)->elecfield;
    // (S_x + kx) -> denergy -= grid.dt * (S_x+kx)-> rhoM * (S_x+kx)-> u1M *
    // (S_x+kx)->elecfield;
  }

  // cout << ", " << -grid.dt * (S_x+1)-> rho * (S_x+1)-> u1 *
  // (S_x+1)->elecfield << ". " <<  (S_x+1)->denergy << endl;

  if (FLAG_SAVEFLUX) {
    for (int kx = 0; kx < Nx; kx++) {
      dm1[kx] = grid.dt * S_x[kx].rho * S_x[kx].elecfield;
      denergy[kx] = grid.dt * S_x[kx].rho * S_x[kx].u1 * S_x[kx].elecfield;
    }
    save_macro<double>(dm1, "dm1_elecfield");
    save_macro<double>(denergy, "denergy_elecfield");

    for (int kx = 0; kx < Nx; kx++) {
      drho[kx] = S_x[kx].drho;
      dm1[kx] = S_x[kx].dm1;
      denergy[kx] = S_x[kx].denergy;
    }
    save_macro<double>(drho, "drho_all");
    save_macro<double>(dm1, "dm1_all");
    save_macro<double>(denergy, "denergy_all");

    // temperary vector names, don't worry
    for (int kx = 0; kx < Nx; kx++) {
      drho[kx] = S_x[kx].rho_o;
      dm1[kx] = S_x[kx].u1_o;
      rho_m[kx] = S_x[kx].rho;
      u1_m[kx] = S_x[kx].u1_o;
    }
    save_macro<double>(drho, "oldrho");
    save_macro<double>(dm1, "oldu1");
    save_macro<double>(rho_m, "newrho");
    save_macro<double>(u1_m, "newu1");
  }
}

// update the maxwellian part: S_x -> rhoM, u1M, TprtM

void update_maxwellian(std::vector<NeParticleGroup>& S_x,
                       const NumericGridClass& grid) {
  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);
  vector<double> m1(Nx);
  vector<double> energy(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = S_x[kx].rhoM;
    u1[kx] = S_x[kx].u1M;
    Tprt[kx] = S_x[kx].TprtM;
  }

  // convert rho, u1, Tprt to momentum and energy
  uT2mE_x1v3(Nx, rho, u1, Tprt, m1, energy);

  /*
  save_macro<double>(rho, "rhot");
  save_macro<double>(m1, "m1t");
  save_macro<double>(energy, "energyt");
  save_macro<double>(u1, "u1t");
  save_macro<double>(Tprt, "Tprtt");
  */

  // compute rho, m1, energy at time n+1

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] -= S_x[kx].drho;
    m1[kx] -= S_x[kx].dm1;
    energy[kx] -= S_x[kx].denergy;
  }

  // update the maxwellian part in particle group S_x

  mE2uT_x1v3(Nx, rho, m1, energy, u1, Tprt);

  /*
  save_macro<double>(rho, "rhot2");
  save_macro<double>(m1, "m1t2");
  save_macro<double>(energy, "energyt2");
  save_macro<double>(u1, "u1t2");
  save_macro<double>(Tprt, "Tprtt2");
  */

  for (int kx = 0; kx < Nx; kx++) {
    S_x[kx].rhoM = rho[kx];
    S_x[kx].u1M = u1[kx];
    S_x[kx].TprtM = Tprt[kx];
    if (S_x[kx].TprtM < 0) {
      cout << "NEGATIVE TEMP " << kx << ' ' << rho[kx] << ' ' << u1[kx] << ' '
           << Tprt[kx] << ' ' << energy[kx] << ' ' << S_x[kx].denergy << endl;
      cout << S_x[kx].TprtM << endl;
      cout << kx << endl;
      exit(0);
    }
  }
}

// =================================================================
// return the electric energy
double compute_elec_energy(NeParticleGroup* S_x, const NumericGridClass& grid) {
  double E0 = 0.;

  for (int kx = 0; kx < grid.Nx; kx++) {
    double Ekx = (S_x + kx)->elecfield;
    E0 += Ekx * Ekx;
  }

  E0 *= grid.dx;

  return E0;
}

double compute_elec_energy_F(NeParticleGroup* S_x,
                             const NumericGridClass& grid) {
  double E0 = 0.;

  for (int kx = 0; kx < grid.Nx; kx++) {
    double Ekx = (S_x + kx)->elecfield_F;
    E0 += Ekx * Ekx;
  }

  E0 *= grid.dx;

  return E0;
}

double compute_total_energy(NeParticleGroup* S_x,
                            const NumericGridClass& grid) {
  double E0 = 0.;

  for (int kx = 0; kx < grid.Nx; kx++) {
    double Ekx = (S_x + kx)->elecfield;
    E0 += Ekx * Ekx;
  }

  E0 *= grid.dx;

  for (int kx = 0; kx < grid.Nx; kx++) {
    E0 += .5 * (S_x + kx)->rhoM *
          ((S_x + kx)->u1M * (S_x + kx)->u1M + 3 * (S_x + kx)->TprtM) * grid.dx;
    E0 += .5 * grid.Neff * ((S_x + kx)->m2P - (S_x + kx)->m2N);
  }
  return E0;
}

double compute_total_energy_F(NeParticleGroup* S_x,
                              const NumericGridClass& grid) {
  double E0 = 0.;

  for (int kx = 0; kx < grid.Nx; kx++) {
    double Ekx = (S_x + kx)->elecfield_F;
    E0 += Ekx * Ekx;
  }

  E0 *= grid.dx;

  for (int kx = 0; kx < grid.Nx; kx++) {
    E0 += .5 * grid.Neff_F * (S_x + kx)->m2F;
  }
  return E0;
}

void compute_change_in_macro_onlykineitc(NeParticleGroup* S_x,
                                         const NumericGridClass& grid) {
  int Nx = grid.Nx;

  vector<double> rho_m(Nx);
  vector<double> u1_m(Nx);
  vector<double> Tprt_m(Nx);
  vector<double> drho(Nx);
  vector<double> dm1(Nx);
  vector<double> denergy(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho_m[kx] = (S_x + kx)->rhoM;
    u1_m[kx] = (S_x + kx)->u1M;
    Tprt_m[kx] = (S_x + kx)->TprtM;
  }

  // change of macro part due to v\cdot\nabla_x M

  Euler_kinetic_x1(rho_m, u1_m, Tprt_m, Nx, grid.dx, grid.dt, grid.bdry_x, drho,
                   dm1, denergy);

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->drho = drho[kx];
    (S_x + kx)->dm1 = dm1[kx];
    (S_x + kx)->denergy = denergy[kx];
    (S_x + kx)->drho_g = 0.;
    (S_x + kx)->dm1_g = 0.;
    (S_x + kx)->denergy_g = 0.;
  }

  // cout << ", " << -grid.dt * (S_x+1)-> rho * (S_x+1)-> u1 *
  // (S_x+1)->elecfield << ". " <<  (S_x+1)->denergy << endl;
}

double myerf(double x) {
  // constants
  double a1 = 0.254829592;
  double a2 = -0.284496736;
  double a3 = 1.421413741;
  double a4 = -1.453152027;
  double a5 = 1.061405429;
  double p = 0.3275911;

  // Save the sign of x
  int sign = 1;
  if (x < 0) sign = -1;
  x = abs(x);

  // A&S formula 7.1.26
  double t = 1.0 / (1.0 + p * x);
  double y =
      1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

  return sign * y;
}
}  // namespace coulomb