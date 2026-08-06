#include "Moments.h"

#include "Numerics.h"

#include <stdexcept>

namespace coulomb {

void uT2mE(const double& rho, double* u, const double& temperature,
           double* momentum, double& energy) {
  for (int component = 0; component < 3; ++component)
    momentum[component] = rho * u[component];

  double velocity_squared = 0.0;
  for (int component = 0; component < 3; ++component)
    velocity_squared += u[component] * u[component];
  energy = 0.5 * rho * (velocity_squared + 3.0 * temperature);
}

void mE2uT(const double& rho, double* momentum, const double& energy,
           double* u, double& temperature) {
  for (int component = 0; component < 3; ++component)
    u[component] = momentum[component] / rho;

  double velocity_squared = 0.0;
  for (int component = 0; component < 3; ++component)
    velocity_squared += u[component] * u[component];
  temperature = (energy * 2.0 / rho - velocity_squared) / 3.0;
}

void uT2mE_x1v3(int grid_size, const std::vector<double>& rho,
                const std::vector<double>& u1,
                const std::vector<double>& temperature,
                std::vector<double>& m1, std::vector<double>& energy) {
  double velocity[3] = {0.0, 0.0, 0.0};
  double momentum[3] = {0.0, 0.0, 0.0};
  double energy_at_cell = 0.0;
  for (int cell = 0; cell < grid_size; ++cell) {
    velocity[0] = u1[cell];
    uT2mE(rho[cell], velocity, temperature[cell], momentum, energy_at_cell);
    m1[cell] = momentum[0];
    energy[cell] = energy_at_cell;
  }
}

void mE2uT_x1v3(int grid_size, const std::vector<double>& rho,
                const std::vector<double>& m1,
                const std::vector<double>& energy,
                std::vector<double>& u1,
                std::vector<double>& temperature) {
  double velocity[3] = {0.0, 0.0, 0.0};
  double momentum[3] = {0.0, 0.0, 0.0};
  double temperature_at_cell = 0.0;
  for (int cell = 0; cell < grid_size; ++cell) {
    momentum[0] = m1[cell];
    mE2uT(rho[cell], momentum, energy[cell], velocity, temperature_at_cell);
    u1[cell] = velocity[0];
    temperature[cell] = temperature_at_cell;
  }
}

void particle2rhomE(const ParticleGroup& group, double effective_particles,
                    double& rho, double* momentum, double& energy) {
  rho = group.m0 * effective_particles;
  momentum[0] = group.m11 * effective_particles;
  momentum[1] = group.m12 * effective_particles;
  momentum[2] = group.m13 * effective_particles;
  energy = 0.5 * group.m2 * effective_particles;
}

void compute_rhouT(int grid_size, const std::vector<ParticleGroup>& groups,
                   double effective_particles, std::vector<double>& rho,
                   std::vector<double>& u1, std::vector<double>& u2,
                   std::vector<double>& u3,
                   std::vector<double>& temperature) {
  for (int cell = 0; cell < grid_size; ++cell) {
    double density = 0.0;
    double momentum[3] = {0.0, 0.0, 0.0};
    double energy = 0.0;
    double velocity[3] = {0.0, 0.0, 0.0};
    double cell_temperature = 0.0;
    particle2rhomE(groups[cell], effective_particles, density, momentum,
                   energy);
    mE2uT(density, momentum, energy, velocity, cell_temperature);
    rho[cell] = density;
    u1[cell] = velocity[0];
    u2[cell] = velocity[1];
    u3[cell] = velocity[2];
    temperature[cell] = cell_temperature;
  }
}

void update_rhouT(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& grid) {
  const int size = grid.Nx;
  const double effective_particles = grid.Neff;
  const double dx = grid.dx;
  std::vector<double> rho_m(size), u1_m(size), temperature_m(size);
  std::vector<double> rho_pn(size), momentum_pn(size), energy_pn(size);
  std::vector<double> rho(size), momentum(size), energy(size);
  std::vector<double> u1(size), temperature(size);

  for (int cell = 0; cell < size; ++cell) {
    rho_m[cell] = groups[cell].rhoM;
    u1_m[cell] = groups[cell].u1M;
    temperature_m[cell] = groups[cell].TprtM;
    rho_pn[cell] = effective_particles / dx *
                   (groups[cell].m0P - groups[cell].m0N);
    momentum_pn[cell] = effective_particles / dx *
                        (groups[cell].m11P - groups[cell].m11N);
    energy_pn[cell] = 0.5 * effective_particles / dx *
                      (groups[cell].m2P - groups[cell].m2N);
    rho[cell] = rho_m[cell];
  }
  uT2mE_x1v3(size, rho_m, u1_m, temperature_m, momentum, energy);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] += rho_pn[cell];
    momentum[cell] += momentum_pn[cell];
    energy[cell] += energy_pn[cell];
  }
  mE2uT_x1v3(size, rho, momentum, energy, u1, temperature);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].rho = rho[cell];
    groups[cell].u1 = u1[cell];
    groups[cell].Tprt = temperature[cell];
  }
}

void update_rhouT_F(std::vector<NeParticleGroup>& groups,
                    const NumericGridClass& grid) {
  const int size = grid.Nx;
  const double effective_particles = grid.Neff_F;
  std::vector<double> rho(size), momentum(size), energy(size);
  std::vector<double> u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = effective_particles / grid.dx * groups[cell].m0F;
    momentum[cell] = effective_particles / grid.dx * groups[cell].m11F;
    energy[cell] = 0.5 * effective_particles / grid.dx * groups[cell].m2F;
  }
  mE2uT_x1v3(size, rho, momentum, energy, u1, temperature);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].rhoF = rho[cell];
    groups[cell].u1F = u1[cell];
    groups[cell].TprtF = temperature[cell];
  }
}

void update_dx_rhouT_M(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid) {
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoM;
    u1[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }
  auto dx_rho = diff_1d_central(rho, size, grid.bdry_x);
  auto dx_u1 = diff_1d_central(u1, size, grid.bdry_x);
  auto dx_temperature = diff_1d_central(temperature, size, grid.bdry_x);
  for (int cell = 0; cell < size; ++cell) {
    dx_rho[cell] /= grid.dx;
    dx_u1[cell] /= grid.dx;
    dx_temperature[cell] /= grid.dx;
    groups[cell].dx_rhoM = dx_rho[cell];
    groups[cell].dx_u1M = dx_u1[cell];
    groups[cell].dx_TprtM = dx_temperature[cell];
  }
}

void update_maxwellian(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid) {
  const int size = grid.Nx;
  std::vector<double> rho(size), velocity(size), temperature(size);
  std::vector<double> momentum(size), energy(size);

  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoM;
    velocity[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }

  uT2mE_x1v3(size, rho, velocity, temperature, momentum, energy);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] -= groups[cell].drho;
    momentum[cell] -= groups[cell].dm1;
    energy[cell] -= groups[cell].denergy;
  }

  mE2uT_x1v3(size, rho, momentum, energy, velocity, temperature);
  for (int cell = 0; cell < size; ++cell) {
    if (!(temperature[cell] >= 0.0))
      throw std::runtime_error("Maxwellian update produced invalid temperature");
    groups[cell].rhoM = rho[cell];
    groups[cell].u1M = velocity[cell];
    groups[cell].TprtM = temperature[cell];
  }
}

void compute_change_in_macro_onlykineitc(
    std::vector<NeParticleGroup>& groups, const NumericGridClass& grid) {
  const int size = grid.Nx;
  std::vector<double> density(size), velocity(size), temperature(size);
  std::vector<double> density_change(size), momentum_change(size),
      energy_change(size);
  for (int cell = 0; cell < size; ++cell) {
    density[cell] = groups[cell].rhoM;
    velocity[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }

  Euler_kinetic_x1(density, velocity, temperature, size, grid.dx, grid.dt,
                   grid.bdry_x, density_change, momentum_change,
                   energy_change);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].drho = density_change[cell];
    groups[cell].dm1 = momentum_change[cell];
    groups[cell].denergy = energy_change[cell];
    groups[cell].drho_g = 0.0;
    groups[cell].dm1_g = 0.0;
    groups[cell].denergy_g = 0.0;
  }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
momentchange_g(NeParticleGroup* groups, const NumericGridClass& grid) {
  const int size = grid.Nx;
  const double effective_particles = grid.Neff;
  const double dx = grid.dx;
  std::vector<double> rho_flux(size), momentum_flux(size), energy_flux(size);
  std::vector<double> rho(size), momentum(size);
  for (int cell = 0; cell < size; ++cell) {
    rho_flux[cell] = effective_particles / dx *
                     (groups[cell].m11P - groups[cell].m11N);
    momentum_flux[cell] = effective_particles / dx *
                          (groups[cell].m21P - groups[cell].m21N);
    energy_flux[cell] = 0.5 * effective_particles / dx *
                        (groups[cell].m31P - groups[cell].m31N);
    rho[cell] = effective_particles / dx *
                (groups[cell].m0P - groups[cell].m0N);
    momentum[cell] = effective_particles / dx *
                     (groups[cell].m11P - groups[cell].m11N);
  }
  auto drho = diff_1d_central(rho_flux, size, grid.bdry_x);
  auto dm1 = diff_1d_central(momentum_flux, size, grid.bdry_x);
  auto denergy = diff_1d_central(energy_flux, size, grid.bdry_x);
  const double time_space_ratio = grid.dt / grid.dx;
  for (int cell = 0; cell < size; ++cell) {
    drho[cell] *= time_space_ratio;
    dm1[cell] *= time_space_ratio;
    denergy[cell] *= time_space_ratio;
    dm1[cell] -= grid.dt * rho[cell] * groups[cell].elecfield;
    denergy[cell] -= grid.dt * momentum[cell] * groups[cell].elecfield;
  }
  return {std::move(drho), std::move(dm1), std::move(denergy)};
}

void momentchange_g_ver2(NeParticleGroup* groups, const NumericGridClass& grid,
                         std::vector<double>& drho,
                         std::vector<double>& dm1,
                         std::vector<double>& denergy) {
  const int size = grid.Nx;
  const double effective_particles = grid.Neff;
  const double dx = grid.dx;
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].computemoments();
    const double rho_old = effective_particles / dx *
                           (groups[cell].m0P_o - groups[cell].m0N_o);
    const double momentum_old = effective_particles / dx *
                                (groups[cell].m11P_o - groups[cell].m11N_o);
    const double energy_old = 0.5 * effective_particles / dx *
                              (groups[cell].m2P_o - groups[cell].m2N_o);
    const double rho_new = effective_particles / dx *
                           (groups[cell].m0P - groups[cell].m0N);
    const double momentum_new = effective_particles / dx *
                                (groups[cell].m11P - groups[cell].m11N);
    const double energy_new = 0.5 * effective_particles / dx *
                              (groups[cell].m2P - groups[cell].m2N);
    drho[cell] = -(rho_new - rho_old);
    dm1[cell] = -(momentum_new - momentum_old);
    denergy[cell] = -(energy_new - energy_old);
  }
}

}  // namespace coulomb
