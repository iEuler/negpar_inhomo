#include "Moments.h"

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationState.h"

#include "MacroOutput.h"
#include "Numerics.h"

#include <stdexcept>

namespace coulomb {

void MomentOperations::update_macro(std::vector<NeParticleGroup> &groups,
                                    const NumericGridClass &grid) {
  for (auto &group : groups)
    group.computemoments();
  MomentOperations::update_primitive(groups, grid);
  MomentOperations::update_maxwellian_derivatives(groups, grid);
}

void MomentOperations::compute_macro_change(
    std::vector<NeParticleGroup> &groups, const NumericGridClass &grid,
    SimulationState &state) {
  const int size = grid.Nx;
  std::vector<double> density(size), velocity(size), temperature(size);
  std::vector<double> density_change(size), momentum_change(size),
      energy_change(size);
  for (int cell = 0; cell < size; ++cell) {
    density[cell] = groups[cell].rhoM;
    velocity[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }

  Numerics::advance_kinetic_euler(density, velocity, temperature, size, grid.dx,
                                  grid.dt, grid.bdry_x, density_change,
                                  momentum_change, energy_change);
  if (state.saveFlux) {
    MacroOutput::save_macro(density_change, "drho_euler", state);
    MacroOutput::save_macro(momentum_change, "dm1_euler", state);
    MacroOutput::save_macro(energy_change, "denergy_euler", state);
  }
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].drho = density_change[cell];
    groups[cell].dm1 = momentum_change[cell];
    groups[cell].denergy = energy_change[cell];
  }

  std::tie(density_change, momentum_change, energy_change) =
      MomentOperations::moment_change(groups.data(), grid);
  if (state.saveFlux) {
    MacroOutput::save_macro(density_change, "drho_g", state);
    MacroOutput::save_macro(momentum_change, "dm1_g", state);
    MacroOutput::save_macro(energy_change, "denergy_g", state);
  }
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].drho_g = density_change[cell];
    groups[cell].dm1_g = momentum_change[cell];
    groups[cell].denergy_g = energy_change[cell];
    groups[cell].drho += density_change[cell];
    groups[cell].dm1 += momentum_change[cell];
    groups[cell].denergy += energy_change[cell];
    groups[cell].dm1 -= grid.dt * groups[cell].rho_o * groups[cell].elecfield;
    groups[cell].denergy -= grid.dt * groups[cell].rho_o * groups[cell].u1_o *
                            groups[cell].elecfield;
  }
}

void MomentOperations::primitive_to_conserved(const double &rho,
                                              const double *u,
                                              const double &temperature,
                                              double *momentum,
                                              double &energy) {
  if (u == nullptr || momentum == nullptr)
    throw std::invalid_argument("MomentOperations::primitive_to_conserved "
                                "velocity or momentum output is null");
  for (int component = 0; component < 3; ++component)
    momentum[component] = rho * u[component];

  double velocity_squared = 0.0;
  for (int component = 0; component < 3; ++component)
    velocity_squared += u[component] * u[component];
  energy = 0.5 * rho * (velocity_squared + 3.0 * temperature);
}

void MomentOperations::conserved_to_primitive(const double &rho,
                                              const double *momentum,
                                              const double &energy, double *u,
                                              double &temperature) {
  if (momentum == nullptr || u == nullptr)
    throw std::invalid_argument("MomentOperations::conserved_to_primitive "
                                "momentum or velocity output is null");
  for (int component = 0; component < 3; ++component)
    u[component] = momentum[component] / rho;

  double velocity_squared = 0.0;
  for (int component = 0; component < 3; ++component)
    velocity_squared += u[component] * u[component];
  temperature = (energy * 2.0 / rho - velocity_squared) / 3.0;
}

void MomentOperations::primitive_to_conserved(
    int grid_size, const std::vector<double> &rho,
    const std::vector<double> &u1, const std::vector<double> &temperature,
    std::vector<double> &m1, std::vector<double> &energy) {
  if (grid_size < 0 || static_cast<std::size_t>(grid_size) != rho.size() ||
      u1.size() != rho.size() || temperature.size() != rho.size() ||
      m1.size() != rho.size() || energy.size() != rho.size())
    throw std::invalid_argument(
        "MomentOperations::primitive_to_conserved input size mismatch");
  double velocity[3] = {0.0, 0.0, 0.0};
  double momentum[3] = {0.0, 0.0, 0.0};
  double energy_at_cell = 0.0;
  for (int cell = 0; cell < grid_size; ++cell) {
    velocity[0] = u1[cell];
    MomentOperations::primitive_to_conserved(
        rho[cell], velocity, temperature[cell], momentum, energy_at_cell);
    m1[cell] = momentum[0];
    energy[cell] = energy_at_cell;
  }
}

void MomentOperations::conserved_to_primitive(
    int grid_size, const std::vector<double> &rho,
    const std::vector<double> &m1, const std::vector<double> &energy,
    std::vector<double> &u1, std::vector<double> &temperature) {
  if (grid_size < 0 || static_cast<std::size_t>(grid_size) != rho.size() ||
      m1.size() != rho.size() || energy.size() != rho.size() ||
      u1.size() != rho.size() || temperature.size() != rho.size())
    throw std::invalid_argument(
        "MomentOperations::conserved_to_primitive input size mismatch");
  double velocity[3] = {0.0, 0.0, 0.0};
  double momentum[3] = {0.0, 0.0, 0.0};
  double temperature_at_cell = 0.0;
  for (int cell = 0; cell < grid_size; ++cell) {
    momentum[0] = m1[cell];
    MomentOperations::conserved_to_primitive(rho[cell], momentum, energy[cell],
                                             velocity, temperature_at_cell);
    u1[cell] = velocity[0];
    temperature[cell] = temperature_at_cell;
  }
}

void MomentOperations::particle_to_conserved(const ParticleGroup &group,
                                             double effective_particles,
                                             double &rho, double *momentum,
                                             double &energy) {
  if (momentum == nullptr)
    throw std::invalid_argument(
        "MomentOperations::particle_to_conserved momentum output is null");
  rho = group.moments.m0 * effective_particles;
  momentum[0] = group.moments.m11 * effective_particles;
  momentum[1] = group.moments.m12 * effective_particles;
  momentum[2] = group.moments.m13 * effective_particles;
  energy = 0.5 * group.moments.m2 * effective_particles;
}

void MomentOperations::compute_primitive(
    int grid_size, const std::vector<ParticleGroup> &groups,
    double effective_particles, std::vector<double> &rho,
    std::vector<double> &u1, std::vector<double> &u2, std::vector<double> &u3,
    std::vector<double> &temperature) {
  if (grid_size < 0 || static_cast<std::size_t>(grid_size) != groups.size() ||
      rho.size() != groups.size() || u1.size() != groups.size() ||
      u2.size() != groups.size() || u3.size() != groups.size() ||
      temperature.size() != groups.size())
    throw std::invalid_argument(
        "MomentOperations::compute_primitive input size mismatch");
  for (int cell = 0; cell < grid_size; ++cell) {
    double density = 0.0;
    double momentum[3] = {0.0, 0.0, 0.0};
    double energy = 0.0;
    double velocity[3] = {0.0, 0.0, 0.0};
    double cell_temperature = 0.0;
    MomentOperations::particle_to_conserved(groups[cell], effective_particles,
                                            density, momentum, energy);
    MomentOperations::conserved_to_primitive(density, momentum, energy,
                                             velocity, cell_temperature);
    rho[cell] = density;
    u1[cell] = velocity[0];
    u2[cell] = velocity[1];
    u3[cell] = velocity[2];
    temperature[cell] = cell_temperature;
  }
}

void MomentOperations::update_primitive(std::vector<NeParticleGroup> &groups,
                                        const NumericGridClass &grid) {
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
    rho_pn[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m0 - groups[cell].negative_moments.m0);
    momentum_pn[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m11 - groups[cell].negative_moments.m11);
    energy_pn[cell] =
        0.5 * effective_particles / dx *
        (groups[cell].positive_moments.m2 - groups[cell].negative_moments.m2);
    rho[cell] = rho_m[cell];
  }
  MomentOperations::primitive_to_conserved(size, rho_m, u1_m, temperature_m,
                                           momentum, energy);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] += rho_pn[cell];
    momentum[cell] += momentum_pn[cell];
    energy[cell] += energy_pn[cell];
  }
  MomentOperations::conserved_to_primitive(size, rho, momentum, energy, u1,
                                           temperature);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].rho = rho[cell];
    groups[cell].u1 = u1[cell];
    groups[cell].Tprt = temperature[cell];
  }
}

void MomentOperations::update_full_primitive(
    std::vector<NeParticleGroup> &groups, const NumericGridClass &grid) {
  const int size = grid.Nx;
  const double effective_particles = grid.Neff_F;
  std::vector<double> rho(size), momentum(size), energy(size);
  std::vector<double> u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = effective_particles / grid.dx * groups[cell].full_moments.m0;
    momentum[cell] =
        effective_particles / grid.dx * groups[cell].full_moments.m11;
    energy[cell] =
        0.5 * effective_particles / grid.dx * groups[cell].full_moments.m2;
  }
  MomentOperations::conserved_to_primitive(size, rho, momentum, energy, u1,
                                           temperature);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].rhoF = rho[cell];
    groups[cell].u1F = u1[cell];
    groups[cell].TprtF = temperature[cell];
  }
}

void MomentOperations::update_maxwellian_derivatives(
    std::vector<NeParticleGroup> &groups, const NumericGridClass &grid) {
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoM;
    u1[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }
  auto dx_rho = Numerics::central_difference(rho, size, grid.bdry_x);
  auto dx_u1 = Numerics::central_difference(u1, size, grid.bdry_x);
  auto dx_temperature =
      Numerics::central_difference(temperature, size, grid.bdry_x);
  for (int cell = 0; cell < size; ++cell) {
    dx_rho[cell] /= grid.dx;
    dx_u1[cell] /= grid.dx;
    dx_temperature[cell] /= grid.dx;
    groups[cell].dx_rhoM = dx_rho[cell];
    groups[cell].dx_u1M = dx_u1[cell];
    groups[cell].dx_TprtM = dx_temperature[cell];
  }
}

void MomentOperations::update_maxwellian(std::vector<NeParticleGroup> &groups,
                                         const NumericGridClass &grid) {
  const int size = grid.Nx;
  std::vector<double> rho(size), velocity(size), temperature(size);
  std::vector<double> momentum(size), energy(size);

  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoM;
    velocity[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }

  MomentOperations::primitive_to_conserved(size, rho, velocity, temperature,
                                           momentum, energy);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] -= groups[cell].drho;
    momentum[cell] -= groups[cell].dm1;
    energy[cell] -= groups[cell].denergy;
  }

  MomentOperations::conserved_to_primitive(size, rho, momentum, energy,
                                           velocity, temperature);
  for (int cell = 0; cell < size; ++cell) {
    if (!(temperature[cell] >= 0.0))
      throw std::runtime_error(
          "Maxwellian update produced invalid temperature");
    groups[cell].rhoM = rho[cell];
    groups[cell].u1M = velocity[cell];
    groups[cell].TprtM = temperature[cell];
  }
}

void MomentOperations::compute_kinetic_macro_change(
    std::vector<NeParticleGroup> &groups, const NumericGridClass &grid) {
  const int size = grid.Nx;
  std::vector<double> density(size), velocity(size), temperature(size);
  std::vector<double> density_change(size), momentum_change(size),
      energy_change(size);
  for (int cell = 0; cell < size; ++cell) {
    density[cell] = groups[cell].rhoM;
    velocity[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }

  Numerics::advance_kinetic_euler(density, velocity, temperature, size, grid.dx,
                                  grid.dt, grid.bdry_x, density_change,
                                  momentum_change, energy_change);
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
MomentOperations::moment_change(const NeParticleGroup *groups,
                                const NumericGridClass &grid) {
  const int size = grid.Nx;
  if (size < 0)
    throw std::invalid_argument(
        "MomentOperations::moment_change grid size is negative");
  if (size > 0 && groups == nullptr)
    throw std::invalid_argument(
        "MomentOperations::moment_change groups pointer is null");
  const double effective_particles = grid.Neff;
  const double dx = grid.dx;
  std::vector<double> rho_flux(size), momentum_flux(size), energy_flux(size);
  std::vector<double> rho(size), momentum(size);
  for (int cell = 0; cell < size; ++cell) {
    rho_flux[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m11 - groups[cell].negative_moments.m11);
    momentum_flux[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m21 - groups[cell].negative_moments.m21);
    energy_flux[cell] =
        0.5 * effective_particles / dx *
        (groups[cell].positive_moments.m31 - groups[cell].negative_moments.m31);
    rho[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m0 - groups[cell].negative_moments.m0);
    momentum[cell] =
        effective_particles / dx *
        (groups[cell].positive_moments.m11 - groups[cell].negative_moments.m11);
  }
  auto drho = Numerics::central_difference(rho_flux, size, grid.bdry_x);
  auto dm1 = Numerics::central_difference(momentum_flux, size, grid.bdry_x);
  auto denergy = Numerics::central_difference(energy_flux, size, grid.bdry_x);
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

} // namespace coulomb
