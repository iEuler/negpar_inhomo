#include "Output.h"

#include "Moments.h"
#include "utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "_global_variables.h"

namespace coulomb {

std::string output_path(const std::string& filename,
                        const SimulationState& state) {
  return (std::filesystem::path(state.outputDirectory) / filename).string();
}

std::string int2str(int value, int digits) {
  unsigned int unsigned_value = value;
  if (value < 0) unsigned_value = static_cast<unsigned int>(-value);

  std::string result;
  while (digits-- > 0) {
    result += static_cast<char>('0' + unsigned_value % 10);
    unsigned_value /= 10;
  }
  if (value < 0) result += '-';
  result = std::string(result.rbegin(), result.rend());
  return '_' + result;
}

void save_complex(int size, std::complex<double>* values,
                  const std::string& filename,
                  const SimulationState& state) {
  if (size < 0 || (size > 0 && values == nullptr))
    throw std::invalid_argument("save_complex received invalid input");

  const auto real_path = output_path(filename + "_r.txt", state);
  const auto imaginary_path = output_path(filename + "_i.txt", state);
  std::ofstream real_file(real_path);
  std::ofstream imaginary_file(imaginary_path);
  if (!real_file || !imaginary_file)
    throw std::runtime_error("Unable to open complex output files");
  for (int index = 0; index < size; ++index) {
    real_file << values[index].real() << '\n';
    imaginary_file << values[index].imag() << '\n';
  }
}

void save_2d(int rows, int columns,
             const std::vector<std::vector<double>>& values,
             const std::string& filename, const SimulationState& state) {
  if (rows < 0 || columns < 0 || static_cast<std::size_t>(rows) != values.size())
    throw std::invalid_argument("save_2d dimensions do not match rows");
  for (const auto& row : values)
    if (static_cast<std::size_t>(columns) != row.size())
      throw std::invalid_argument("save_2d dimensions do not match columns");

  const auto numbered_name =
      filename + (state.filenameWithNumber ? int2str(state.saveIndex) : "") +
      ".txt";
  const auto path = output_path(numbered_name, state);
  std::ofstream file(path);
  if (!file) throw std::runtime_error("Unable to open output file: " + path);
  for (const auto& row : values) {
    for (const auto value : row) file << value << ' ';
    file << '\n';
  }
}

void save_rhouT(int size, const std::vector<double>& rho,
                const std::vector<double>& u1, const std::vector<double>& u2,
                const std::vector<double>& u3,
                const std::vector<double>& temperature,
                const SimulationState& state) {
  if (size < 0 || static_cast<std::size_t>(size) != rho.size() ||
      u1.size() != rho.size() || u2.size() != rho.size() ||
      u3.size() != rho.size() || temperature.size() != rho.size())
    throw std::invalid_argument("save_rhouT input size mismatch");
  save_macro(rho, "rho", state);
  save_macro(u1, "u1", state);
  save_macro(u2, "u2", state);
  save_macro(u3, "u3", state);
  save_macro(temperature, "Tprt", state);
}

void save_rhouT(const std::vector<ParticleGroup>& groups,
                const NumericGridClass& grid, const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), u2(size), u3(size), temperature(size);
  compute_rhouT(size, groups, grid.Neff, rho, u1, u2, u3, temperature);
  save_rhouT(size, rho, u1, u2, u3, temperature, state);
}

void save_rhouT(std::vector<NeParticleGroup>& groups,
                const NumericGridClass& grid, const SimulationState& state) {
  update_rhouT(groups, grid);
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rho;
    u1[cell] = groups[cell].u1;
    temperature[cell] = groups[cell].Tprt;
  }
  save_macro(rho, "rho", state);
  save_macro(u1, "u1", state);
  save_macro(temperature, "Tprt", state);
}

void save_rhouT_F(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& grid,
                  const SimulationState& state) {
  update_rhouT_F(groups, grid);
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoF;
    u1[cell] = groups[cell].u1F;
    temperature[cell] = groups[cell].TprtF;
  }
  save_macro(rho, "rhoF", state);
  save_macro(u1, "u1F", state);
  save_macro(temperature, "TprtF", state);
}

void save_rhouT_maxwellian(const std::vector<NeParticleGroup>& groups,
                           const NumericGridClass& grid,
                           const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<double> rho(size), u1(size), temperature(size);
  for (int cell = 0; cell < size; ++cell) {
    rho[cell] = groups[cell].rhoM;
    u1[cell] = groups[cell].u1M;
    temperature[cell] = groups[cell].TprtM;
  }
  save_macro(rho, "rhoM", state);
  save_macro(u1, "u1M", state);
  save_macro(temperature, "TprtM", state);
}

void save_E(std::vector<NeParticleGroup>& groups,
            const NumericGridClass& grid, const SimulationState& state) {
  update_rhouT(groups, grid);
  const int size = grid.Nx;
  std::vector<double> electric_field(size), electric_field_full(size);
  for (int cell = 0; cell < size; ++cell) {
    electric_field[cell] = groups[cell].elecfield;
    electric_field_full[cell] = groups[cell].elecfield_F;
  }
  save_macro(electric_field, "elecfield", state);
  save_macro(electric_field_full, "elecfield_F", state);
}

void save_E(std::vector<ParticleGroup>& groups, const NumericGridClass& grid,
            const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<double> electric_field(size);
  for (int cell = 0; cell < size; ++cell)
    electric_field[cell] = groups[cell].elecfield;
  save_macro(electric_field, "E", state);
}

void save_NpNn(std::vector<NeParticleGroup>& groups,
               const NumericGridClass& grid, const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<int> positive(size), negative(size);
  for (int cell = 0; cell < size; ++cell) {
    positive[cell] = groups[cell].size(ParticleKind::Positive);
    negative[cell] = groups[cell].size(ParticleKind::Negative);
  }
  save_macro(positive, "Np", state);
  save_macro(negative, "Nn", state);
}

void save_d_rhouT(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& grid,
                  const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<double> density_change(size), momentum_change(size),
      energy_change(size);
  for (int cell = 0; cell < size; ++cell) {
    density_change[cell] = groups[cell].drho;
    momentum_change[cell] = groups[cell].dm1;
    energy_change[cell] = groups[cell].denergy;
  }
  save_macro(density_change, "drho", state);
  save_macro(momentum_change, "dm1", state);
  save_macro(energy_change, "denergy", state);
}

void save_dx_rhouT_M(std::vector<NeParticleGroup>& groups,
                     const NumericGridClass& grid,
                     const SimulationState& state) {
  const int size = grid.Nx;
  std::vector<double> density_gradient(size), velocity_gradient(size),
      temperature_gradient(size);
  for (int cell = 0; cell < size; ++cell) {
    density_gradient[cell] = groups[cell].dx_rhoM;
    velocity_gradient[cell] = groups[cell].dx_u1M;
    temperature_gradient[cell] = groups[cell].dx_TprtM;
  }
  save_macro(density_gradient, "dxrho", state);
  save_macro(velocity_gradient, "dxu1", state);
  save_macro(temperature_gradient, "dxTprt", state);
}

void save_m012_PN(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& /*grid*/,
                  const SimulationState& state) {
  const int size = static_cast<int>(groups.size());
  std::vector<double> zeroth(size), first(size), second(size);
  for (int cell = 0; cell < size; ++cell) {
    groups[cell].computemoments();
    zeroth[cell] = groups[cell].m0P - groups[cell].m0N;
    first[cell] = groups[cell].m11P - groups[cell].m11N;
    second[cell] = groups[cell].m2P - groups[cell].m2N;
  }
  save_macro(zeroth, "m0", state);
  save_macro(first, "m1", state);
  save_macro(second, "m2", state);
}

void save_grids(const NumericGridClass& grid, const SimulationState& state) {
  save_macro(grid.x, "x", state);
  save_macro(grid.vx, "v", state);
}

void save_dist(std::vector<ParticleGroup>& groups,
               const NumericGridClass& grid, const SimulationState& state) {
  const int size = grid.Nx;
  const int bins = grid.Nv;
  std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
  for (int cell = 0; cell < size; ++cell) {
    const auto& particles = groups[cell].list();
    std::vector<double> velocities;
    velocities.reserve(particles.size());
    for (const auto& particle : particles) velocities.push_back(particle.velocity(0));
    histinfo_fixbar(velocities, counts[cell], grid.vmin, grid.vmax);
  }

  std::vector<std::vector<double>> distribution(
      size, std::vector<double>(bins));
  for (int cell = 0; cell < size; ++cell)
    for (int bin = 0; bin < bins; ++bin)
      distribution[cell][bin] = counts[cell][bin] * grid.Neff;
  save_2d(size, bins, distribution, "dist", state);
}

void save_dist(std::vector<NeParticleGroup>& groups,
               const NumericGridClass& grid, char particle_type,
               const SimulationState& state) {
  const int size = grid.Nx;
  const int bins = grid.Nv;
  std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
  for (int cell = 0; cell < size; ++cell) {
    const auto& particles = groups[cell].list(particle_type);
    std::vector<double> velocities;
    velocities.reserve(particles.size());
    for (const auto& particle : particles) velocities.push_back(particle.velocity(0));
    histinfo_fixbar(velocities, counts[cell], grid.vmin, grid.vmax);
  }

  const double effective_particles =
      particle_type == 'f' ? grid.Neff_F : grid.Neff;
  const double coefficient = effective_particles / grid.dx / grid.dv;
  std::vector<std::vector<double>> distribution(
      size, std::vector<double>(bins));
  for (int cell = 0; cell < size; ++cell)
    for (int bin = 0; bin < bins; ++bin)
      distribution[cell][bin] = counts[cell][bin] * coefficient;

  const std::string name = particle_type == 'p'
                               ? "distp"
                               : (particle_type == 'n' ? "distn" : "distf");
  save_2d(size, bins, distribution, name, state);
}

void save_dist(std::vector<NeParticleGroup>& groups,
               const NumericGridClass& grid, const SimulationState& state) {
  save_dist(groups, grid, 'p', state);
  save_dist(groups, grid, 'n', state);
  save_dist(groups, grid, 'f', state);
}

void save_particle1d1d(std::vector<ParticleGroup>& groups,
                       const NumericGridClass& grid,
                       const SimulationState& state) {
  const std::string name =
      "particle" + (state.filenameWithNumber ? int2str(state.saveIndex) : "") +
      ".txt";
  std::ofstream file(output_path(name, state));
  if (!file) throw std::runtime_error("Unable to open particle output file");
  for (int cell = 0; cell < grid.Nx; ++cell)
    for (const auto& particle : groups[cell].list())
      file << particle.position() << ' ' << particle.velocity(0) << ' ';
}

void save_particle1d1d(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid, char particle_type,
                       int quantity, const SimulationState& state) {
  if (quantity != 1 && quantity != 2)
    throw std::invalid_argument("particle output quantity must be 1 or 2");
  const std::string prefix = particle_type == 'p' ? "particleP" : "particleN";
  const std::string name =
      prefix + (state.filenameWithNumber ? int2str(state.saveIndex) : "") + ".txt";
  std::ofstream file(output_path(name, state));
  if (!file) throw std::runtime_error("Unable to open particle output file");
  for (int cell = 0; cell < grid.Nx; ++cell) {
    for (const auto& particle : groups[cell].list(particle_type)) {
      double value = particle.velocity(0);
      if (quantity == 2) {
        const auto velocity = particle.velocity();
        value = 0.0;
        for (const double component : velocity) value += component * component;
        value = std::sqrt(value);
      }
      file << value << '\n';
    }
  }
}

void save_particle1d1d(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid,
                       const SimulationState& state) {
  save_particle1d1d(groups, grid, 'p', 1, state);
  save_particle1d1d(groups, grid, 'n', 1, state);
}

void save_particleenergy(std::vector<NeParticleGroup>& groups,
                         const NumericGridClass& grid,
                         const SimulationState& state) {
  save_particle1d1d(groups, grid, 'p', 2, state);
  save_particle1d1d(groups, grid, 'n', 2, state);
}

void save_homo_rdist(const SimulationState& state) {
  constexpr int bin_count = 100;
  constexpr double rmax = 5.0;
  const double dr = rmax / bin_count;
  std::vector<double> radii(bin_count);
  for (int bin = 0; bin < bin_count; ++bin)
    radii[bin] = (bin + 0.5) * dr;
  save_macro(radii, "rdist", state);
}

void save_homo_rdist(int bin_count, const SimulationState& state) {
  if (bin_count <= 0)
    throw std::invalid_argument("homogeneous radial distribution needs bins");
  constexpr double rmax = 10.0;
  const double dr = rmax / bin_count;
  std::vector<double> radii(bin_count);
  for (int bin = 0; bin < bin_count; ++bin)
    radii[bin] = (bin + 0.5) * dr;
  save_macro(radii, "rdist", state);
}

void save_homo_dist(const NeParticleGroup& group, int bin_count,
                    int case_index, const SimulationState& state) {
  if (bin_count <= 0 || case_index < 0 || case_index > 2)
    throw std::invalid_argument("invalid homogeneous distribution parameters");
  constexpr double rmax = 10.0;
  const char* suffixes[] = {"", "_before", "_after"};
  const std::string suffix = suffixes[case_index];

  const auto save_species = [&](char particle_type, const char* name) {
    const auto& particles = group.list(particle_type);
    std::vector<double> speeds;
    speeds.reserve(particles.size());
    for (const auto& particle : particles) {
      const auto velocity = particle.velocity();
      double speed_squared = 0.0;
      for (const double component : velocity)
        speed_squared += component * component;
      speeds.push_back(std::sqrt(speed_squared));
    }
    std::vector<int> histogram(bin_count);
    histinfo_fixbar(speeds, histogram, 0.0, rmax);
    save_macro(histogram, std::string(name) + suffix, state);
  };

  save_species('p', "pdist");
  save_species('n', "ndist");
  save_species('f', "fdist");
}

void saveparameter(const ParaClass& parameters, const NumericGridClass& grid,
                   const SimulationState& state) {
  std::ofstream file(output_path("parameter.txt", state));
  if (!file) throw std::runtime_error("Unable to open parameter output file");
  file << std::setprecision(15)
       << "coeff_binarycoll " << parameters.coeff_binarycoll << '\n'
       << "resample_ratio " << parameters.resample_ratio << '\n'
       << "Npickup_neg " << parameters.Npickup_neg << '\n'
       << "Nfreq " << parameters.Nfreq << '\n'
       << "xmax " << grid.xmax << '\n'
       << "xmin " << grid.xmin << '\n'
       << "vmax " << grid.vmax << '\n'
       << "vmin " << grid.vmin << '\n'
       << "tmax " << grid.tmax << '\n'
       << "Nx " << grid.Nx << '\n'
       << "Nt " << grid.Nt << '\n'
       << "Nv " << grid.Nv << '\n'
       << "dx " << grid.dx << '\n'
       << "dv " << grid.dv << '\n'
       << "dt " << grid.dt << '\n'
       << "Neff " << grid.Neff << '\n'
       << "Neff_F " << grid.Neff_F << '\n'
       << "collision_kernel " << collision_name(parameters.collisionType)
       << '\n'
       << "lambda_Poisson " << parameters.lambda_Poisson << '\n'
       << "resample_spatial_ratio " << parameters.resample_spatial_ratio
       << '\n'
       << "sync_time_interval " << parameters.sync_time_interval << '\n'
       << "resample_sync_ratio " << parameters.resample_sync_ratio << '\n';

  std::ofstream second_file(output_path("parameter2.txt", state));
  if (!second_file)
    throw std::runtime_error("Unable to open secondary parameter output file");
  second_file << std::setprecision(15)
              << "method_binarycoll "
              << binary_collision_name(parameters.method_binarycoll) << '\n'
              << "bdry_x " << grid.bdry_x << '\n'
              << "bdry_v " << grid.bdry_v << '\n'
              << "method " << method_name(parameters.method) << '\n';
}

void save_initial(IniValClass& initial_data, const SimulationState& state) {
  std::ofstream file(output_path("parameter.txt", state), std::ios::app);
  if (!file) throw std::runtime_error("Unable to append parameter output");
  file << std::setprecision(15);
  if (initial_data.probname == "LandauDamping") {
    file << "LD_alpha " << initial_data.LD_alpha << '\n'
         << "totalmass " << initial_data.totalmass << '\n';
  } else if (initial_data.probname == "BumpOnTail") {
    file << "BOT_beta " << initial_data.BOT_beta << '\n'
         << "BOT_rho0 " << initial_data.BOT_rho0 << '\n'
         << "BOT_Tprt " << initial_data.BOT_Tprt << '\n'
         << "BOT_dTprt " << initial_data.BOT_dTprt << '\n'
         << "BOT_Tx " << initial_data.BOT_Tx << '\n'
         << "BOT_ub " << initial_data.BOT_ub << '\n';
  }
  std::ofstream second_file(output_path("parameter2.txt", state), std::ios::app);
  if (!second_file)
    throw std::runtime_error("Unable to append secondary parameter output");
  second_file << "problem_name " << initial_data.probname << '\n';
}

}  // namespace coulomb
