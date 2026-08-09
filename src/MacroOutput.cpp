#include "MacroOutput.h"

#include <fstream>
#include <stdexcept>
#include <vector>

#include "Grid.h"
#include "Moments.h"
#include "ParticleGroup.h"
#include "ParticleOutput.h"

namespace coulomb {
void save_macro_evolution(std::vector<NeParticleGroup>& groups,
                          const NumericGridClass& grid,
                          const SimulationState& state) {
  for (int cell = 0; cell < grid.Nx; ++cell) groups[cell].computemoments();
  update_macro(groups, grid);
  save_rhouT(groups, grid, state);
  save_rhouT_F(groups, grid, state);
  save_E(groups, grid, state);
  save_dist(groups, grid, state);
  save_rhouT_maxwellian(groups, grid, state);

  std::ofstream file(output_path("numOfDist.txt", state));
  if (!file)
    throw std::runtime_error("Unable to open distribution count output");
  file << state.saveIndex;
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
    zeroth[cell] = groups[cell].positive_moments.m0 - groups[cell].negative_moments.m0;
    first[cell] = groups[cell].positive_moments.m11 - groups[cell].negative_moments.m11;
    second[cell] = groups[cell].positive_moments.m2 - groups[cell].negative_moments.m2;
  }
  save_macro(zeroth, "m0", state);
  save_macro(first, "m1", state);
  save_macro(second, "m2", state);
}

}  // namespace coulomb
