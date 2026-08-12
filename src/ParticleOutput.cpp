#include "ParticleOutput.h"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Grid.h"
#include "Histogram.h"
#include "MacroOutput.h"
#include "OutputPaths.h"
#include "ParticleGroup.h"

namespace coulomb {
void ParticleOutput::save_distribution(std::vector<ParticleGroup> &groups,
                                       const NumericGridClass &grid,
                                       const SimulationState &state) {
  const int size = grid.Nx;
  const int bins = grid.Nv;
  std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
  for (int cell = 0; cell < size; ++cell) {
    const auto &particles = groups[cell].list();
    std::vector<double> velocities;
    velocities.reserve(particles.size());
    for (const auto &particle : particles)
      velocities.push_back(particle.velocity(0));
    Histogram{}.fixed_bins(velocities, counts[cell], grid.vmin, grid.vmax);
  }

  std::vector<std::vector<double>> distribution(size,
                                                std::vector<double>(bins));
  for (int cell = 0; cell < size; ++cell)
    for (int bin = 0; bin < bins; ++bin)
      distribution[cell][bin] = counts[cell][bin] * grid.Neff;
  MacroOutput{}.save_2d(size, bins, distribution, "dist", state);
}

void ParticleOutput::save_distribution(std::vector<NeParticleGroup> &groups,
                                       const NumericGridClass &grid,
                                       ParticleKind kind,
                                       const SimulationState &state) {
  const int size = grid.Nx;
  const int bins = grid.Nv;
  std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
  for (int cell = 0; cell < size; ++cell) {
    const auto &particles = groups[cell].list(kind);
    std::vector<double> velocities;
    velocities.reserve(particles.size());
    for (const auto &particle : particles)
      velocities.push_back(particle.velocity(0));
    Histogram{}.fixed_bins(velocities, counts[cell], grid.vmin, grid.vmax);
  }

  const double effective_particles =
      kind == ParticleKind::Full ? grid.Neff_F : grid.Neff;
  const double coefficient = effective_particles / grid.dx / grid.dv;
  std::vector<std::vector<double>> distribution(size,
                                                std::vector<double>(bins));
  for (int cell = 0; cell < size; ++cell)
    for (int bin = 0; bin < bins; ++bin)
      distribution[cell][bin] = counts[cell][bin] * coefficient;

  const std::string name =
      kind == ParticleKind::Positive
          ? "distp"
          : (kind == ParticleKind::Negative ? "distn" : "distf");
  MacroOutput{}.save_2d(size, bins, distribution, name, state);
}

void ParticleOutput::save_distribution(std::vector<NeParticleGroup> &groups,
                                       const NumericGridClass &grid,
                                       const SimulationState &state) {
  ParticleOutput::save_distribution(groups, grid, ParticleKind::Positive,
                                    state);
  ParticleOutput::save_distribution(groups, grid, ParticleKind::Negative,
                                    state);
  ParticleOutput::save_distribution(groups, grid, ParticleKind::Full, state);
}

void ParticleOutput::save_phase_space(std::vector<ParticleGroup> &groups,
                                      const NumericGridClass &grid,
                                      const SimulationState &state) {
  const std::string name =
      "particle" +
      (state.filenameWithNumber
           ? OutputPaths(state).format_index(state.saveIndex)
           : "") +
      ".txt";
  std::ofstream file(OutputPaths(state).resolve(name));
  if (!file)
    throw std::runtime_error("Unable to open particle output file");
  for (int cell = 0; cell < grid.Nx; ++cell)
    for (const auto &particle : groups[cell].list())
      file << particle.position() << ' ' << particle.velocity(0) << ' ';
}

void ParticleOutput::save_phase_space(std::vector<NeParticleGroup> &groups,
                                      const NumericGridClass &grid,
                                      ParticleKind kind, int quantity,
                                      const SimulationState &state) {
  if (quantity != 1 && quantity != 2)
    throw std::invalid_argument("particle output quantity must be 1 or 2");
  if (kind == ParticleKind::Full)
    throw std::invalid_argument(
        "particle output requires a signed particle kind");
  const std::string prefix =
      kind == ParticleKind::Positive ? "particleP" : "particleN";
  const std::string name =
      prefix +
      (state.filenameWithNumber
           ? OutputPaths(state).format_index(state.saveIndex)
           : "") +
      ".txt";
  std::ofstream file(OutputPaths(state).resolve(name));
  if (!file)
    throw std::runtime_error("Unable to open particle output file");
  for (int cell = 0; cell < grid.Nx; ++cell) {
    for (const auto &particle : groups[cell].list(kind)) {
      double value = particle.velocity(0);
      if (quantity == 2) {
        const auto velocity = particle.velocity();
        value = 0.0;
        for (const double component : velocity)
          value += component * component;
        value = std::sqrt(value);
      }
      file << value << '\n';
    }
  }
}

void ParticleOutput::save_phase_space(std::vector<NeParticleGroup> &groups,
                                      const NumericGridClass &grid,
                                      const SimulationState &state) {
  ParticleOutput::save_phase_space(groups, grid, ParticleKind::Positive, 1,
                                   state);
  ParticleOutput::save_phase_space(groups, grid, ParticleKind::Negative, 1,
                                   state);
}

void ParticleOutput::save_energy(std::vector<NeParticleGroup> &groups,
                                 const NumericGridClass &grid,
                                 const SimulationState &state) {
  ParticleOutput::save_phase_space(groups, grid, ParticleKind::Positive, 2,
                                   state);
  ParticleOutput::save_phase_space(groups, grid, ParticleKind::Negative, 2,
                                   state);
}

void ParticleOutput::save_homogeneous_radial_distribution(
    const SimulationState &state) {
  constexpr int bin_count = 100;
  constexpr double rmax = 5.0;
  const double dr = rmax / bin_count;
  std::vector<double> radii(bin_count);
  for (int bin = 0; bin < bin_count; ++bin)
    radii[bin] = (bin + 0.5) * dr;
  MacroOutput{}.save_macro(radii, "rdist", state);
}

void ParticleOutput::save_homogeneous_radial_distribution(
    int bin_count, const SimulationState &state) {
  if (bin_count <= 0)
    throw std::invalid_argument("homogeneous radial distribution needs bins");
  constexpr double rmax = 10.0;
  const double dr = rmax / bin_count;
  std::vector<double> radii(bin_count);
  for (int bin = 0; bin < bin_count; ++bin)
    radii[bin] = (bin + 0.5) * dr;
  MacroOutput{}.save_macro(radii, "rdist", state);
}

void ParticleOutput::save_homogeneous_distribution(
    const NeParticleGroup &group, int bin_count, int case_index,
    const SimulationState &state) {
  if (bin_count <= 0 || case_index < 0 || case_index > 2)
    throw std::invalid_argument("invalid homogeneous distribution parameters");
  constexpr double rmax = 10.0;
  const char *suffixes[] = {"", "_before", "_after"};
  const std::string suffix = suffixes[case_index];

  const auto save_species = [&](ParticleKind kind, const char *name) {
    const auto &particles = group.list(kind);
    std::vector<double> speeds;
    speeds.reserve(particles.size());
    for (const auto &particle : particles) {
      const auto velocity = particle.velocity();
      double speed_squared = 0.0;
      for (const double component : velocity)
        speed_squared += component * component;
      speeds.push_back(std::sqrt(speed_squared));
    }
    std::vector<int> histogram(bin_count);
    Histogram{}.fixed_bins(speeds, histogram, 0.0, rmax);
    MacroOutput{}.save_macro(histogram, std::string(name) + suffix, state);
  };

  save_species(ParticleKind::Positive, "pdist");
  save_species(ParticleKind::Negative, "ndist");
  save_species(ParticleKind::Full, "fdist");
}

} // namespace coulomb
