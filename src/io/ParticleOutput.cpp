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
void ParticleOutput::saveDistribution(std::vector<ParticleGroup>& groups,
									  const NumericGridClass& grid,
									  const SimulationState& state) {
	const int size = grid.nx;
	const int bins = grid.nv;
	std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
	for (int cell = 0; cell < size; ++cell) {
		const auto& particles = groups[cell].list();
		std::vector<double> velocities;
		velocities.reserve(particles.size());
		for (const auto& particle : particles)
			velocities.push_back(particle.velocity(0));
		Histogram{}.fixedBins(velocities, counts[cell], grid.vmin, grid.vmax);
	}

	std::vector<std::vector<double>> distribution(size,
												  std::vector<double>(bins));
	for (int cell = 0; cell < size; ++cell)
		for (int bin = 0; bin < bins; ++bin)
			distribution[cell][bin] = counts[cell][bin] * grid.neff;
	MacroOutput{}.save2D(size, bins, distribution, "dist", state);
}

void ParticleOutput::saveDistribution(std::vector<NeParticleGroup>& groups,
									  const NumericGridClass& grid,
									  ParticleKind kind,
									  const SimulationState& state) {
	const int size = grid.nx;
	const int bins = grid.nv;
	std::vector<std::vector<int>> counts(size, std::vector<int>(bins));
	for (int cell = 0; cell < size; ++cell) {
		const auto& particles = groups[cell].list(kind);
		std::vector<double> velocities;
		velocities.reserve(particles.size());
		for (const auto& particle : particles)
			velocities.push_back(particle.velocity(0));
		Histogram{}.fixedBins(velocities, counts[cell], grid.vmin, grid.vmax);
	}

	const double effectiveParticles =
		kind == ParticleKind::Full ? grid.neffF : grid.neff;
	const double coefficient = effectiveParticles / grid.dx / grid.dv;
	std::vector<std::vector<double>> distribution(size,
												  std::vector<double>(bins));
	for (int cell = 0; cell < size; ++cell)
		for (int bin = 0; bin < bins; ++bin)
			distribution[cell][bin] = counts[cell][bin] * coefficient;

	const std::string name =
		kind == ParticleKind::Positive
			? "distp"
			: (kind == ParticleKind::Negative ? "distn" : "distf");
	MacroOutput{}.save2D(size, bins, distribution, name, state);
}

void ParticleOutput::saveDistribution(std::vector<NeParticleGroup>& groups,
									  const NumericGridClass& grid,
									  const SimulationState& state) {
	ParticleOutput::saveDistribution(groups, grid, ParticleKind::Positive,
									 state);
	ParticleOutput::saveDistribution(groups, grid, ParticleKind::Negative,
									 state);
	ParticleOutput::saveDistribution(groups, grid, ParticleKind::Full, state);
}

void ParticleOutput::savePhaseSpace(std::vector<ParticleGroup>& groups,
									const NumericGridClass& grid,
									const SimulationState& state) {
	const std::string name =
		"particle" +
		(state.filenameWithNumber
			 ? OutputPaths(state).formatIndex(state.saveIndex)
			 : "") +
		".txt";
	std::ofstream file(OutputPaths(state).resolve(name));
	if (!file)
		throw std::runtime_error("Unable to open particle output file");
	for (int cell = 0; cell < grid.nx; ++cell)
		for (const auto& particle : groups[cell].list())
			file << particle.position() << ' ' << particle.velocity(0) << ' ';
}

void ParticleOutput::savePhaseSpace(std::vector<NeParticleGroup>& groups,
									const NumericGridClass& grid,
									ParticleKind kind, int quantity,
									const SimulationState& state) {
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
			 ? OutputPaths(state).formatIndex(state.saveIndex)
			 : "") +
		".txt";
	std::ofstream file(OutputPaths(state).resolve(name));
	if (!file)
		throw std::runtime_error("Unable to open particle output file");
	for (int cell = 0; cell < grid.nx; ++cell) {
		for (const auto& particle : groups[cell].list(kind)) {
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

void ParticleOutput::savePhaseSpace(std::vector<NeParticleGroup>& groups,
									const NumericGridClass& grid,
									const SimulationState& state) {
	ParticleOutput::savePhaseSpace(groups, grid, ParticleKind::Positive, 1,
								   state);
	ParticleOutput::savePhaseSpace(groups, grid, ParticleKind::Negative, 1,
								   state);
}

void ParticleOutput::saveEnergy(std::vector<NeParticleGroup>& groups,
								const NumericGridClass& grid,
								const SimulationState& state) {
	ParticleOutput::savePhaseSpace(groups, grid, ParticleKind::Positive, 2,
								   state);
	ParticleOutput::savePhaseSpace(groups, grid, ParticleKind::Negative, 2,
								   state);
}

void ParticleOutput::saveHomogeneousRadialDistribution(
	const SimulationState& state) {
	constexpr int binCount = 100;
	constexpr double rMax = 5.0;
	const double dr = rMax / binCount;
	std::vector<double> radii(binCount);
	for (int bin = 0; bin < binCount; ++bin)
		radii[bin] = (bin + 0.5) * dr;
	MacroOutput{}.saveMacro(radii, "rdist", state);
}

void ParticleOutput::saveHomogeneousRadialDistribution(
	int binCount, const SimulationState& state) {
	if (binCount <= 0)
		throw std::invalid_argument(
			"homogeneous radial distribution needs bins");
	constexpr double rMax = 10.0;
	const double dr = rMax / binCount;
	std::vector<double> radii(binCount);
	for (int bin = 0; bin < binCount; ++bin)
		radii[bin] = (bin + 0.5) * dr;
	MacroOutput{}.saveMacro(radii, "rdist", state);
}

void ParticleOutput::saveHomogeneousDistribution(const NeParticleGroup& group,
												 int binCount, int caseIndex,
												 const SimulationState& state) {
	if (binCount <= 0 || caseIndex < 0 || caseIndex > 2)
		throw std::invalid_argument(
			"invalid homogeneous distribution parameters");
	constexpr double rMax = 10.0;
	const char* suffixes[] = {"", "_before", "_after"};
	const std::string suffix = suffixes[caseIndex];

	const auto saveSpecies = [&](ParticleKind kind, const char* name) {
		const auto& particles = group.list(kind);
		std::vector<double> speeds;
		speeds.reserve(particles.size());
		for (const auto& particle : particles) {
			const auto velocity = particle.velocity();
			double speedSquared = 0.0;
			for (const double component : velocity)
				speedSquared += component * component;
			speeds.push_back(std::sqrt(speedSquared));
		}
		std::vector<int> histogram(binCount);
		Histogram{}.fixedBins(speeds, histogram, 0.0, rMax);
		MacroOutput{}.saveMacro(histogram, std::string(name) + suffix, state);
	};

	saveSpecies(ParticleKind::Positive, "pdist");
	saveSpecies(ParticleKind::Negative, "ndist");
	saveSpecies(ParticleKind::Full, "fdist");
}

} // namespace coulomb
