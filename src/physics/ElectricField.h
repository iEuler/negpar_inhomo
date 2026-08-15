#pragma once

#include <vector>

#include "SimulationTypes.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParticleGroup;

class ElectricFieldSolver {
  public:
	explicit ElectricFieldSolver(
		const NumericGridClass& grid,
		HdpCouplingMode couplingMode = HdpCouplingMode::Decoupled)
		: gridRef(grid), couplingMode(couplingMode) {}

	std::vector<double> solvePoisson(const std::vector<double>& rho,
									 double lambda);
	void update(std::vector<ParticleGroup>& groups);
	void update(std::vector<NeParticleGroup>& groups);
	void updatePic(std::vector<NeParticleGroup>& groups);
	void updateFromCoarse(std::vector<NeParticleGroup>& groups);
	void clear(std::vector<NeParticleGroup>& groups);
	void updateFromDensity(std::vector<NeParticleGroup>& groups);

  private:
	const NumericGridClass& gridRef;
	HdpCouplingMode couplingMode;
};

} // namespace coulomb
