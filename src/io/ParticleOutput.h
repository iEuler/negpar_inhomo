#pragma once

#include <vector>

#include "Particle.h"

namespace coulomb {

class NumericGridClass;
class ParticleGroup;
class NeParticleGroup;
struct SimulationState;

class ParticleOutput {
  public:
	void saveDistribution(std::vector<ParticleGroup>& groups,
						  const NumericGridClass& grid,
						  const SimulationState& state);
	void saveDistribution(std::vector<NeParticleGroup>& groups,
						  const NumericGridClass& grid, ParticleKind kind,
						  const SimulationState& state);
	void saveDistribution(std::vector<NeParticleGroup>& groups,
						  const NumericGridClass& grid,
						  const SimulationState& state);
	void savePhaseSpace(std::vector<ParticleGroup>& groups,
						const NumericGridClass& grid,
						const SimulationState& state);
	void savePhaseSpace(std::vector<NeParticleGroup>& groups,
						const NumericGridClass& grid, ParticleKind kind,
						int quantity, const SimulationState& state);
	void savePhaseSpace(std::vector<NeParticleGroup>& groups,
						const NumericGridClass& grid,
						const SimulationState& state);
	void saveEnergy(std::vector<NeParticleGroup>& groups,
					const NumericGridClass& grid, const SimulationState& state);
	void saveHomogeneousRadialDistribution(const SimulationState& state);
	void saveHomogeneousRadialDistribution(int binCount,
										   const SimulationState& state);
	void saveHomogeneousDistribution(const NeParticleGroup& group, int binCount,
									 int caseIndex,
									 const SimulationState& state);
};

} // namespace coulomb
