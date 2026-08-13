#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class Particle1D3D;
class ParticleGroup;
enum class ParticleKind;

class Advection {
  public:
	Advection(const NumericGridClass& grid, SimulationState& state)
		: gridRef(grid), stateRef(state) {}

	Particle1D3D moveParticle(const Particle1D3D& particle,
							  double electricField);
	int findParticleGroup(Particle1D3D& particle) const;

	void relocateParticle(std::vector<ParticleGroup>& groups, int groupBefore,
						  int particleIndex, int groupAfter);
	void resetMovedFlags(std::vector<ParticleGroup>& groups);
	void advance(std::vector<ParticleGroup>& groups);

	void relocateParticle(std::vector<NeParticleGroup>& groups,
						  ParticleKind kind, int groupBefore, int particleIndex,
						  int groupAfter);
	void resetMovedFlags(std::vector<NeParticleGroup>& groups,
						 ParticleKind kind);

	void advance(std::vector<NeParticleGroup>& groups, ParticleKind kind);
	void advance(std::vector<NeParticleGroup>& groups);

  private:
	const NumericGridClass& gridRef;
	SimulationState& stateRef;
};

} // namespace coulomb
