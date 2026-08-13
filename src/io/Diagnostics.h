#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
enum class ParticleKind;

// Diagnostics are pure reductions over particle and grid data. They do not
// mutate run counters, output policy, or timing, so they intentionally remain
// state-free while stateful callers pass SimulationState to the APIs that need
// it.

class Diagnostics {
  public:
	explicit Diagnostics(const NumericGridClass& grid) : gridRef(grid) {}

	double electricEnergy(const std::vector<NeParticleGroup>& groups) const;
	double fullElectricEnergy(const std::vector<NeParticleGroup>& groups) const;
	double totalEnergy(const std::vector<NeParticleGroup>& groups) const;
	double fullTotalEnergy(const std::vector<NeParticleGroup>& groups) const;
	int particleCount(const std::vector<NeParticleGroup>& groups, int size,
					  ParticleKind kind) const;

  private:
	const NumericGridClass& gridRef;
};

} // namespace coulomb
