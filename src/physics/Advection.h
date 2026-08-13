#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class Particle1d3d;
class ParticleGroup;
enum class ParticleKind;

class Advection {
  public:
	Advection(const NumericGridClass& grid, SimulationState& state)
		: grid_(grid), state_(state) {}

	Particle1d3d move_particle(const Particle1d3d& particle,
							   double electric_field);
	int find_particle_group(Particle1d3d& particle) const;

	void relocate_particle(std::vector<ParticleGroup>& groups, int group_before,
						   int particle_index, int group_after);
	void reset_moved_flags(std::vector<ParticleGroup>& groups);
	void advance(std::vector<ParticleGroup>& groups);

	void relocate_particle(std::vector<NeParticleGroup>& groups,
						   ParticleKind kind, int group_before,
						   int particle_index, int group_after);
	void reset_moved_flags(std::vector<NeParticleGroup>& groups,
						   ParticleKind kind);

	void advance(std::vector<NeParticleGroup>& groups, ParticleKind kind);
	void advance(std::vector<NeParticleGroup>& groups);

  private:
	const NumericGridClass& grid_;
	SimulationState& state_;
};

} // namespace coulomb
