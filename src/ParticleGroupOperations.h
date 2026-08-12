#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
enum class ParticleKind;
struct RandomContext;

class ParticleGroupOperations {
  public:
	void merge_signed(NeParticleGroup& groups,
					  const NeParticleGroup& new_groups);
	void merge_full(NeParticleGroup& groups, const NeParticleGroup& new_groups);
	void merge(NeParticleGroup& groups, const NeParticleGroup& new_groups,
			   const std::vector<ParticleKind>& particle_types);
	void assign_positions(NeParticleGroup& groups, double xmin, double xmax,
						  RandomContext& random);
};

} // namespace coulomb
