#pragma once

#include "ParticleGroup.h"

namespace coulomb {

struct ParticlePartition {
	NeParticleGroup core;
	NeParticleGroup tail;
};

class ParticlePartitioning {
  public:
	static ParticlePartition split(const NeParticleGroup& source,
								   double standardDeviationCutoff);
};

} // namespace coulomb
