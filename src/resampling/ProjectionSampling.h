#pragma once

#include <vector>

#include "SimulationTypes.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
struct RandomContext;

class ProjectionSampling {
  public:
	void sampleHomogeneous(NeParticleGroup& groups,
						   const NumericGridClass& grid, RandomContext& random,
						   ProjectionMode mode = ProjectionMode::FullMicroMacro);
	void sample(std::vector<NeParticleGroup>& groups,
				const NumericGridClass& grid, RandomContext& random);
	void sample(std::vector<NeParticleGroup>& groups,
				const NumericGridClass& grid, RandomContext& random,
				ProjectionMode mode);
};

} // namespace coulomb
