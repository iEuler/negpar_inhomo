#pragma once

#include <vector>

#include "SimulationState.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;

class ParticleResampling {
  public:
	ParticleResampling(NumericGridClass& grid, ParaClass& parameters,
					   SimulationState& state)
		: gridRef(grid), parametersRef(parameters), stateRef(state) {}

	bool resampleHomogeneous(NeParticleGroup& groups,
							 double previousSignedWeight = -1.0,
							 double previousFullWeight = -1.0);
	void resample(std::vector<NeParticleGroup>& groups);
	void resampleFullHomogeneous(NeParticleGroup& groups,
								 double newEffectiveParticles,
								 double effectiveParticles, int frequency,
								 double dxSpace);
	void resampleFull(std::vector<NeParticleGroup>& groups,
					  double newEffectiveParticles, int frequency);
	void resampleFullPreservingMass(std::vector<NeParticleGroup>& groups,
									int oldCount);
	void synchronizeCoarse(std::vector<NeParticleGroup>& groups);

  private:
	NumericGridClass& gridRef;
	ParaClass& parametersRef;
	SimulationState& stateRef;
};

} // namespace coulomb
