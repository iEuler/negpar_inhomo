#pragma once

#include <initializer_list>
#include <vector>

#include "Particle.h"

namespace coulomb {
class NeParticleGroup;
class Particle1D3D;

namespace resampling {

class ResamplingVelocity {
  public:
	NeParticleGroup normalizeSigned(const NeParticleGroup& particles);
	NeParticleGroup normalizeKinds(
		const NeParticleGroup& particles,
		std::initializer_list<ParticleKind> kinds);
	void restore(std::vector<Particle1D3D>& particles,
				 const std::vector<double>& velocityBounds);
};

} // namespace resampling
} // namespace coulomb
