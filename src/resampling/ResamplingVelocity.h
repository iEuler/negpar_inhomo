#pragma once

#include <vector>

namespace coulomb {
class NeParticleGroup;
class Particle1D3D;

namespace resampling {

class ResamplingVelocity {
  public:
	NeParticleGroup normalizeSigned(const NeParticleGroup& particles);
	void restore(std::vector<Particle1D3D>& particles,
				 const std::vector<double>& velocityBounds);
};

} // namespace resampling
} // namespace coulomb
