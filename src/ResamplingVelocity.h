#pragma once

#include <vector>

namespace coulomb {
class NeParticleGroup;
class Particle1d3d;

namespace resampling {

NeParticleGroup normalize_signed_velocities(const NeParticleGroup& particles);
void restore_velocities(std::vector<Particle1d3d>& particles,
                        const std::vector<double>& velocity_bounds);

}  // namespace resampling
}  // namespace coulomb
