#include "ResamplingVelocity.h"

#include <array>
#include <stdexcept>

#include "Constants.h"
#include "Particle.h"
#include "ParticleGroup.h"

namespace coulomb::resampling {
namespace {

std::array<double, 3> velocity_spans(
    const std::vector<double>& velocity_bounds) {
  if (velocity_bounds.size() != 6)
    throw std::invalid_argument("velocity bounds must contain three ranges");

  std::array<double, 3> spans{};
  for (std::size_t component = 0; component < spans.size(); ++component) {
    spans[component] = velocity_bounds[2 * component + 1] -
                       velocity_bounds[2 * component];
    if (!(spans[component] > 0.0))
      throw std::invalid_argument("velocity bounds must have positive spans");
  }
  return spans;
}

}  // namespace

NeParticleGroup normalize_signed_velocities(
    const NeParticleGroup& particles) {
  const auto spans = velocity_spans(particles.xyz_minmax);
  NeParticleGroup normalized;

  for (const auto kind :
       {ParticleKind::Positive, ParticleKind::Negative}) {
    for (const auto& particle : particles.list(kind)) {
      std::array<double, 3> velocity{};
      for (std::size_t component = 0; component < velocity.size(); ++component)
        velocity[component] =
            (particle.velocity(static_cast<int>(component)) -
             particles.xyz_minmax[2 * component]) *
            2.0 * pi / spans[component];
      Particle1d3d normalized_particle;
      normalized_particle.set_velocity(velocity);
      normalized.push_back(normalized_particle, kind);
    }
  }

  normalized.rhoM = particles.rhoM;
  normalized.u1M =
      (particles.u1M - particles.xyz_minmax[0]) * 2.0 * pi / spans[0];
  normalized.u2M =
      (particles.u2M - particles.xyz_minmax[2]) * 2.0 * pi / spans[1];
  normalized.u3M =
      (particles.u3M - particles.xyz_minmax[4]) * 2.0 * pi / spans[2];
  normalized.T1M = particles.TprtM * 4.0 * pi * pi /
                   (spans[0] * spans[0]);
  normalized.T2M = particles.TprtM * 4.0 * pi * pi /
                   (spans[1] * spans[1]);
  normalized.T3M = particles.TprtM * 4.0 * pi * pi /
                   (spans[2] * spans[2]);
  return normalized;
}

void restore_velocities(std::vector<Particle1d3d>& particles,
                        const std::vector<double>& velocity_bounds) {
  const auto spans = velocity_spans(velocity_bounds);
  for (auto& particle : particles) {
    std::array<double, 3> velocity{};
    for (std::size_t component = 0; component < velocity.size(); ++component)
      velocity[component] = velocity_bounds[2 * component] +
                            particle.velocity(static_cast<int>(component)) *
                                spans[component] / (2.0 * pi);
    particle.set_velocity(velocity);
  }
}

}  // namespace coulomb::resampling
