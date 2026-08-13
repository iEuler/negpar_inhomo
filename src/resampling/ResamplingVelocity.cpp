#include "ResamplingVelocity.h"

#include <array>
#include <stdexcept>

#include "Constants.h"
#include "Particle.h"
#include "ParticleGroup.h"

namespace coulomb::resampling {
namespace {

std::array<double, 3> velocitySpans(const std::vector<double>& velocityBounds) {
	if (velocityBounds.size() != 6)
		throw std::invalid_argument(
			"velocity bounds must contain three ranges");

	std::array<double, 3> spans{};
	for (std::size_t component = 0; component < spans.size(); ++component) {
		spans[component] =
			velocityBounds[2 * component + 1] - velocityBounds[2 * component];
		if (!(spans[component] > 0.0))
			throw std::invalid_argument(
				"velocity bounds must have positive spans");
	}
	return spans;
}

} // namespace

NeParticleGroup
ResamplingVelocity::normalizeSigned(const NeParticleGroup& particles) {
	const auto spans = velocitySpans(particles.xyzMinMax);
	NeParticleGroup normalized;

	for (const auto kind : {ParticleKind::Positive, ParticleKind::Negative}) {
		for (const auto& particle : particles.list(kind)) {
			std::array<double, 3> velocity{};
			for (std::size_t component = 0; component < velocity.size();
				 ++component)
				velocity[component] =
					(particle.velocity(static_cast<int>(component)) -
					 particles.xyzMinMax[2 * component]) *
					2.0 * pi / spans[component];
			Particle1D3D normalizedParticle;
			normalizedParticle.setVelocity(velocity);
			normalized.pushBack(normalizedParticle, kind);
		}
	}

	normalized.rhoM = particles.rhoM;
	normalized.u1M =
		(particles.u1M - particles.xyzMinMax[0]) * 2.0 * pi / spans[0];
	normalized.u2M =
		(particles.u2M - particles.xyzMinMax[2]) * 2.0 * pi / spans[1];
	normalized.u3M =
		(particles.u3M - particles.xyzMinMax[4]) * 2.0 * pi / spans[2];
	normalized.t1M = particles.tprtM * 4.0 * pi * pi / (spans[0] * spans[0]);
	normalized.t2M = particles.tprtM * 4.0 * pi * pi / (spans[1] * spans[1]);
	normalized.t3M = particles.tprtM * 4.0 * pi * pi / (spans[2] * spans[2]);
	return normalized;
}

void ResamplingVelocity::restore(std::vector<Particle1D3D>& particles,
								 const std::vector<double>& velocityBounds) {
	const auto spans = velocitySpans(velocityBounds);
	for (auto& particle : particles) {
		std::array<double, 3> velocity{};
		for (std::size_t component = 0; component < velocity.size();
			 ++component)
			velocity[component] =
				velocityBounds[2 * component] +
				particle.velocity(static_cast<int>(component)) *
					spans[component] / (2.0 * pi);
		particle.setVelocity(velocity);
	}
}

} // namespace coulomb::resampling
