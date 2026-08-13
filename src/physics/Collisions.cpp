#include "Collisions.h"

#include <cmath>

#include "Constants.h"
#include "Particle.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {

std::pair<std::vector<double>, std::vector<double>>
CollisionOperator::collidePair(const std::vector<double>& velocity1,
							   const std::vector<double>& velocity2) {
	const auto& parameters = parametersRef;
	auto& random = randomContext;
	std::vector<double> velocity1After(3), velocity2After(3);

	if (parameters.methodBinaryColl == BinaryCollisionMethod::TA) {
		std::vector<double> relativeVelocity(3), velocityChange(3);
		for (int component = 0; component < 3; ++component)
			relativeVelocity[component] =
				velocity1[component] - velocity2[component];

		double relativeSpeed = 0.0;
		for (int component = 0; component < 3; ++component)
			relativeSpeed +=
				relativeVelocity[component] * relativeVelocity[component];
		relativeSpeed = std::sqrt(relativeSpeed);

		const double variance = parameters.coeffBinaryColl * parameters.dt /
								(relativeSpeed * relativeSpeed * relativeSpeed);
		const double delta =
			std::sqrt(variance) * RandomSampling(random).normal();
		const double phi = 2.0 * pi * RandomSampling(random).uniform();
		const double sine = 2.0 * delta / (1.0 + delta * delta);
		const double cosine = 1.0 - 2.0 * delta * delta / (1.0 + delta * delta);

		const double perpendicularSpeed =
			std::sqrt(relativeVelocity[0] * relativeVelocity[0] +
					  relativeVelocity[1] * relativeVelocity[1]) +
			1e-10;
		velocityChange[0] = (relativeVelocity[0] / perpendicularSpeed) *
								relativeVelocity[2] * sine * std::cos(phi) -
							(relativeVelocity[1] / perpendicularSpeed) *
								relativeSpeed * sine * std::sin(phi) -
							relativeVelocity[0] * (1.0 - cosine);
		velocityChange[1] = (relativeVelocity[1] / perpendicularSpeed) *
								relativeVelocity[2] * sine * std::cos(phi) +
							(relativeVelocity[0] / perpendicularSpeed) *
								relativeSpeed * sine * std::sin(phi) -
							relativeVelocity[1] * (1.0 - cosine);
		velocityChange[2] = -perpendicularSpeed * sine * std::cos(phi) -
							relativeVelocity[2] * (1.0 - cosine);

		for (int component = 0; component < 3; ++component) {
			velocity1After[component] =
				velocity1[component] + 0.5 * velocityChange[component];
			velocity2After[component] =
				velocity2[component] - 0.5 * velocityChange[component];
		}
	}

	return {std::move(velocity1After), std::move(velocity2After)};
}

void CollisionOperator::collideHomogeneous(std::vector<Particle1D3D>& particles,
										   int particleCount) {
	auto& random = randomContext;
	const auto permutation =
		RandomSampling(random).permutation(particleCount, particleCount);
	for (int pair = 0; pair < particleCount / 2; ++pair) {
		const int first = permutation[2 * pair] - 1;
		const int second = permutation[2 * pair + 1] - 1;
		const auto velocities = collidePair(particles[first].velocity(),
											particles[second].velocity());
		particles[first].setVelocity(velocities.first);
		particles[second].setVelocity(velocities.second);
	}
}

} // namespace coulomb
