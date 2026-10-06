#include "Collisions.h"

#include <cmath>
#include <stdexcept>

#include "Constants.h"
#include "Particle.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {

std::pair<std::vector<double>, std::vector<double>>
CollisionOperator::collidePair(const std::vector<double>& velocity1,
							   const std::vector<double>& velocity2, double backgroundDensity) {
	if (velocity1.size() != 3 || velocity2.size() != 3 ||
		!std::isfinite(backgroundDensity) || backgroundDensity < 0.0)
		throw std::invalid_argument("Invalid collision velocities or background density");
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
		if (relativeSpeed == 0.0 || backgroundDensity == 0.0 || parameters.coeffBinaryColl == 0.0)
			return {velocity1, velocity2};

		const double variance = backgroundDensity * parameters.coeffBinaryColl * parameters.dt /
								(relativeSpeed * relativeSpeed * relativeSpeed);
		const double delta =
			std::sqrt(variance) * RandomSampling(random).normal();
		const double phi = 2.0 * pi * RandomSampling(random).uniform();
		const double sine = 2.0 * delta / (1.0 + delta * delta);
		const double cosine = 1.0 - 2.0 * delta * delta / (1.0 + delta * delta);

		const double perpendicularSpeed =
			std::sqrt(relativeVelocity[0] * relativeVelocity[0] +
					  relativeVelocity[1] * relativeVelocity[1]);
		if (perpendicularSpeed <= 1e-14 * relativeSpeed) {
			velocityChange[0] = relativeSpeed * sine * std::cos(phi);
			velocityChange[1] = relativeSpeed * sine * std::sin(phi);
			velocityChange[2] = -relativeVelocity[2] * (1.0 - cosine);
		} else {
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
		}

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
										   int particleCount, double backgroundDensity) {
	auto& random = randomContext;
	const auto permutation =
		RandomSampling(random).permutation(particleCount, particleCount);
	for (int pair = 0; pair < particleCount / 2; ++pair) {
		const int first = permutation[2 * pair] - 1;
		const int second = permutation[2 * pair + 1] - 1;
		const auto velocities = collidePair(particles[first].velocity(),
											particles[second].velocity(), backgroundDensity);
		particles[first].setVelocity(velocities.first);
		particles[second].setVelocity(velocities.second);
	}
}

} // namespace coulomb
