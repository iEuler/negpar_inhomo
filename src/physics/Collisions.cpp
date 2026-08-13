#include "Collisions.h"

#include <cmath>

#include "Constants.h"
#include "Particle.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {

std::pair<std::vector<double>, std::vector<double>>
CollisionOperator::collide_pair(const std::vector<double>& velocity1,
								const std::vector<double>& velocity2) {
	const auto& parameters = parameters_;
	auto& random = random_;
	std::vector<double> velocity1_after(3), velocity2_after(3);

	if (parameters.method_binarycoll == BinaryCollisionMethod::TA) {
		std::vector<double> relative_velocity(3), velocity_change(3);
		for (int component = 0; component < 3; ++component)
			relative_velocity[component] =
				velocity1[component] - velocity2[component];

		double relative_speed = 0.0;
		for (int component = 0; component < 3; ++component)
			relative_speed +=
				relative_velocity[component] * relative_velocity[component];
		relative_speed = std::sqrt(relative_speed);

		const double variance =
			parameters.coeff_binarycoll * parameters.dt /
			(relative_speed * relative_speed * relative_speed);
		const double delta =
			std::sqrt(variance) * RandomSampling(random).normal();
		const double phi = 2.0 * pi * RandomSampling(random).uniform();
		const double sine = 2.0 * delta / (1.0 + delta * delta);
		const double cosine = 1.0 - 2.0 * delta * delta / (1.0 + delta * delta);

		const double perpendicular_speed =
			std::sqrt(relative_velocity[0] * relative_velocity[0] +
					  relative_velocity[1] * relative_velocity[1]) +
			1e-10;
		velocity_change[0] = (relative_velocity[0] / perpendicular_speed) *
								 relative_velocity[2] * sine * std::cos(phi) -
							 (relative_velocity[1] / perpendicular_speed) *
								 relative_speed * sine * std::sin(phi) -
							 relative_velocity[0] * (1.0 - cosine);
		velocity_change[1] = (relative_velocity[1] / perpendicular_speed) *
								 relative_velocity[2] * sine * std::cos(phi) +
							 (relative_velocity[0] / perpendicular_speed) *
								 relative_speed * sine * std::sin(phi) -
							 relative_velocity[1] * (1.0 - cosine);
		velocity_change[2] = -perpendicular_speed * sine * std::cos(phi) -
							 relative_velocity[2] * (1.0 - cosine);

		for (int component = 0; component < 3; ++component) {
			velocity1_after[component] =
				velocity1[component] + 0.5 * velocity_change[component];
			velocity2_after[component] =
				velocity2[component] - 0.5 * velocity_change[component];
		}
	}

	return {std::move(velocity1_after), std::move(velocity2_after)};
}

void CollisionOperator::collide_homogeneous(
	std::vector<Particle1d3d>& particles, int particle_count) {
	auto& random = random_;
	const auto permutation =
		RandomSampling(random).permutation(particle_count, particle_count);
	for (int pair = 0; pair < particle_count / 2; ++pair) {
		const int first = permutation[2 * pair] - 1;
		const int second = permutation[2 * pair + 1] - 1;
		const auto velocities = collide_pair(particles[first].velocity(),
											 particles[second].velocity());
		particles[first].set_velocity(velocities.first);
		particles[second].set_velocity(velocities.second);
	}
}

} // namespace coulomb
