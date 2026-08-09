#include "Collisions.h"

#include <cmath>

#include "Constants.h"
#include "Particle.h"
#include "SimulationConfig.h"
#include "RandomSampling.h"

namespace coulomb {

std::pair<std::vector<double>, std::vector<double>> coulombBinary3d(
    const std::vector<double>& velocity1,
    const std::vector<double>& velocity2, const ParaClass& parameters,
    RandomContext& random) {
  std::vector<double> velocity1_after(3), velocity2_after(3);

  if (parameters.method_binarycoll == BinaryCollisionMethod::TA) {
    std::vector<double> relative_velocity(3), velocity_change(3);
    for (int component = 0; component < 3; ++component)
      relative_velocity[component] =
          velocity1[component] - velocity2[component];

    double relative_speed = 0.0;
    for (int component = 0; component < 3; ++component)
      relative_speed += relative_velocity[component] *
                        relative_velocity[component];
    relative_speed = std::sqrt(relative_speed);

    const double variance = parameters.coeff_binarycoll * parameters.dt /
                            (relative_speed * relative_speed * relative_speed);
    const double delta = std::sqrt(variance) * myrandn(random);
    const double phi = 2.0 * pi * myrand(random);
    const double sine = 2.0 * delta / (1.0 + delta * delta);
    const double cosine =
        1.0 - 2.0 * delta * delta / (1.0 + delta * delta);

    const double perpendicular_speed =
        std::sqrt(relative_velocity[0] * relative_velocity[0] +
                  relative_velocity[1] * relative_velocity[1]) +
        1e-10;
    velocity_change[0] =
        (relative_velocity[0] / perpendicular_speed) * relative_velocity[2] *
            sine * std::cos(phi) -
        (relative_velocity[1] / perpendicular_speed) * relative_speed * sine *
            std::sin(phi) -
        relative_velocity[0] * (1.0 - cosine);
    velocity_change[1] =
        (relative_velocity[1] / perpendicular_speed) * relative_velocity[2] *
            sine * std::cos(phi) +
        (relative_velocity[0] / perpendicular_speed) * relative_speed * sine *
            std::sin(phi) -
        relative_velocity[1] * (1.0 - cosine);
    velocity_change[2] =
        -perpendicular_speed * sine * std::cos(phi) -
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

void coulomb_collision_homo(std::vector<Particle1d3d>& particles,
                            int particle_count,
                            const ParaClass& parameters,
                            RandomContext& random) {
  const auto permutation = myrandperm(particle_count, particle_count, random);
  for (int pair = 0; pair < particle_count / 2; ++pair) {
    const int first = permutation[2 * pair] - 1;
    const int second = permutation[2 * pair + 1] - 1;
    const auto velocities = coulombBinary3d(
        particles[first].velocity(), particles[second].velocity(), parameters,
        random);
    particles[first].set_velocity(velocities.first);
    particles[second].set_velocity(velocities.second);
  }
}

}  // namespace coulomb
