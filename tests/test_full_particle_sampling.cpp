#include <cmath>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "FullParticleSampling.h"
#include "ParticleGroup.h"
#include "RandomContext.h"
#include "RandomSampling.h"

namespace {

coulomb::NeParticleGroup full_particle_fixture() {
  coulomb::NeParticleGroup particles;
  int position_index = 0;
  for (const double vx : {-1.0, 1.0}) {
    for (const double vy : {-1.0, 1.0}) {
      for (const double vz : {-1.0, 1.0}) {
        const double position = -0.9 + 0.2 * position_index++;
        const coulomb::Particle1d3d particle(position, {vx, vy, vz});
        particles.push_back(particle, coulomb::ParticleKind::Positive);
        particles.push_back(particle, coulomb::ParticleKind::Negative);
      }
    }
  }
  particles.rhoM = 1.0;
  particles.u1M = 0.0;
  particles.u2M = 0.0;
  particles.u3M = 0.0;
  particles.TprtM = 0.25;
  return particles;
}

void require_same_list(const coulomb::NeParticleGroup &first,
                       const coulomb::NeParticleGroup &second,
                       coulomb::ParticleKind kind) {
  REQUIRE(first.size(kind) == second.size(kind));
  for (int index = 0; index < first.size(kind); ++index) {
    const auto &lhs = first.list(index, kind);
    const auto &rhs = second.list(index, kind);
    REQUIRE(lhs.position() == rhs.position());
    REQUIRE(lhs.velocity() == rhs.velocity());
  }
}

} // namespace

TEST_CASE("negpar.unit.resampling.full-particle reconstruction replays exactly "
          "for a fixed seed",
          "[resampling][full-particle][sampling]") {
  auto first_input = full_particle_fixture();
  auto second_input = first_input;
  const auto original = first_input;

  coulomb::RandomContext first_random;
  coulomb::RandomContext second_random;
  first_random.reseed(20260809);
  second_random.reseed(20260809);

  auto first = coulomb::FullParticleSampling{}.resample(
      first_input, 2, 0.1, 0.05, 1.0, first_random);
  auto second = coulomb::FullParticleSampling{}.resample(
      second_input, 2, 0.1, 0.05, 1.0, second_random);

  REQUIRE(first.size(coulomb::ParticleKind::Positive) == 0);
  REQUIRE(first.size(coulomb::ParticleKind::Negative) == 0);
  REQUIRE(first.size(coulomb::ParticleKind::Full) == 12);
  REQUIRE(first.size(coulomb::ParticleKind::Full) ==
          second.size(coulomb::ParticleKind::Full));
  require_same_list(first, second, coulomb::ParticleKind::Full);

  require_same_list(first_input, original, coulomb::ParticleKind::Positive);
  require_same_list(first_input, original, coulomb::ParticleKind::Negative);
  require_same_list(first_input, original, coulomb::ParticleKind::Full);

  for (const auto &particle : first.list(coulomb::ParticleKind::Full)) {
    REQUIRE(std::isfinite(particle.position()));
    REQUIRE(particle.position() == 0.0);
    for (int component = 0; component < 3; ++component) {
      const double velocity = particle.velocity(component);
      REQUIRE(std::isfinite(velocity));
      REQUIRE(velocity >= first_input.xyz_minmax[2 * component]);
      REQUIRE(velocity <= first_input.xyz_minmax[2 * component + 1]);
    }
  }

  first.computemoments();
  const double sampled_mass = 0.05 * first.full_moments.m0;
  const double sampled_momentum = 0.05 * first.full_moments.m11;
  REQUIRE(sampled_mass == Catch::Approx(1.0).margin(0.5));
  REQUIRE(sampled_momentum == Catch::Approx(0.0).margin(0.5));
}
