#include <type_traits>

#include <catch2/catch_test_macros.hpp>

#include "Resampler.h"
#include "ResamplerHelper.h"
#include "_global_variables.h"
#include "utils.h"

TEST_CASE("experimental resampler owns isolated symbols and explicit RNG",
          "[experimental][resampler]") {
  static_assert(std::is_constructible_v<coulomb::experimental::Resampler,
                                        double, double, size_t, bool, double>);

  coulomb::NeParticleGroup first;
  coulomb::NeParticleGroup second;
  coulomb::RandomContext first_random;
  coulomb::RandomContext second_random;
  coulomb::reseed_random(first_random, 9876);
  coulomb::reseed_random(second_random, 9876);

  double first_bound = 1.0;
  double second_bound = 1.0;
  const std::vector<double> sample{coulomb::pi, coulomb::pi, coulomb::pi};

  coulomb::experimental::acceptSampled(sample, first, 0.25, first_bound,
                                        false, first_random);
  coulomb::experimental::acceptSampled(sample, second, 0.25, second_bound,
                                        false, second_random);

  REQUIRE(first_bound == second_bound);
  REQUIRE(first.size(coulomb::ParticleKind::Positive) ==
          second.size(coulomb::ParticleKind::Positive));
  REQUIRE(first.size(coulomb::ParticleKind::Negative) ==
          second.size(coulomb::ParticleKind::Negative));
}
