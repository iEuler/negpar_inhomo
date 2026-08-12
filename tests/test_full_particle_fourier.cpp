#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Constants.h"
#include "FullParticleFourier.h"
#include "ParticleGroup.h"

TEST_CASE("negpar.unit.resampling.full-particle Fourier interpolation "
          "preserves dimensions and order",
          "[resampling][full-particle][fourier]") {
  const std::vector<std::complex<double>> coefficients(8, {1.0, 0.0});

  const auto first = coulomb::FullParticleFourier::interpolate_derivatives(
      coefficients, 2, 2, 2, 2);
  const auto second = coulomb::FullParticleFourier::interpolate_derivatives(
      coefficients, 2, 2, 2, 2);

  REQUIRE(first == second);
  REQUIRE(first.size() == 10);
  for (const auto &component : first)
    REQUIRE(component.size() == 64);

  const auto at_index = coulomb::FullParticleFourier::values_at(first, 17);
  REQUIRE(at_index.size() == 10);
  for (std::size_t component = 0; component < first.size(); ++component)
    REQUIRE(at_index[component] == first[component][17]);
}

TEST_CASE("negpar.unit.resampling.full-particle Maxwellian derivatives are "
          "finite and ordered",
          "[resampling][full-particle][fourier]") {
  constexpr int frequency = 2;
  constexpr int augmentation = 2;
  constexpr int side = frequency * augmentation;
  std::vector<std::vector<double>> derivatives(
      10, std::vector<double>(side * side * side, 0.0));

  const std::vector<double> velocity{coulomb::pi, coulomb::pi, coulomb::pi};
  const std::vector<double> temperature{1.0, 2.0, 3.0};
  coulomb::FullParticleFourier::add_maxwellian(
      1.0, velocity, temperature, 0.0, derivatives, frequency, augmentation);

  for (const auto &component : derivatives)
    REQUIRE(std::all_of(component.begin(), component.end(),
                        [](double value) { return std::isfinite(value); }));

  const int center = 2 + side * (2 + side * 2);
  const double value = derivatives[0][center];
  REQUIRE(value > 0.0);
  REQUIRE(derivatives[1][center] == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(derivatives[2][center] == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(derivatives[3][center] == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(derivatives[4][center] == Catch::Approx(-value));
  REQUIRE(derivatives[5][center] == Catch::Approx(-value / 2.0));
  REQUIRE(derivatives[6][center] == Catch::Approx(-value / 3.0));
  REQUIRE(derivatives[7][center] == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(derivatives[8][center] == Catch::Approx(0.0).margin(1e-15));
  REQUIRE(derivatives[9][center] == Catch::Approx(0.0).margin(1e-15));

  const int left = 2 + side * (2 + side * 1);
  const int right = 2 + side * (2 + side * 3);
  REQUIRE(derivatives[0][left] == Catch::Approx(derivatives[0][right]));
  REQUIRE(derivatives[1][left] == Catch::Approx(-derivatives[1][right]));
}

TEST_CASE("negpar.unit.resampling.full-particle Fourier coefficients and "
          "filtering are deterministic",
          "[resampling][full-particle][fourier]") {
  coulomb::NeParticleGroup particles;
  particles.push_back(
      coulomb::Particle1d3d({coulomb::pi, coulomb::pi, coulomb::pi}),
      coulomb::ParticleKind::Positive);
  particles.push_back(coulomb::Particle1d3d(
                          {0.5 * coulomb::pi, coulomb::pi, 1.5 * coulomb::pi}),
                      coulomb::ParticleKind::Negative);

  auto first =
      coulomb::FullParticleFourier::approximate_transform(particles, 2, 2, 2);
  auto second =
      coulomb::FullParticleFourier::approximate_transform(particles, 2, 2, 2);
  REQUIRE(first == second);
  REQUIRE(first.size() == 8);

  const auto original = first;
  std::vector<int> flags(first.size(), 0);
  coulomb::FullParticleFourier::filter(first, flags,
                                       static_cast<int>(first.size()));
  REQUIRE(first == original);
  REQUIRE(std::all_of(flags.begin(), flags.end(),
                      [](int flag) { return flag == 1; }));

  const std::vector<double> fixture{7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const auto bounds = coulomb::FullParticleFourier::upper_bound(2, fixture);
  REQUIRE(bounds.size() == fixture.size());
  for (const double bound : bounds)
    REQUIRE(bound == Catch::Approx(7.0));
}
