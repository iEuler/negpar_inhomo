#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Grid.h"
#include "InitialConditions.h"
#include "ParticleGroup.h"
#include "ParticleInitialization.h"
#include "RandomContext.h"
#include "SimulationConfig.h"
#include "RandomSampling.h"

TEST_CASE("negpar.unit.initialization.selected initial conditions preserve Landau defaults",
          "[initialization][conditions]") {
  coulomb::NumericGridClass grid;
  auto initial_data = coulomb::make_initial_conditions(grid);

  REQUIRE(initial_data.probname == "LandauDamping");
  REQUIRE(initial_data.LD_alpha == Catch::Approx(0.4));

  coulomb::configure_initial_condition(initial_data, grid, 0);
  REQUIRE(initial_data.rho ==
          Catch::Approx(1.0 + initial_data.LD_alpha * std::sin(grid.x[0])));
  REQUIRE(initial_data.velocity[0] == 0.0);
  REQUIRE(initial_data.velocity[1] == 0.0);
  REQUIRE(initial_data.velocity[2] == 0.0);
  REQUIRE(initial_data.Tprt == 1.0);
}

TEST_CASE("negpar.unit.initialization.particle initialization replays a small fixed-seed fixture",
          "[initialization][particles]") {
  coulomb::IniValClass initial_data;
  initial_data.probname = "Delta";
  initial_data.rho = 1.0;
  initial_data.velocity[0] = 1.0;
  initial_data.velocity[1] = -0.5;
  initial_data.velocity[2] = 0.25;
  initial_data.Tprt = 1.0;

  coulomb::NeParticleGroup first;
  coulomb::NeParticleGroup second;
  first.set_xrange(-0.5, 0.5);
  second.set_xrange(-0.5, 0.5);

  coulomb::RandomContext first_random;
  coulomb::RandomContext second_random;
  coulomb::reseed_random(first_random, 1401);
  coulomb::reseed_random(second_random, 1401);

  coulomb::initialize_Negpar(first, initial_data, 0.1, 0.1, 1.0,
                             first_random);
  coulomb::initialize_Negpar(second, initial_data, 0.1, 0.1, 1.0,
                             second_random);

  REQUIRE(first.size(coulomb::ParticleKind::Positive) == 0);
  REQUIRE(first.size(coulomb::ParticleKind::Negative) == 0);
  REQUIRE(first.size(coulomb::ParticleKind::Full) == 10);
  REQUIRE(second.size(coulomb::ParticleKind::Full) == 10);

  for (int index = 0; index < first.size(coulomb::ParticleKind::Full);
       ++index) {
    const auto& lhs = first.list(index, coulomb::ParticleKind::Full);
    const auto& rhs = second.list(index, coulomb::ParticleKind::Full);
    REQUIRE(lhs.position() == rhs.position());
    REQUIRE(lhs.velocity() == rhs.velocity());
    REQUIRE(lhs.position() >= -0.5);
    REQUIRE(lhs.position() <= 0.5);
    REQUIRE((lhs.velocity() ==
             std::vector<double>{1.0, -0.5, 0.25}));
  }

  REQUIRE(first.full_moments.m0 == Catch::Approx(10.0));
  REQUIRE(first.full_moments.m11 == Catch::Approx(10.0));
  REQUIRE(first.full_moments.m12 == Catch::Approx(-5.0));
  REQUIRE(first.full_moments.m13 == Catch::Approx(2.5));
}

TEST_CASE("negpar.unit.initialization.Two-Stream preprocessing produces finite positive coefficients",
          "[initialization][conditions][two-stream]") {
  coulomb::IniValClass initial_data;
  coulomb::initialize_TwoStreamInstab(initial_data);

  REQUIRE(initial_data.TSI_rhof > 0.0);
  REQUIRE(initial_data.TSI_rhop > 0.0);
  REQUIRE(initial_data.TSI_Tprt > 0.0);
  REQUIRE(initial_data.TSI_max_f_over_M > 0.0);
  REQUIRE(std::isfinite(initial_data.TSI_rhof));
  REQUIRE(std::isfinite(initial_data.TSI_rhop));
  REQUIRE(std::isfinite(initial_data.TSI_Tprt));
  REQUIRE(std::isfinite(initial_data.TSI_max_f_over_M));
  REQUIRE(initial_data.TSI_m21 + initial_data.TSI_m22 +
              initial_data.TSI_m23 ==
          Catch::Approx(3.0 * initial_data.TSI_rhof *
                        initial_data.TSI_Tprt));
}
