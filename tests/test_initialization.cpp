#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Grid.h"
#include "InitialConditions.h"
#include "ParticleGroup.h"
#include "ParticleInitialization.h"
#include "RandomContext.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

TEST_CASE("negpar.unit.initialization.selected initial conditions preserve "
		  "Landau defaults",
		  "[initialization][conditions]") {
	coulomb::NumericGridClass grid;
	auto initialData = coulomb::InitialConditions{}.create(grid);

	REQUIRE(initialData.problemName == "LandauDamping");
	REQUIRE(initialData.ldAlpha == Catch::Approx(0.4));

	coulomb::InitialConditions{}.configure(initialData, grid, 0);
	REQUIRE(initialData.rho ==
			Catch::Approx(1.0 + initialData.ldAlpha * std::sin(grid.x[0])));
	REQUIRE(initialData.velocity[0] == 0.0);
	REQUIRE(initialData.velocity[1] == 0.0);
	REQUIRE(initialData.velocity[2] == 0.0);
	REQUIRE(initialData.tprt == 1.0);
}

TEST_CASE("negpar.unit.initialization.particle initialization replays a small "
		  "fixed-seed fixture",
		  "[initialization][particles]") {
	coulomb::IniValClass initialData;
	initialData.problemName = "Delta";
	initialData.rho = 1.0;
	initialData.velocity[0] = 1.0;
	initialData.velocity[1] = -0.5;
	initialData.velocity[2] = 0.25;
	initialData.tprt = 1.0;

	coulomb::NeParticleGroup first;
	coulomb::NeParticleGroup second;
	first.setXRange(-0.5, 0.5);
	second.setXRange(-0.5, 0.5);

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(1401);
	secondRandom.reseed(1401);

	coulomb::ParticleInitialization{}.initialize(first, initialData, 0.1, 0.1,
												 1.0, firstRandom);
	coulomb::ParticleInitialization{}.initialize(second, initialData, 0.1, 0.1,
												 1.0, secondRandom);

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
		REQUIRE((lhs.velocity() == std::vector<double>{1.0, -0.5, 0.25}));
	}

	REQUIRE(first.fullMoments.m0 == Catch::Approx(10.0));
	REQUIRE(first.fullMoments.m11 == Catch::Approx(10.0));
	REQUIRE(first.fullMoments.m12 == Catch::Approx(-5.0));
	REQUIRE(first.fullMoments.m13 == Catch::Approx(2.5));
}

TEST_CASE("negpar.unit.initialization.Two-Stream preprocessing produces finite "
		  "positive coefficients",
		  "[initialization][conditions][two-stream]") {
	coulomb::IniValClass initialData;
	coulomb::InitialConditions{}.configureTwoStream(initialData);

	REQUIRE(initialData.tsiRhof > 0.0);
	REQUIRE(initialData.tsiRhop > 0.0);
	REQUIRE(initialData.tsiTprt > 0.0);
	REQUIRE(initialData.tsiMaxFOverM > 0.0);
	REQUIRE(std::isfinite(initialData.tsiRhof));
	REQUIRE(std::isfinite(initialData.tsiRhop));
	REQUIRE(std::isfinite(initialData.tsiTprt));
	REQUIRE(std::isfinite(initialData.tsiMaxFOverM));
	REQUIRE(initialData.tsiM21 + initialData.tsiM22 + initialData.tsiM23 ==
			Catch::Approx(3.0 * initialData.tsiRhof * initialData.tsiTprt));
}
