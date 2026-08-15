#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationTypes.h"
#include "WeightedHdpCoupling.h"

TEST_CASE("weighted HDP uses the deviational variance coefficient",
		  "[weighted-hdp]") {
	coulomb::NumericGridClass grid(1, coulomb::SimulationMethod::HDP);
	grid.dx = 1.0;
	grid.neff = 2.0;
	grid.neffF = 3.0;

	coulomb::NeParticleGroup group;
	group.rhoM = 4.0;
	group.u1M = 0.0;
	group.u2M = 0.0;
	group.u3M = 0.0;
	group.tprtM = 1.0;
	group.positiveMoments.m0 = 2.0;
	group.negativeMoments.m0 = 1.0;
	group.fullMoments.m0 = 1.0;
	group.pushBack(coulomb::Particle1D3D{}, coulomb::ParticleKind::Positive);
	group.pushBack(coulomb::Particle1D3D{}, coulomb::ParticleKind::Negative);
	group.pushBack(coulomb::Particle1D3D{}, coulomb::ParticleKind::Full);

	const double expectedWeight = 8.0 / 17.0;
	const double expectedSignedDensity = 6.0;
	const double expectedFullDensity = 3.0;
	const double expectedDensity =
		(1.0 - expectedWeight) * expectedSignedDensity +
		expectedWeight * expectedFullDensity;

	REQUIRE(coulomb::WeightedHdpCoupling::blendWeight(group, grid) ==
			Catch::Approx(expectedWeight));
	REQUIRE(coulomb::WeightedHdpCoupling::density(group, grid) ==
			Catch::Approx(expectedDensity));
}

TEST_CASE("weighted HDP falls back to the signed estimate for an empty cell",
		  "[weighted-hdp]") {
	coulomb::NumericGridClass grid(1, coulomb::SimulationMethod::HDP);
	coulomb::NeParticleGroup group;
	group.rhoM = 2.0;
	group.tprtM = 1.0;

	REQUIRE(coulomb::WeightedHdpCoupling::blendWeight(group, grid) == 0.0);
	REQUIRE(coulomb::WeightedHdpCoupling::density(group, grid) ==
			Catch::Approx(2.0));
}
