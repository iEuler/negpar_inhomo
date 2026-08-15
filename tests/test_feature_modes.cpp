#include <cmath>

#include <catch2/catch_test_macros.hpp>

#include "Grid.h"
#include "NegativeParticleCollisions.h"
#include "ParticleGroup.h"
#include "ProjectionSampling.h"
#include "RandomContext.h"
#include "SimulationConfig.h"

namespace {

coulomb::Particle1D3D particle(double x, double y, double z) {
	return coulomb::Particle1D3D({x, y, z});
}

coulomb::NeParticleGroup collisionFixture() {
	coulomb::NeParticleGroup group;
	group.u1M = 0.2;
	group.u2M = -0.1;
	group.u3M = 0.3;
	group.tprtM = 1.1;
	group.rhoM = 1.0;
	group.rho = 1.0;
	group.alphaNeg = 0.0;
	group.alphaPos = 0.0;
	group.rMax = 0.0;
	group.pushBack(particle(1.0, 0.2, -0.5), coulomb::ParticleKind::Positive);
	group.pushBack(particle(-0.5, 1.2, 0.7), coulomb::ParticleKind::Negative);
	group.pushBack(particle(0.0, 0.0, 0.0), coulomb::ParticleKind::Full);
	group.pushBack(particle(0.5, -1.0, 1.0), coulomb::ParticleKind::Full);
	return group;
}

} // namespace

TEST_CASE("linearized Coulomb uses Maxwellian partners without changing full "
		  "particles",
		  "[collisions][linearized]") {
	auto group = collisionFixture();
	const auto fullBefore = group.list(coulomb::ParticleKind::Full);
	coulomb::ParaClass parameters;
	parameters.collisionCoupling = coulomb::CollisionCoupling::Linearized;
	coulomb::NumericGridClass grid(1);
	coulomb::RandomContext random;
	random.reseed(1234);
	coulomb::NegativeParticleCollisions(grid, parameters, random)
		.collideWithFull(group);
	REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(group.size(coulomb::ParticleKind::Negative) == 1);
	for (std::size_t index = 0; index < fullBefore.size(); ++index)
		REQUIRE(group.list(static_cast<int>(index), coulomb::ParticleKind::Full)
					.velocity() == fullBefore[index].velocity());
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative})
		for (const auto& one : group.list(kind))
			for (int component = 0; component < 3; ++component)
				REQUIRE(std::isfinite(one.velocity(component)));
}

TEST_CASE("disabling Delta-M leaves collision-source counts unchanged",
		  "[collisions][delta-m]") {
	auto group = collisionFixture();
	coulomb::ParaClass parameters;
	parameters.deltaMMode = coulomb::DeltaMMode::Disabled;
	coulomb::NumericGridClass grid(1);
	grid.neff = 0.1;
	coulomb::RandomContext random;
	random.reseed(5678);
	coulomb::NegativeParticleCollisions(grid, parameters, random)
		.collideHomogeneous(group);
	REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(group.size(coulomb::ParticleKind::Negative) == 1);
}

TEST_CASE("Maxwellian-only projection is an explicit selectable mode",
		  "[projection][mode]") {
	coulomb::NeParticleGroup group;
	group.rhoM = 1.0;
	group.tprtM = 1.0;
	group.setXRange(-1.0, 1.0);
	coulomb::NumericGridClass grid(1);
	coulomb::RandomContext random;
	random.reseed(9012);
	coulomb::ProjectionSampling{}.sampleHomogeneous(
		group, grid, random, coulomb::ProjectionMode::MaxwellianOnly);
	REQUIRE(group.size(coulomb::ParticleKind::Positive) == 0);
	REQUIRE(group.size(coulomb::ParticleKind::Negative) == 0);
}
