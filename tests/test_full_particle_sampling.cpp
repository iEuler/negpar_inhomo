#include <cmath>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "FullParticleSampling.h"
#include "ParticleGroup.h"
#include "RandomContext.h"
#include "RandomSampling.h"

namespace {

coulomb::NeParticleGroup fullParticleFixture() {
	coulomb::NeParticleGroup particles;
	int positionIndex = 0;
	for (const double vx : {-1.0, 1.0}) {
		for (const double vy : {-1.0, 1.0}) {
			for (const double vz : {-1.0, 1.0}) {
				const double position = -0.9 + 0.2 * positionIndex++;
				const coulomb::Particle1D3D particle(position, {vx, vy, vz});
				particles.pushBack(particle, coulomb::ParticleKind::Positive);
				particles.pushBack(particle, coulomb::ParticleKind::Negative);
			}
		}
	}
	particles.rhoM = 1.0;
	particles.u1M = 0.0;
	particles.u2M = 0.0;
	particles.u3M = 0.0;
	particles.tprtM = 0.25;
	return particles;
}

void requireSameList(const coulomb::NeParticleGroup& first,
					 const coulomb::NeParticleGroup& second,
					 coulomb::ParticleKind kind) {
	REQUIRE(first.size(kind) == second.size(kind));
	for (int index = 0; index < first.size(kind); ++index) {
		const auto& lhs = first.list(index, kind);
		const auto& rhs = second.list(index, kind);
		REQUIRE(lhs.position() == rhs.position());
		REQUIRE(lhs.velocity() == rhs.velocity());
	}
}

} // namespace

TEST_CASE("negpar.unit.resampling.full-particle reconstruction replays exactly "
		  "for a fixed seed",
		  "[resampling][full-particle][sampling]") {
	auto firstInput = fullParticleFixture();
	auto secondInput = firstInput;
	const auto original = firstInput;

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(20260809);
	secondRandom.reseed(20260809);

	auto first = coulomb::FullParticleSampling{}.resample(
		firstInput, 8, 0.1, 0.005, 1.0, firstRandom);
	auto second = coulomb::FullParticleSampling{}.resample(
		secondInput, 8, 0.1, 0.005, 1.0, secondRandom);

	REQUIRE(first.size(coulomb::ParticleKind::Positive) == 0);
	REQUIRE(first.size(coulomb::ParticleKind::Negative) == 0);
	REQUIRE(first.size(coulomb::ParticleKind::Full) > 0);
	REQUIRE(first.size(coulomb::ParticleKind::Full) ==
			second.size(coulomb::ParticleKind::Full));
	requireSameList(first, second, coulomb::ParticleKind::Full);

	requireSameList(firstInput, original, coulomb::ParticleKind::Positive);
	requireSameList(firstInput, original, coulomb::ParticleKind::Negative);
	requireSameList(firstInput, original, coulomb::ParticleKind::Full);

	for (const auto& particle : first.list(coulomb::ParticleKind::Full)) {
		REQUIRE(std::isfinite(particle.position()));
		REQUIRE(particle.position() == 0.0);
		for (int component = 0; component < 3; ++component) {
			const double velocity = particle.velocity(component);
			REQUIRE(std::isfinite(velocity));
			REQUIRE(velocity >= firstInput.xyzMinMax[2 * component]);
			REQUIRE(velocity <= firstInput.xyzMinMax[2 * component + 1]);
		}
	}

	first.computeMoments();
	const double sampledMass = 0.005 * first.fullMoments.m0;
	const double sampledMomentum = 0.005 * first.fullMoments.m11;
	REQUIRE(sampledMass == Catch::Approx(1.0).margin(0.25));
	REQUIRE(sampledMomentum == Catch::Approx(0.0).margin(0.15));
}

TEST_CASE("full reconstruction retains Maxwellian mass with no signed particles", "[resampling][full-particle][research]") {
	coulomb::NeParticleGroup source;
	source.rhoM = 1.0;
	source.tprtM = 1.0;
	source.u1M = 0.0;
	source.u2M = 0.0;
	source.u3M = 0.0;
	coulomb::RandomContext random;
	random.reseed(20000);
	auto result = coulomb::FullParticleSampling{}.resample(source, 8, 0.001, 0.001, 1.0, random);
	result.computeMoments();
	REQUIRE(0.001 * result.fullMoments.m0 == Catch::Approx(1.0).margin(0.15));
	REQUIRE(source.xyzMinMax[0] <= -6.0);
	REQUIRE(source.xyzMinMax[1] >= 6.0);
}
