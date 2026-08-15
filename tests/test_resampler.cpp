#include <cmath>
#include <stdexcept>
#include <type_traits>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Constants.h"
#include "EffectiveWeightSelector.h"
#include "ParticlePartition.h"
#include "RandomContext.h"
#include "RandomSampling.h"
#include "Resampler.h"
#include "ResamplerHelper.h"
#include "ResamplingNumerics.h"
#include "WeightedFourierCoupling.h"

namespace { // Fourier resampler fixtures

coulomb::NeParticleGroup signedFixture() {
	coulomb::NeParticleGroup particles;
	for (int component = 0; component < 3; ++component) {
		for (const double sign : {-1.0, 1.0}) {
			std::vector<double> velocity(3, 0.0);
			velocity[component] = sign;
			const coulomb::Particle1D3D anchor(velocity);
			particles.pushBack(anchor, coulomb::ParticleKind::Positive);
			particles.pushBack(anchor, coulomb::ParticleKind::Negative);
		}
	}

	for (int index = 0; index < 24; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 24.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.5 * std::cos(angle), 0.5 * std::sin(angle),
								   0.25 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Positive);
	}
	for (int index = 0; index < 8; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 8.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.2 * std::cos(angle), 0.2 * std::sin(angle),
								   0.1 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Negative);
	}
	return particles;
}

void requireSameParticles(const coulomb::NeParticleGroup& first,
						  const coulomb::NeParticleGroup& second,
						  coulomb::ParticleKind kind) {
	REQUIRE(first.size(kind) == second.size(kind));
	for (int index = 0; index < first.size(kind); ++index) {
		for (int component = 0; component < 3; ++component) {
			REQUIRE(first.list(index, kind).velocity(component) ==
					second.list(index, kind).velocity(component));
		}
	}
}

TEST_CASE("negpar.unit.resampling.weighted Fourier reconstruction replays",
		  "[resampling][fourier][weighted][reproducibility]") {
	auto particles = signedFixture();
	for (int index = 0; index < 12; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 12.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.35 * std::cos(angle),
								   0.35 * std::sin(angle), 0.1}),
			coulomb::ParticleKind::Full);
	}
	coulomb::resampling::FourierResamplerConfig config;
	config.effectiveParticleWeight = 0.1;
	config.fullParticleWeight = 0.2;
	config.weightedCoupling = true;
	config.frequencyCount = 4;
	config.maxSamplingAttempts = 10'000;

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(440011);
	secondRandom.reseed(440011);
	const auto originalBounds = particles.xyzMinMax;
	const auto first =
		coulomb::resampling::FourierResampler(particles, config).resample(
			firstRandom);
	const auto second =
		coulomb::resampling::FourierResampler(particles, config).resample(
			secondRandom);
	requireSameParticles(first, second, coulomb::ParticleKind::Positive);
	requireSameParticles(first, second, coulomb::ParticleKind::Negative);
	REQUIRE(particles.xyzMinMax == originalBounds);
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		for (const auto& particle : first.list(kind))
			for (int component = 0; component < 3; ++component)
				REQUIRE(std::isfinite(particle.velocity(component)));
	}
}

} // namespace

TEST_CASE("negpar.unit.resampling.Fourier resampler uses explicit RNG",
		  "[resampling][fourier]") {
	static_assert(
		std::is_constructible_v<coulomb::resampling::FourierResampler,
								const coulomb::NeParticleGroup&,
								coulomb::resampling::FourierResamplerConfig>);

	coulomb::NeParticleGroup first;
	coulomb::NeParticleGroup second;
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(9876);
	secondRandom.reseed(9876);

	double firstBound = 1.0;
	double secondBound = 1.0;
	const std::vector<double> sample{coulomb::pi, coulomb::pi, coulomb::pi};

	coulomb::resampling::ResamplerHelper(firstRandom)
		.acceptSample(sample, first, 0.25, firstBound);
	coulomb::resampling::ResamplerHelper(secondRandom)
		.acceptSample(sample, second, 0.25, secondBound);

	REQUIRE(firstBound == secondBound);
	REQUIRE(first.size(coulomb::ParticleKind::Positive) ==
			second.size(coulomb::ParticleKind::Positive));
	REQUIRE(first.size(coulomb::ParticleKind::Negative) ==
			second.size(coulomb::ParticleKind::Negative));
}

TEST_CASE("negpar.unit.resampling.Fourier resampler configuration rejects "
		  "invalid grids",
		  "[resampling][fourier][validation]") {
	coulomb::NeParticleGroup particles;
	coulomb::resampling::FourierResamplerConfig config;

	config.effectiveParticleWeight = 0.0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.effectiveParticleWeight = 1.0;
	config.frequencyCount = 3;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.frequencyCount = 4;
	config.maxSamplingAttempts = 0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.resampling.weighted Fourier coupling clamps finite "
		  "frequency weights",
		  "[resampling][fourier][weighted]") {
	const std::complex<double> full{0.1, 0.2};
	const std::complex<double> positive{0.4, -0.1};
	const std::complex<double> negative{-0.2, 0.3};
	const double weight = coulomb::resampling::WeightedFourierCoupling::
		optimalWeight(full, positive, negative, 0.2, 0.1, 12, 20, 8);
	REQUIRE(std::isfinite(weight));
	REQUIRE(weight >= 0.0);
	REQUIRE(weight <= 1.0);
	const auto blended = coulomb::resampling::WeightedFourierCoupling::blend(
		full, positive - negative, weight);
	REQUIRE(std::isfinite(blended.real()));
	REQUIRE(std::isfinite(blended.imag()));
	REQUIRE(coulomb::resampling::WeightedFourierCoupling::blend(
				full, positive - negative, 0.0) == positive - negative);
	REQUIRE(coulomb::resampling::WeightedFourierCoupling::blend(
				full, positive - negative, 1.0) == full);
}

TEST_CASE("negpar.unit.resampling.partial partition keeps typed core and tail",
		  "[resampling][partial]") {
	coulomb::NeParticleGroup source;
	source.u1M = 0.0;
	source.u2M = 0.0;
	source.u3M = 0.0;
	source.tprtM = 1.0;
	source.pushBack(coulomb::Particle1D3D(0.25, {0.5, 0.0, 0.0}),
				   coulomb::ParticleKind::Positive);
	source.pushBack(coulomb::Particle1D3D(0.75, {4.0, 0.0, 0.0}),
				   coulomb::ParticleKind::Positive);
	source.pushBack(coulomb::Particle1D3D(0.5, {0.0, 0.0, 0.0}),
				   coulomb::ParticleKind::Full);

	const auto partition = coulomb::ParticlePartitioning::split(source, 3.0);
	REQUIRE(partition.core.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(partition.tail.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(partition.core.size(coulomb::ParticleKind::Full) == 1);
	REQUIRE(partition.tail.size(coulomb::ParticleKind::Full) == 0);
	REQUIRE(partition.tail.list(0, coulomb::ParticleKind::Positive).position() ==
			Catch::Approx(0.75));
}

TEST_CASE("negpar.unit.resampling.adaptive effective weights stay bounded "
		  "for degenerate cells",
		  "[resampling][adaptive]") {
	const std::vector<coulomb::EffectiveWeightCell> cells{
		{0.0, 1.0, 0, 0}, {1.0, 1.0, 20, 1}};
	const auto selection = coulomb::EffectiveWeightSelector{}.select(
		0.001, 0.002, 0.0001, 0.005, 0.0001, 0.005, 0.205, 3.277, 0.01,
		10.0, cells);
	REQUIRE(std::isfinite(selection.signedWeight));
	REQUIRE(std::isfinite(selection.fullWeight));
	REQUIRE(selection.signedWeight >= 0.0001);
	REQUIRE(selection.signedWeight <= 0.005);
	REQUIRE(selection.fullWeight >= 0.0001);
	REQUIRE(selection.fullWeight <= 0.005);
}

TEST_CASE(
	"negpar.unit.resampling.Fourier interpolation uses every first derivative",
	"[resampling][fourier][interpolation]") {
	const std::vector<double> derivatives{1.0, 2.0, 3.0, 5.0, 0.0,
										  0.0, 0.0, 0.0, 0.0, 0.0};
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.evaluateQuadraticTaylor(
				0.1, 0.2, 0.4, derivatives) == Catch::Approx(3.8));
}

TEST_CASE("negpar.unit.resampling.Fourier upper bounds cover their own "
		  "interpolation cell",
		  "[resampling][fourier][bounds]") {
	coulomb::Vector3D values(2, std::vector(2, std::vector<double>(2, 0.0)));
	values[0][0][0] = 7.0;

	coulomb::RandomContext random;
	const auto bounds =
		coulomb::resampling::ResamplerHelper(random).upperBound(values);
	for (const auto& plane : bounds)
		for (const auto& row : plane)
			for (const double bound : row)
				REQUIRE(bound == Catch::Approx(7.0));
}

TEST_CASE(
	"negpar.unit.resampling.Fourier signed resampling conserves low moments",
	"[resampling][fourier][conservation]") {
	auto particles = signedFixture();
	auto original = particles;
	original.computeMoments();

	coulomb::resampling::FourierResamplerConfig config;
	config.effectiveParticleWeight = 0.1;
	config.frequencyCount = 4;
	config.maxSamplingAttempts = 10'000;

	SECTION("exact Fourier transform") { config.useApproximation = false; }
	SECTION("approximate Fourier transform") { config.useApproximation = true; }

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(20260808);
	secondRandom.reseed(20260808);

	const coulomb::resampling::FourierResampler firstResampler(particles,
															   config);
	const coulomb::resampling::FourierResampler secondResampler(particles,
																config);
	auto first = firstResampler.resample(firstRandom);
	const auto second = secondResampler.resample(secondRandom);

	requireSameParticles(first, second, coulomb::ParticleKind::Positive);
	requireSameParticles(first, second, coulomb::ParticleKind::Negative);
	REQUIRE(particles.xyzMinMax ==
			std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

	first.computeMoments();
	const double originalMass =
		config.effectiveParticleWeight *
		(original.positiveMoments.m0 - original.negativeMoments.m0);
	const double sampledMass =
		config.effectiveParticleWeight *
		(first.positiveMoments.m0 - first.negativeMoments.m0);
	const double originalEnergy =
		config.effectiveParticleWeight *
		(original.positiveMoments.m2 - original.negativeMoments.m2);
	const double sampledEnergy =
		config.effectiveParticleWeight *
		(first.positiveMoments.m2 - first.negativeMoments.m2);

	CAPTURE(originalMass, sampledMass, originalEnergy, sampledEnergy,
			first.positiveMoments.m0, first.negativeMoments.m0);
	REQUIRE(sampledMass == Catch::Approx(originalMass).margin(0.8));
	REQUIRE(sampledEnergy == Catch::Approx(originalEnergy).margin(1.0));

	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		for (const auto& particle : first.list(kind)) {
			for (int component = 0; component < 3; ++component) {
				REQUIRE(std::isfinite(particle.velocity(component)));
				REQUIRE(std::abs(particle.velocity(component)) <= 1.000001);
			}
		}
	}
}
