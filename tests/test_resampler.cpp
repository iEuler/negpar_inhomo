#include <cmath>
#include <stdexcept>
#include <type_traits>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Constants.h"
#include "RandomContext.h"
#include "RandomSampling.h"
#include "Resampler.h"
#include "ResamplerHelper.h"
#include "ResamplingNumerics.h"

namespace { // Fourier resampler fixtures

coulomb::NeParticleGroup signed_fixture() {
	coulomb::NeParticleGroup particles;
	for (int component = 0; component < 3; ++component) {
		for (const double sign : {-1.0, 1.0}) {
			std::vector<double> velocity(3, 0.0);
			velocity[component] = sign;
			const coulomb::Particle1d3d anchor(velocity);
			particles.push_back(anchor, coulomb::ParticleKind::Positive);
			particles.push_back(anchor, coulomb::ParticleKind::Negative);
		}
	}

	for (int index = 0; index < 24; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 24.0;
		particles.push_back(
			coulomb::Particle1d3d({0.5 * std::cos(angle), 0.5 * std::sin(angle),
								   0.25 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Positive);
	}
	for (int index = 0; index < 8; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 8.0;
		particles.push_back(
			coulomb::Particle1d3d({0.2 * std::cos(angle), 0.2 * std::sin(angle),
								   0.1 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Negative);
	}
	return particles;
}

void require_same_particles(const coulomb::NeParticleGroup& first,
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

} // namespace

TEST_CASE("negpar.unit.resampling.Fourier resampler uses explicit RNG",
		  "[resampling][fourier]") {
	static_assert(
		std::is_constructible_v<coulomb::resampling::FourierResampler,
								const coulomb::NeParticleGroup&,
								coulomb::resampling::FourierResamplerConfig>);

	coulomb::NeParticleGroup first;
	coulomb::NeParticleGroup second;
	coulomb::RandomContext first_random;
	coulomb::RandomContext second_random;
	first_random.reseed(9876);
	second_random.reseed(9876);

	double first_bound = 1.0;
	double second_bound = 1.0;
	const std::vector<double> sample{coulomb::pi, coulomb::pi, coulomb::pi};

	coulomb::resampling::ResamplerHelper(first_random)
		.accept_sample(sample, first, 0.25, first_bound);
	coulomb::resampling::ResamplerHelper(second_random)
		.accept_sample(sample, second, 0.25, second_bound);

	REQUIRE(first_bound == second_bound);
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

	config.effective_particle_weight = 0.0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.effective_particle_weight = 1.0;
	config.frequency_count = 3;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.frequency_count = 4;
	config.max_sampling_attempts = 0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);
}

TEST_CASE(
	"negpar.unit.resampling.Fourier interpolation uses every first derivative",
	"[resampling][fourier][interpolation]") {
	const std::vector<double> derivatives{1.0, 2.0, 3.0, 5.0, 0.0,
										  0.0, 0.0, 0.0, 0.0, 0.0};
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.evaluate_quadratic_taylor(
				0.1, 0.2, 0.4, derivatives) == Catch::Approx(3.8));
}

TEST_CASE("negpar.unit.resampling.Fourier upper bounds cover their own "
		  "interpolation cell",
		  "[resampling][fourier][bounds]") {
	coulomb::Vector3D values(2, std::vector(2, std::vector<double>(2, 0.0)));
	values[0][0][0] = 7.0;

	coulomb::RandomContext random;
	const auto bounds =
		coulomb::resampling::ResamplerHelper(random).upper_bound(values);
	for (const auto& plane : bounds)
		for (const auto& row : plane)
			for (const double bound : row)
				REQUIRE(bound == Catch::Approx(7.0));
}

TEST_CASE(
	"negpar.unit.resampling.Fourier signed resampling conserves low moments",
	"[resampling][fourier][conservation]") {
	auto particles = signed_fixture();
	auto original = particles;
	original.computemoments();

	coulomb::resampling::FourierResamplerConfig config;
	config.effective_particle_weight = 0.1;
	config.frequency_count = 4;
	config.max_sampling_attempts = 10'000;

	SECTION("exact Fourier transform") { config.use_approximation = false; }
	SECTION("approximate Fourier transform") {
		config.use_approximation = true;
	}

	coulomb::RandomContext first_random;
	coulomb::RandomContext second_random;
	first_random.reseed(20260808);
	second_random.reseed(20260808);

	const coulomb::resampling::FourierResampler first_resampler(particles,
																config);
	const coulomb::resampling::FourierResampler second_resampler(particles,
																 config);
	auto first = first_resampler.resample(first_random);
	const auto second = second_resampler.resample(second_random);

	require_same_particles(first, second, coulomb::ParticleKind::Positive);
	require_same_particles(first, second, coulomb::ParticleKind::Negative);
	REQUIRE(particles.xyz_minmax ==
			std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

	first.computemoments();
	const double original_mass =
		config.effective_particle_weight *
		(original.positive_moments.m0 - original.negative_moments.m0);
	const double sampled_mass =
		config.effective_particle_weight *
		(first.positive_moments.m0 - first.negative_moments.m0);
	const double original_energy =
		config.effective_particle_weight *
		(original.positive_moments.m2 - original.negative_moments.m2);
	const double sampled_energy =
		config.effective_particle_weight *
		(first.positive_moments.m2 - first.negative_moments.m2);

	CAPTURE(original_mass, sampled_mass, original_energy, sampled_energy,
			first.positive_moments.m0, first.negative_moments.m0);
	REQUIRE(sampled_mass == Catch::Approx(original_mass).margin(0.8));
	REQUIRE(sampled_energy == Catch::Approx(original_energy).margin(1.0));

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
