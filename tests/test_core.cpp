#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "Advection.h"
#include "Collisions.h"
#include "Constants.h"
#include "Diagnostics.h"
#include "FFT.h"
#include "Grid.h"
#include "Histogram.h"
#include "MacroOutput.h"
#include "Moments.h"
#include "NegativeParticleCollisions.h"
#include "NegativeParticleSampling.h"
#include "Numerics.h"
#include "OutputPaths.h"
#include "Particle.h"
#include "ParticleConservation.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "ParticleOutput.h"
#include "ParticleResampling.h"
#include "RandomContext.h"
#include "RandomSampling.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"
#include "RunOptions.h"
#include "SimulationConfig.h"
#include "SimulationState.h"
#include "SimulationTypes.h"
#include "TensorTypes.h"

namespace {

using coulomb::Particle1D3D;

Particle1D3D particle(double x, double vx, double vy, double vz) {
	return Particle1D3D(x, {vx, vy, vz});
}

} // namespace

TEST_CASE("negpar.unit.grid.default grid is initialized", "[grid]") {
	coulomb::NumericGridClass grid;

	REQUIRE(grid.nx == 100);
	REQUIRE(grid.nv == 200);
	REQUIRE(grid.x.size() == static_cast<size_t>(grid.nx));
	REQUIRE(grid.vx.size() == static_cast<size_t>(grid.nv));
	REQUIRE(grid.dx > 0.0);
	REQUIRE(grid.dv > 0.0);
	REQUIRE(grid.dt > 0.0);
}

TEST_CASE("negpar.unit.grid.grid rejects invalid configuration",
		  "[grid][validation]") {
	REQUIRE_THROWS_AS(coulomb::NumericGridClass(0), std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::NumericGridClass(-4), std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::NumericGridClass(
						  10, static_cast<coulomb::SimulationMethod>(99)),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.methodName(
						  static_cast<coulomb::SimulationMethod>(99)),
					  std::invalid_argument);
}

TEST_CASE(
	"negpar.unit.parameters.typed collision mode preserves the legacy name",
	"[parameters]") {
	const coulomb::ParaClass parameters;
	const coulomb::IniValClass initialData;
	coulomb::RandomContext random;
	random.reseed(2468);
	REQUIRE(parameters.methodBinaryColl == coulomb::BinaryCollisionMethod::TA);
	REQUIRE(initialData.rho == 0.0);
	REQUIRE(initialData.velocity[0] == 0.0);
	REQUIRE(std::string(coulomb::SimulationTypes{}.binaryCollisionName(
				parameters.methodBinaryColl)) == "TA");
	REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.collisionName(
						  static_cast<coulomb::CollisionType>(99)),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.binaryCollisionName(
						  static_cast<coulomb::BinaryCollisionMethod>(99)),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.advection.particle advection wraps periodic boundaries",
		  "[advection][boundary]") {
	coulomb::NumericGridClass grid(4);
	coulomb::SimulationState state;
	grid.bdryX = coulomb::BoundaryCondition::Periodic;
	grid.dt = 1.0;
	const double start = grid.xmax - 0.25;
	const auto moved = coulomb::Advection(grid, state)
						   .moveParticle(particle(start, 1.0, 2.0, 3.0), 0.0);

	REQUIRE(moved.position() == Catch::Approx(grid.xmin + 0.75));
	REQUIRE(moved.velocity(0) == Catch::Approx(1.0));
	REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("negpar.unit.advection.particle advection reflects at nonperiodic "
		  "boundaries",
		  "[advection][boundary]") {
	coulomb::NumericGridClass grid(4);
	coulomb::SimulationState state;
	grid.bdryX = coulomb::BoundaryCondition::Reflective;
	grid.dt = 1.0;
	const auto moved =
		coulomb::Advection(grid, state)
			.moveParticle(particle(grid.xmax - 0.25, 1.0, 2.0, 3.0), 0.5);

	REQUIRE(moved.position() == Catch::Approx(grid.xmax - 0.75));
	REQUIRE(moved.velocity(0) == Catch::Approx(-1.5));
	REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("negpar.unit.advection.particle exactly at the upper boundary "
		  "belongs to the last cell",
		  "[advection][boundary]") {
	coulomb::NumericGridClass grid(4);
	coulomb::SimulationState state;
	grid.bdryX = coulomb::BoundaryCondition::Reflective;
	grid.dt = grid.xmax - (grid.xmax - 0.5);
	coulomb::Advection advection(grid, state);
	auto moved =
		advection.moveParticle(particle(grid.xmax - 0.5, 1.0, 0.0, 0.0), 0.0);

	REQUIRE(moved.position() == Catch::Approx(grid.xmax));
	REQUIRE(advection.findParticleGroup(moved) == grid.nx - 1);
}

TEST_CASE("negpar.unit.advection.particle group lookup rejects positions "
		  "outside the grid",
		  "[advection][validation]") {
	coulomb::NumericGridClass grid(4);
	coulomb::SimulationState state;
	coulomb::Advection advection(grid, state);
	auto above = particle(grid.xmax + 1.0, 0.0, 0.0, 0.0);
	auto below = particle(grid.xmin - 1.0, 0.0, 0.0, 0.0);

	REQUIRE_THROWS_AS(advection.findParticleGroup(above), std::out_of_range);
	REQUIRE_THROWS_AS(advection.findParticleGroup(below), std::out_of_range);
}

TEST_CASE(
	"negpar.unit.numerics.numerics error function preserves the approximation",
	"[numerics]") {
	REQUIRE(std::abs(coulomb::Numerics{}.errorFunction(0.0)) < 2e-9);
	REQUIRE(std::abs(coulomb::Numerics{}.errorFunction(1.0) - 0.84270079) <
			2e-7);
	REQUIRE(std::abs(coulomb::Numerics{}.errorFunction(-1.0) + 0.84270079) <
			2e-7);
}

TEST_CASE("negpar.unit.numerics.one-dimensional shifts preserve legacy "
		  "boundary conventions",
		  "[numerics][boundaries]") {
	const std::vector<double> values = {1.0, 2.0, 3.0};
	const std::vector<double> shifted = {2.0, 3.0, 1.0};
	const std::vector<double> extended = {-1.0, 1.0, 2.0};
	const std::vector<double> difference = {-0.5, 1.0, -0.5};
	REQUIRE(coulomb::SimulationTypes{}.decodeBoundary('p') ==
			coulomb::BoundaryCondition::Periodic);
	REQUIRE(coulomb::SimulationTypes{}.decodeBoundary('n') ==
			coulomb::BoundaryCondition::Reflective);
	REQUIRE(coulomb::SimulationTypes{}.encodeBoundary(
				coulomb::BoundaryCondition::Reflective) == 'n');
	REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.decodeBoundary('x'),
					  std::invalid_argument);
	REQUIRE(coulomb::Numerics{}.circularShift(values, 3, 1) == shifted);
	REQUIRE(coulomb::Numerics{}.circularShift(values, 3, 4) == shifted);
	REQUIRE(coulomb::Numerics{}.boundaryShift(values, 3, -1, -1.0) == extended);
	REQUIRE(coulomb::Numerics{}.centralDifference(
				values, 3, coulomb::BoundaryCondition::Periodic) == difference);
	REQUIRE(
		coulomb::Numerics{}.circularShift(std::vector<double>{}, 0, 8).empty());
	REQUIRE_THROWS_AS(coulomb::Numerics{}.circularShift(values, 2, 1),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(
		coulomb::Numerics{}.centralDifference(
			values, 3, static_cast<coulomb::BoundaryCondition>(99)),
		std::invalid_argument);
}

TEST_CASE("negpar.unit.numerics.kinetic Euler fluxes handle periodic and "
		  "constant states",
		  "[numerics][euler]") {
	const std::vector<double> values = {1.0, 2.0, 4.0};
	REQUIRE(
		coulomb::Numerics{}.limitedFlux(values, 3, 2.0, 1.0, 1, values,
										coulomb::BoundaryCondition::Periodic) ==
		std::vector<double>{-1.5, 0.5, 1.0});
	REQUIRE(
		coulomb::Numerics{}.limitedFlux(values, 3, 2.0, 1.0, -1, values,
										coulomb::BoundaryCondition::Periodic) ==
		std::vector<double>{0.5, 1.0, -1.5});

	std::vector<double> densityChange(3), momentumChange(3), energyChange(3);
	const std::vector<double> constant(3, 1.0);
	coulomb::Numerics{}.advanceKineticEuler(
		constant, constant, constant, 3, 1.0, 0.1,
		coulomb::BoundaryCondition::Periodic, densityChange, momentumChange,
		energyChange);
	REQUIRE(densityChange == std::vector<double>{0.0, 0.0, 0.0});
	REQUIRE(momentumChange == std::vector<double>{0.0, 0.0, 0.0});
	REQUIRE(energyChange == std::vector<double>{0.0, 0.0, 0.0});
	REQUIRE_THROWS_AS(coulomb::Numerics{}.advanceKineticEuler(
						  constant, constant, std::vector<double>(3, 0.0), 3,
						  1.0, 0.1, coulomb::BoundaryCondition::Periodic,
						  densityChange, momentumChange, energyChange),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFT reshape helpers preserve row-major ordering",
		  "[fft][numerics]") {
	const std::vector<int> flat = {0, 1, 2, 3, 4, 5, 6, 7};
	const auto shaped = coulomb::TensorReshape{}.to3D(flat, 2, 2, 2);

	REQUIRE(shaped[0][0][0] == 0);
	REQUIRE(shaped[0][1][1] == 3);
	REQUIRE(shaped[1][0][0] == 4);
	REQUIRE(coulomb::TensorReshape{}.to1D(shaped) == flat);
	REQUIRE_THROWS_AS(coulomb::TensorReshape{}.to3D(flat, 3, 2, 2),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFT reshape helper rejects ragged input",
		  "[fft][validation]") {
	const std::vector<std::vector<std::vector<int>>> ragged = {
		{{1, 2}, {3}}, {{4, 5}, {6, 7}}};
	REQUIRE_THROWS_AS(coulomb::TensorReshape{}.to1D(ragged),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFTW wrappers have explicit ownership and validate "
		  "dimensions",
		  "[fft][resource]") {
	static_assert(!std::is_copy_constructible_v<coulomb::Fft1D>);
	static_assert(!std::is_copy_assignable_v<coulomb::Fft1D>);
	static_assert(!std::is_move_constructible_v<coulomb::Fft1D>);
	static_assert(!std::is_move_assignable_v<coulomb::Fft1D>);
	static_assert(!std::is_copy_constructible_v<coulomb::Fft3D>);
	static_assert(!std::is_copy_assignable_v<coulomb::Fft3D>);
	static_assert(!std::is_move_constructible_v<coulomb::Fft3D>);
	static_assert(!std::is_move_assignable_v<coulomb::Fft3D>);

	REQUIRE_THROWS_AS(coulomb::Fft1D(0), std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::Fft3D(0, 2, 2), std::invalid_argument);

	// Repeated construction and destruction exercises plan-before-buffer
	// cleanup under the sanitizer preset. FFTW's inverse is intentionally
	// unnormalized.
	for (int iteration = 0; iteration < 8; ++iteration) {
		coulomb::Fft1D fft(4);
		const std::vector<double> input{1.0, -2.0, 0.5, 3.0};
		const auto transformed = fft.fft(input);
		const auto restored = fft.ifft(transformed);
		for (std::size_t index = 0; index < input.size(); ++index) {
			REQUIRE(restored[index] ==
					Catch::Approx(4.0 * input[index]).margin(1e-12));
		}
		REQUIRE_THROWS_AS(fft.fft(std::vector<double>(3)),
						  std::invalid_argument);
	}

	const coulomb::Vector3D input3D = {{{1.0, 2.0}, {3.0, 4.0}},
									   {{-1.0, -2.0}, {-3.0, -4.0}}};
	coulomb::Fft3D fft3D(2, 2, 2);
	const auto restored3D = fft3D.ifft(fft3D.fft(input3D));
	for (std::size_t first = 0; first < 2; ++first) {
		for (std::size_t second = 0; second < 2; ++second) {
			for (std::size_t third = 0; third < 2; ++third) {
				REQUIRE(restored3D[first][second][third] ==
						Catch::Approx(8.0 * input3D[first][second][third])
							.margin(1e-12));
			}
		}
	}
}

TEST_CASE("negpar.unit.moments.particle group computes moments",
		  "[particles][moments]") {
	REQUIRE_THROWS_AS(coulomb::Particle1D3D(std::vector<double>{1.0, 2.0}),
					  std::invalid_argument);
	coulomb::Particle1D3D particleWithVelocity;
	REQUIRE_THROWS_AS(
		particleWithVelocity.setVelocity(std::vector<double>{1.0, 2.0}),
		std::invalid_argument);
	REQUIRE_THROWS_AS(particleWithVelocity.velocity(-1), std::out_of_range);
	REQUIRE_THROWS_AS(particleWithVelocity.velocity(3), std::out_of_range);

	coulomb::ParticleGroup group;
	group.pushBack(particle(0.0, 1.0, 2.0, 3.0));
	group.pushBack(particle(0.0, -1.0, 4.0, 0.0));

	group.computeMoments();

	REQUIRE(group.size() == 2);
	REQUIRE(group.moments.m0 == 2.0);
	REQUIRE(group.moments.m11 == 0.0);
	REQUIRE(group.moments.m12 == 6.0);
	REQUIRE(group.moments.m13 == 3.0);
	REQUIRE(group.moments.m2 == 31.0);

	group.erase(0);
	REQUIRE(group.size() == 1);
	REQUIRE(group.list(0).velocity(0) == Catch::Approx(-1.0));
	REQUIRE_THROWS_AS(group.list(-1), std::out_of_range);
	REQUIRE_THROWS_AS(group.list(1), std::out_of_range);
	REQUIRE_THROWS_AS(group.erase(-1), std::out_of_range);
	REQUIRE_THROWS_AS(group.erase(1), std::out_of_range);
}

TEST_CASE("negpar.unit.moments.moment conversions round trip", "[moments]") {
	double velocity[3] = {1.0, -2.0, 0.5};
	double momentum[3] = {};
	double energy = 0.0;
	coulomb::MomentOperations{}.primitiveToConserved(2.0, velocity, 3.0,
													 momentum, energy);

	REQUIRE(momentum[0] == 2.0);
	REQUIRE(momentum[1] == -4.0);
	REQUIRE(momentum[2] == 1.0);
	REQUIRE(energy == Catch::Approx(2.0 * (1.0 + 4.0 + 0.25 + 9.0) / 2.0));

	double recoveredVelocity[3] = {};
	double recoveredTemperature = 0.0;
	coulomb::MomentOperations{}.conservedToPrimitive(
		2.0, momentum, energy, recoveredVelocity, recoveredTemperature);
	for (int component = 0; component < 3; ++component)
		REQUIRE(recoveredVelocity[component] ==
				Catch::Approx(velocity[component]));
	REQUIRE(recoveredTemperature == Catch::Approx(3.0));

	REQUIRE_THROWS_AS(coulomb::MomentOperations{}.primitiveToConserved(
						  2.0, nullptr, 3.0, momentum, energy),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(
		coulomb::MomentOperations{}.conservedToPrimitive(
			2.0, nullptr, energy, recoveredVelocity, recoveredTemperature),
		std::invalid_argument);

	coulomb::ParticleGroup emptyGroup;
	double density = 0.0;
	double ignoredEnergy = 0.0;
	REQUIRE_THROWS_AS(coulomb::MomentOperations{}.particleToConserved(
						  emptyGroup, 1.0, density, nullptr, ignoredEnergy),
					  std::invalid_argument);

	std::vector<coulomb::ParticleGroup> groups(1);
	std::vector<double> output(1);
	REQUIRE_THROWS_AS(
		coulomb::MomentOperations{}.computePrimitive(
			2, groups, 1.0, output, output, output, output, output),
		std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::MomentOperations{}.primitiveToConserved(
						  2, output, output, output, output, output),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::MomentOperations{}.conservedToPrimitive(
						  2, output, output, output, output, output),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.moments.macroscopic moment updates reconstruct particle "
		  "fields",
		  "[moments][reconstruction]") {
	coulomb::NumericGridClass grid(3);
	std::vector<coulomb::NeParticleGroup> groups(3);
	REQUIRE_THROWS_AS(coulomb::MomentOperations{}.momentChange(nullptr, grid),
					  std::invalid_argument);
	for (auto& group : groups) {
		group.rhoM = 1.0;
		group.u1M = 0.0;
		group.tprtM = 1.0;
		group.computeMoments();
	}

	coulomb::MomentOperations{}.updatePrimitive(groups, grid);
	for (const auto& group : groups) {
		REQUIRE(group.rho == Catch::Approx(1.0));
		REQUIRE(group.u1 == Catch::Approx(0.0));
		REQUIRE(group.tprt == Catch::Approx(1.0));
	}

	groups[0].pushBack(particle(0.0, 2.0, 0.0, 0.0),
					   coulomb::ParticleKind::Full);
	for (auto& group : groups)
		group.computeMoments();
	coulomb::MomentOperations{}.updateFullPrimitive(groups, grid);
	REQUIRE(groups[0].rhoF == Catch::Approx(grid.neffF / grid.dx));
	REQUIRE(groups[0].u1F == Catch::Approx(2.0));
}

TEST_CASE(
	"negpar.unit.moments.Maxwellian update applies conservative moment changes",
	"[moments][maxwellian]") {
	coulomb::NumericGridClass grid(1);
	std::vector<coulomb::NeParticleGroup> groups(1);
	groups[0].rhoM = 2.0;
	groups[0].u1M = 1.0;
	groups[0].tprtM = 3.0;
	groups[0].drho = 0.5;
	groups[0].dm1 = 0.25;
	groups[0].dEnergy = 1.0;

	coulomb::MomentOperations{}.updateMaxwellian(groups, grid);

	REQUIRE(groups[0].rhoM == Catch::Approx(1.5));
	REQUIRE(groups[0].u1M == Catch::Approx(1.75 / 1.5));
	REQUIRE(groups[0].tprtM == Catch::Approx(3.5462962963));
}

TEST_CASE("negpar.unit.resampling.full-particle count rescaling preserves "
		  "effective mass",
		  "[resampling][conservation]") {
	coulomb::NumericGridClass grid(2);
	grid.neffF = 0.25;
	std::vector<coulomb::NeParticleGroup> groups(2);

	for (auto& group : groups) {
		group.rhoM = 1.0;
		for (int particleIndex = 0; particleIndex < 30; ++particleIndex)
			group.pushBack(particle(0.0, particleIndex * 0.01, 0.0, 0.0),
						   coulomb::ParticleKind::Full);
	}

	const double massBefore = groups[0].rhoM * grid.dx * groups.size();
	coulomb::RandomContext random;
	random.reseed(2468);

	// Passing zero forces the routine to perform the count-reduction path.
	coulomb::ParaClass parameters;
	coulomb::SimulationState state;
	state.random.reseed(2468);
	coulomb::ParticleResampling(grid, parameters, state)
		.resampleFullPreservingMass(groups, 0);

	const int fullCount = coulomb::Diagnostics(grid).particleCount(
		groups, grid.nx, coulomb::ParticleKind::Full);
	REQUIRE(fullCount == 50);
	REQUIRE(grid.neffF * fullCount == Catch::Approx(massBefore).margin(1e-12));
}

TEST_CASE("negpar.unit.resampling.shared resampling numerics map modes and "
		  "evaluate all derivatives",
		  "[resampling][numerics]") {
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.frequencies(4) ==
			std::vector<double>{0.0, 1.0, 2.0, -1.0});
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.augmentedLocations(
				4, 2) == std::vector<std::size_t>{0, 1, 2, 7});
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.imaginaryFrequencies(4) ==
			std::vector<std::complex<double>>{
				{0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}, {0.0, -1.0}});

	const std::vector<double> derivatives{1.0, 2.0, 3.0, 5.0, 0.0,
										  0.0, 0.0, 0.0, 0.0, 0.0};
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.evaluateQuadraticTaylor(
				0.1, 0.2, 0.4, derivatives) == Catch::Approx(3.8));
	REQUIRE_THROWS_AS(coulomb::resampling::ResamplingNumerics{}.frequencies(0),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(
		coulomb::resampling::ResamplingNumerics{}.augmentedLocations(4, 0),
		std::invalid_argument);
	REQUIRE_THROWS_AS(
		coulomb::resampling::ResamplingNumerics{}.evaluateQuadraticTaylor(
			0.0, 0.0, 0.0, std::vector<double>(9)),
		std::invalid_argument);
}

TEST_CASE(
	"negpar.unit.resampling.shared resampling velocity transforms round trip",
	"[resampling][velocity]") {
	coulomb::NeParticleGroup particles;
	particles.rhoM = 2.0;
	particles.u1M = 0.25;
	particles.u2M = -0.5;
	particles.u3M = 0.75;
	particles.tprtM = 1.5;
	particles.xyzMinMax = {-2.0, 2.0, -3.0, 1.0, 0.0, 4.0};
	particles.pushBack(coulomb::Particle1D3D({-1.0, -2.0, 1.0}),
					   coulomb::ParticleKind::Positive);
	particles.pushBack(coulomb::Particle1D3D({1.0, 0.0, 3.0}),
					   coulomb::ParticleKind::Negative);

	auto normalized =
		coulomb::resampling::ResamplingVelocity{}.normalizeSigned(particles);
	REQUIRE(normalized.rhoM == Catch::Approx(particles.rhoM));
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		auto restored = normalized.list(kind);
		coulomb::resampling::ResamplingVelocity{}.restore(restored,
														  particles.xyzMinMax);
		for (int component = 0; component < 3; ++component)
			REQUIRE(restored.front().velocity(component) ==
					Catch::Approx(particles.list(0, kind).velocity(component)));
	}

	REQUIRE_THROWS_AS(coulomb::resampling::ResamplingVelocity{}.restore(
						  normalized.list(coulomb::ParticleKind::Positive),
						  std::vector<double>(5)),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.collisions.binary collision preserves pair momentum and "
		  "kinetic energy",
		  "[collisions][invariants]") {
	// A deterministic invariant check; particle ordering is not involved.
	coulomb::ParaClass parameters;
	parameters.dt = 0.01;
	coulomb::RandomContext random;
	random.reseed(2468);

	const std::vector<double> first = {1.0, -0.5, 0.75};
	const std::vector<double> second = {-0.25, 1.5, -1.0};
	const auto collided = coulomb::CollisionOperator(parameters, random)
							  .collidePair(first, second);

	double energyBefore = 0.0;
	double energyAfter = 0.0;
	for (int component = 0; component < 3; ++component) {
		REQUIRE(
			collided.first[component] + collided.second[component] ==
			Catch::Approx(first[component] + second[component]).margin(1e-14));
		energyBefore += first[component] * first[component] +
						second[component] * second[component];
		energyAfter += collided.first[component] * collided.first[component] +
					   collided.second[component] * collided.second[component];
	}
	REQUIRE(energyAfter == Catch::Approx(energyBefore).margin(1e-10));
}

TEST_CASE("negpar.unit.collisions.binary collision is reproducible for an "
		  "explicit seed",
		  "[collisions][reproducibility]") {
	const coulomb::ParaClass parameters;
	const std::vector<double> first = {0.5, 1.0, -0.75};
	const std::vector<double> second = {-1.0, 0.25, 1.5};
	coulomb::RandomContext random;

	random.reseed(13579);
	const auto expected = coulomb::CollisionOperator(parameters, random)
							  .collidePair(first, second);
	random.reseed(13579);
	const auto repeated = coulomb::CollisionOperator(parameters, random)
							  .collidePair(first, second);

	REQUIRE(repeated.first == expected.first);
	REQUIRE(repeated.second == expected.second);
}

TEST_CASE("negpar.unit.collisions.negative-particle collision kernels replay "
		  "for an explicit seed",
		  "[collisions][negative-particle][reproducibility]") {
	const auto makeGroup = [] {
		coulomb::NeParticleGroup group;
		group.u1M = 0.25;
		group.u2M = -0.5;
		group.u3M = 0.75;
		group.tprtM = 1.25;
		group.pushBack(particle(0.0, 1.0, 0.0, 0.0),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, -1.0, 0.5, 0.0),
					   coulomb::ParticleKind::Negative);
		group.pushBack(particle(0.0, 0.0, 1.0, 0.5),
					   coulomb::ParticleKind::Full);
		group.pushBack(particle(0.0, 0.5, -1.0, 1.0),
					   coulomb::ParticleKind::Full);
		return group;
	};
	const auto requireSameParticles =
		[](const coulomb::NeParticleGroup& first,
		   const coulomb::NeParticleGroup& second) {
			for (const auto kind : {coulomb::ParticleKind::Positive,
									coulomb::ParticleKind::Negative,
									coulomb::ParticleKind::Full}) {
				REQUIRE(first.size(kind) == second.size(kind));
				for (int index = 0; index < first.size(kind); ++index)
					REQUIRE(first.list(index, kind).velocity() ==
							second.list(index, kind).velocity());
			}
		};

	coulomb::ParaClass parameters;
	coulomb::NumericGridClass grid(1);
	auto first = makeGroup();
	auto repeated = makeGroup();
	const auto fullBefore = first.list(coulomb::ParticleKind::Full);
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext repeatedRandom;
	firstRandom.reseed(112233);
	repeatedRandom.reseed(112233);
	coulomb::NegativeParticleCollisions(grid, parameters, firstRandom)
		.collideWithFull(first);
	coulomb::NegativeParticleCollisions(grid, parameters, repeatedRandom)
		.collideWithFull(repeated);
	requireSameParticles(first, repeated);
	REQUIRE(first.size(coulomb::ParticleKind::Full) ==
			static_cast<int>(fullBefore.size()));
	for (std::size_t index = 0; index < fullBefore.size(); ++index)
		REQUIRE(first.list(static_cast<int>(index), coulomb::ParticleKind::Full)
					.velocity() == fullBefore[index].velocity());

	parameters.dt = 0.1;
	parameters.coeffBinaryColl = 10.0;
	first = makeGroup();
	repeated = makeGroup();
	firstRandom.reseed(445566);
	repeatedRandom.reseed(445566);
	coulomb::NegativeParticleCollisions(grid, parameters, firstRandom)
		.collideBgkHomogeneous(first);
	coulomb::NegativeParticleCollisions(grid, parameters, repeatedRandom)
		.collideBgkHomogeneous(repeated);
	requireSameParticles(first, repeated);
	REQUIRE(first.size(coulomb::ParticleKind::Positive) == 0);
	REQUIRE(first.size(coulomb::ParticleKind::Negative) == 0);
	for (const auto& full : first.list(coulomb::ParticleKind::Full)) {
		for (int component = 0; component < 3; ++component)
			REQUIRE(std::isfinite(full.velocity(component)));
	}
}

TEST_CASE("negpar.unit.collisions.negative-particle collision pipeline replays "
		  "for an explicit seed",
		  "[collisions][negative-particle][pipeline][reproducibility]") {
	const auto makeGroup = [] {
		coulomb::NeParticleGroup group;
		group.setXRange(-1.0, 1.0);
		group.rhoM = 1.0;
		group.u1M = 0.2;
		group.u2M = -0.1;
		group.u3M = 0.3;
		group.tprtM = 1.1;
		group.rho = 1.2;
		group.alphaNeg = 0.0;
		group.alphaPos = 0.0;
		group.rMax = 0.0;
		group.pushBack(particle(-0.8, 1.0, 0.2, -0.5),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(-0.4, -0.5, 1.2, 0.7),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, 0.3, -0.4, 0.2),
					   coulomb::ParticleKind::Negative);
		group.pushBack(particle(0.2, -1.0, 0.1, 0.4),
					   coulomb::ParticleKind::Full);
		group.pushBack(particle(0.4, 0.6, -0.8, 1.1),
					   coulomb::ParticleKind::Full);
		group.pushBack(particle(0.6, 1.3, 0.5, -0.7),
					   coulomb::ParticleKind::Full);
		group.pushBack(particle(0.8, -0.2, 1.4, 0.9),
					   coulomb::ParticleKind::Full);
		group.computeMoments();
		return group;
	};

	coulomb::ParaClass parameters;
	coulomb::NumericGridClass grid(1);
	grid.neff = 0.2;
	const auto before = makeGroup();
	auto first = before;
	auto repeated = before;
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext repeatedRandom;
	firstRandom.reseed(246810);
	repeatedRandom.reseed(246810);
	coulomb::NegativeParticleCollisions(grid, parameters, firstRandom)
		.collideHomogeneous(first);
	coulomb::NegativeParticleCollisions(grid, parameters, repeatedRandom)
		.collideHomogeneous(repeated);

	bool velocityChanged = false;
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative,
		  coulomb::ParticleKind::Full}) {
		REQUIRE(first.size(kind) == before.size(kind));
		REQUIRE(first.size(kind) == repeated.size(kind));
		for (int index = 0; index < first.size(kind); ++index) {
			REQUIRE(first.list(index, kind).position() ==
					before.list(index, kind).position());
			REQUIRE(first.list(index, kind).velocity() ==
					repeated.list(index, kind).velocity());
			velocityChanged =
				velocityChanged || first.list(index, kind).velocity() !=
									   before.list(index, kind).velocity();
			for (int component = 0; component < 3; ++component)
				REQUIRE(
					std::isfinite(first.list(index, kind).velocity(component)));
		}
	}
	REQUIRE(velocityChanged);
}

TEST_CASE("negpar.unit.negative_particles.negative-particle source sampling "
		  "replays for an explicit seed",
		  "[negative-particle][sampling][reproducibility]") {
	const auto makeGroup = [] {
		coulomb::NeParticleGroup group;
		group.rhoM = 1.0;
		group.u1M = 0.2;
		group.u2M = -0.1;
		group.u3M = 0.3;
		group.tprtM = 1.1;
		group.rho = 1.2;
		for (int repeat = 0; repeat < 100; ++repeat) {
			group.pushBack(particle(0.0, 1.0, 0.2, -0.5),
						   coulomb::ParticleKind::Positive);
			group.pushBack(particle(0.0, -0.5, 1.2, 0.7),
						   coulomb::ParticleKind::Positive);
			group.pushBack(particle(0.0, 0.8, -1.0, 1.4),
						   coulomb::ParticleKind::Positive);
			group.pushBack(particle(0.0, 0.3, -0.4, 0.2),
						   coulomb::ParticleKind::Negative);
			group.pushBack(particle(0.0, -0.7, 0.9, -1.1),
						   coulomb::ParticleKind::Negative);
		}
		group.computeMoments();
		return group;
	};

	coulomb::ParaClass parameters;
	auto source = makeGroup();
	coulomb::NegativeParticleSampling{}.updateBounds(source, parameters);
	REQUIRE(std::isfinite(source.alphaNeg));
	REQUIRE(std::isfinite(source.alphaPos));
	REQUIRE(source.alphaNeg > 0.0);
	REQUIRE(source.alphaPos > 0.0);
	REQUIRE(source.rMax > 0.0);

	const std::vector<double> mean{source.u1M, source.u2M, source.u3M};
	const double expectedAtMean =
		source.rhoM / std::pow(std::sqrt(2.0 * coulomb::pi * source.tprtM), 3);
	REQUIRE(coulomb::NegativeParticleSampling{}.evaluateMaxwellian(
				mean, source) == Catch::Approx(expectedAtMean).margin(1e-15));

	constexpr double effectiveParticles = 0.002;
	coulomb::RandomContext firstCountRandom;
	coulomb::RandomContext repeatedCountRandom;
	firstCountRandom.reseed(991122);
	repeatedCountRandom.reseed(991122);
	const int firstVirtualCount =
		coulomb::NegativeParticleSampling{}.estimateVirtualCount(
			source, effectiveParticles, firstCountRandom);
	const int repeatedVirtualCount =
		coulomb::NegativeParticleSampling{}.estimateVirtualCount(
			source, effectiveParticles, repeatedCountRandom);
	REQUIRE(firstVirtualCount == repeatedVirtualCount);
	REQUIRE(firstVirtualCount > 0);

	auto repeatedSource = source;
	coulomb::NeParticleGroup firstSample;
	coulomb::NeParticleGroup repeatedSample;
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext repeatedRandom;
	firstRandom.reseed(334455);
	repeatedRandom.reseed(334455);
	coulomb::NegativeParticleSampling{}.sampleDelta(
		source, firstSample, parameters, effectiveParticles, firstRandom);
	coulomb::NegativeParticleSampling{}.sampleDelta(
		repeatedSource, repeatedSample, parameters, effectiveParticles,
		repeatedRandom);

	int sampledParticles = 0;
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		REQUIRE(firstSample.size(kind) == repeatedSample.size(kind));
		sampledParticles += firstSample.size(kind);
		for (int index = 0; index < firstSample.size(kind); ++index) {
			REQUIRE(firstSample.list(index, kind).velocity() ==
					repeatedSample.list(index, kind).velocity());
			for (int component = 0; component < 3; ++component)
				REQUIRE(std::isfinite(
					firstSample.list(index, kind).velocity(component)));
		}
	}
	REQUIRE(sampledParticles > 0);
}

TEST_CASE("negpar.unit.negative_particles.signed particle conservation "
		  "restores target moments reproducibly",
		  "[negative-particle][conservation][reproducibility]") {
	const auto makeGroup = [] {
		coulomb::NeParticleGroup group;
		group.pushBack(particle(0.0, 1.0, 0.2, -0.5),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, -0.5, 1.2, 0.7),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, 0.8, -1.0, 1.4),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, -1.2, 0.6, -0.9),
					   coulomb::ParticleKind::Positive);
		group.pushBack(particle(0.0, 0.3, -0.4, 0.2),
					   coulomb::ParticleKind::Negative);
		group.pushBack(particle(0.0, -0.7, 0.9, -1.1),
					   coulomb::ParticleKind::Negative);
		return group;
	};
	const auto signedMoments = [](coulomb::NeParticleGroup& group,
								  double effectiveParticles) {
		group.computeMoments();
		return std::array<double, 7>{
			effectiveParticles *
				(group.positiveMoments.m0 - group.negativeMoments.m0),
			effectiveParticles *
				(group.positiveMoments.m11 - group.negativeMoments.m11),
			effectiveParticles *
				(group.positiveMoments.m12 - group.negativeMoments.m12),
			effectiveParticles *
				(group.positiveMoments.m13 - group.negativeMoments.m13),
			effectiveParticles *
				(group.positiveMoments.m21 - group.negativeMoments.m21),
			effectiveParticles *
				(group.positiveMoments.m22 - group.negativeMoments.m22),
			effectiveParticles *
				(group.positiveMoments.m23 - group.negativeMoments.m23)};
	};

	constexpr double effectiveParticles = 0.75;
	auto baseline = makeGroup();
	const auto target = signedMoments(baseline, effectiveParticles);
	auto first = makeGroup();
	auto repeated = makeGroup();
	const std::array<double, 3> perturbedVelocity{1.1, 0.1, -0.3};
	first.list(0, coulomb::ParticleKind::Positive)
		.setVelocity(perturbedVelocity);
	repeated.list(0, coulomb::ParticleKind::Positive)
		.setVelocity(perturbedVelocity);

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext repeatedRandom;
	firstRandom.reseed(778899);
	repeatedRandom.reseed(778899);
	coulomb::ParticleConservation{}.enforce(
		target[0], target[1], target[2], target[3], target[4], target[5],
		target[6], first, effectiveParticles, true, firstRandom);
	coulomb::ParticleConservation{}.enforce(
		target[0], target[1], target[2], target[3], target[4], target[5],
		target[6], repeated, effectiveParticles, true, repeatedRandom);

	const auto actual = signedMoments(first, effectiveParticles);
	for (std::size_t index = 0; index < target.size(); ++index)
		REQUIRE(actual[index] == Catch::Approx(target[index]).margin(1e-12));
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		REQUIRE(first.size(kind) == repeated.size(kind));
		for (int index = 0; index < first.size(kind); ++index)
			REQUIRE(first.list(index, kind).velocity() ==
					repeated.list(index, kind).velocity());
	}
}

TEST_CASE("negpar.unit.particles.positive, negative, and full particle groups "
		  "are independent",
		  "[particles]") {
	coulomb::NeParticleGroup group;
	REQUIRE(group.isResampled == false);
	REQUIRE(group.rhoM == 0.0);
	REQUIRE(group.positiveMoments.m0 == 0.0);
	group.pushBack(particle(0.0, 1.0, 0.0, 0.0),
				   coulomb::ParticleKind::Positive);
	group.pushBack(particle(0.0, 2.0, 0.0, 0.0),
				   coulomb::ParticleKind::Negative);
	group.pushBack(particle(0.0, 3.0, 0.0, 0.0), coulomb::ParticleKind::Full);

	REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(group.size(coulomb::ParticleKind::Negative) == 1);
	REQUIRE(group.size(coulomb::ParticleKind::Full) == 1);

	group.computeMoments();
	REQUIRE(group.positiveMoments.m0 == 1.0);
	REQUIRE(group.positiveMoments.m11 == 1.0);
	REQUIRE(group.negativeMoments.m2 == 4.0);
	REQUIRE(group.fullMoments.m31 == 27.0);

	group.copyMoments();
	REQUIRE(group.previousPositiveMoments.m11 == 1.0);
	REQUIRE(group.previousNegativeMoments.m2 == 4.0);

	group.clear(coulomb::ParticleKind::Negative);
	REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(group.size(coulomb::ParticleKind::Negative) == 0);
	REQUIRE(group.size(coulomb::ParticleKind::Full) == 1);

	REQUIRE(coulomb::ParticleKindCodec{}.decode('p') ==
			coulomb::ParticleKind::Positive);
	REQUIRE(coulomb::ParticleKindCodec{}.encode(coulomb::ParticleKind::Full) ==
			'f');
	REQUIRE_THROWS_AS(coulomb::ParticleKindCodec{}.decode('x'),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.particles.particle-group merge and position operations "
		  "preserve typed data",
		  "[particles][operations][reproducibility]") {
	coulomb::NeParticleGroup source;
	source.pushBack(particle(1.0, 1.0, 0.0, 0.0),
					coulomb::ParticleKind::Positive);
	source.pushBack(particle(2.0, 0.0, 2.0, 0.0),
					coulomb::ParticleKind::Negative);
	source.pushBack(particle(3.0, 0.0, 0.0, 3.0), coulomb::ParticleKind::Full);

	coulomb::NeParticleGroup merged;
	merged.pushBack(particle(0.0, -1.0, -1.0, -1.0),
					coulomb::ParticleKind::Positive);
	coulomb::ParticleGroupOperations{}.mergeSigned(merged, source);
	REQUIRE(merged.size(coulomb::ParticleKind::Positive) == 2);
	REQUIRE(merged.size(coulomb::ParticleKind::Negative) == 1);
	REQUIRE(merged.size(coulomb::ParticleKind::Full) == 0);
	REQUIRE(merged.list(1, coulomb::ParticleKind::Positive).velocity() ==
			source.list(0, coulomb::ParticleKind::Positive).velocity());
	REQUIRE(merged.list(0, coulomb::ParticleKind::Negative).velocity() ==
			source.list(0, coulomb::ParticleKind::Negative).velocity());

	coulomb::ParticleGroupOperations{}.mergeFull(merged, source);
	REQUIRE(merged.size(coulomb::ParticleKind::Positive) == 2);
	REQUIRE(merged.size(coulomb::ParticleKind::Negative) == 1);
	REQUIRE(merged.size(coulomb::ParticleKind::Full) == 1);
	REQUIRE(merged.list(0, coulomb::ParticleKind::Full).velocity() ==
			source.list(0, coulomb::ParticleKind::Full).velocity());

	coulomb::NeParticleGroup selected;
	coulomb::ParticleGroupOperations{}.merge(
		selected, source,
		{coulomb::ParticleKind::Negative, coulomb::ParticleKind::Full});
	REQUIRE(selected.size(coulomb::ParticleKind::Positive) == 0);
	REQUIRE(selected.size(coulomb::ParticleKind::Negative) == 1);
	REQUIRE(selected.size(coulomb::ParticleKind::Full) == 1);

	auto first = source;
	auto repeated = source;
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext repeatedRandom;
	firstRandom.reseed(8675309);
	repeatedRandom.reseed(8675309);
	coulomb::ParticleGroupOperations{}.assignPositions(first, -2.0, 3.0,
													   firstRandom);
	coulomb::ParticleGroupOperations{}.assignPositions(repeated, -2.0, 3.0,
													   repeatedRandom);
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative,
		  coulomb::ParticleKind::Full}) {
		for (int index = 0; index < first.size(kind); ++index) {
			REQUIRE(first.list(index, kind).position() ==
					repeated.list(index, kind).position());
			REQUIRE(first.list(index, kind).position() >= -2.0);
			REQUIRE(first.list(index, kind).position() < 3.0);
			REQUIRE(first.list(index, kind).velocity() ==
					source.list(index, kind).velocity());
		}
	}
}

TEST_CASE("negpar.unit.diagnostics.particle count diagnostics validate the "
		  "requested grid size",
		  "[diagnostics][particles]") {
	std::vector<coulomb::NeParticleGroup> groups(2);
	coulomb::NumericGridClass grid(2);
	coulomb::Diagnostics diagnostics(grid);
	groups[0].pushBack(particle(0.0, 1.0, 0.0, 0.0),
					   coulomb::ParticleKind::Positive);
	groups[1].pushBack(particle(0.0, 2.0, 0.0, 0.0),
					   coulomb::ParticleKind::Negative);

	REQUIRE(diagnostics.particleCount(groups, 2,
									  coulomb::ParticleKind::Positive) == 1);
	REQUIRE(diagnostics.particleCount(groups, 2,
									  coulomb::ParticleKind::Negative) == 1);
	REQUIRE_THROWS_AS(
		diagnostics.particleCount(groups, 3, coulomb::ParticleKind::Positive),
		std::invalid_argument);
}

TEST_CASE("negpar.unit.random.random permutation validates its requested size",
		  "[utils]") {
	coulomb::RandomContext random;
	random.reseed(1234);
	REQUIRE_THROWS(coulomb::RandomSampling(random).permutation(2, 3));

	const auto permutation = coulomb::RandomSampling(random).permutation(10, 5);
	REQUIRE(permutation.size() == 5);
	for (const int value : permutation)
		REQUIRE(value >= 1);
	for (const int value : permutation)
		REQUIRE(value <= 10);
}

TEST_CASE("negpar.unit.random.random generator can be explicitly seeded",
		  "[utils][reproducibility]") {
	coulomb::RandomContext random;
	random.reseed(12345);
	const double first = coulomb::RandomSampling(random).uniform();
	const double normal = coulomb::RandomSampling(random).normal();

	random.reseed(12345);
	REQUIRE(coulomb::RandomSampling(random).uniform() == first);
	REQUIRE(coulomb::RandomSampling(random).normal() == normal);
}

TEST_CASE(
	"negpar.unit.random.random context owns an independently reproducible seed",
	"[utils][reproducibility]") {
	STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<coulomb::RandomContext>);
	STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<coulomb::RandomContext>);
	STATIC_REQUIRE(std::is_move_constructible_v<coulomb::RandomContext>);
	STATIC_REQUIRE(std::is_move_assignable_v<coulomb::RandomContext>);

	coulomb::RandomContext firstContext;
	coulomb::RandomContext secondContext;
	firstContext.reseed(24680);
	secondContext.reseed(24680);

	REQUIRE(coulomb::RandomSampling(firstContext).uniform() ==
			coulomb::RandomSampling(secondContext).uniform());
	REQUIRE(coulomb::RandomSampling(firstContext).normal() ==
			coulomb::RandomSampling(secondContext).normal());
}

TEST_CASE(
	"negpar.unit.random.random generator is reproducible per OpenMP thread",
	"[utils][openmp][reproducibility]") {
	const int threadCount = std::max(1, std::min(4, omp_get_max_threads()));
	constexpr int drawsPerThread = 4096;
	const auto sampleCount =
		static_cast<std::size_t>(threadCount) * drawsPerThread;
	std::vector<double> first(sampleCount), second(sampleCount);
	coulomb::RandomContext random;

	const auto draw = [&](std::vector<double>& uniform,
						  std::vector<double>& normal) {
#pragma omp parallel num_threads(threadCount)
		{
			const int thread = omp_get_thread_num();
			for (int drawIndex = 0; drawIndex < drawsPerThread; ++drawIndex) {
				const auto index =
					static_cast<std::size_t>(thread) * drawsPerThread +
					drawIndex;
				uniform[index] = coulomb::RandomSampling(random).uniform();
				normal[index] = coulomb::RandomSampling(random).normal();
			}
		}
	};

	random.reseed(9876);
	draw(first, second);

	std::vector<double> firstAgain(sampleCount), secondAgain(sampleCount);
	random.reseed(9876);
	draw(firstAgain, secondAgain);

	REQUIRE(firstAgain == first);
	REQUIRE(secondAgain == second);
	REQUIRE(std::all_of(first.begin(), first.end(), [](double value) {
		return std::isfinite(value) && value > 0.0 && value < 1.0;
	}));
	REQUIRE(std::all_of(second.begin(), second.end(),
						[](double value) { return std::isfinite(value); }));
}

TEST_CASE("negpar.unit.collisions.parallel collisions preserve per-cell "
		  "invariants and replay",
		  "[collisions][openmp][invariants][reproducibility]") {
	const int threadCount = std::max(1, std::min(4, omp_get_max_threads()));
	if (threadCount < 2)
		SKIP("OpenMP runtime exposes only one thread");

	const auto makeCells = [threadCount] {
		std::vector<std::vector<coulomb::Particle1D3D>> cells(threadCount);
		for (int cell = 0; cell < threadCount; ++cell) {
			for (int index = 0; index < 8; ++index) {
				const double offset = static_cast<double>(cell * 8 + index);
				cells[cell].push_back(particle(cell + 0.25, 0.5 + 0.07 * offset,
											   -1.0 + 0.11 * offset,
											   0.25 - 0.05 * offset));
			}
		}
		return cells;
	};
	const auto summarize = [](const auto& cells) {
		std::vector<std::array<double, 4>> summaries(cells.size());
		for (std::size_t cell = 0; cell < cells.size(); ++cell) {
			for (const auto& value : cells[cell]) {
				for (int component = 0; component < 3; ++component) {
					const double velocity = value.velocity(component);
					summaries[cell][component] += velocity;
					summaries[cell][3] += velocity * velocity;
				}
			}
		}
		return summaries;
	};

	auto first = makeCells();
	auto repeated = first;
	const auto before = summarize(first);
	coulomb::ParaClass parameters;
	parameters.dt = 0.01;
	coulomb::RandomContext random;

	const auto collide = [&](auto& cells) {
#pragma omp parallel for num_threads(threadCount) schedule(static, 1)
		for (int cell = 0; cell < threadCount; ++cell) {
			coulomb::CollisionOperator(parameters, random)
				.collideHomogeneous(cells[cell],
									static_cast<int>(cells[cell].size()));
		}
	};

	random.reseed(424242);
	collide(first);
	random.reseed(424242);
	collide(repeated);

	const auto after = summarize(first);
	for (std::size_t cell = 0; cell < before.size(); ++cell) {
		for (int component = 0; component < 3; ++component) {
			REQUIRE(after[cell][component] ==
					Catch::Approx(before[cell][component]).margin(1e-12));
		}
		REQUIRE(after[cell][3] == Catch::Approx(before[cell][3]).margin(1e-10));
		REQUIRE(repeated[cell].size() == first[cell].size());
		for (std::size_t index = 0; index < first[cell].size(); ++index) {
			REQUIRE(repeated[cell][index].velocity() ==
					first[cell][index].velocity());
		}
	}
}

TEST_CASE(
	"negpar.unit.cli.run options parse seed, threads, and output directory",
	"[cli]") {
	char arg0[] = "negpar_inhomo";
	char arg1[] = "--seed";
	char arg2[] = "12345";
	char arg3[] = "--threads";
	char arg4[] = "1";
	char arg5[] = "--steps";
	char arg6[] = "10";
	char arg7[] = "--output-dir";
	char arg8[] = "run-test";
	char* argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8};

	const auto options = coulomb::RunOptions::parse(9, argv);

	REQUIRE(options.seed == 12345);
	REQUIRE(options.steps == 10);
	REQUIRE(options.threads == 1);
	REQUIRE(options.outputDirectory == "run-test");
}

TEST_CASE("negpar.unit.cli.run options reject invalid values",
		  "[cli][validation]") {
	char arg0[] = "negpar_inhomo";
	char arg1[] = "--threads";
	char arg2[] = "0";
	char* argv[] = {arg0, arg1, arg2};

	REQUIRE_THROWS_AS(coulomb::RunOptions::parse(3, argv),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.cli.run options reject malformed numeric values",
		  "[cli][validation]") {
	char arg0[] = "negpar_inhomo";
	char arg1[] = "--steps";
	char arg2[] = "10steps";
	char* argv[] = {arg0, arg1, arg2};

	REQUIRE_THROWS_AS(coulomb::RunOptions::parse(3, argv),
					  std::invalid_argument);
}

TEST_CASE("negpar.unit.cli.runtime state reset clears cross-run counters",
		  "[cli][state]") {
	coulomb::SimulationState state;
	state.saveIndex = 9;
	state.filenameWithNumber = true;
	state.saveFlux = true;
	state.movedCount = 4;
	state.resampleCount = 5;
	state.syncTime = 1.25;

	coulomb::RunOptions::resetRuntimeState(state);

	REQUIRE(state.saveIndex == 0);
	REQUIRE_FALSE(state.filenameWithNumber);
	REQUIRE_FALSE(state.saveFlux);
	REQUIRE(state.movedCount == 0);
	REQUIRE(state.resampleCount == 0);
	REQUIRE(state.syncTime == 0.0);
}

TEST_CASE(
	"negpar.unit.output.output filename numbering preserves legacy format",
	"[output]") {
	coulomb::SimulationState state;
	coulomb::OutputPaths paths(state);
	REQUIRE(paths.formatIndex(7) == "_007");
	REQUIRE(paths.formatIndex(-4) == "_-004");
	REQUIRE(paths.formatIndex(12, 2) == "_12");
}

TEST_CASE(
	"negpar.unit.output.fixed-bin histogram preserves legacy edge clamping",
	"[output][histogram]") {
	const std::vector<double> values{-1.0, 0.0, 0.49, 0.5, 0.99, 1.0, 2.0};
	std::vector<int> counts(2);

	coulomb::Histogram{}.fixedBins(values, counts, 0.0, 1.0);

	REQUIRE(counts == std::vector<int>{3, 4});
}

TEST_CASE("negpar.unit.output.macro output writer preserves numbered precision",
		  "[output]") {
	const auto directory =
		std::filesystem::temp_directory_path() / "negpar_save_macro_test";
	std::filesystem::remove_all(directory);
	std::filesystem::create_directories(directory);

	coulomb::SimulationState state;
	state.outputDirectory = directory.string();
	state.filenameWithNumber = true;
	state.saveIndex = 7;

	coulomb::MacroOutput{}.saveMacro(std::vector<double>{1.0 / 3.0, -2.0},
									 "values", state);
	coulomb::MacroOutput{}.save2D(2, 2, {{1.0, 2.0}, {3.0, 4.0}}, "matrix",
								  state);
	std::complex<double> spectrum[] = {{1.0, -2.0}, {3.0, 4.0}};
	coulomb::MacroOutput{}.saveComplex(2, spectrum, "spectrum", state);

	std::ifstream file(directory / "values_007.txt");
	REQUIRE(file.good());
	std::string contents((std::istreambuf_iterator<char>(file)),
						 std::istreambuf_iterator<char>());
	file.close();
	const bool matrixExists =
		std::filesystem::exists(directory / "matrix_007.txt");
	const bool realExists =
		std::filesystem::exists(directory / "spectrum_r.txt");
	const bool imaginaryExists =
		std::filesystem::exists(directory / "spectrum_i.txt");
	std::vector<std::string> filenames;
	for (const auto& entry : std::filesystem::directory_iterator(directory))
		filenames.push_back(entry.path().filename().string());
	std::sort(filenames.begin(), filenames.end());

	std::ifstream matrixFile(directory / "matrix_007.txt");
	std::string matrixContents((std::istreambuf_iterator<char>(matrixFile)),
							   std::istreambuf_iterator<char>());
	std::ifstream realFile(directory / "spectrum_r.txt");
	std::string realContents((std::istreambuf_iterator<char>(realFile)),
							 std::istreambuf_iterator<char>());
	std::ifstream imaginaryFile(directory / "spectrum_i.txt");
	std::string imaginaryContents(
		(std::istreambuf_iterator<char>(imaginaryFile)),
		std::istreambuf_iterator<char>());
	matrixFile.close();
	realFile.close();
	imaginaryFile.close();
	bool invalidDimensions = false;
	try {
		coulomb::MacroOutput{}.save2D(1, 2, {{1.0}}, "invalid", state);
	}
	catch (const std::invalid_argument&) {
		invalidDimensions = true;
	}
	std::filesystem::remove_all(directory);

	REQUIRE(contents.find("0.333333333333333") != std::string::npos);
	REQUIRE(contents.find("-2") != std::string::npos);
	REQUIRE(matrixExists);
	REQUIRE(realExists);
	REQUIRE(imaginaryExists);
	REQUIRE((filenames ==
			 std::vector<std::string>{"matrix_007.txt", "spectrum_i.txt",
									  "spectrum_r.txt", "values_007.txt"}));
	REQUIRE(matrixContents == "1 2 \n3 4 \n");
	REQUIRE(realContents == "1\n3\n");
	REQUIRE(imaginaryContents == "-2\n4\n");
	REQUIRE(invalidDimensions);
	REQUIRE_THROWS_AS(
		coulomb::ParticleOutput{}.saveHomogeneousRadialDistribution(0, state),
		std::invalid_argument);
	REQUIRE_THROWS_AS(coulomb::OutputPaths(state).resolve("../escape.txt"),
					  std::invalid_argument);
	REQUIRE_THROWS_AS(
		coulomb::OutputPaths(state).resolve(
			(std::filesystem::temp_directory_path() / "escape.txt").string()),
		std::invalid_argument);
}
