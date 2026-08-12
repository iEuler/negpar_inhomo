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

using coulomb::Particle1d3d;

Particle1d3d particle(double x, double vx, double vy, double vz) {
  return Particle1d3d(x, {vx, vy, vz});
}

} // namespace

TEST_CASE("negpar.unit.grid.default grid is initialized", "[grid]") {
  coulomb::NumericGridClass grid;

  REQUIRE(grid.Nx == 100);
  REQUIRE(grid.Nv == 200);
  REQUIRE(grid.x.size() == static_cast<size_t>(grid.Nx));
  REQUIRE(grid.vx.size() == static_cast<size_t>(grid.Nv));
  REQUIRE(grid.dx > 0.0);
  REQUIRE(grid.dv > 0.0);
  REQUIRE(grid.dt > 0.0);
}

TEST_CASE("negpar.unit.grid.grid rejects invalid configuration",
          "[grid][validation]") {
  REQUIRE_THROWS_AS(coulomb::NumericGridClass(0), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::NumericGridClass(-4), std::invalid_argument);
  REQUIRE_THROWS_AS(
      coulomb::NumericGridClass(10, static_cast<coulomb::SimulationMethod>(99)),
      std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.method_name(
                        static_cast<coulomb::SimulationMethod>(99)),
                    std::invalid_argument);
}

TEST_CASE(
    "negpar.unit.parameters.typed collision mode preserves the legacy name",
    "[parameters]") {
  const coulomb::ParaClass parameters;
  const coulomb::IniValClass initial_data;
  coulomb::RandomContext random;
  random.reseed(2468);
  REQUIRE(parameters.method_binarycoll == coulomb::BinaryCollisionMethod::TA);
  REQUIRE(initial_data.rho == 0.0);
  REQUIRE(initial_data.velocity[0] == 0.0);
  REQUIRE(std::string(coulomb::SimulationTypes{}.binary_collision_name(
              parameters.method_binarycoll)) == "TA");
  REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.collision_name(
                        static_cast<coulomb::CollisionType>(99)),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.binary_collision_name(
                        static_cast<coulomb::BinaryCollisionMethod>(99)),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.advection.particle advection wraps periodic boundaries",
          "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = coulomb::BoundaryCondition::Periodic;
  grid.dt = 1.0;
  const double start = grid.xmax - 0.25;
  const auto moved = coulomb::Advection(grid, state)
                         .move_particle(particle(start, 1.0, 2.0, 3.0), 0.0);

  REQUIRE(moved.position() == Catch::Approx(grid.xmin + 0.75));
  REQUIRE(moved.velocity(0) == Catch::Approx(1.0));
  REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("negpar.unit.advection.particle advection reflects at nonperiodic "
          "boundaries",
          "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = coulomb::BoundaryCondition::Reflective;
  grid.dt = 1.0;
  const auto moved =
      coulomb::Advection(grid, state)
          .move_particle(particle(grid.xmax - 0.25, 1.0, 2.0, 3.0), 0.5);

  REQUIRE(moved.position() == Catch::Approx(grid.xmax - 0.75));
  REQUIRE(moved.velocity(0) == Catch::Approx(-1.5));
  REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("negpar.unit.advection.particle exactly at the upper boundary "
          "belongs to the last cell",
          "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = coulomb::BoundaryCondition::Reflective;
  grid.dt = grid.xmax - (grid.xmax - 0.5);
  coulomb::Advection advection(grid, state);
  auto moved =
      advection.move_particle(particle(grid.xmax - 0.5, 1.0, 0.0, 0.0), 0.0);

  REQUIRE(moved.position() == Catch::Approx(grid.xmax));
  REQUIRE(advection.find_particle_group(moved) == grid.Nx - 1);
}

TEST_CASE("negpar.unit.advection.particle group lookup rejects positions "
          "outside the grid",
          "[advection][validation]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  coulomb::Advection advection(grid, state);
  auto above = particle(grid.xmax + 1.0, 0.0, 0.0, 0.0);
  auto below = particle(grid.xmin - 1.0, 0.0, 0.0, 0.0);

  REQUIRE_THROWS_AS(advection.find_particle_group(above), std::out_of_range);
  REQUIRE_THROWS_AS(advection.find_particle_group(below), std::out_of_range);
}

TEST_CASE(
    "negpar.unit.numerics.numerics error function preserves the approximation",
    "[numerics]") {
  REQUIRE(std::abs(coulomb::Numerics{}.error_function(0.0)) < 2e-9);
  REQUIRE(std::abs(coulomb::Numerics{}.error_function(1.0) - 0.84270079) <
          2e-7);
  REQUIRE(std::abs(coulomb::Numerics{}.error_function(-1.0) + 0.84270079) <
          2e-7);
}

TEST_CASE("negpar.unit.numerics.one-dimensional shifts preserve legacy "
          "boundary conventions",
          "[numerics][boundaries]") {
  const std::vector<double> values = {1.0, 2.0, 3.0};
  const std::vector<double> shifted = {2.0, 3.0, 1.0};
  const std::vector<double> extended = {-1.0, 1.0, 2.0};
  const std::vector<double> difference = {-0.5, 1.0, -0.5};
  REQUIRE(coulomb::SimulationTypes{}.decode_boundary('p') ==
          coulomb::BoundaryCondition::Periodic);
  REQUIRE(coulomb::SimulationTypes{}.decode_boundary('n') ==
          coulomb::BoundaryCondition::Reflective);
  REQUIRE(coulomb::SimulationTypes{}.encode_boundary(
              coulomb::BoundaryCondition::Reflective) == 'n');
  REQUIRE_THROWS_AS(coulomb::SimulationTypes{}.decode_boundary('x'),
                    std::invalid_argument);
  REQUIRE(coulomb::Numerics{}.circular_shift(values, 3, 1) == shifted);
  REQUIRE(coulomb::Numerics{}.circular_shift(values, 3, 4) == shifted);
  REQUIRE(coulomb::Numerics{}.boundary_shift(values, 3, -1, -1.0) == extended);
  REQUIRE(coulomb::Numerics{}.central_difference(
              values, 3, coulomb::BoundaryCondition::Periodic) == difference);
  REQUIRE(
      coulomb::Numerics{}.circular_shift(std::vector<double>{}, 0, 8).empty());
  REQUIRE_THROWS_AS(coulomb::Numerics{}.circular_shift(values, 2, 1),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::Numerics{}.central_difference(
                        values, 3, static_cast<coulomb::BoundaryCondition>(99)),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.numerics.kinetic Euler fluxes handle periodic and "
          "constant states",
          "[numerics][euler]") {
  const std::vector<double> values = {1.0, 2.0, 4.0};
  REQUIRE(
      coulomb::Numerics{}.limited_flux(values, 3, 2.0, 1.0, 1, values,
                                       coulomb::BoundaryCondition::Periodic) ==
      std::vector<double>{-1.5, 0.5, 1.0});
  REQUIRE(
      coulomb::Numerics{}.limited_flux(values, 3, 2.0, 1.0, -1, values,
                                       coulomb::BoundaryCondition::Periodic) ==
      std::vector<double>{0.5, 1.0, -1.5});

  std::vector<double> density_change(3), momentum_change(3), energy_change(3);
  const std::vector<double> constant(3, 1.0);
  coulomb::Numerics{}.advance_kinetic_euler(
      constant, constant, constant, 3, 1.0, 0.1,
      coulomb::BoundaryCondition::Periodic, density_change, momentum_change,
      energy_change);
  REQUIRE(density_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE(momentum_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE(energy_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE_THROWS_AS(coulomb::Numerics{}.advance_kinetic_euler(
                        constant, constant, std::vector<double>(3, 0.0), 3, 1.0,
                        0.1, coulomb::BoundaryCondition::Periodic,
                        density_change, momentum_change, energy_change),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFT reshape helpers preserve row-major ordering",
          "[fft][numerics]") {
  const std::vector<int> flat = {0, 1, 2, 3, 4, 5, 6, 7};
  const auto shaped = coulomb::TensorReshape{}.to_3d(flat, 2, 2, 2);

  REQUIRE(shaped[0][0][0] == 0);
  REQUIRE(shaped[0][1][1] == 3);
  REQUIRE(shaped[1][0][0] == 4);
  REQUIRE(coulomb::TensorReshape{}.to_1d(shaped) == flat);
  REQUIRE_THROWS_AS(coulomb::TensorReshape{}.to_3d(flat, 3, 2, 2),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFT reshape helper rejects ragged input",
          "[fft][validation]") {
  const std::vector<std::vector<std::vector<int>>> ragged = {{{1, 2}, {3}},
                                                             {{4, 5}, {6, 7}}};
  REQUIRE_THROWS_AS(coulomb::TensorReshape{}.to_1d(ragged),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.fft.FFTW wrappers have explicit ownership and validate "
          "dimensions",
          "[fft][resource]") {
  static_assert(!std::is_copy_constructible_v<coulomb::FFT1D>);
  static_assert(!std::is_copy_assignable_v<coulomb::FFT1D>);
  static_assert(!std::is_move_constructible_v<coulomb::FFT1D>);
  static_assert(!std::is_move_assignable_v<coulomb::FFT1D>);
  static_assert(!std::is_copy_constructible_v<coulomb::FFT3D>);
  static_assert(!std::is_copy_assignable_v<coulomb::FFT3D>);
  static_assert(!std::is_move_constructible_v<coulomb::FFT3D>);
  static_assert(!std::is_move_assignable_v<coulomb::FFT3D>);

  REQUIRE_THROWS_AS(coulomb::FFT1D(0), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::FFT3D(0, 2, 2), std::invalid_argument);

  // Repeated construction and destruction exercises plan-before-buffer cleanup
  // under the sanitizer preset. FFTW's inverse is intentionally unnormalized.
  for (int iteration = 0; iteration < 8; ++iteration) {
    coulomb::FFT1D fft(4);
    const std::vector<double> input{1.0, -2.0, 0.5, 3.0};
    const auto transformed = fft.fft(input);
    const auto restored = fft.ifft(transformed);
    for (std::size_t index = 0; index < input.size(); ++index) {
      REQUIRE(restored[index] ==
              Catch::Approx(4.0 * input[index]).margin(1e-12));
    }
    REQUIRE_THROWS_AS(fft.fft(std::vector<double>(3)), std::invalid_argument);
  }

  const coulomb::Vector3D input3d = {{{1.0, 2.0}, {3.0, 4.0}},
                                     {{-1.0, -2.0}, {-3.0, -4.0}}};
  coulomb::FFT3D fft3d(2, 2, 2);
  const auto restored3d = fft3d.ifft(fft3d.fft(input3d));
  for (std::size_t first = 0; first < 2; ++first) {
    for (std::size_t second = 0; second < 2; ++second) {
      for (std::size_t third = 0; third < 2; ++third) {
        REQUIRE(
            restored3d[first][second][third] ==
            Catch::Approx(8.0 * input3d[first][second][third]).margin(1e-12));
      }
    }
  }
}

TEST_CASE("negpar.unit.moments.particle group computes moments",
          "[particles][moments]") {
  REQUIRE_THROWS_AS(coulomb::Particle1d3d(std::vector<double>{1.0, 2.0}),
                    std::invalid_argument);
  coulomb::Particle1d3d particle_with_velocity;
  REQUIRE_THROWS_AS(
      particle_with_velocity.set_velocity(std::vector<double>{1.0, 2.0}),
      std::invalid_argument);
  REQUIRE_THROWS_AS(particle_with_velocity.velocity(-1), std::out_of_range);
  REQUIRE_THROWS_AS(particle_with_velocity.velocity(3), std::out_of_range);

  coulomb::ParticleGroup group;
  group.push_back(particle(0.0, 1.0, 2.0, 3.0));
  group.push_back(particle(0.0, -1.0, 4.0, 0.0));

  group.computemoments();

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
  coulomb::MomentOperations{}.primitive_to_conserved(2.0, velocity, 3.0,
                                                     momentum, energy);

  REQUIRE(momentum[0] == 2.0);
  REQUIRE(momentum[1] == -4.0);
  REQUIRE(momentum[2] == 1.0);
  REQUIRE(energy == Catch::Approx(2.0 * (1.0 + 4.0 + 0.25 + 9.0) / 2.0));

  double recovered_velocity[3] = {};
  double recovered_temperature = 0.0;
  coulomb::MomentOperations{}.conserved_to_primitive(
      2.0, momentum, energy, recovered_velocity, recovered_temperature);
  for (int component = 0; component < 3; ++component)
    REQUIRE(recovered_velocity[component] ==
            Catch::Approx(velocity[component]));
  REQUIRE(recovered_temperature == Catch::Approx(3.0));

  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.primitive_to_conserved(
                        2.0, nullptr, 3.0, momentum, energy),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      coulomb::MomentOperations{}.conserved_to_primitive(
          2.0, nullptr, energy, recovered_velocity, recovered_temperature),
      std::invalid_argument);

  coulomb::ParticleGroup empty_group;
  double density = 0.0;
  double ignored_energy = 0.0;
  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.particle_to_conserved(
                        empty_group, 1.0, density, nullptr, ignored_energy),
                    std::invalid_argument);

  std::vector<coulomb::ParticleGroup> groups(1);
  std::vector<double> output(1);
  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.compute_primitive(
                        2, groups, 1.0, output, output, output, output, output),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.primitive_to_conserved(
                        2, output, output, output, output, output),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.conserved_to_primitive(
                        2, output, output, output, output, output),
                    std::invalid_argument);
}

TEST_CASE("negpar.unit.moments.macroscopic moment updates reconstruct particle "
          "fields",
          "[moments][reconstruction]") {
  coulomb::NumericGridClass grid(3);
  std::vector<coulomb::NeParticleGroup> groups(3);
  REQUIRE_THROWS_AS(coulomb::MomentOperations{}.moment_change(nullptr, grid),
                    std::invalid_argument);
  for (auto &group : groups) {
    group.rhoM = 1.0;
    group.u1M = 0.0;
    group.TprtM = 1.0;
    group.computemoments();
  }

  coulomb::MomentOperations{}.update_primitive(groups, grid);
  for (const auto &group : groups) {
    REQUIRE(group.rho == Catch::Approx(1.0));
    REQUIRE(group.u1 == Catch::Approx(0.0));
    REQUIRE(group.Tprt == Catch::Approx(1.0));
  }

  groups[0].push_back(particle(0.0, 2.0, 0.0, 0.0),
                      coulomb::ParticleKind::Full);
  for (auto &group : groups)
    group.computemoments();
  coulomb::MomentOperations{}.update_full_primitive(groups, grid);
  REQUIRE(groups[0].rhoF == Catch::Approx(grid.Neff_F / grid.dx));
  REQUIRE(groups[0].u1F == Catch::Approx(2.0));
}

TEST_CASE(
    "negpar.unit.moments.Maxwellian update applies conservative moment changes",
    "[moments][maxwellian]") {
  coulomb::NumericGridClass grid(1);
  std::vector<coulomb::NeParticleGroup> groups(1);
  groups[0].rhoM = 2.0;
  groups[0].u1M = 1.0;
  groups[0].TprtM = 3.0;
  groups[0].drho = 0.5;
  groups[0].dm1 = 0.25;
  groups[0].denergy = 1.0;

  coulomb::MomentOperations{}.update_maxwellian(groups, grid);

  REQUIRE(groups[0].rhoM == Catch::Approx(1.5));
  REQUIRE(groups[0].u1M == Catch::Approx(1.75 / 1.5));
  REQUIRE(groups[0].TprtM == Catch::Approx(3.5462962963));
}

TEST_CASE("negpar.unit.resampling.full-particle count rescaling preserves "
          "effective mass",
          "[resampling][conservation]") {
  coulomb::NumericGridClass grid(2);
  grid.Neff_F = 0.25;
  std::vector<coulomb::NeParticleGroup> groups(2);

  for (auto &group : groups) {
    group.rhoM = 1.0;
    for (int particle_index = 0; particle_index < 30; ++particle_index)
      group.push_back(particle(0.0, particle_index * 0.01, 0.0, 0.0),
                      coulomb::ParticleKind::Full);
  }

  const double mass_before = groups[0].rhoM * grid.dx * groups.size();
  coulomb::RandomContext random;
  random.reseed(2468);

  // Passing zero forces the routine to perform the count-reduction path.
  coulomb::ParaClass parameters;
  coulomb::SimulationState state;
  state.random.reseed(2468);
  coulomb::ParticleResampling(grid, parameters, state)
      .resample_full_preserving_mass(groups, 0);

  const int full_count = coulomb::Diagnostics(grid).particle_count(
      groups, grid.Nx, coulomb::ParticleKind::Full);
  REQUIRE(full_count == 50);
  REQUIRE(grid.Neff_F * full_count == Catch::Approx(mass_before).margin(1e-12));
}

TEST_CASE("negpar.unit.resampling.shared resampling numerics map modes and "
          "evaluate all derivatives",
          "[resampling][numerics]") {
  REQUIRE(coulomb::resampling::ResamplingNumerics{}.frequencies(4) ==
          std::vector<double>{0.0, 1.0, 2.0, -1.0});
  REQUIRE(coulomb::resampling::ResamplingNumerics{}.augmented_locations(4, 2) ==
          std::vector<std::size_t>{0, 1, 2, 7});
  REQUIRE(coulomb::resampling::ResamplingNumerics{}.imaginary_frequencies(4) ==
          std::vector<std::complex<double>>{
              {0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}, {0.0, -1.0}});

  const std::vector<double> derivatives{1.0, 2.0, 3.0, 5.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0};
  REQUIRE(coulomb::resampling::ResamplingNumerics{}.evaluate_quadratic_taylor(
              0.1, 0.2, 0.4, derivatives) == Catch::Approx(3.8));
  REQUIRE_THROWS_AS(coulomb::resampling::ResamplingNumerics{}.frequencies(0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      coulomb::resampling::ResamplingNumerics{}.augmented_locations(4, 0),
      std::invalid_argument);
  REQUIRE_THROWS_AS(
      coulomb::resampling::ResamplingNumerics{}.evaluate_quadratic_taylor(
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
  particles.TprtM = 1.5;
  particles.xyz_minmax = {-2.0, 2.0, -3.0, 1.0, 0.0, 4.0};
  particles.push_back(coulomb::Particle1d3d({-1.0, -2.0, 1.0}),
                      coulomb::ParticleKind::Positive);
  particles.push_back(coulomb::Particle1d3d({1.0, 0.0, 3.0}),
                      coulomb::ParticleKind::Negative);

  auto normalized =
      coulomb::resampling::ResamplingVelocity{}.normalize_signed(particles);
  REQUIRE(normalized.rhoM == Catch::Approx(particles.rhoM));
  for (const auto kind :
       {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
    auto restored = normalized.list(kind);
    coulomb::resampling::ResamplingVelocity{}.restore(restored,
                                                      particles.xyz_minmax);
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
                            .collide_pair(first, second);

  double energy_before = 0.0;
  double energy_after = 0.0;
  for (int component = 0; component < 3; ++component) {
    REQUIRE(collided.first[component] + collided.second[component] ==
            Catch::Approx(first[component] + second[component]).margin(1e-14));
    energy_before += first[component] * first[component] +
                     second[component] * second[component];
    energy_after += collided.first[component] * collided.first[component] +
                    collided.second[component] * collided.second[component];
  }
  REQUIRE(energy_after == Catch::Approx(energy_before).margin(1e-10));
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
                            .collide_pair(first, second);
  random.reseed(13579);
  const auto repeated = coulomb::CollisionOperator(parameters, random)
                            .collide_pair(first, second);

  REQUIRE(repeated.first == expected.first);
  REQUIRE(repeated.second == expected.second);
}

TEST_CASE("negpar.unit.collisions.negative-particle collision kernels replay "
          "for an explicit seed",
          "[collisions][negative-particle][reproducibility]") {
  const auto make_group = [] {
    coulomb::NeParticleGroup group;
    group.u1M = 0.25;
    group.u2M = -0.5;
    group.u3M = 0.75;
    group.TprtM = 1.25;
    group.push_back(particle(0.0, 1.0, 0.0, 0.0),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, -1.0, 0.5, 0.0),
                    coulomb::ParticleKind::Negative);
    group.push_back(particle(0.0, 0.0, 1.0, 0.5), coulomb::ParticleKind::Full);
    group.push_back(particle(0.0, 0.5, -1.0, 1.0), coulomb::ParticleKind::Full);
    return group;
  };
  const auto require_same_particles =
      [](const coulomb::NeParticleGroup &first,
         const coulomb::NeParticleGroup &second) {
        for (const auto kind :
             {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative,
              coulomb::ParticleKind::Full}) {
          REQUIRE(first.size(kind) == second.size(kind));
          for (int index = 0; index < first.size(kind); ++index)
            REQUIRE(first.list(index, kind).velocity() ==
                    second.list(index, kind).velocity());
        }
      };

  coulomb::ParaClass parameters;
  coulomb::NumericGridClass grid(1);
  auto first = make_group();
  auto repeated = make_group();
  const auto full_before = first.list(coulomb::ParticleKind::Full);
  coulomb::RandomContext first_random;
  coulomb::RandomContext repeated_random;
  first_random.reseed(112233);
  repeated_random.reseed(112233);
  coulomb::NegativeParticleCollisions(grid, parameters, first_random)
      .collide_with_full(first);
  coulomb::NegativeParticleCollisions(grid, parameters, repeated_random)
      .collide_with_full(repeated);
  require_same_particles(first, repeated);
  REQUIRE(first.size(coulomb::ParticleKind::Full) ==
          static_cast<int>(full_before.size()));
  for (std::size_t index = 0; index < full_before.size(); ++index)
    REQUIRE(first.list(static_cast<int>(index), coulomb::ParticleKind::Full)
                .velocity() == full_before[index].velocity());

  parameters.dt = 0.1;
  parameters.coeff_binarycoll = 10.0;
  first = make_group();
  repeated = make_group();
  first_random.reseed(445566);
  repeated_random.reseed(445566);
  coulomb::NegativeParticleCollisions(grid, parameters, first_random)
      .collide_bgk_homogeneous(first);
  coulomb::NegativeParticleCollisions(grid, parameters, repeated_random)
      .collide_bgk_homogeneous(repeated);
  require_same_particles(first, repeated);
  REQUIRE(first.size(coulomb::ParticleKind::Positive) == 0);
  REQUIRE(first.size(coulomb::ParticleKind::Negative) == 0);
  for (const auto &full : first.list(coulomb::ParticleKind::Full)) {
    for (int component = 0; component < 3; ++component)
      REQUIRE(std::isfinite(full.velocity(component)));
  }
}

TEST_CASE("negpar.unit.collisions.negative-particle collision pipeline replays "
          "for an explicit seed",
          "[collisions][negative-particle][pipeline][reproducibility]") {
  const auto make_group = [] {
    coulomb::NeParticleGroup group;
    group.set_xrange(-1.0, 1.0);
    group.rhoM = 1.0;
    group.u1M = 0.2;
    group.u2M = -0.1;
    group.u3M = 0.3;
    group.TprtM = 1.1;
    group.rho = 1.2;
    group.alpha_neg = 0.0;
    group.alpha_pos = 0.0;
    group.rmax = 0.0;
    group.push_back(particle(-0.8, 1.0, 0.2, -0.5),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(-0.4, -0.5, 1.2, 0.7),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, 0.3, -0.4, 0.2),
                    coulomb::ParticleKind::Negative);
    group.push_back(particle(0.2, -1.0, 0.1, 0.4), coulomb::ParticleKind::Full);
    group.push_back(particle(0.4, 0.6, -0.8, 1.1), coulomb::ParticleKind::Full);
    group.push_back(particle(0.6, 1.3, 0.5, -0.7), coulomb::ParticleKind::Full);
    group.push_back(particle(0.8, -0.2, 1.4, 0.9), coulomb::ParticleKind::Full);
    group.computemoments();
    return group;
  };

  coulomb::ParaClass parameters;
  coulomb::NumericGridClass grid(1);
  grid.Neff = 0.2;
  const auto before = make_group();
  auto first = before;
  auto repeated = before;
  coulomb::RandomContext first_random;
  coulomb::RandomContext repeated_random;
  first_random.reseed(246810);
  repeated_random.reseed(246810);
  coulomb::NegativeParticleCollisions(grid, parameters, first_random)
      .collide_homogeneous(first);
  coulomb::NegativeParticleCollisions(grid, parameters, repeated_random)
      .collide_homogeneous(repeated);

  bool velocity_changed = false;
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
      velocity_changed =
          velocity_changed || first.list(index, kind).velocity() !=
                                  before.list(index, kind).velocity();
      for (int component = 0; component < 3; ++component)
        REQUIRE(std::isfinite(first.list(index, kind).velocity(component)));
    }
  }
  REQUIRE(velocity_changed);
}

TEST_CASE("negpar.unit.negative_particles.negative-particle source sampling "
          "replays for an explicit seed",
          "[negative-particle][sampling][reproducibility]") {
  const auto make_group = [] {
    coulomb::NeParticleGroup group;
    group.rhoM = 1.0;
    group.u1M = 0.2;
    group.u2M = -0.1;
    group.u3M = 0.3;
    group.TprtM = 1.1;
    group.rho = 1.2;
    for (int repeat = 0; repeat < 100; ++repeat) {
      group.push_back(particle(0.0, 1.0, 0.2, -0.5),
                      coulomb::ParticleKind::Positive);
      group.push_back(particle(0.0, -0.5, 1.2, 0.7),
                      coulomb::ParticleKind::Positive);
      group.push_back(particle(0.0, 0.8, -1.0, 1.4),
                      coulomb::ParticleKind::Positive);
      group.push_back(particle(0.0, 0.3, -0.4, 0.2),
                      coulomb::ParticleKind::Negative);
      group.push_back(particle(0.0, -0.7, 0.9, -1.1),
                      coulomb::ParticleKind::Negative);
    }
    group.computemoments();
    return group;
  };

  coulomb::ParaClass parameters;
  auto source = make_group();
  coulomb::NegativeParticleSampling{}.update_bounds(source, parameters);
  REQUIRE(std::isfinite(source.alpha_neg));
  REQUIRE(std::isfinite(source.alpha_pos));
  REQUIRE(source.alpha_neg > 0.0);
  REQUIRE(source.alpha_pos > 0.0);
  REQUIRE(source.rmax > 0.0);

  const std::vector<double> mean{source.u1M, source.u2M, source.u3M};
  const double expected_at_mean =
      source.rhoM / std::pow(std::sqrt(2.0 * coulomb::pi * source.TprtM), 3);
  REQUIRE(coulomb::NegativeParticleSampling{}.evaluate_maxwellian(
              mean, source) == Catch::Approx(expected_at_mean).margin(1e-15));

  constexpr double effective_particles = 0.002;
  coulomb::RandomContext first_count_random;
  coulomb::RandomContext repeated_count_random;
  first_count_random.reseed(991122);
  repeated_count_random.reseed(991122);
  const int first_virtual_count =
      coulomb::NegativeParticleSampling{}.estimate_virtual_count(
          source, effective_particles, first_count_random);
  const int repeated_virtual_count =
      coulomb::NegativeParticleSampling{}.estimate_virtual_count(
          source, effective_particles, repeated_count_random);
  REQUIRE(first_virtual_count == repeated_virtual_count);
  REQUIRE(first_virtual_count > 0);

  auto repeated_source = source;
  coulomb::NeParticleGroup first_sample;
  coulomb::NeParticleGroup repeated_sample;
  coulomb::RandomContext first_random;
  coulomb::RandomContext repeated_random;
  first_random.reseed(334455);
  repeated_random.reseed(334455);
  coulomb::NegativeParticleSampling{}.sample_delta(
      source, first_sample, parameters, effective_particles, first_random);
  coulomb::NegativeParticleSampling{}.sample_delta(
      repeated_source, repeated_sample, parameters, effective_particles,
      repeated_random);

  int sampled_particles = 0;
  for (const auto kind :
       {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
    REQUIRE(first_sample.size(kind) == repeated_sample.size(kind));
    sampled_particles += first_sample.size(kind);
    for (int index = 0; index < first_sample.size(kind); ++index) {
      REQUIRE(first_sample.list(index, kind).velocity() ==
              repeated_sample.list(index, kind).velocity());
      for (int component = 0; component < 3; ++component)
        REQUIRE(
            std::isfinite(first_sample.list(index, kind).velocity(component)));
    }
  }
  REQUIRE(sampled_particles > 0);
}

TEST_CASE("negpar.unit.negative_particles.signed particle conservation "
          "restores target moments reproducibly",
          "[negative-particle][conservation][reproducibility]") {
  const auto make_group = [] {
    coulomb::NeParticleGroup group;
    group.push_back(particle(0.0, 1.0, 0.2, -0.5),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, -0.5, 1.2, 0.7),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, 0.8, -1.0, 1.4),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, -1.2, 0.6, -0.9),
                    coulomb::ParticleKind::Positive);
    group.push_back(particle(0.0, 0.3, -0.4, 0.2),
                    coulomb::ParticleKind::Negative);
    group.push_back(particle(0.0, -0.7, 0.9, -1.1),
                    coulomb::ParticleKind::Negative);
    return group;
  };
  const auto signed_moments = [](coulomb::NeParticleGroup &group,
                                 double effective_particles) {
    group.computemoments();
    return std::array<double, 7>{
        effective_particles *
            (group.positive_moments.m0 - group.negative_moments.m0),
        effective_particles *
            (group.positive_moments.m11 - group.negative_moments.m11),
        effective_particles *
            (group.positive_moments.m12 - group.negative_moments.m12),
        effective_particles *
            (group.positive_moments.m13 - group.negative_moments.m13),
        effective_particles *
            (group.positive_moments.m21 - group.negative_moments.m21),
        effective_particles *
            (group.positive_moments.m22 - group.negative_moments.m22),
        effective_particles *
            (group.positive_moments.m23 - group.negative_moments.m23)};
  };

  constexpr double effective_particles = 0.75;
  auto baseline = make_group();
  const auto target = signed_moments(baseline, effective_particles);
  auto first = make_group();
  auto repeated = make_group();
  const std::array<double, 3> perturbed_velocity{1.1, 0.1, -0.3};
  first.list(0, coulomb::ParticleKind::Positive)
      .set_velocity(perturbed_velocity);
  repeated.list(0, coulomb::ParticleKind::Positive)
      .set_velocity(perturbed_velocity);

  coulomb::RandomContext first_random;
  coulomb::RandomContext repeated_random;
  first_random.reseed(778899);
  repeated_random.reseed(778899);
  coulomb::ParticleConservation{}.enforce(
      target[0], target[1], target[2], target[3], target[4], target[5],
      target[6], first, effective_particles, true, first_random);
  coulomb::ParticleConservation{}.enforce(
      target[0], target[1], target[2], target[3], target[4], target[5],
      target[6], repeated, effective_particles, true, repeated_random);

  const auto actual = signed_moments(first, effective_particles);
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
  REQUIRE(group.positive_moments.m0 == 0.0);
  group.push_back(particle(0.0, 1.0, 0.0, 0.0),
                  coulomb::ParticleKind::Positive);
  group.push_back(particle(0.0, 2.0, 0.0, 0.0),
                  coulomb::ParticleKind::Negative);
  group.push_back(particle(0.0, 3.0, 0.0, 0.0), coulomb::ParticleKind::Full);

  REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
  REQUIRE(group.size(coulomb::ParticleKind::Negative) == 1);
  REQUIRE(group.size(coulomb::ParticleKind::Full) == 1);

  group.computemoments();
  REQUIRE(group.positive_moments.m0 == 1.0);
  REQUIRE(group.positive_moments.m11 == 1.0);
  REQUIRE(group.negative_moments.m2 == 4.0);
  REQUIRE(group.full_moments.m31 == 27.0);

  group.copymoments();
  REQUIRE(group.previous_positive_moments.m11 == 1.0);
  REQUIRE(group.previous_negative_moments.m2 == 4.0);

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
  source.push_back(particle(1.0, 1.0, 0.0, 0.0),
                   coulomb::ParticleKind::Positive);
  source.push_back(particle(2.0, 0.0, 2.0, 0.0),
                   coulomb::ParticleKind::Negative);
  source.push_back(particle(3.0, 0.0, 0.0, 3.0), coulomb::ParticleKind::Full);

  coulomb::NeParticleGroup merged;
  merged.push_back(particle(0.0, -1.0, -1.0, -1.0),
                   coulomb::ParticleKind::Positive);
  coulomb::ParticleGroupOperations{}.merge_signed(merged, source);
  REQUIRE(merged.size(coulomb::ParticleKind::Positive) == 2);
  REQUIRE(merged.size(coulomb::ParticleKind::Negative) == 1);
  REQUIRE(merged.size(coulomb::ParticleKind::Full) == 0);
  REQUIRE(merged.list(1, coulomb::ParticleKind::Positive).velocity() ==
          source.list(0, coulomb::ParticleKind::Positive).velocity());
  REQUIRE(merged.list(0, coulomb::ParticleKind::Negative).velocity() ==
          source.list(0, coulomb::ParticleKind::Negative).velocity());

  coulomb::ParticleGroupOperations{}.merge_full(merged, source);
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
  coulomb::RandomContext first_random;
  coulomb::RandomContext repeated_random;
  first_random.reseed(8675309);
  repeated_random.reseed(8675309);
  coulomb::ParticleGroupOperations{}.assign_positions(first, -2.0, 3.0,
                                                      first_random);
  coulomb::ParticleGroupOperations{}.assign_positions(repeated, -2.0, 3.0,
                                                      repeated_random);
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
  groups[0].push_back(particle(0.0, 1.0, 0.0, 0.0),
                      coulomb::ParticleKind::Positive);
  groups[1].push_back(particle(0.0, 2.0, 0.0, 0.0),
                      coulomb::ParticleKind::Negative);

  REQUIRE(diagnostics.particle_count(groups, 2,
                                     coulomb::ParticleKind::Positive) == 1);
  REQUIRE(diagnostics.particle_count(groups, 2,
                                     coulomb::ParticleKind::Negative) == 1);
  REQUIRE_THROWS_AS(
      diagnostics.particle_count(groups, 3, coulomb::ParticleKind::Positive),
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

  coulomb::RandomContext first_context;
  coulomb::RandomContext second_context;
  first_context.reseed(24680);
  second_context.reseed(24680);

  REQUIRE(coulomb::RandomSampling(first_context).uniform() ==
          coulomb::RandomSampling(second_context).uniform());
  REQUIRE(coulomb::RandomSampling(first_context).normal() ==
          coulomb::RandomSampling(second_context).normal());
}

TEST_CASE(
    "negpar.unit.random.random generator is reproducible per OpenMP thread",
    "[utils][openmp][reproducibility]") {
  const int thread_count = std::max(1, std::min(4, omp_get_max_threads()));
  constexpr int draws_per_thread = 4096;
  const auto sample_count =
      static_cast<std::size_t>(thread_count) * draws_per_thread;
  std::vector<double> first(sample_count), second(sample_count);
  coulomb::RandomContext random;

  const auto draw = [&](std::vector<double> &uniform,
                        std::vector<double> &normal) {
#pragma omp parallel num_threads(thread_count)
    {
      const int thread = omp_get_thread_num();
      for (int draw_index = 0; draw_index < draws_per_thread; ++draw_index) {
        const auto index =
            static_cast<std::size_t>(thread) * draws_per_thread + draw_index;
        uniform[index] = coulomb::RandomSampling(random).uniform();
        normal[index] = coulomb::RandomSampling(random).normal();
      }
    }
  };

  random.reseed(9876);
  draw(first, second);

  std::vector<double> first_again(sample_count), second_again(sample_count);
  random.reseed(9876);
  draw(first_again, second_again);

  REQUIRE(first_again == first);
  REQUIRE(second_again == second);
  REQUIRE(std::all_of(first.begin(), first.end(), [](double value) {
    return std::isfinite(value) && value > 0.0 && value < 1.0;
  }));
  REQUIRE(std::all_of(second.begin(), second.end(),
                      [](double value) { return std::isfinite(value); }));
}

TEST_CASE("negpar.unit.collisions.parallel collisions preserve per-cell "
          "invariants and replay",
          "[collisions][openmp][invariants][reproducibility]") {
  const int thread_count = std::max(1, std::min(4, omp_get_max_threads()));
  if (thread_count < 2)
    SKIP("OpenMP runtime exposes only one thread");

  const auto make_cells = [thread_count] {
    std::vector<std::vector<coulomb::Particle1d3d>> cells(thread_count);
    for (int cell = 0; cell < thread_count; ++cell) {
      for (int index = 0; index < 8; ++index) {
        const double offset = static_cast<double>(cell * 8 + index);
        cells[cell].push_back(particle(cell + 0.25, 0.5 + 0.07 * offset,
                                       -1.0 + 0.11 * offset,
                                       0.25 - 0.05 * offset));
      }
    }
    return cells;
  };
  const auto summarize = [](const auto &cells) {
    std::vector<std::array<double, 4>> summaries(cells.size());
    for (std::size_t cell = 0; cell < cells.size(); ++cell) {
      for (const auto &value : cells[cell]) {
        for (int component = 0; component < 3; ++component) {
          const double velocity = value.velocity(component);
          summaries[cell][component] += velocity;
          summaries[cell][3] += velocity * velocity;
        }
      }
    }
    return summaries;
  };

  auto first = make_cells();
  auto repeated = first;
  const auto before = summarize(first);
  coulomb::ParaClass parameters;
  parameters.dt = 0.01;
  coulomb::RandomContext random;

  const auto collide = [&](auto &cells) {
#pragma omp parallel for num_threads(thread_count) schedule(static, 1)
    for (int cell = 0; cell < thread_count; ++cell) {
      coulomb::CollisionOperator(parameters, random)
          .collide_homogeneous(cells[cell],
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
  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8};

  const auto options = coulomb::RunOptions::parse(9, argv);

  REQUIRE(options.seed == 12345);
  REQUIRE(options.steps == 10);
  REQUIRE(options.threads == 1);
  REQUIRE(options.output_directory == "run-test");
}

TEST_CASE("negpar.unit.cli.run options reject invalid values",
          "[cli][validation]") {
  char arg0[] = "negpar_inhomo";
  char arg1[] = "--threads";
  char arg2[] = "0";
  char *argv[] = {arg0, arg1, arg2};

  REQUIRE_THROWS_AS(coulomb::RunOptions::parse(3, argv), std::invalid_argument);
}

TEST_CASE("negpar.unit.cli.run options reject malformed numeric values",
          "[cli][validation]") {
  char arg0[] = "negpar_inhomo";
  char arg1[] = "--steps";
  char arg2[] = "10steps";
  char *argv[] = {arg0, arg1, arg2};

  REQUIRE_THROWS_AS(coulomb::RunOptions::parse(3, argv), std::invalid_argument);
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

  coulomb::RunOptions::reset_runtime_state(state);

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
  REQUIRE(paths.format_index(7) == "_007");
  REQUIRE(paths.format_index(-4) == "_-004");
  REQUIRE(paths.format_index(12, 2) == "_12");
}

TEST_CASE(
    "negpar.unit.output.fixed-bin histogram preserves legacy edge clamping",
    "[output][histogram]") {
  const std::vector<double> values{-1.0, 0.0, 0.49, 0.5, 0.99, 1.0, 2.0};
  std::vector<int> counts(2);

  coulomb::Histogram{}.fixed_bins(values, counts, 0.0, 1.0);

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

  coulomb::MacroOutput{}.save_macro(std::vector<double>{1.0 / 3.0, -2.0},
                                    "values", state);
  coulomb::MacroOutput{}.save_2d(2, 2, {{1.0, 2.0}, {3.0, 4.0}}, "matrix",
                                 state);
  std::complex<double> spectrum[] = {{1.0, -2.0}, {3.0, 4.0}};
  coulomb::MacroOutput{}.save_complex(2, spectrum, "spectrum", state);

  std::ifstream file(directory / "values_007.txt");
  REQUIRE(file.good());
  std::string contents((std::istreambuf_iterator<char>(file)),
                       std::istreambuf_iterator<char>());
  file.close();
  const bool matrix_exists =
      std::filesystem::exists(directory / "matrix_007.txt");
  const bool real_exists =
      std::filesystem::exists(directory / "spectrum_r.txt");
  const bool imaginary_exists =
      std::filesystem::exists(directory / "spectrum_i.txt");
  std::vector<std::string> filenames;
  for (const auto &entry : std::filesystem::directory_iterator(directory))
    filenames.push_back(entry.path().filename().string());
  std::sort(filenames.begin(), filenames.end());

  std::ifstream matrix_file(directory / "matrix_007.txt");
  std::string matrix_contents((std::istreambuf_iterator<char>(matrix_file)),
                              std::istreambuf_iterator<char>());
  std::ifstream real_file(directory / "spectrum_r.txt");
  std::string real_contents((std::istreambuf_iterator<char>(real_file)),
                            std::istreambuf_iterator<char>());
  std::ifstream imaginary_file(directory / "spectrum_i.txt");
  std::string imaginary_contents(
      (std::istreambuf_iterator<char>(imaginary_file)),
      std::istreambuf_iterator<char>());
  matrix_file.close();
  real_file.close();
  imaginary_file.close();
  bool invalid_dimensions = false;
  try {
    coulomb::MacroOutput{}.save_2d(1, 2, {{1.0}}, "invalid", state);
  } catch (const std::invalid_argument &) {
    invalid_dimensions = true;
  }
  std::filesystem::remove_all(directory);

  REQUIRE(contents.find("0.333333333333333") != std::string::npos);
  REQUIRE(contents.find("-2") != std::string::npos);
  REQUIRE(matrix_exists);
  REQUIRE(real_exists);
  REQUIRE(imaginary_exists);
  REQUIRE((filenames ==
           std::vector<std::string>{"matrix_007.txt", "spectrum_i.txt",
                                    "spectrum_r.txt", "values_007.txt"}));
  REQUIRE(matrix_contents == "1 2 \n3 4 \n");
  REQUIRE(real_contents == "1\n3\n");
  REQUIRE(imaginary_contents == "-2\n4\n");
  REQUIRE(invalid_dimensions);
  REQUIRE_THROWS_AS(
      coulomb::ParticleOutput{}.save_homogeneous_radial_distribution(0, state),
      std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::OutputPaths(state).resolve("../escape.txt"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      coulomb::OutputPaths(state).resolve(
          (std::filesystem::temp_directory_path() / "escape.txt").string()),
      std::invalid_argument);
}
