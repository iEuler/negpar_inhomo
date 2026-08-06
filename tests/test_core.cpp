#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <filesystem>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <omp.h>

#include "Classes.h"
#include "Advection.h"
#include "Collisions.h"
#include "Diagnostics.h"
#include "FFT.h"
#include "Moments.h"
#include "Numerics.h"
#include "Output.h"
#include "RunOptions.h"
#include "_global_variables.h"
#include "utils.h"

namespace {

using coulomb::Particle1d3d;

Particle1d3d particle(double x, double vx, double vy, double vz) {
  return Particle1d3d(x, {vx, vy, vz});
}

}  // namespace

TEST_CASE("default grid is initialized", "[grid]") {
  coulomb::NumericGridClass grid;

  REQUIRE(grid.Nx == 100);
  REQUIRE(grid.Nv == 200);
  REQUIRE(grid.x.size() == static_cast<size_t>(grid.Nx));
  REQUIRE(grid.vx.size() == static_cast<size_t>(grid.Nv));
  REQUIRE(grid.dx > 0.0);
  REQUIRE(grid.dv > 0.0);
  REQUIRE(grid.dt > 0.0);
}

TEST_CASE("grid rejects invalid configuration", "[grid][validation]") {
  REQUIRE_THROWS_AS(coulomb::NumericGridClass(0), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::NumericGridClass(-4), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::NumericGridClass(10, "invalid"),
                    std::invalid_argument);
}

TEST_CASE("typed collision mode preserves the legacy name", "[parameters]") {
  const coulomb::ParaClass parameters;
  coulomb::RandomContext random;
  coulomb::reseed_random(random, 2468);
  REQUIRE(parameters.method_binarycoll ==
          coulomb::BinaryCollisionMethod::TA);
  REQUIRE(std::string(coulomb::binary_collision_name(
              parameters.method_binarycoll)) == "TA");
}

TEST_CASE("particle advection wraps periodic boundaries", "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = 'p';
  grid.dt = 1.0;
  const double start = grid.xmax - 0.25;
  const auto moved = coulomb::moveparticle(particle(start, 1.0, 2.0, 3.0),
                                            0.0, grid, state);

  REQUIRE(moved.position() == Catch::Approx(grid.xmin + 0.75));
  REQUIRE(moved.velocity(0) == Catch::Approx(1.0));
  REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("particle advection reflects at nonperiodic boundaries",
          "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = 'n';
  grid.dt = 1.0;
  const auto moved = coulomb::moveparticle(particle(grid.xmax - 0.25, 1.0,
                                                     2.0, 3.0),
                                            0.5, grid, state);

  REQUIRE(moved.position() == Catch::Approx(grid.xmax - 0.75));
  REQUIRE(moved.velocity(0) == Catch::Approx(-1.5));
  REQUIRE(moved.velocity(1) == Catch::Approx(2.0));
}

TEST_CASE("particle exactly at the upper boundary belongs to the last cell",
          "[advection][boundary]") {
  coulomb::NumericGridClass grid(4);
  coulomb::SimulationState state;
  grid.bdry_x = 'n';
  grid.dt = grid.xmax - (grid.xmax - 0.5);
  auto moved = coulomb::moveparticle(
      particle(grid.xmax - 0.5, 1.0, 0.0, 0.0), 0.0, grid, state);

  REQUIRE(moved.position() == Catch::Approx(grid.xmax));
  REQUIRE(coulomb::findparticlegroup(moved, grid) == grid.Nx - 1);
}

TEST_CASE("particle group lookup rejects positions outside the grid",
          "[advection][validation]") {
  coulomb::NumericGridClass grid(4);
  auto above = particle(grid.xmax + 1.0, 0.0, 0.0, 0.0);
  auto below = particle(grid.xmin - 1.0, 0.0, 0.0, 0.0);

  REQUIRE_THROWS_AS(coulomb::findparticlegroup(above, grid),
                    std::out_of_range);
  REQUIRE_THROWS_AS(coulomb::findparticlegroup(below, grid),
                    std::out_of_range);
}

TEST_CASE("numerics error function preserves the approximation", "[numerics]") {
  REQUIRE(std::abs(coulomb::myerf(0.0)) < 2e-9);
  REQUIRE(std::abs(coulomb::myerf(1.0) - 0.84270079) < 2e-7);
  REQUIRE(std::abs(coulomb::myerf(-1.0) + 0.84270079) < 2e-7);
}

TEST_CASE("one-dimensional shifts preserve legacy boundary conventions",
          "[numerics][boundaries]") {
  const std::vector<double> values = {1.0, 2.0, 3.0};
  const std::vector<double> shifted = {2.0, 3.0, 1.0};
  const std::vector<double> extended = {-1.0, 1.0, 2.0};
  const std::vector<double> difference = {-0.5, 1.0, -0.5};
  REQUIRE(coulomb::cshift_1d(values, 3, 1) == shifted);
  REQUIRE(coulomb::cshift_1d(values, 3, 4) == shifted);
  REQUIRE(coulomb::eoshift_1d(values, 3, -1, -1.0) == extended);
  REQUIRE(coulomb::diff_1d_central(values, 3, 'p') == difference);
  REQUIRE(coulomb::cshift_1d(std::vector<double>{}, 0, 8).empty());
  REQUIRE_THROWS_AS(coulomb::cshift_1d(values, 2, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::diff_1d_central(values, 3, 'x'),
                    std::invalid_argument);
}

TEST_CASE("kinetic Euler fluxes handle periodic and constant states",
          "[numerics][euler]") {
  const std::vector<double> values = {1.0, 2.0, 4.0};
  REQUIRE(coulomb::limiter_x1_o2(values, 3, 2.0, 1.0, 1, values, 'p') ==
          std::vector<double>{-1.5, 0.5, 1.0});
  REQUIRE(coulomb::limiter_x1_o2(values, 3, 2.0, 1.0, -1, values, 'p') ==
          std::vector<double>{0.5, 1.0, -1.5});

  std::vector<double> density_change(3), momentum_change(3), energy_change(3);
  const std::vector<double> constant(3, 1.0);
  coulomb::Euler_kinetic_x1(constant, constant, constant, 3, 1.0, 0.1, 'p',
                            density_change, momentum_change, energy_change);
  REQUIRE(density_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE(momentum_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE(energy_change == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE_THROWS_AS(coulomb::Euler_kinetic_x1(
                        constant, constant, std::vector<double>(3, 0.0), 3,
                        1.0, 0.1, 'p', density_change, momentum_change,
                        energy_change),
                    std::invalid_argument);
}

TEST_CASE("FFT reshape helpers preserve row-major ordering", "[fft][numerics]") {
  const std::vector<int> flat = {0, 1, 2, 3, 4, 5, 6, 7};
  const auto shaped = coulomb::reshape1dTo3d(flat, 2, 2, 2);

  REQUIRE(shaped[0][0][0] == 0);
  REQUIRE(shaped[0][1][1] == 3);
  REQUIRE(shaped[1][0][0] == 4);
  REQUIRE(coulomb::reshape3dTo1d(shaped) == flat);
  REQUIRE_THROWS_AS(coulomb::reshape1dTo3d(flat, 3, 2, 2),
                    std::invalid_argument);
}

TEST_CASE("FFT reshape helper rejects ragged input", "[fft][validation]") {
  const std::vector<std::vector<std::vector<int>>> ragged = {
      {{1, 2}, {3}}, {{4, 5}, {6, 7}}};
  REQUIRE_THROWS_AS(coulomb::reshape3dTo1d(ragged), std::invalid_argument);
}

TEST_CASE("FFTW wrappers have explicit ownership and validate dimensions",
          "[fft][resource]") {
  static_assert(!std::is_copy_constructible_v<coulomb::FFT1D>);
  static_assert(!std::is_copy_assignable_v<coulomb::FFT1D>);
  static_assert(!std::is_copy_constructible_v<coulomb::FFT3D>);
  static_assert(!std::is_copy_assignable_v<coulomb::FFT3D>);

  REQUIRE_THROWS_AS(coulomb::FFT1D(0), std::invalid_argument);
  REQUIRE_THROWS_AS(coulomb::FFT3D(0, 2, 2), std::invalid_argument);
}

TEST_CASE("particle group computes moments", "[particles][moments]") {
  coulomb::ParticleGroup group;
  group.push_back(particle(0.0, 1.0, 2.0, 3.0));
  group.push_back(particle(0.0, -1.0, 4.0, 0.0));

  group.computemoments();

  REQUIRE(group.size() == 2);
  REQUIRE(group.m0 == 2.0);
  REQUIRE(group.m11 == 0.0);
  REQUIRE(group.m12 == 6.0);
  REQUIRE(group.m13 == 3.0);
  REQUIRE(group.m2 == 31.0);

  group.erase(0);
  REQUIRE(group.size() == 1);
}

TEST_CASE("moment conversions round trip", "[moments]") {
  double velocity[3] = {1.0, -2.0, 0.5};
  double momentum[3] = {};
  double energy = 0.0;
  coulomb::uT2mE(2.0, velocity, 3.0, momentum, energy);

  REQUIRE(momentum[0] == 2.0);
  REQUIRE(momentum[1] == -4.0);
  REQUIRE(momentum[2] == 1.0);
  REQUIRE(energy == Catch::Approx(2.0 * (1.0 + 4.0 + 0.25 + 9.0) / 2.0));

  double recovered_velocity[3] = {};
  double recovered_temperature = 0.0;
  coulomb::mE2uT(2.0, momentum, energy, recovered_velocity,
                 recovered_temperature);
  for (int component = 0; component < 3; ++component)
    REQUIRE(recovered_velocity[component] == Catch::Approx(velocity[component]));
  REQUIRE(recovered_temperature == Catch::Approx(3.0));
}

TEST_CASE("macroscopic moment updates reconstruct particle fields",
          "[moments][reconstruction]") {
  coulomb::NumericGridClass grid(3);
  std::vector<coulomb::NeParticleGroup> groups(3);
  for (auto& group : groups) {
    group.rhoM = 1.0;
    group.u1M = 0.0;
    group.TprtM = 1.0;
    group.computemoments();
  }

  coulomb::update_rhouT(groups, grid);
  for (const auto& group : groups) {
    REQUIRE(group.rho == Catch::Approx(1.0));
    REQUIRE(group.u1 == Catch::Approx(0.0));
    REQUIRE(group.Tprt == Catch::Approx(1.0));
  }

  groups[0].push_back(particle(0.0, 2.0, 0.0, 0.0),
                      coulomb::ParticleKind::Full);
  for (auto& group : groups) group.computemoments();
  coulomb::update_rhouT_F(groups, grid);
  REQUIRE(groups[0].rhoF == Catch::Approx(grid.Neff_F / grid.dx));
  REQUIRE(groups[0].u1F == Catch::Approx(2.0));
}

TEST_CASE("Maxwellian update applies conservative moment changes",
          "[moments][maxwellian]") {
  coulomb::NumericGridClass grid(1);
  std::vector<coulomb::NeParticleGroup> groups(1);
  groups[0].rhoM = 2.0;
  groups[0].u1M = 1.0;
  groups[0].TprtM = 3.0;
  groups[0].drho = 0.5;
  groups[0].dm1 = 0.25;
  groups[0].denergy = 1.0;

  coulomb::update_maxwellian(groups, grid);

  REQUIRE(groups[0].rhoM == Catch::Approx(1.5));
  REQUIRE(groups[0].u1M == Catch::Approx(1.75 / 1.5));
  REQUIRE(groups[0].TprtM == Catch::Approx(3.5462962963));
}

TEST_CASE("binary collision preserves pair momentum and kinetic energy",
          "[collisions][invariants]") {
  // A deterministic invariant check; particle ordering is not involved.
  coulomb::ParaClass parameters;
  parameters.dt = 0.01;
  coulomb::RandomContext random;
  coulomb::reseed_random(random, 2468);

  const std::vector<double> first = {1.0, -0.5, 0.75};
  const std::vector<double> second = {-0.25, 1.5, -1.0};
  const auto collided =
      coulomb::coulombBinary3d(first, second, parameters, random);

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

TEST_CASE("binary collision is reproducible for an explicit seed",
          "[collisions][reproducibility]") {
  const coulomb::ParaClass parameters;
  const std::vector<double> first = {0.5, 1.0, -0.75};
  const std::vector<double> second = {-1.0, 0.25, 1.5};
  coulomb::RandomContext random;

  coulomb::reseed_random(random, 13579);
  const auto expected =
      coulomb::coulombBinary3d(first, second, parameters, random);
  coulomb::reseed_random(random, 13579);
  const auto repeated =
      coulomb::coulombBinary3d(first, second, parameters, random);

  REQUIRE(repeated.first == expected.first);
  REQUIRE(repeated.second == expected.second);
}

TEST_CASE("positive, negative, and full particle groups are independent",
          "[particles]") {
  coulomb::NeParticleGroup group;
  group.push_back(particle(0.0, 1.0, 0.0, 0.0),
                  coulomb::ParticleKind::Positive);
  group.push_back(particle(0.0, 2.0, 0.0, 0.0),
                  coulomb::ParticleKind::Negative);
  group.push_back(particle(0.0, 3.0, 0.0, 0.0), coulomb::ParticleKind::Full);

  REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
  REQUIRE(group.size(coulomb::ParticleKind::Negative) == 1);
  REQUIRE(group.size(coulomb::ParticleKind::Full) == 1);

  group.clear(coulomb::ParticleKind::Negative);
  REQUIRE(group.size(coulomb::ParticleKind::Positive) == 1);
  REQUIRE(group.size(coulomb::ParticleKind::Negative) == 0);
  REQUIRE(group.size(coulomb::ParticleKind::Full) == 1);

  REQUIRE(coulomb::particle_kind_from_code('p') ==
          coulomb::ParticleKind::Positive);
  REQUIRE(coulomb::particle_kind_code(coulomb::ParticleKind::Full) == 'f');
  REQUIRE_THROWS_AS(group.size('x'), std::invalid_argument);
  REQUIRE_THROWS_AS(group.clear('x'), std::invalid_argument);
}

TEST_CASE("particle count diagnostics validate the requested grid size",
          "[diagnostics][particles]") {
  std::vector<coulomb::NeParticleGroup> groups(2);
  groups[0].push_back(particle(0.0, 1.0, 0.0, 0.0),
                      coulomb::ParticleKind::Positive);
  groups[1].push_back(particle(0.0, 2.0, 0.0, 0.0),
                      coulomb::ParticleKind::Negative);

  REQUIRE(coulomb::count_particle_number(groups, 2, 'p') == 1);
  REQUIRE(coulomb::count_particle_number(groups, 2,
                                          coulomb::ParticleKind::Negative) ==
          1);
  REQUIRE_THROWS_AS(coulomb::count_particle_number(groups, 3, 'p'),
                    std::invalid_argument);
}

TEST_CASE("random permutation validates its requested size", "[utils]") {
  coulomb::RandomContext random;
  coulomb::reseed_random(random, 1234);
  REQUIRE_THROWS(coulomb::myrandperm(2, 3, random));

  const auto permutation = coulomb::myrandperm(10, 5, random);
  REQUIRE(permutation.size() == 5);
  for (const int value : permutation) REQUIRE(value >= 1);
  for (const int value : permutation) REQUIRE(value <= 10);
}

TEST_CASE("random generator can be explicitly seeded", "[utils][reproducibility]") {
  coulomb::RandomContext random;
  coulomb::reseed_random(random, 12345);
  const double first = coulomb::myrand(random);
  const double normal = coulomb::myrandn(random);

  coulomb::reseed_random(random, 12345);
  REQUIRE(coulomb::myrand(random) == first);
  REQUIRE(coulomb::myrandn(random) == normal);
}

TEST_CASE("random context owns an independently reproducible seed",
          "[utils][reproducibility]") {
  coulomb::RandomContext first_context;
  coulomb::RandomContext second_context;
  coulomb::reseed_random(first_context, 24680);
  coulomb::reseed_random(second_context, 24680);

  REQUIRE(coulomb::myrand(first_context) == coulomb::myrand(second_context));
  REQUIRE(coulomb::myrandn(first_context) == coulomb::myrandn(second_context));
}

TEST_CASE("random generator is reproducible per OpenMP thread",
          "[utils][openmp][reproducibility]") {
  const int thread_count = omp_get_max_threads();
  std::vector<double> first(thread_count), second(thread_count);
  coulomb::RandomContext random;

  coulomb::reseed_random(random, 9876);
#pragma omp parallel for
  for (int thread = 0; thread < thread_count; ++thread) {
    first[thread] = coulomb::myrand(random);
    second[thread] = coulomb::myrandn(random);
  }

  std::vector<double> first_again(thread_count), second_again(thread_count);
  coulomb::reseed_random(random, 9876);
#pragma omp parallel for
  for (int thread = 0; thread < thread_count; ++thread) {
    first_again[thread] = coulomb::myrand(random);
    second_again[thread] = coulomb::myrandn(random);
  }

  REQUIRE(first_again == first);
  REQUIRE(second_again == second);
}

TEST_CASE("run options parse seed, threads, and output directory", "[cli]") {
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

  const auto options = coulomb::parse_run_options(9, argv);

  REQUIRE(options.seed == 12345);
  REQUIRE(options.steps == 10);
  REQUIRE(options.threads == 1);
  REQUIRE(options.output_directory == "run-test");
}

TEST_CASE("run options reject invalid values", "[cli][validation]") {
  char arg0[] = "negpar_inhomo";
  char arg1[] = "--threads";
  char arg2[] = "0";
  char* argv[] = {arg0, arg1, arg2};

  REQUIRE_THROWS_AS(coulomb::parse_run_options(3, argv), std::invalid_argument);
}

TEST_CASE("run options reject malformed numeric values", "[cli][validation]") {
  char arg0[] = "negpar_inhomo";
  char arg1[] = "--steps";
  char arg2[] = "10steps";
  char* argv[] = {arg0, arg1, arg2};

  REQUIRE_THROWS_AS(coulomb::parse_run_options(3, argv), std::invalid_argument);
}

TEST_CASE("runtime state reset clears cross-run counters", "[cli][state]") {
  coulomb::SimulationState state;
  state.saveIndex = 9;
  state.filenameWithNumber = true;
  state.saveFlux = true;
  state.movedCount = 4;
  state.resampleCount = 5;
  state.syncTime = 1.25;

  coulomb::reset_runtime_state(state);

  REQUIRE(state.saveIndex == 0);
  REQUIRE_FALSE(state.filenameWithNumber);
  REQUIRE_FALSE(state.saveFlux);
  REQUIRE(state.movedCount == 0);
  REQUIRE(state.resampleCount == 0);
  REQUIRE(state.syncTime == 0.0);
}

TEST_CASE("output filename numbering preserves legacy format", "[output]") {
  REQUIRE(coulomb::int2str(7) == "_007");
  REQUIRE(coulomb::int2str(-4) == "_-004");
  REQUIRE(coulomb::int2str(12, 2) == "_12");
}

TEST_CASE("macro output writer preserves numbered precision", "[output]") {
  const auto directory =
      std::filesystem::temp_directory_path() / "negpar_save_macro_test";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  coulomb::SimulationState state;
  state.outputDirectory = directory.string();
  state.filenameWithNumber = true;
  state.saveIndex = 7;

  coulomb::save_macro(std::vector<double>{1.0 / 3.0, -2.0}, "values", state);
  coulomb::save_2d(2, 2, {{1.0, 2.0}, {3.0, 4.0}}, "matrix", state);
  std::complex<double> spectrum[] = {{1.0, -2.0}, {3.0, 4.0}};
  coulomb::save_complex(2, spectrum, "spectrum", state);

  std::ifstream file(directory / "values_007.txt");
  REQUIRE(file.good());
  std::string contents((std::istreambuf_iterator<char>(file)),
                       std::istreambuf_iterator<char>());
  file.close();
  const bool matrix_exists = std::filesystem::exists(directory / "matrix_007.txt");
  const bool real_exists = std::filesystem::exists(directory / "spectrum_r.txt");
  const bool imaginary_exists =
      std::filesystem::exists(directory / "spectrum_i.txt");
  bool invalid_dimensions = false;
  try {
    coulomb::save_2d(1, 2, {{1.0}}, "invalid", state);
  } catch (const std::invalid_argument&) {
    invalid_dimensions = true;
  }
  std::filesystem::remove_all(directory);

  REQUIRE(contents.find("0.333333333333333") != std::string::npos);
  REQUIRE(contents.find("-2") != std::string::npos);
  REQUIRE(matrix_exists);
  REQUIRE(real_exists);
  REQUIRE(imaginary_exists);
  REQUIRE(invalid_dimensions);
  REQUIRE_THROWS_AS(coulomb::save_homo_rdist(0, state), std::invalid_argument);
}
