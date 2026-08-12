#include "NegativeParticleCollisions.h"

#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>

#include "Collisions.h"
#include "Grid.h"
#include "NegativeParticleSampling.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {

void NegativeParticleCollisions::collide_with_full(NeParticleGroup &groups) {
  const auto &parameters = parameters_;
  auto &random = random_;
  const int full_count = groups.size(ParticleKind::Full);
  const int positive_count = groups.size(ParticleKind::Positive);
  const int negative_count = groups.size(ParticleKind::Negative);
  if (full_count < positive_count + negative_count) {
    std::cout << "Too few F particles." << std::endl;
    std::cout << "(" << positive_count << ", " << negative_count << ", "
              << full_count << ") " << std::endl;
  }

  auto &positive = groups.list(ParticleKind::Positive);
  auto &negative = groups.list(ParticleKind::Negative);
  auto &full = groups.list(ParticleKind::Full);
  const auto permutation = RandomSampling(random).permutation(
      full_count, positive_count + negative_count);

  for (int index = 0; index < positive_count; ++index) {
    const int full_index = permutation[index] - 1;
    const auto velocities = CollisionOperator(parameters, random)
                                .collide_pair(positive[index].velocity(),
                                              full[full_index].velocity());
    positive[index].set_velocity(velocities.first);
  }
  for (int index = 0; index < negative_count; ++index) {
    const int full_index = permutation[index + positive_count] - 1;
    const auto velocities = CollisionOperator(parameters, random)
                                .collide_pair(negative[index].velocity(),
                                              full[full_index].velocity());
    negative[index].set_velocity(velocities.first);
  }
}

void NegativeParticleCollisions::collide_homogeneous(NeParticleGroup &S_x) {
  const auto &para = parameters_;
  const double Neff = grid_.Neff;
  auto &random = random_;
  NeParticleGroup S_x_new;

  NegativeParticleSampling{}.sample_delta(S_x, S_x_new, para, Neff, random);
  ParticleGroupOperations{}.assign_positions(S_x_new, S_x.get_xmin(),
                                             S_x.get_xmax(), random);
  collide_with_full(S_x);
  ParticleGroupOperations{}.merge_signed(S_x, S_x_new);

  auto &Sf = S_x.list(ParticleKind::Full);
  CollisionOperator(para, random)
      .collide_homogeneous(Sf, S_x.size(ParticleKind::Full));
}

void NegativeParticleCollisions::collide(std::vector<NeParticleGroup> &S_x) {
  const auto &grid = grid_;
  const auto &para = parameters_;
  NegativeParticleSampling{}.update_bounds(S_x, grid, para);
  for (int kx = 0; kx < grid.Nx; kx++)
    collide_homogeneous(S_x[kx]);
}

void NegativeParticleCollisions::collide_parallel(
    std::vector<NeParticleGroup> &S_x) {
  const auto &grid = grid_;
  const auto &para = parameters_;
  NegativeParticleSampling{}.update_bounds(S_x, grid, para);
#pragma omp parallel if (para.FLAG_USE_OPENMP)
  {
#pragma omp for
    for (int kx = 0; kx < grid.Nx; kx++)
      collide_homogeneous(S_x[kx]);
  }
}

void NegativeParticleCollisions::collide_bgk_homogeneous(
    NeParticleGroup &groups) {
  auto &parameters = parameters_;
  auto &random = random_;
  const int positive_count = groups.size(ParticleKind::Positive);
  const int negative_count = groups.size(ParticleKind::Negative);
  const int full_count = groups.size(ParticleKind::Full);
  const int positive_remove = RandomSampling(random).stochastic_floor(
      positive_count * (parameters.dt * parameters.coeff_binarycoll));
  const int negative_remove = RandomSampling(random).stochastic_floor(
      negative_count * (parameters.dt * parameters.coeff_binarycoll));

  for (int index = 0; index < positive_remove; ++index) {
    const int remove_index = static_cast<int>(
        RandomSampling(random).uniform() * groups.size(ParticleKind::Positive));
    groups.erase(remove_index, ParticleKind::Positive);
  }
  for (int index = 0; index < negative_remove; ++index) {
    const int remove_index = static_cast<int>(
        RandomSampling(random).uniform() * groups.size(ParticleKind::Negative));
    groups.erase(remove_index, ParticleKind::Negative);
  }

  const double change_rate = parameters.dt * parameters.coeff_binarycoll;
  const double sqrt_temperature = std::sqrt(groups.TprtM);
  std::array<double, 3> velocity{};
  auto &full = groups.list(ParticleKind::Full);
  for (int index = 0; index < full_count; ++index) {
    if (RandomSampling(random).uniform() < change_rate) {
      velocity[0] =
          groups.u1M + sqrt_temperature * RandomSampling(random).normal();
      velocity[1] =
          groups.u2M + sqrt_temperature * RandomSampling(random).normal();
      velocity[2] =
          groups.u3M + sqrt_temperature * RandomSampling(random).normal();
      full[index].set_velocity(velocity);
    }
  }
}

void NegativeParticleCollisions::collide_bgk(
    std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell)
    collide_bgk_homogeneous(groups[cell]);
}

} // namespace coulomb
