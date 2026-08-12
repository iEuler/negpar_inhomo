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

void NegativeParticleCollisions::collide_with_full(NeParticleGroup &groups,
                                                   const ParaClass &parameters,
                                                   RandomContext &random) {
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
  const auto permutation = RandomSampling::permutation(
      full_count, positive_count + negative_count, random);

  for (int index = 0; index < positive_count; ++index) {
    const int full_index = permutation[index] - 1;
    const auto velocities = CollisionOperator::collide_pair(
        positive[index].velocity(), full[full_index].velocity(), parameters,
        random);
    positive[index].set_velocity(velocities.first);
  }
  for (int index = 0; index < negative_count; ++index) {
    const int full_index = permutation[index + positive_count] - 1;
    const auto velocities = CollisionOperator::collide_pair(
        negative[index].velocity(), full[full_index].velocity(), parameters,
        random);
    negative[index].set_velocity(velocities.first);
  }
}

void NegativeParticleCollisions::collide_homogeneous(NeParticleGroup &S_x,
                                                     const ParaClass &para,
                                                     double Neff,
                                                     RandomContext &random) {
  NeParticleGroup S_x_new;

  NegativeParticleSampling::sample_delta(S_x, S_x_new, para, Neff, random);
  ParticleGroupOperations::assign_positions(S_x_new, S_x.get_xmin(),
                                            S_x.get_xmax(), random);
  NegativeParticleCollisions::collide_with_full(S_x, para, random);
  ParticleGroupOperations::merge_signed(S_x, S_x_new);

  auto &Sf = S_x.list(ParticleKind::Full);
  CollisionOperator::collide_homogeneous(Sf, S_x.size(ParticleKind::Full), para,
                                         random);
}

void NegativeParticleCollisions::collide(std::vector<NeParticleGroup> &S_x,
                                         const NumericGridClass &grid,
                                         const ParaClass &para,
                                         RandomContext &random) {
  NegativeParticleSampling::update_bounds(S_x, grid, para);
  for (int kx = 0; kx < grid.Nx; kx++)
    NegativeParticleCollisions::collide_homogeneous(S_x[kx], para, grid.Neff,
                                                    random);
}

void NegativeParticleCollisions::collide_parallel(
    std::vector<NeParticleGroup> &S_x, const NumericGridClass &grid,
    const ParaClass &para, RandomContext &random) {
  NegativeParticleSampling::update_bounds(S_x, grid, para);
#pragma omp parallel if (para.FLAG_USE_OPENMP)
  {
#pragma omp for
    for (int kx = 0; kx < grid.Nx; kx++)
      NegativeParticleCollisions::collide_homogeneous(S_x[kx], para, grid.Neff,
                                                      random);
  }
}

void NegativeParticleCollisions::collide_bgk_homogeneous(
    NeParticleGroup &groups, ParaClass &parameters, RandomContext &random) {
  const int positive_count = groups.size(ParticleKind::Positive);
  const int negative_count = groups.size(ParticleKind::Negative);
  const int full_count = groups.size(ParticleKind::Full);
  const int positive_remove = RandomSampling::stochastic_floor(
      positive_count * (parameters.dt * parameters.coeff_binarycoll), random);
  const int negative_remove = RandomSampling::stochastic_floor(
      negative_count * (parameters.dt * parameters.coeff_binarycoll), random);

  for (int index = 0; index < positive_remove; ++index) {
    const int remove_index = static_cast<int>(
        RandomSampling::uniform(random) * groups.size(ParticleKind::Positive));
    groups.erase(remove_index, ParticleKind::Positive);
  }
  for (int index = 0; index < negative_remove; ++index) {
    const int remove_index = static_cast<int>(
        RandomSampling::uniform(random) * groups.size(ParticleKind::Negative));
    groups.erase(remove_index, ParticleKind::Negative);
  }

  const double change_rate = parameters.dt * parameters.coeff_binarycoll;
  const double sqrt_temperature = std::sqrt(groups.TprtM);
  std::array<double, 3> velocity{};
  auto &full = groups.list(ParticleKind::Full);
  for (int index = 0; index < full_count; ++index) {
    if (RandomSampling::uniform(random) < change_rate) {
      velocity[0] =
          groups.u1M + sqrt_temperature * RandomSampling::normal(random);
      velocity[1] =
          groups.u2M + sqrt_temperature * RandomSampling::normal(random);
      velocity[2] =
          groups.u3M + sqrt_temperature * RandomSampling::normal(random);
      full[index].set_velocity(velocity);
    }
  }
}

void NegativeParticleCollisions::collide_bgk(
    std::vector<NeParticleGroup> &groups, NumericGridClass &grid,
    ParaClass &parameters, RandomContext &random) {
  for (int cell = 0; cell < grid.Nx; ++cell)
    NegativeParticleCollisions::collide_bgk_homogeneous(groups[cell],
                                                        parameters, random);
}

} // namespace coulomb
