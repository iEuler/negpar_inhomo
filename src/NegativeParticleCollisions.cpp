#include "NegativeParticleCollisions.h"

#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>

#include "Collisions.h"
#include "Grid.h"
#include "NegativeParticleSampling.h"
#include "ParticleGroupOperations.h"
#include "ParticleGroup.h"
#include "SimulationConfig.h"
#include "RandomSampling.h"

namespace coulomb {

void coulomb_collision_homo_PFNF(NeParticleGroup& groups,
                                 const ParaClass& parameters,
                                 RandomContext& random) {
  const int full_count = groups.size(ParticleKind::Full);
  const int positive_count = groups.size(ParticleKind::Positive);
  const int negative_count = groups.size(ParticleKind::Negative);
  if (full_count < positive_count + negative_count) {
    std::cout << "Too few F particles." << std::endl;
    std::cout << "(" << positive_count << ", " << negative_count << ", "
              << full_count << ") " << std::endl;
  }

  auto& positive = groups.list(ParticleKind::Positive);
  auto& negative = groups.list(ParticleKind::Negative);
  auto& full = groups.list(ParticleKind::Full);
  const auto permutation =
      myrandperm(full_count, positive_count + negative_count, random);

  for (int index = 0; index < positive_count; ++index) {
    const int full_index = permutation[index] - 1;
    const auto velocities = coulombBinary3d(
        positive[index].velocity(), full[full_index].velocity(), parameters,
        random);
    positive[index].set_velocity(velocities.first);
  }
  for (int index = 0; index < negative_count; ++index) {
    const int full_index = permutation[index + positive_count] - 1;
    const auto velocities = coulombBinary3d(
        negative[index].velocity(), full[full_index].velocity(), parameters,
        random);
    negative[index].set_velocity(velocities.first);
  }
}

void NegPar_collision_homo(NeParticleGroup &S_x, const ParaClass &para,
                           double Neff, RandomContext& random) {
  NeParticleGroup S_x_new;

  samplefromDeltam(S_x, S_x_new, para, Neff, random);
  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax(), random);
  coulomb_collision_homo_PFNF(S_x, para, random);
  merge_NeParticleGroup(S_x, S_x_new);

  auto &Sf = S_x.list(ParticleKind::Full);
  coulomb_collision_homo(Sf, S_x.size(ParticleKind::Full), para, random);
}

void NegPar_collision(std::vector<NeParticleGroup> &S_x,
                      const NumericGridClass &grid, const ParaClass &para,
                      RandomContext& random) {
  finddeltambound_inhomo(S_x, grid, para);
  for (int kx = 0; kx < grid.Nx; kx++)
    NegPar_collision_homo(S_x[kx], para, grid.Neff, random);
}

void NegPar_collision_openmp(std::vector<NeParticleGroup> &S_x,
                             const NumericGridClass &grid,
                             const ParaClass &para, RandomContext& random) {
  finddeltambound_inhomo(S_x, grid, para);
#pragma omp parallel if (para.FLAG_USE_OPENMP)
  {
#pragma omp for
    for (int kx = 0; kx < grid.Nx; kx++)
      NegPar_collision_homo(S_x[kx], para, grid.Neff, random);
  }
}

void NegPar_BGK_collision_homo(NeParticleGroup& groups,
                               ParaClass& parameters,
                               RandomContext& random) {
  const int positive_count = groups.size(ParticleKind::Positive);
  const int negative_count = groups.size(ParticleKind::Negative);
  const int full_count = groups.size(ParticleKind::Full);
  const int positive_remove = myfloor(
      positive_count * (parameters.dt * parameters.coeff_binarycoll), random);
  const int negative_remove = myfloor(
      negative_count * (parameters.dt * parameters.coeff_binarycoll), random);

  for (int index = 0; index < positive_remove; ++index) {
    const int remove_index = static_cast<int>(
        myrand(random) * groups.size(ParticleKind::Positive));
    groups.erase(remove_index, ParticleKind::Positive);
  }
  for (int index = 0; index < negative_remove; ++index) {
    const int remove_index = static_cast<int>(
        myrand(random) * groups.size(ParticleKind::Negative));
    groups.erase(remove_index, ParticleKind::Negative);
  }

  const double change_rate = parameters.dt * parameters.coeff_binarycoll;
  const double sqrt_temperature = std::sqrt(groups.TprtM);
  std::array<double, 3> velocity{};
  auto& full = groups.list(ParticleKind::Full);
  for (int index = 0; index < full_count; ++index) {
    if (myrand(random) < change_rate) {
      velocity[0] = groups.u1M + sqrt_temperature * myrandn(random);
      velocity[1] = groups.u2M + sqrt_temperature * myrandn(random);
      velocity[2] = groups.u3M + sqrt_temperature * myrandn(random);
      full[index].set_velocity(velocity);
    }
  }
}

void NegPar_BGK_collision(std::vector<NeParticleGroup>& groups,
                          NumericGridClass& grid, ParaClass& parameters,
                          RandomContext& random) {
  for (int cell = 0; cell < grid.Nx; ++cell)
    NegPar_BGK_collision_homo(groups[cell], parameters, random);
}

}  // namespace coulomb
