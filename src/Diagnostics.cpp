#include "Diagnostics.h"

#include "Grid.h"
#include "ParticleGroup.h"

#include <stdexcept>

namespace coulomb {

double Diagnostics::electric_energy(const std::vector<NeParticleGroup> &S_x,
                                    const NumericGridClass &grid) {
  double energy = 0.0;
  for (int kx = 0; kx < grid.Nx; ++kx) {
    const double field = S_x[kx].elecfield;
    energy += field * field;
  }
  return energy * grid.dx;
}

double
Diagnostics::full_electric_energy(const std::vector<NeParticleGroup> &S_x,
                                  const NumericGridClass &grid) {
  double energy = 0.0;
  for (int kx = 0; kx < grid.Nx; ++kx) {
    const double field = S_x[kx].elecfield_F;
    energy += field * field;
  }
  return energy * grid.dx;
}

double Diagnostics::total_energy(const std::vector<NeParticleGroup> &S_x,
                                 const NumericGridClass &grid) {
  double energy = Diagnostics::electric_energy(S_x, grid);
  for (int kx = 0; kx < grid.Nx; ++kx) {
    const auto &group = S_x[kx];
    energy += 0.5 * group.rhoM * (group.u1M * group.u1M + 3.0 * group.TprtM) *
              grid.dx;
    energy += 0.5 * grid.Neff *
              (group.positive_moments.m2 - group.negative_moments.m2);
  }
  return energy;
}

double Diagnostics::full_total_energy(const std::vector<NeParticleGroup> &S_x,
                                      const NumericGridClass &grid) {
  double energy = Diagnostics::full_electric_energy(S_x, grid);
  for (int kx = 0; kx < grid.Nx; ++kx)
    energy += 0.5 * grid.Neff_F * S_x[kx].full_moments.m2;
  return energy;
}

int Diagnostics::particle_count(const std::vector<NeParticleGroup> &groups,
                                int size, ParticleKind kind) {
  if (size < 0 || static_cast<std::size_t>(size) > groups.size())
    throw std::invalid_argument("particle count size exceeds group count");
  int count = 0;
  for (int cell = 0; cell < size; ++cell)
    count += groups[cell].size(kind);
  return count;
}

} // namespace coulomb
