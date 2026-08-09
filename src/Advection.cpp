#include "Advection.h"

#include "Grid.h"
#include "Particle.h"
#include "ParticleGroup.h"
#include "ElectricField.h"


#include <iostream>
#include <stdexcept>

namespace coulomb {

Particle1d3d moveparticle(const Particle1d3d& particle, double elecfield,
                          const NumericGridClass& grid, SimulationState& state) {
  Particle1d3d moved;
  double xnew = particle.position() + grid.dt * particle.velocity(0);
  auto velocity = particle.velocity();
  double vxnew = velocity[0] + grid.dt * elecfield;

  if (grid.bdry_x == BoundaryCondition::Periodic) {
    const double xperiod = grid.xmax - grid.xmin;
    while (xnew >= grid.xmax) xnew -= xperiod;
    while (xnew < grid.xmin) xnew += xperiod;
  } else if (grid.bdry_x == BoundaryCondition::Reflective) {
    while ((xnew >= grid.xmax) || (xnew < grid.xmin)) {
      if (xnew >= grid.xmax) {
        xnew = 2 * grid.xmax - xnew;
        vxnew = -vxnew;
        // An exact hit on the upper wall is already a valid reflected
        // position; repeating the loop would otherwise reflect it forever.
        if (xnew == grid.xmax) break;
      }
      if (xnew < grid.xmin) {
        xnew = 2 * grid.xmin - xnew;
        vxnew = -vxnew;
      }
    }
  }

  velocity[0] = vxnew;
  state.movedCount++;
  return Particle1d3d(xnew, velocity, true);
}

int findparticlegroup(Particle1d3d& particle, const NumericGridClass& grid) {
  if (particle.position() < grid.xmin)
    throw std::out_of_range("Particle position is below the spatial grid");
  if (particle.position() == grid.xmax)
    return grid.Nx - 1;

  int group = static_cast<int>((particle.position() - grid.xmin) / grid.dx);
  if (group >= grid.Nx) {
    throw std::out_of_range("Particle position is above the spatial grid");
  }
  if (group < 0)
    throw std::out_of_range("Particle position is below the spatial grid");
  return group;
}

void relocateparticle(std::vector<ParticleGroup>& groups, int group_before,
                      int particle_index, int group_after) {
  if (group_before != group_after) {
    groups[group_after].push_back(
        groups[group_before].list().at(particle_index));
    groups[group_before].erase(particle_index);
  }
}

void reset_flag_moved(std::vector<ParticleGroup>& groups, int grid_size) {
  for (int group = 0; group < grid_size; ++group) {
    auto& particles = groups[group].list();
    for (int index = 0; index < groups[group].size(); ++index)
      particles[index].flag_moved = false;
  }
}

void particleadvection(std::vector<ParticleGroup>& groups,
                       const NumericGridClass& grid, SimulationState& state) {
  updateelecfiled(groups, grid);

  for (int group = 0; group < grid.Nx; ++group) {
    auto& particles = groups[group].list();
    const double elecfield = groups[group].elecfield;
    int index = 0;
    while (index < groups[group].size()) {
      if (!particles[index].flag_moved) {
        particles[index] = moveparticle(particles[index], elecfield, grid, state);
        const int destination = findparticlegroup(particles[index], grid);
        relocateparticle(groups, group, index, destination);
      } else {
        ++index;
      }
    }
  }
  reset_flag_moved(groups, grid.Nx);
}

void relocateparticle(std::vector<NeParticleGroup>& groups, ParticleKind kind,
                      int group_before, int particle_index, int group_after) {
  if (group_before != group_after) {
    groups[group_after].push_back(
        groups[group_before].list(kind).at(particle_index), kind);
    groups[group_before].erase(particle_index, kind);
  }
}

void reset_flag_moved(std::vector<NeParticleGroup>& groups, ParticleKind kind,
                      int grid_size) {
  for (int group = 0; group < grid_size; ++group) {
    auto& particles = groups[group].list(kind);
    for (int index = 0; index < groups[group].size(kind); ++index) {
      if (!particles[index].flag_moved) std::cout << "NOT MOVED\n";
      particles[index].flag_moved = false;
    }
  }
}

void particleadvection(std::vector<NeParticleGroup>& groups, ParticleKind kind,
                       const NumericGridClass& grid, SimulationState& state) {
  for (int group = 0; group < grid.Nx; ++group) {
    auto& particles = groups[group].list(kind);
    double elecfield = groups[group].elecfield;
    if (kind == ParticleKind::Full) elecfield = groups[group].elecfield_F;
    int index = 0;
    while (index < groups[group].size(kind)) {
      if (!particles[index].flag_moved) {
        particles[index] = moveparticle(particles[index], elecfield, grid, state);
        const int destination = findparticlegroup(particles[index], grid);
        relocateparticle(groups, kind, group, index, destination);
      } else {
        ++index;
      }
    }
  }
  reset_flag_moved(groups, kind, grid.Nx);
  state.movedCount = 0;
}

void particleadvection(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid, SimulationState& state) {
  particleadvection(groups, ParticleKind::Positive, grid, state);
  particleadvection(groups, ParticleKind::Negative, grid, state);
  particleadvection(groups, ParticleKind::Full, grid, state);
}

}  // namespace coulomb
