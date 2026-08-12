#include "Advection.h"

#include "ElectricField.h"
#include "Grid.h"
#include "Particle.h"
#include "ParticleGroup.h"

#include <iostream>
#include <stdexcept>

namespace coulomb {

Particle1d3d Advection::move_particle(const Particle1d3d &particle,
                                      double elecfield,
                                      const NumericGridClass &grid,
                                      SimulationState &state) {
  Particle1d3d moved;
  double xnew = particle.position() + grid.dt * particle.velocity(0);
  auto velocity = particle.velocity();
  double vxnew = velocity[0] + grid.dt * elecfield;

  if (grid.bdry_x == BoundaryCondition::Periodic) {
    const double xperiod = grid.xmax - grid.xmin;
    while (xnew >= grid.xmax)
      xnew -= xperiod;
    while (xnew < grid.xmin)
      xnew += xperiod;
  } else if (grid.bdry_x == BoundaryCondition::Reflective) {
    while ((xnew >= grid.xmax) || (xnew < grid.xmin)) {
      if (xnew >= grid.xmax) {
        xnew = 2 * grid.xmax - xnew;
        vxnew = -vxnew;
        // An exact hit on the upper wall is already a valid reflected
        // position; repeating the loop would otherwise reflect it forever.
        if (xnew == grid.xmax)
          break;
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

int Advection::find_particle_group(Particle1d3d &particle,
                                   const NumericGridClass &grid) {
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

void Advection::relocate_particle(std::vector<ParticleGroup> &groups,
                                  int group_before, int particle_index,
                                  int group_after) {
  if (group_before != group_after) {
    groups[group_after].push_back(
        groups[group_before].list().at(particle_index));
    groups[group_before].erase(particle_index);
  }
}

void Advection::reset_moved_flags(std::vector<ParticleGroup> &groups,
                                  int grid_size) {
  for (int group = 0; group < grid_size; ++group) {
    auto &particles = groups[group].list();
    for (int index = 0; index < groups[group].size(); ++index)
      particles[index].flag_moved = false;
  }
}

void Advection::advance(std::vector<ParticleGroup> &groups,
                        const NumericGridClass &grid, SimulationState &state) {
  ElectricFieldSolver::update(groups, grid);

  for (int group = 0; group < grid.Nx; ++group) {
    auto &particles = groups[group].list();
    const double elecfield = groups[group].elecfield;
    int index = 0;
    while (index < groups[group].size()) {
      if (!particles[index].flag_moved) {
        particles[index] =
            Advection::move_particle(particles[index], elecfield, grid, state);
        const int destination =
            Advection::find_particle_group(particles[index], grid);
        Advection::relocate_particle(groups, group, index, destination);
      } else {
        ++index;
      }
    }
  }
  Advection::reset_moved_flags(groups, grid.Nx);
}

void Advection::relocate_particle(std::vector<NeParticleGroup> &groups,
                                  ParticleKind kind, int group_before,
                                  int particle_index, int group_after) {
  if (group_before != group_after) {
    groups[group_after].push_back(
        groups[group_before].list(kind).at(particle_index), kind);
    groups[group_before].erase(particle_index, kind);
  }
}

void Advection::reset_moved_flags(std::vector<NeParticleGroup> &groups,
                                  ParticleKind kind, int grid_size) {
  for (int group = 0; group < grid_size; ++group) {
    auto &particles = groups[group].list(kind);
    for (int index = 0; index < groups[group].size(kind); ++index) {
      if (!particles[index].flag_moved)
        std::cout << "NOT MOVED\n";
      particles[index].flag_moved = false;
    }
  }
}

void Advection::advance(std::vector<NeParticleGroup> &groups, ParticleKind kind,
                        const NumericGridClass &grid, SimulationState &state) {
  for (int group = 0; group < grid.Nx; ++group) {
    auto &particles = groups[group].list(kind);
    double elecfield = groups[group].elecfield;
    if (kind == ParticleKind::Full)
      elecfield = groups[group].elecfield_F;
    int index = 0;
    while (index < groups[group].size(kind)) {
      if (!particles[index].flag_moved) {
        particles[index] =
            Advection::move_particle(particles[index], elecfield, grid, state);
        const int destination =
            Advection::find_particle_group(particles[index], grid);
        Advection::relocate_particle(groups, kind, group, index, destination);
      } else {
        ++index;
      }
    }
  }
  Advection::reset_moved_flags(groups, kind, grid.Nx);
  state.movedCount = 0;
}

void Advection::advance(std::vector<NeParticleGroup> &groups,
                        const NumericGridClass &grid, SimulationState &state) {
  Advection::advance(groups, ParticleKind::Positive, grid, state);
  Advection::advance(groups, ParticleKind::Negative, grid, state);
  Advection::advance(groups, ParticleKind::Full, grid, state);
}

} // namespace coulomb
