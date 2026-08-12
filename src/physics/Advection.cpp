#include "Advection.h"

#include "ElectricField.h"
#include "Grid.h"
#include "Particle.h"
#include "ParticleGroup.h"

#include <iostream>
#include <stdexcept>

namespace coulomb {

Particle1d3d Advection::move_particle(const Particle1d3d& particle,
									  double electric_field) {
	double xnew = particle.position() + grid_.dt * particle.velocity(0);
	auto velocity = particle.velocity();
	double vxnew = velocity[0] + grid_.dt * electric_field;
	if (grid_.bdry_x == BoundaryCondition::Periodic) {
		const double period = grid_.xmax - grid_.xmin;
		while (xnew >= grid_.xmax)
			xnew -= period;
		while (xnew < grid_.xmin)
			xnew += period;
	} else if (grid_.bdry_x == BoundaryCondition::Reflective) {
		while (xnew >= grid_.xmax || xnew < grid_.xmin) {
			if (xnew >= grid_.xmax) {
				xnew = 2 * grid_.xmax - xnew;
				vxnew = -vxnew;
				if (xnew == grid_.xmax)
					break;
			}
			if (xnew < grid_.xmin) {
				xnew = 2 * grid_.xmin - xnew;
				vxnew = -vxnew;
			}
		}
	}
	velocity[0] = vxnew;
	++state_.movedCount;
	return Particle1d3d(xnew, velocity, true);
}

int Advection::find_particle_group(Particle1d3d& particle) const {
	if (particle.position() < grid_.xmin)
		throw std::out_of_range("Particle position is below the spatial grid");
	if (particle.position() == grid_.xmax)
		return grid_.Nx - 1;
	const int group =
		static_cast<int>((particle.position() - grid_.xmin) / grid_.dx);
	if (group >= grid_.Nx)
		throw std::out_of_range("Particle position is above the spatial grid");
	if (group < 0)
		throw std::out_of_range("Particle position is below the spatial grid");
	return group;
}

void Advection::relocate_particle(std::vector<ParticleGroup>& groups,
								  int group_before, int particle_index,
								  int group_after) {
	if (group_before != group_after) {
		groups[group_after].push_back(
			groups[group_before].list().at(particle_index));
		groups[group_before].erase(particle_index);
	}
}

void Advection::reset_moved_flags(std::vector<ParticleGroup>& groups) {
	for (int group = 0; group < grid_.Nx; ++group)
		for (auto& particle : groups[group].list())
			particle.flag_moved = false;
}

void Advection::advance(std::vector<ParticleGroup>& groups) {
	ElectricFieldSolver(grid_).update(groups);
	for (int group = 0; group < grid_.Nx; ++group) {
		auto& particles = groups[group].list();
		const double electric_field = groups[group].elecfield;
		int index = 0;
		while (index < groups[group].size()) {
			if (!particles[index].flag_moved) {
				particles[index] =
					move_particle(particles[index], electric_field);
				const int destination = find_particle_group(particles[index]);
				relocate_particle(groups, group, index, destination);
			} else
				++index;
		}
	}
	reset_moved_flags(groups);
}

void Advection::relocate_particle(std::vector<NeParticleGroup>& groups,
								  ParticleKind kind, int group_before,
								  int particle_index, int group_after) {
	if (group_before != group_after) {
		groups[group_after].push_back(
			groups[group_before].list(kind).at(particle_index), kind);
		groups[group_before].erase(particle_index, kind);
	}
}

void Advection::reset_moved_flags(std::vector<NeParticleGroup>& groups,
								  ParticleKind kind) {
	for (int group = 0; group < grid_.Nx; ++group)
		for (auto& particle : groups[group].list(kind)) {
			if (!particle.flag_moved)
				std::cout << "NOT MOVED\n";
			particle.flag_moved = false;
		}
}

void Advection::advance(std::vector<NeParticleGroup>& groups,
						ParticleKind kind) {
	for (int group = 0; group < grid_.Nx; ++group) {
		auto& particles = groups[group].list(kind);
		const double electric_field = kind == ParticleKind::Full
										  ? groups[group].elecfield_F
										  : groups[group].elecfield;
		int index = 0;
		while (index < groups[group].size(kind)) {
			if (!particles[index].flag_moved) {
				particles[index] =
					move_particle(particles[index], electric_field);
				const int destination = find_particle_group(particles[index]);
				relocate_particle(groups, kind, group, index, destination);
			} else
				++index;
		}
	}
	reset_moved_flags(groups, kind);
	state_.movedCount = 0;
}

void Advection::advance(std::vector<NeParticleGroup>& groups) {
	advance(groups, ParticleKind::Positive);
	advance(groups, ParticleKind::Negative);
	advance(groups, ParticleKind::Full);
}

} // namespace coulomb
