#include "Advection.h"

#include "ElectricField.h"
#include "Grid.h"
#include "Particle.h"
#include "ParticleGroup.h"

#include <iostream>
#include <stdexcept>

namespace coulomb {

Particle1D3D Advection::moveParticle(const Particle1D3D& particle,
									 double electricField) {
	double xnew = particle.position() + gridRef.dt * particle.velocity(0);
	auto velocity = particle.velocity();
	double vxnew = velocity[0] + gridRef.dt * electricField;
	if (gridRef.bdryX == BoundaryCondition::Periodic) {
		const double period = gridRef.xmax - gridRef.xmin;
		while (xnew >= gridRef.xmax)
			xnew -= period;
		while (xnew < gridRef.xmin)
			xnew += period;
	} else if (gridRef.bdryX == BoundaryCondition::Reflective) {
		while (xnew >= gridRef.xmax || xnew < gridRef.xmin) {
			if (xnew >= gridRef.xmax) {
				xnew = 2 * gridRef.xmax - xnew;
				vxnew = -vxnew;
				if (xnew == gridRef.xmax)
					break;
			}
			if (xnew < gridRef.xmin) {
				xnew = 2 * gridRef.xmin - xnew;
				vxnew = -vxnew;
			}
		}
	}
	velocity[0] = vxnew;
	++stateRef.movedCount;
	return Particle1D3D(xnew, velocity, true);
}

int Advection::findParticleGroup(Particle1D3D& particle) const {
	if (particle.position() < gridRef.xmin)
		throw std::out_of_range("Particle position is below the spatial grid");
	if (particle.position() == gridRef.xmax)
		return gridRef.nx - 1;
	const int group =
		static_cast<int>((particle.position() - gridRef.xmin) / gridRef.dx);
	if (group >= gridRef.nx)
		throw std::out_of_range("Particle position is above the spatial grid");
	if (group < 0)
		throw std::out_of_range("Particle position is below the spatial grid");
	return group;
}

void Advection::relocateParticle(std::vector<ParticleGroup>& groups,
								 int groupBefore, int particleIndex,
								 int groupAfter) {
	if (groupBefore != groupAfter) {
		groups[groupAfter].pushBack(
			groups[groupBefore].list().at(particleIndex));
		groups[groupBefore].erase(particleIndex);
	}
}

void Advection::resetMovedFlags(std::vector<ParticleGroup>& groups) {
	for (int group = 0; group < gridRef.nx; ++group)
		for (auto& particle : groups[group].list())
			particle.flagMoved = false;
}

void Advection::advance(std::vector<ParticleGroup>& groups) {
	ElectricFieldSolver(gridRef).update(groups);
	for (int group = 0; group < gridRef.nx; ++group) {
		auto& particles = groups[group].list();
		const double electricField = groups[group].elecField;
		int index = 0;
		while (index < groups[group].size()) {
			if (!particles[index].flagMoved) {
				particles[index] =
					moveParticle(particles[index], electricField);
				const int destination = findParticleGroup(particles[index]);
				relocateParticle(groups, group, index, destination);
			} else
				++index;
		}
	}
	resetMovedFlags(groups);
}

void Advection::relocateParticle(std::vector<NeParticleGroup>& groups,
								 ParticleKind kind, int groupBefore,
								 int particleIndex, int groupAfter) {
	if (groupBefore != groupAfter) {
		groups[groupAfter].pushBack(
			groups[groupBefore].list(kind).at(particleIndex), kind);
		groups[groupBefore].erase(particleIndex, kind);
	}
}

void Advection::resetMovedFlags(std::vector<NeParticleGroup>& groups,
								ParticleKind kind) {
	for (int group = 0; group < gridRef.nx; ++group)
		for (auto& particle : groups[group].list(kind)) {
			if (!particle.flagMoved)
				std::cout << "NOT MOVED\n";
			particle.flagMoved = false;
		}
}

void Advection::advance(std::vector<NeParticleGroup>& groups,
						ParticleKind kind) {
	for (int group = 0; group < gridRef.nx; ++group) {
		auto& particles = groups[group].list(kind);
		const double electricField = kind == ParticleKind::Full
										 ? groups[group].elecFieldF
										 : groups[group].elecField;
		int index = 0;
		while (index < groups[group].size(kind)) {
			if (!particles[index].flagMoved) {
				particles[index] =
					moveParticle(particles[index], electricField);
				const int destination = findParticleGroup(particles[index]);
				relocateParticle(groups, kind, group, index, destination);
			} else
				++index;
		}
	}
	resetMovedFlags(groups, kind);
	stateRef.movedCount = 0;
}

void Advection::advance(std::vector<NeParticleGroup>& groups) {
	advance(groups, ParticleKind::Positive);
	advance(groups, ParticleKind::Negative);
	advance(groups, ParticleKind::Full);
}

} // namespace coulomb
