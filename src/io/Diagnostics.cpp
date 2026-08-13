#include "Diagnostics.h"

#include "Grid.h"
#include "ParticleGroup.h"

#include <stdexcept>

namespace coulomb {

double
Diagnostics::electricEnergy(const std::vector<NeParticleGroup>& sX) const {
	const auto& grid = gridRef;
	double energy = 0.0;
	for (int kx = 0; kx < grid.nx; ++kx) {
		const double field = sX[kx].elecField;
		energy += field * field;
	}
	return energy * grid.dx;
}

double
Diagnostics::fullElectricEnergy(const std::vector<NeParticleGroup>& sX) const {
	const auto& grid = gridRef;
	double energy = 0.0;
	for (int kx = 0; kx < grid.nx; ++kx) {
		const double field = sX[kx].elecFieldF;
		energy += field * field;
	}
	return energy * grid.dx;
}

double Diagnostics::totalEnergy(const std::vector<NeParticleGroup>& sX) const {
	const auto& grid = gridRef;
	double energy = electricEnergy(sX);
	for (int kx = 0; kx < grid.nx; ++kx) {
		const auto& group = sX[kx];
		energy += 0.5 * group.rhoM *
				  (group.u1M * group.u1M + 3.0 * group.tprtM) * grid.dx;
		energy += 0.5 * grid.neff *
				  (group.positiveMoments.m2 - group.negativeMoments.m2);
	}
	return energy;
}

double
Diagnostics::fullTotalEnergy(const std::vector<NeParticleGroup>& sX) const {
	const auto& grid = gridRef;
	double energy = fullElectricEnergy(sX);
	for (int kx = 0; kx < grid.nx; ++kx)
		energy += 0.5 * grid.neffF * sX[kx].fullMoments.m2;
	return energy;
}

int Diagnostics::particleCount(const std::vector<NeParticleGroup>& groups,
							   int size, ParticleKind kind) const {
	if (size < 0 || static_cast<std::size_t>(size) > groups.size())
		throw std::invalid_argument("particle count size exceeds group count");
	int count = 0;
	for (int cell = 0; cell < size; ++cell)
		count += groups[cell].size(kind);
	return count;
}

} // namespace coulomb
