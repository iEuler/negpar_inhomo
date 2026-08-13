#include "Initialization.h"

#include <iostream>

#include "Grid.h"
#include "InitialConditions.h"
#include "Moments.h"
#include "ParticleGroup.h"
#include "ParticleInitialization.h"
#include "RunMetadataOutput.h"
#include "SimulationConfig.h"
#include "SimulationState.h"

namespace coulomb {

void Initialization::initialize(NumericGridClass& grid,
								std::vector<NeParticleGroup>& groups,
								SimulationState& state) {
	auto initialData = InitialConditions{}.create(grid);

	for (int cell = 0; cell < grid.nx; ++cell) {
		InitialConditions{}.configure(initialData, grid, cell);
		const double center = grid.x[cell];
		const double dx = grid.x[1] - grid.x[0];
		groups[cell].setXRange(center - dx / 2.0, center + dx / 2.0);
		groups[cell].resetFlagResampled();
		ParticleInitialization{}.initialize(groups[cell], initialData,
											grid.neff, grid.neffF, grid.dx,
											state.random);
	}

	RunMetadataOutput{}.saveInitialConditions(initialData, state);
	MomentOperations{}.updateMacro(groups, grid);
}

} // namespace coulomb
