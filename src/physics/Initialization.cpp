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
	auto initial_data = InitialConditions{}.create(grid);

	for (int cell = 0; cell < grid.Nx; ++cell) {
		InitialConditions{}.configure(initial_data, grid, cell);
		const double center = grid.x[cell];
		const double dx = grid.x[1] - grid.x[0];
		groups[cell].set_xrange(center - dx / 2.0, center + dx / 2.0);
		groups[cell].reset_flag_resampled();
		ParticleInitialization{}.initialize(groups[cell], initial_data,
											grid.Neff, grid.Neff_F, grid.dx,
											state.random);
	}

	RunMetadataOutput{}.save_initial_conditions(initial_data, state);
	MomentOperations{}.update_macro(groups, grid);
}

} // namespace coulomb
