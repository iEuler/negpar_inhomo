#include "Initialization.h"

#include <iostream>

#include "Grid.h"
#include "InitialConditions.h"
#include "Moments.h"
#include "Output.h"
#include "ParticleGroup.h"
#include "ParticleInitialization.h"
#include "SimulationConfig.h"
#include "SimulationState.h"

namespace coulomb {

void initialize_distri_Negpar(NumericGridClass& grid,
                              std::vector<NeParticleGroup>& groups,
                              SimulationState& state) {
  auto initial_data = make_initial_conditions(grid);

  for (int cell = 0; cell < grid.Nx; ++cell) {
    configure_initial_condition(initial_data, grid, cell);
    const double center = grid.x[cell];
    const double dx = grid.x[1] - grid.x[0];
    groups[cell].set_xrange(center - dx / 2.0, center + dx / 2.0);
    groups[cell].reset_flag_resampled();
    initialize_Negpar(groups[cell], initial_data, grid.Neff, grid.Neff_F,
                      grid.dx, state.random);
  }

  save_initial(initial_data, state);
  update_macro(groups, grid);
}

}  // namespace coulomb
