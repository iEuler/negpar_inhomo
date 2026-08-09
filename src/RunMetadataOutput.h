#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;
class ParaClass;
struct SimulationState;

void save_grids(const NumericGridClass& grid, const SimulationState& state);
void saveparameter(const ParaClass& parameters, const NumericGridClass& grid,
                   const SimulationState& state);
void save_initial(IniValClass& initial_data, const SimulationState& state);

}  // namespace coulomb
