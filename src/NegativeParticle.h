#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct SimulationState;

void Negpar_inhomo_onestep(std::vector<NeParticleGroup>& groups,
                           NumericGridClass& grid, ParaClass& parameters,
                           SimulationState& state);
void Negpar_inhomo_onestep_PIC(std::vector<NeParticleGroup>& groups,
                               NumericGridClass& grid, ParaClass& parameters,
                               SimulationState& state);

}  // namespace coulomb
