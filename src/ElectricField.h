#pragma once

#include <vector>

#include "Classes.h"

namespace coulomb {

std::vector<double> PoissonSolver(const std::vector<double>& rho, int grid_size,
                                  double domain_size, double lambda);

void updateelecfiled(std::vector<ParticleGroup>& groups,
                     const NumericGridClass& grid);
void updateelecfiled(std::vector<NeParticleGroup>& groups,
                     const NumericGridClass& grid);
void updateelecfiled_PIC(std::vector<NeParticleGroup>& groups,
                         const NumericGridClass& grid);
void updateelecfiled_fromcoarse(std::vector<NeParticleGroup>& groups,
                                const NumericGridClass& grid);
void updateelecfiled_zero(std::vector<NeParticleGroup>& groups,
                          const NumericGridClass& grid);
void updateelecfiled_rho(std::vector<NeParticleGroup>& groups,
                         const NumericGridClass& grid);

}  // namespace coulomb
