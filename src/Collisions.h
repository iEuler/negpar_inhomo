#pragma once

#include <utility>
#include <vector>

#include "Classes.h"

namespace coulomb {

std::pair<std::vector<double>, std::vector<double>> coulombBinary3d(
    const std::vector<double>& velocity1,
    const std::vector<double>& velocity2, const ParaClass& parameters,
    RandomContext& random);

void coulomb_collision_homo(std::vector<Particle1d3d>& particles,
                            int particle_count,
                            const ParaClass& parameters,
                            RandomContext& random);

}  // namespace coulomb
