#pragma once

#include <string>

namespace coulomb {

struct SimulationState;

std::string output_path(const std::string& filename,
                        const SimulationState& state);
std::string int2str(int value, int digits = 3);

}  // namespace coulomb
