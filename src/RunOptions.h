#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "SimulationState.h"

namespace coulomb {

struct RunOptions {
  std::optional<std::uint32_t> seed;
  std::optional<int> steps;
  int threads{};
  std::string output_directory;
};

RunOptions parse_run_options(int argc, char** argv);
void reset_runtime_state(SimulationState& state);
void apply_run_options(const RunOptions& options, SimulationState& state);

}  // namespace coulomb
