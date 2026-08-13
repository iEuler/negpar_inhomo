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

	static RunOptions parse(int argc, char** argv);
	static void reset_runtime_state(SimulationState& state);
	void apply(SimulationState& state) const;
};

} // namespace coulomb
