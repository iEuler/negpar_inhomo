#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

#include "ConfigFile.h"
#include "SimulationState.h"

namespace coulomb {

struct RunOptions {
	std::optional<std::uint32_t> seed;
	std::optional<int> steps;
	int threads{};
	std::string outputDirectory;
	ParaClass parameters;
	std::optional<std::filesystem::path> configPath;
	bool validateConfigOnly{};
	bool printEffectiveConfig{};

	static RunOptions parse(int argc, char** argv);
	static void resetRuntimeState(SimulationState& state);
	void apply(SimulationState& state) const;
	std::string effectiveConfigJson() const;
};

} // namespace coulomb
