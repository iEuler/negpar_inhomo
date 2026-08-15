#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

#include "SimulationConfig.h"

namespace coulomb {

struct ConfigRuntimeValues {
	std::optional<std::uint32_t> seed;
	std::optional<int> steps;
	std::optional<int> threads;
	std::optional<std::string> outputDirectory;
};

struct ConfigFileValues {
	ParaClass parameters;
	ConfigRuntimeValues runtime;
};

class ConfigFile {
  public:
	static ConfigFileValues load(const std::filesystem::path& path);
	static std::string serialize(const ConfigFileValues& values);
};

} // namespace coulomb
