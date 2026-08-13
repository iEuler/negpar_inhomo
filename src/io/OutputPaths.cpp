#include "OutputPaths.h"

#include <filesystem>
#include <stdexcept>

#include "SimulationState.h"

namespace coulomb {
std::string OutputPaths::resolve(const std::string& filename) const {
	const std::filesystem::path relative(filename);
	if (relative.empty() || relative.is_absolute())
		throw std::invalid_argument("output filename must be relative");
	for (const auto& component : relative)
		if (component == "..")
			throw std::invalid_argument(
				"output filename escapes output directory");
	return (std::filesystem::path(stateRef.outputDirectory) / relative)
		.string();
}

std::string OutputPaths::formatIndex(int value, int digits) const {
	unsigned int unsignedValue = value;
	if (value < 0)
		unsignedValue = static_cast<unsigned int>(-value);

	std::string result;
	while (digits-- > 0) {
		result += static_cast<char>('0' + unsignedValue % 10);
		unsignedValue /= 10;
	}
	if (value < 0)
		result += '-';
	result = std::string(result.rbegin(), result.rend());
	return '_' + result;
}

} // namespace coulomb
