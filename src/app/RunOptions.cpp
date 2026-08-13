#include "RunOptions.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <omp.h>
#include <stdexcept>

#include "RandomContext.h"

namespace coulomb {

namespace {

const char* usage() {
	return "Usage: negpar_inhomo [--seed <uint32>] [--threads <count>] "
		   "[--steps <count>] [--output-dir <path>]";
}

std::uint64_t parseInteger(const std::string& option, const char* value) {
	try {
		std::size_t parsedLength = 0;
		const std::string text(value);
		const auto parsed = std::stoull(text, &parsedLength, 10);
		if (parsedLength != text.size())
			throw std::invalid_argument("contains trailing characters");
		return parsed;
	}
	catch (const std::exception& error) {
		throw std::invalid_argument("Invalid value for " + option + ": " +
									error.what());
	}
}

} // namespace

void RunOptions::resetRuntimeState(SimulationState& state) {
	state = SimulationState{};
}

RunOptions RunOptions::parse(int argc, char** argv) {
	if (argc > 1 && ((argc - 1) % 2 != 0))
		throw std::invalid_argument(usage());

	RunOptions options;
	options.threads = std::max(1, omp_get_max_threads() - 2);
	options.outputDirectory = "result";

	for (int optionIndex = 1; optionIndex < argc; optionIndex += 2) {
		const std::string option(argv[optionIndex]);
		const char* value = argv[optionIndex + 1];

		if (option == "--output-dir") {
			options.outputDirectory = value;
			if (options.outputDirectory.empty())
				throw std::invalid_argument(
					"Output directory must not be empty");
			continue;
		}

		const auto parsed = parseInteger(option, value);
		if (option == "--seed") {
			if (parsed > std::numeric_limits<std::uint32_t>::max())
				throw std::invalid_argument("Seed is too large");
			options.seed = static_cast<std::uint32_t>(parsed);
		} else if (option == "--threads") {
			if (parsed == 0 || parsed > std::numeric_limits<int>::max())
				throw std::invalid_argument("Thread count is out of range");
			options.threads = static_cast<int>(parsed);
		} else if (option == "--steps") {
			if (parsed == 0 || parsed > std::numeric_limits<int>::max())
				throw std::invalid_argument("Step count is out of range");
			options.steps = static_cast<int>(parsed);
		} else {
			throw std::invalid_argument("Unknown option: " + option + "\n" +
										usage());
		}
	}

	return options;
}

void RunOptions::apply(SimulationState& state) const {
	RunOptions::resetRuntimeState(state);
	const auto selectedSeed = seed.value_or(RandomContext::generateSeed());
	state.random.reseed(selectedSeed);
	omp_set_num_threads(threads);

	state.outputDirectory = outputDirectory;
	std::error_code outputError;
	std::filesystem::create_directories(state.outputDirectory, outputError);
	if (outputError)
		throw std::runtime_error("Unable to create output directory '" +
								 state.outputDirectory +
								 "': " + outputError.message());

	std::ofstream metadata(std::filesystem::path(state.outputDirectory) /
						   "run_metadata.txt");
	if (!metadata)
		throw std::runtime_error("Unable to write run metadata in '" +
								 state.outputDirectory + "'");
	metadata << "seed " << state.random.seed << '\n'
			 << "threads " << threads << '\n'
			 << "steps " << (steps ? std::to_string(*steps) : "default") << '\n'
			 << "output_directory " << state.outputDirectory << '\n'
			 << "rng_engine std::mt19937\n"
			 << "uniform_distribution std::uniform_real_distribution<double>\n"
			 << "normal_distribution std::normal_distribution<double>\n"
			 << "thread_seed_derivation "
				"std_seed_seq_base_seed_and_openmp_thread_id\n"
			 << "cross_platform_bitwise_identity false\n"
#if defined(_MSC_VER)
			 << "compiler msvc " << _MSC_VER << '\n'
#if defined(_MSVC_STL_VERSION)
			 << "standard_library msvc_stl " << _MSVC_STL_VERSION << '\n'
#else
			 << "standard_library msvc_stl unknown\n"
#endif
#elif defined(__clang__)
			 << "compiler clang " << __clang_major__ << '.' << __clang_minor__
			 << '.' << __clang_patchlevel__ << '\n'
			 << "standard_library implementation-dependent\n"
#elif defined(__GNUC__)
           << "compiler gcc " << __GNUC__ << '.' << __GNUC_MINOR__ << '.'
           << __GNUC_PATCHLEVEL__ << '\n'
           << "standard_library libstdc++ " << __GLIBCXX__ << '\n'
#else
           << "compiler unknown\n"
           << "standard_library unknown\n"
#endif
#ifdef NDEBUG
			 << "build_type release\n";
#else
			 << "build_type debug\n";
#endif

	std::cout << "seed = " << state.random.seed << std::endl;
	std::cout << "threads = " << threads << std::endl;
	std::cout << "output_dir = " << state.outputDirectory << std::endl;
}

} // namespace coulomb
