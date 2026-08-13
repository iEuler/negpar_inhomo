#include "RandomSampling.h"

#include <cstdint>
#include <limits>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>

namespace coulomb {
namespace {

struct ThreadRandomState {
	std::mt19937 engine;
	std::uniform_real_distribution<> uniform{0, 1};
	std::normal_distribution<> normal{0, 1};
	std::uint64_t ownerId = 0;
	std::uint32_t seed = 0;
	std::uint32_t threadId = std::numeric_limits<std::uint32_t>::max();
	std::uint64_t generation = 0;
};

ThreadRandomState& randomState(RandomContext& context) {
	thread_local ThreadRandomState state;
	const auto threadId = static_cast<std::uint32_t>(omp_get_thread_num());
	if (state.ownerId != context.instanceId() || state.seed != context.seed ||
		state.generation != context.generation || state.threadId != threadId) {
		std::seed_seq sequence{context.seed, threadId};
		state.engine.seed(sequence);
		state.uniform.reset();
		state.normal.reset();
		state.ownerId = context.instanceId();
		state.seed = context.seed;
		state.threadId = threadId;
		state.generation = context.generation;
	}
	return state;
}

} // namespace

double RandomSampling::uniform() {
	double value = -1.0;
	auto& state = randomState(context);
	while (value <= 0.0 || value >= 1.0)
		value = state.uniform(state.engine);
	return value;
}

double RandomSampling::normal() {
	auto& state = randomState(context);
	return state.normal(state.engine);
}

std::vector<int> RandomSampling::permutation(int inputSize, int outputSize) {
	if (outputSize > inputSize) {
		throw std::runtime_error("Nout [" + std::to_string(outputSize) +
								 "] must not be larger than Nin [" +
								 std::to_string(inputSize) + "].");
	}

	std::vector<int> permutation(outputSize);
	std::vector<int> candidates(inputSize);
	for (int index = 0; index < inputSize; ++index) {
		candidates[index] = index + 1;
	}

	for (int index = 0; index < outputSize; ++index) {
		const int remaining = inputSize - index;
		const int selected = static_cast<int>(remaining * uniform());
		permutation[index] = candidates[selected];
		candidates[selected] = candidates[remaining - 1];
	}
	return permutation;
}

int RandomSampling::stochasticFloor(double value) {
	int result = static_cast<int>(value);
	const double fraction = value - result;
	if (uniform() < fraction)
		++result;
	return result;
}

} // namespace coulomb
