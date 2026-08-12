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
  std::uint64_t owner_id = 0;
  std::uint32_t seed = 0;
  std::uint32_t thread_id = std::numeric_limits<std::uint32_t>::max();
  std::uint64_t generation = 0;
};

ThreadRandomState &random_state(RandomContext &context) {
  thread_local ThreadRandomState state;
  const auto thread_id = static_cast<std::uint32_t>(omp_get_thread_num());
  if (state.owner_id != context.instance_id() || state.seed != context.seed ||
      state.generation != context.generation || state.thread_id != thread_id) {
    std::seed_seq sequence{context.seed, thread_id};
    state.engine.seed(sequence);
    state.uniform.reset();
    state.normal.reset();
    state.owner_id = context.instance_id();
    state.seed = context.seed;
    state.thread_id = thread_id;
    state.generation = context.generation;
  }
  return state;
}

} // namespace

double RandomSampling::uniform(RandomContext &context) {
  double value = -1.0;
  auto &state = random_state(context);
  while (value <= 0.0 || value >= 1.0)
    value = state.uniform(state.engine);
  return value;
}

double RandomSampling::normal(RandomContext &context) {
  auto &state = random_state(context);
  return state.normal(state.engine);
}

std::vector<int> RandomSampling::permutation(int input_size, int output_size,
                                             RandomContext &context) {
  if (output_size > input_size) {
    throw std::runtime_error("Nout [" + std::to_string(output_size) +
                             "] must not be larger than Nin [" +
                             std::to_string(input_size) + "].");
  }

  std::vector<int> permutation(output_size);
  std::vector<int> candidates(input_size);
  for (int index = 0; index < input_size; ++index) {
    candidates[index] = index + 1;
  }

  for (int index = 0; index < output_size; ++index) {
    const int remaining = input_size - index;
    const int selected =
        static_cast<int>(remaining * RandomSampling::uniform(context));
    permutation[index] = candidates[selected];
    candidates[selected] = candidates[remaining - 1];
  }
  return permutation;
}

int RandomSampling::stochastic_floor(double value, RandomContext &context) {
  int result = static_cast<int>(value);
  const double fraction = value - result;
  if (RandomSampling::uniform(context) < fraction)
    ++result;
  return result;
}

} // namespace coulomb
