#include "RandomContext.h"

#include <atomic>
#include <random>

namespace coulomb {

RandomContext::RandomContext() {
  static std::atomic<std::uint64_t> next_id{1};
  instance_id_ = next_id.fetch_add(1, std::memory_order_relaxed);
}

RandomContext::RandomContext(RandomContext&& other) noexcept : RandomContext() {
  seed = other.seed;
  generation = other.generation;
}

RandomContext& RandomContext::operator=(RandomContext&& other) noexcept {
  if (this != &other) {
    seed = other.seed;
    ++generation;
  }
  return *this;
}

void reseed_random(RandomContext& context, std::uint32_t seed) {
  context.seed = seed;
  ++context.generation;
}

std::uint32_t generate_random_seed() {
  std::random_device random_device;
  return random_device();
}

}  // namespace coulomb
