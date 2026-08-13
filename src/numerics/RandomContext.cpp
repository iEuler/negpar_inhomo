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

void RandomContext::reseed(std::uint32_t new_seed) {
	seed = new_seed;
	++generation;
}

std::uint32_t RandomContext::generate_seed() {
	std::random_device random_device;
	return random_device();
}

} // namespace coulomb
