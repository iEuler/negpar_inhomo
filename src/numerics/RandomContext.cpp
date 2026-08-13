#include "RandomContext.h"

#include <atomic>
#include <random>

namespace coulomb {

RandomContext::RandomContext() {
	static std::atomic<std::uint64_t> nextId{1};
	instanceIdValue = nextId.fetch_add(1, std::memory_order_relaxed);
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

void RandomContext::reseed(std::uint32_t newSeed) {
	seed = newSeed;
	++generation;
}

std::uint32_t RandomContext::generateSeed() {
	std::random_device randomDevice;
	return randomDevice();
}

} // namespace coulomb
