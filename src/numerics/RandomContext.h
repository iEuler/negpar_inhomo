#pragma once

#include <cstdint>

namespace coulomb {

// Owns the base seed from which each OpenMP thread derives its local engine.
// The engine itself remains thread-local in RandomSampling.cpp, while its seed
// source is part of the simulation's explicit state.
struct RandomContext {
	RandomContext();
	RandomContext(const RandomContext&) = delete;
	RandomContext& operator=(const RandomContext&) = delete;
	RandomContext(RandomContext&& other) noexcept;
	RandomContext& operator=(RandomContext&& other) noexcept;

	std::uint64_t instanceId() const noexcept { return instanceIdValue; }
	void reseed(std::uint32_t newSeed);
	static std::uint32_t generateSeed();

	std::uint32_t seed{};
	std::uint64_t generation{};

  private:
	std::uint64_t instanceIdValue{};
};

} // namespace coulomb
