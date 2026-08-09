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

  std::uint64_t instance_id() const noexcept { return instance_id_; }

  std::uint32_t seed{};
  std::uint64_t generation{};

 private:
  std::uint64_t instance_id_{};
};

void reseed_random(RandomContext& context, std::uint32_t seed);
std::uint32_t generate_random_seed();

}  // namespace coulomb
