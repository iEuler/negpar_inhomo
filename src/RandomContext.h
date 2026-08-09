#pragma once

#include <cstdint>

namespace coulomb {

// Owns the base seed from which each OpenMP thread derives its local engine.
// The engine itself remains thread-local in utils.cpp, but its seed source is
// part of the simulation's explicit state.
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

}  // namespace coulomb
