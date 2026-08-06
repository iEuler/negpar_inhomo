#pragma once
#include <ctime>
#include <cstdint>
#include <string>

namespace coulomb {

inline constexpr double pi = 3.1415926535897932;

// Owns the base seed from which each OpenMP thread derives its local engine.
// The engine itself remains thread-local in utils.cpp, but its seed source is
// now part of the simulation's explicit state.
struct RandomContext {
  std::uint32_t seed{};
  std::uint64_t generation{};
};

struct SimulationState {
  int saveIndex{};
  bool filenameWithNumber{};
  bool saveFlux{};
  int movedCount{};
  int resampleCount{};
  double syncTime{};
  std::string outputDirectory{"result"};
  RandomContext random{};

  std::clock_t t0All{};
  std::clock_t t1All{};
  std::clock_t t0Collision{};
  std::clock_t t1Collision{};
  std::clock_t t0Advection{};
  std::clock_t t1Advection{};
  std::clock_t t0Resampling{};
  std::clock_t t1Resampling{};
};

}  // namespace coulomb
