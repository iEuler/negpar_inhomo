#pragma once

#include "RunOptions.h"

namespace coulomb {

class Simulation {
 public:
  Simulation(RunOptions options, SimulationState& state);

  int run();

 private:
  RunOptions options_;
  SimulationState& state_;
};

}  // namespace coulomb
