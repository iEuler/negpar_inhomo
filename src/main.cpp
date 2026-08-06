#include <exception>
#include <iostream>

#include "RunOptions.h"
#include "Simulation.h"

int main(int argc, char** argv) {
  try {
    const auto options = coulomb::parse_run_options(argc, argv);
    coulomb::SimulationState state;
    return coulomb::run_simulation(options, state);
  } catch (const std::exception& error) {
    std::cerr << error.what() << std::endl;
    return 2;
  }
}
