#include <exception>
#include <iostream>

#include "RunOptions.h"
#include "Simulation.h"

int main(int argc, char** argv) {
	try {
		const auto options = coulomb::RunOptions::parse(argc, argv);
		coulomb::SimulationState state;
		options.apply(state);
		coulomb::Simulation simulation(options, state);
		return simulation.run();
	}
	catch (const std::exception& error) {
		std::cerr << error.what() << std::endl;
		return 2;
	}
}
