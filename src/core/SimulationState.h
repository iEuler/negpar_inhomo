#pragma once

#include "RandomContext.h"

#include <ctime>
#include <string>

namespace coulomb {

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

} // namespace coulomb
