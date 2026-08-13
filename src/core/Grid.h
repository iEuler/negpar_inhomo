#pragma once

#include <vector>

#include "SimulationTypes.h"

namespace coulomb {

class NumericGridClass {
  public:
	int nx, nt, nv;
	double xmax, xmin, vmax, vmin, tmax, dx, dv, dt, neff, neffF;
	std::vector<double> x, vx;
	BoundaryCondition bdryX, bdryV;
	double lambdaPoisson;

	NumericGridClass(int nX, SimulationMethod method);
	NumericGridClass(int nX) : NumericGridClass(nX, SimulationMethod::HDP) {}
	NumericGridClass() : NumericGridClass(100) {}
};

} // namespace coulomb
