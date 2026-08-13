#include "Grid.h"

#include <stdexcept>

#include "Constants.h"

namespace coulomb {

NumericGridClass::NumericGridClass(int nX, SimulationMethod method) {
	if (nX <= 0)
		throw std::invalid_argument("NumericGridClass requires nX > 0");
	if (method != SimulationMethod::HDP && method != SimulationMethod::PIC)
		throw std::invalid_argument("NumericGridClass method is invalid");

	lambdaPoisson = 10.0;
	xmax = 2 * pi / 0.5;
	xmin = 0;
	vmax = 6.0;
	vmin = -vmax;
	tmax = 10;
	nx = nX;
	dx = (xmax - xmin) / nx;
	dt = dx / 2 / vmax;
	if (nx == 1)
		dt = 0.001;
	nt = static_cast<int>(tmax / dt);

	if (method == SimulationMethod::HDP) {
		neff = 1e-5;
		neffF = 1e-5;
	} else {
		neff = 1e-4;
		neffF = 5e-7;
	}

	x.resize(nx);
	for (int kx = 0; kx < nx; ++kx)
		x[kx] = xmin + (kx + 0.5) * dx;
	nv = 200;
	dv = (vmax - vmin) / nv;
	vx.resize(nv);
	for (int kv = 0; kv < nv; ++kv)
		vx[kv] = vmin + (kv + 0.5) * dv;
	bdryX = BoundaryCondition::Periodic;
	bdryV = BoundaryCondition::Periodic;
}

} // namespace coulomb
