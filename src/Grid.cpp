#include "Grid.h"

#include <stdexcept>

#include "Constants.h"

namespace coulomb {

NumericGridClass::NumericGridClass(int n_x, SimulationMethod method) {
	if (n_x <= 0)
		throw std::invalid_argument("NumericGridClass requires n_x > 0");
	if (method != SimulationMethod::HDP && method != SimulationMethod::PIC)
		throw std::invalid_argument("NumericGridClass method is invalid");

	lambda_Poisson = 10.0;
	xmax = 2 * pi / 0.5;
	xmin = 0;
	vmax = 6.0;
	vmin = -vmax;
	tmax = 10;
	Nx = n_x;
	dx = (xmax - xmin) / Nx;
	dt = dx / 2 / vmax;
	if (Nx == 1)
		dt = 0.001;
	Nt = static_cast<int>(tmax / dt);

	if (method == SimulationMethod::HDP) {
		Neff = 1e-5;
		Neff_F = 1e-5;
	} else {
		Neff = 1e-4;
		Neff_F = 5e-7;
	}

	x.resize(Nx);
	for (int kx = 0; kx < Nx; ++kx)
		x[kx] = xmin + (kx + 0.5) * dx;
	Nv = 200;
	dv = (vmax - vmin) / Nv;
	vx.resize(Nv);
	for (int kv = 0; kv < Nv; ++kv)
		vx[kv] = vmin + (kv + 0.5) * dv;
	bdry_x = BoundaryCondition::Periodic;
	bdry_v = BoundaryCondition::Periodic;
}

} // namespace coulomb
