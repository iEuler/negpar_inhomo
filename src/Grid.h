#pragma once

#include <vector>

#include "SimulationTypes.h"

namespace coulomb {

class NumericGridClass {
 public:
  int Nx, Nt, Nv;
  double xmax, xmin, vmax, vmin, tmax, dx, dv, dt, Neff, Neff_F;
  std::vector<double> x, vx;
  BoundaryCondition bdry_x, bdry_v;
  double lambda_Poisson;

  NumericGridClass(int n_x, SimulationMethod method);
  NumericGridClass(int n_x) : NumericGridClass(n_x, SimulationMethod::HDP) {}
  NumericGridClass() : NumericGridClass(100) {}
};

}  // namespace coulomb
