#pragma once

#include <vector>

namespace coulomb {

double myerf(double x);

std::vector<double> cshift_1d(const std::vector<double>& values, int size,
                              int shift_distance);
std::vector<double> eoshift_1d(const std::vector<double>& values, int size,
                               int shift_distance, double boundary);
std::vector<double> diff_1d_central(const std::vector<double>& values, int size,
                                    char boundary);

// First-order kinetic fluxes used by the macroscopic Euler update.
std::vector<double> limiter_x1_o2(const std::vector<double>& values, int size,
                                  double dx, double dt, int velocity_sign,
                                  const std::vector<double>& boundary_values,
                                  char boundary);
void Euler_kinetic_x1(const std::vector<double>& density,
                      const std::vector<double>& velocity,
                      const std::vector<double>& temperature, int size,
                      double dx, double dt, char boundary,
                      std::vector<double>& density_change,
                      std::vector<double>& momentum_change,
                      std::vector<double>& energy_change);

}  // namespace coulomb
