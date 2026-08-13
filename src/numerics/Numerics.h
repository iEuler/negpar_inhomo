#pragma once

#include <vector>

namespace coulomb {

enum class BoundaryCondition;

class Numerics {
  public:
	double error_function(double x) const;
	std::vector<double> circular_shift(const std::vector<double>& values,
									   int size, int shift_distance) const;
	std::vector<double> boundary_shift(const std::vector<double>& values,
									   int size, int shift_distance,
									   double boundary) const;
	std::vector<double> central_difference(const std::vector<double>& values,
										   int size,
										   BoundaryCondition boundary) const;
	std::vector<double> limited_flux(const std::vector<double>& values,
									 int size, double dx, double dt,
									 int velocity_sign,
									 const std::vector<double>& boundary_values,
									 BoundaryCondition boundary) const;
	void advance_kinetic_euler(const std::vector<double>& density,
							   const std::vector<double>& velocity,
							   const std::vector<double>& temperature, int size,
							   double dx, double dt, BoundaryCondition boundary,
							   std::vector<double>& density_change,
							   std::vector<double>& momentum_change,
							   std::vector<double>& energy_change) const;
};

} // namespace coulomb
