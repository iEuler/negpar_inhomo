#pragma once

#include <vector>

namespace coulomb {

enum class BoundaryCondition;

class Numerics {
  public:
	double errorFunction(double x) const;
	std::vector<double> circularShift(const std::vector<double>& values,
									  int size, int shiftDistance) const;
	std::vector<double> boundaryShift(const std::vector<double>& values,
									  int size, int shiftDistance,
									  double boundary) const;
	std::vector<double> centralDifference(const std::vector<double>& values,
										  int size,
										  BoundaryCondition boundary) const;
	std::vector<double> limitedFlux(const std::vector<double>& values, int size,
									double dx, double dt, int velocitySign,
									const std::vector<double>& boundaryValues,
									BoundaryCondition boundary) const;
	void advanceKineticEuler(const std::vector<double>& density,
							 const std::vector<double>& velocity,
							 const std::vector<double>& temperature, int size,
							 double dx, double dt, BoundaryCondition boundary,
							 std::vector<double>& densityChange,
							 std::vector<double>& momentumChange,
							 std::vector<double>& energyChange) const;
};

} // namespace coulomb
