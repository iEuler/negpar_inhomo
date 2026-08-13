#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace coulomb::resampling {

class ResamplingNumerics {
  public:
	std::vector<double> frequencies(std::size_t count);
	std::vector<std::size_t> augmentedLocations(std::size_t count,
												std::size_t augmentationFactor);
	std::vector<std::complex<double>> imaginaryFrequencies(std::size_t count);
	double evaluateQuadraticTaylor(double deltaX, double deltaY, double deltaZ,
								   const std::vector<double>& derivatives);
};

} // namespace coulomb::resampling
