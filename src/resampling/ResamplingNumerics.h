#pragma once

#include <complex>
#include <cstddef>
#include <vector>
#include "RandomContext.h"

namespace coulomb::resampling {

// One independent uniform point per unit interval of cumulative expected
// proposals. Interior strata contribute deterministically; only boundary
// strata require a draw. Cells follow the resampler's lexicographic grid order.
class StratifiedProposalAllocator {
  public:
	int allocate(double expected, RandomContext& random);
  private:
	double cumulative{0.0};
	double boundaryUniform{0.0};
	bool initialized{false};
	bool boundaryIncluded{false};
};


class ResamplingNumerics {
  public:
	std::vector<double> frequencies(std::size_t count);
	std::vector<std::size_t> augmentedLocations(std::size_t count,
												std::size_t augmentationFactor);
	std::vector<std::complex<double>> imaginaryFrequencies(std::size_t count);
	double quadraticEnvelope(double halfWidth, const std::vector<double>& derivatives);
	std::vector<double> wrapPeriodicSample(const std::vector<double>& sample);
	double evaluateQuadraticTaylor(double deltaX, double deltaY, double deltaZ,
								   const std::vector<double>& derivatives);
};

} // namespace coulomb::resampling
