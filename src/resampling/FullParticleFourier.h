#pragma once

// Deterministic Fourier numerics used by full-particle reconstruction.

#include <complex>
#include <vector>

namespace coulomb {

class NeParticleGroup;

class FullParticleFourier {
  public:
	std::vector<std::complex<double>>
	approximateTransform(NeParticleGroup& groups, int frequencyX,
						 int frequencyY, int frequencyZ);
	void filter(std::vector<std::complex<double>>& coefficients,
				std::vector<int>& flags, int size);
	std::vector<double>
	interpolateCoarse(const std::vector<std::complex<double>>& coefficients,
					  int frequencyX, int frequencyY, int frequencyZ);
	std::vector<std::vector<double>> interpolateDerivatives(
		const std::vector<std::complex<double>>& coefficients, int frequencyX,
		int frequencyY, int frequencyZ, int augmentationFactor);
	std::vector<double> upperBound(int count,
								   const std::vector<double>& values);
	std::vector<double> valuesAt(const std::vector<std::vector<double>>& values,
								 int index);
	void addMaxwellian(double density, std::vector<double> velocity,
					   std::vector<double> temperature,
					   double effectiveParticles,
					   std::vector<std::vector<double>>& derivatives,
					   int frequency, int augmentationFactor);
};
} // namespace coulomb
