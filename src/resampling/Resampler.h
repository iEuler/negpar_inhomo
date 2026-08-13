#pragma once
#include <complex>
#include <vector>

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

struct FourierResamplerConfig {
	double effectiveParticleWeight{1.0};
	size_t frequencyCount{30};
	bool useApproximation{true};
	size_t maxSamplingAttempts{1'000'000};
};

class FourierResampler {
  public:
	explicit FourierResampler(const NeParticleGroup& particles,
							  FourierResamplerConfig config = {});

	void reinit(const NeParticleGroup& particles) {
		particlesValue = particles;
	}

	NeParticleGroup resample(RandomContext& random) const;

  private:
	NeParticleGroup particlesValue;
	double neff;
	size_t nfreq;
	bool useApproximation;
	size_t augFactor = 2;
	size_t maxSamplingAttempts;

	VectorComplex3D fft3D(NeParticleGroup& sX) const;
	VectorComplex3D fft3DApprox(NeParticleGroup& sX) const;
	std::vector<Vector3D>
	derivativesFromFft(const VectorComplex3D& fourierCoeff) const;

	Vector3D funcOnAugGrid(const VectorComplex3D& fourierCoeff) const;
	Vector3D derivativesFromFftOneTerm(const VectorComplex3D& fourierCoeff,
									   int orderx, int ordery,
									   int orderz) const;

	VectorComplex3D fft3DApproxOneterm(const Vector3D& f, int orderx,
									   int ordery, int orderz) const;
};
} // namespace coulomb::resampling
