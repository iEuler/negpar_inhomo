#pragma once
#include <complex>
#include <vector>

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

enum class ResamplingEnvelope { LegacyAdaptive, CertifiedQuadratic };
enum class ResamplingCellGeometry { LegacyShifted, PeriodicWrapped };
enum class ResamplingProposalAllocation { IndependentRounding, Stratified };
struct FourierResamplerDiagnostics {
	size_t attempts{0};
	size_t envelopeIncreases{0};
};

struct FourierResamplerConfig {
	double effectiveParticleWeight{1.0};
	double sourceSignedParticleWeight{0.0};
	double fullParticleWeight{1.0};
	size_t frequencyCount{30};
	bool useApproximation{true};
	bool weightedCoupling{false};
	ResamplingEnvelope envelope{ResamplingEnvelope::LegacyAdaptive};
	ResamplingCellGeometry cellGeometry{ResamplingCellGeometry::LegacyShifted};
	ResamplingProposalAllocation proposalAllocation{ResamplingProposalAllocation::IndependentRounding};
	// Optional physical bounds [xmin,xmax,ymin,ymax,zmin,zmax]. Empty uses extrema.
	std::vector<double> fixedVelocityBounds{};
	size_t maxSamplingAttempts{1'000'000};
};

class FourierResampler {
  public:
	explicit FourierResampler(const NeParticleGroup& particles,
							  FourierResamplerConfig config = {});

	void reinit(const NeParticleGroup& particles) {
		particlesValue = particles;
	}

	NeParticleGroup resample(RandomContext& random, FourierResamplerDiagnostics* diagnostics = nullptr) const;

  private:
	NeParticleGroup particlesValue;
	double neff;
	double outputNeff;
	double fullNeff;
	size_t nfreq;
	bool useApproximation;
	bool weightedCoupling;
	ResamplingEnvelope envelope;
	ResamplingCellGeometry cellGeometry;
	ResamplingProposalAllocation proposalAllocation;
	std::vector<double> fixedVelocityBounds;
	size_t augFactor = 2;
	size_t maxSamplingAttempts;

	VectorComplex3D fft3D(NeParticleGroup& sX) const;
	VectorComplex3D fft3DApprox(NeParticleGroup& sX) const;
	VectorComplex3D fft3DForKind(const NeParticleGroup& source,
								 ParticleKind kind,
								 double effectiveWeight) const;
	std::complex<double> maxwellianCoefficient(
		const NeParticleGroup& normalized, int kx, int ky, int kz) const;
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
