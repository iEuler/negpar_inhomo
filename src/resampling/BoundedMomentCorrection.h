#pragma once
#include <array>
#include <cstddef>
#include "ParticleGroup.h"

namespace coulomb::resampling {
using SignedLowMoments = std::array<double,7>; // mass, three momenta, three diagonal second moments

enum class MomentCorrectionStatus {
    Success=0, IncompatibleMass, InsufficientParticles, Singular,
    SupportLimited, IterationLimit, DisplacementLimit
};
struct MomentCorrectionConfig {
    std::size_t maxIterations{40};
    double relativeTolerance{1e-10};
    double maxRmsDisplacementFraction{0.1}; // fraction of core radius
};
struct MomentCorrectionResult {
    MomentCorrectionStatus status{MomentCorrectionStatus::IterationLimit};
    std::size_t iterations{0};
    int removed{0};
    double rmsDisplacement{0.0}, maxDisplacement{0.0};
    double normalizedResidual{0.0};
    bool success() const { return status==MomentCorrectionStatus::Success; }
};
class BoundedMomentCorrection {
  public:
    static SignedLowMoments moments(const NeParticleGroup& source,double weight);
    // Transactional: failure leaves every particle and metadata unchanged.
    // Mass adjustment deletes uniformly selected particles of the excess sign.
    // Each velocity step minimizes its Euclidean norm subject to linearized
    // momentum/diagonal-second constraints on the currently free particles.
    // This is a local iterative correction, not a global optimum certificate.
    MomentCorrectionResult apply(NeParticleGroup& core,const SignedLowMoments& target,
        double weight,const std::array<double,3>& center,double radius,
        RandomContext& random,MomentCorrectionConfig config={}) const;
};
} // namespace coulomb::resampling
