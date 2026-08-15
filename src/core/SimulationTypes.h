#pragma once

#include <stdexcept>

namespace coulomb {

enum class CollisionType { NoCollision, Coulomb, BGK };
enum class BinaryCollisionMethod { TA };
enum class BoundaryCondition { Periodic, Reflective };
enum class SimulationMethod { HDP, PIC };
enum class HdpCouplingMode { Decoupled, VarianceWeighted };
enum class EffectiveWeightPolicy { Fixed, QuadraticAdaptive };
enum class CollisionCoupling { Standard, Linearized };
enum class ProjectionMode { FullMicroMacro, MaxwellianOnly };
enum class DeltaMMode { Enabled, Disabled };

class SimulationTypes {
  public:
	BoundaryCondition decodeBoundary(char code) const;
	char encodeBoundary(BoundaryCondition condition) const;

	const char* methodName(SimulationMethod method) const {
		switch (method) {
		case SimulationMethod::HDP:
			return "HDP";
		case SimulationMethod::PIC:
			return "PIC";
		}
		throw std::invalid_argument("unknown simulation method");
	}

	const char* collisionName(CollisionType type) const {
		switch (type) {
		case CollisionType::NoCollision:
			return "NO_COLLISION";
		case CollisionType::Coulomb:
			return "COULOMB_COLLISION";
		case CollisionType::BGK:
			return "BGK_COLLISION";
		}
		throw std::invalid_argument("unknown collision type");
	}

	const char* binaryCollisionName(BinaryCollisionMethod method) const {
		switch (method) {
		case BinaryCollisionMethod::TA:
			return "TA";
		}
		throw std::invalid_argument("unknown binary collision method");
	}

	const char* hdpCouplingName(HdpCouplingMode mode) const {
		switch (mode) {
		case HdpCouplingMode::Decoupled:
			return "DECOUPLED";
		case HdpCouplingMode::VarianceWeighted:
			return "VARIANCE_WEIGHTED";
		}
		throw std::invalid_argument("unknown HDP coupling mode");
	}

	const char* effectiveWeightPolicyName(EffectiveWeightPolicy policy) const {
		switch (policy) {
		case EffectiveWeightPolicy::Fixed:
			return "FIXED";
		case EffectiveWeightPolicy::QuadraticAdaptive:
			return "QUADRATIC_ADAPTIVE";
		}
		throw std::invalid_argument("unknown effective-weight policy");
	}

	const char* collisionCouplingName(CollisionCoupling coupling) const {
		switch (coupling) {
		case CollisionCoupling::Standard:
			return "STANDARD";
		case CollisionCoupling::Linearized:
			return "LINEARIZED";
		}
		throw std::invalid_argument("unknown collision coupling");
	}

	const char* projectionModeName(ProjectionMode mode) const {
		switch (mode) {
		case ProjectionMode::FullMicroMacro:
			return "FULL_MICRO_MACRO";
		case ProjectionMode::MaxwellianOnly:
			return "MAXWELLIAN_ONLY";
		}
		throw std::invalid_argument("unknown projection mode");
	}

	const char* deltaMModeName(DeltaMMode mode) const {
		switch (mode) {
		case DeltaMMode::Enabled:
			return "ENABLED";
		case DeltaMMode::Disabled:
			return "DISABLED";
		}
		throw std::invalid_argument("unknown Delta-M mode");
	}
};

} // namespace coulomb
