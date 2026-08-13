#pragma once

#include <stdexcept>

namespace coulomb {

enum class CollisionType { NoCollision, Coulomb, BGK };
enum class BinaryCollisionMethod { TA };
enum class BoundaryCondition { Periodic, Reflective };
enum class SimulationMethod { HDP, PIC };

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
};

} // namespace coulomb
