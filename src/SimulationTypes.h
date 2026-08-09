#pragma once

#include <stdexcept>

namespace coulomb {

enum class CollisionType { NoCollision, Coulomb, BGK };
enum class BinaryCollisionMethod { TA };
enum class BoundaryCondition { Periodic, Reflective };
enum class SimulationMethod { HDP, PIC };

BoundaryCondition boundary_condition_from_code(char code);
char boundary_condition_code(BoundaryCondition condition);

inline const char* method_name(SimulationMethod method) {
  switch (method) {
    case SimulationMethod::HDP: return "HDP";
    case SimulationMethod::PIC: return "PIC";
  }
  throw std::invalid_argument("unknown simulation method");
}

inline const char* collision_name(CollisionType type) {
  switch (type) {
    case CollisionType::NoCollision: return "NO_COLLISION";
    case CollisionType::Coulomb: return "COULOMB_COLLISION";
    case CollisionType::BGK: return "BGK_COLLISION";
  }
  throw std::invalid_argument("unknown collision type");
}

inline const char* binary_collision_name(BinaryCollisionMethod method) {
  switch (method) {
    case BinaryCollisionMethod::TA: return "TA";
  }
  throw std::invalid_argument("unknown binary collision method");
}

}  // namespace coulomb
