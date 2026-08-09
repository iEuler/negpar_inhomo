#include "SimulationTypes.h"

namespace coulomb {

BoundaryCondition boundary_condition_from_code(char code) {
  switch (code) {
    case 'p': return BoundaryCondition::Periodic;
    case 'n': return BoundaryCondition::Reflective;
    default:
      throw std::invalid_argument("boundary condition must be 'p' or 'n'");
  }
}

char boundary_condition_code(BoundaryCondition condition) {
  switch (condition) {
    case BoundaryCondition::Periodic: return 'p';
    case BoundaryCondition::Reflective: return 'n';
  }
  throw std::invalid_argument("unknown boundary condition");
}

}  // namespace coulomb
