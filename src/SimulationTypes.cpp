#include "SimulationTypes.h"

namespace coulomb {

BoundaryCondition SimulationTypes::decode_boundary(char code) const {
  switch (code) {
  case 'p':
    return BoundaryCondition::Periodic;
  case 'n':
    return BoundaryCondition::Reflective;
  default:
    throw std::invalid_argument("boundary condition must be 'p' or 'n'");
  }
}

char SimulationTypes::encode_boundary(BoundaryCondition condition) const {
  switch (condition) {
  case BoundaryCondition::Periodic:
    return 'p';
  case BoundaryCondition::Reflective:
    return 'n';
  }
  throw std::invalid_argument("unknown boundary condition");
}

} // namespace coulomb
