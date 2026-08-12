#pragma once

#include <vector>

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParaClass;
struct SimulationState;

class SimulationSteps {
public:
  SimulationSteps(NumericGridClass &grid, ParaClass &parameters,
                  SimulationState &state)
      : grid_(grid), parameters_(parameters), state_(state) {}

  void advance_hdp(std::vector<NeParticleGroup> &groups);
  void advance_pic(std::vector<NeParticleGroup> &groups);

private:
  NumericGridClass &grid_;
  ParaClass &parameters_;
  SimulationState &state_;
};

} // namespace coulomb
