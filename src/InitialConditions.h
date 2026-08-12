#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;

class InitialConditions {
public:
  static IniValClass create(NumericGridClass &grid);
  static void configure(IniValClass &initial_data, const NumericGridClass &grid,
                        int cell);
  static void configure_two_stream(IniValClass &initial_data);
};

} // namespace coulomb
