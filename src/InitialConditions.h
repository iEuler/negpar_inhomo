#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;

class InitialConditions {
  public:
	IniValClass create(NumericGridClass& grid);
	void configure(IniValClass& initial_data, const NumericGridClass& grid,
				   int cell);
	void configure_two_stream(IniValClass& initial_data);
};

} // namespace coulomb
