#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;

class InitialConditions {
  public:
	IniValClass create(NumericGridClass& grid);
	void configure(IniValClass& initialData, const NumericGridClass& grid,
				   int cell);
	void configureTwoStream(IniValClass& initialData);
};

} // namespace coulomb
