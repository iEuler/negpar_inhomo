#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;

class InitialConditions {
  public:
	IniValClass create(NumericGridClass& grid, double landauAmplitude = 0.4);
	void configure(IniValClass& initialData, const NumericGridClass& grid,
				   int cell);
	void configureTwoStream(IniValClass& initialData);
};

} // namespace coulomb
