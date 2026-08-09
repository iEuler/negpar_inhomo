#pragma once

namespace coulomb {

class IniValClass;
class NumericGridClass;

IniValClass make_initial_conditions(NumericGridClass& grid);
void configure_initial_condition(IniValClass& initial_data,
                                 const NumericGridClass& grid, int cell);
void initialize_TwoStreamInstab(IniValClass& initial_data);

}  // namespace coulomb
