#pragma once

#include <vector>

#include "RandomContext.h"

namespace coulomb {

double myrand(RandomContext& context);
double myrandn(RandomContext& context);
std::vector<int> myrandperm(int input_size, int output_size,
                            RandomContext& context);
int myfloor(double value, RandomContext& context);

}  // namespace coulomb
