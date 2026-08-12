#pragma once

#include <vector>

#include "RandomContext.h"

namespace coulomb {

class RandomSampling {
public:
  static double uniform(RandomContext &context);
  static double normal(RandomContext &context);
  static std::vector<int> permutation(int input_size, int output_size,
                                      RandomContext &context);
  static int stochastic_floor(double value, RandomContext &context);
};

} // namespace coulomb
