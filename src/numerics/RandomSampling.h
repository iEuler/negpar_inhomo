#pragma once

#include <vector>

#include "RandomContext.h"

namespace coulomb {

class RandomSampling {
  public:
	explicit RandomSampling(RandomContext& context) : context_(context) {}

	double uniform();
	double normal();
	std::vector<int> permutation(int input_size, int output_size);
	int stochastic_floor(double value);

  private:
	RandomContext& context_;
};

} // namespace coulomb
