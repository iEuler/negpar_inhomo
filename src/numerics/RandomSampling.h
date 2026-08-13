#pragma once

#include <vector>

#include "RandomContext.h"

namespace coulomb {

class RandomSampling {
  public:
	explicit RandomSampling(RandomContext& context) : context(context) {}

	double uniform();
	double normal();
	std::vector<int> permutation(int inputSize, int outputSize);
	int stochasticFloor(double value);

  private:
	RandomContext& context;
};

} // namespace coulomb
