#pragma once

#include <cstdint>
#include <vector>

#include "_global_variables.h"

namespace coulomb {

double myrand(RandomContext& context);
double myrandn(RandomContext& context);
std::vector<int> myrandperm(int Nin, int Nout, RandomContext& context);
int myfloor(double x, RandomContext& context);
void reseed_random(RandomContext& context, std::uint32_t seed);
std::uint32_t generate_random_seed();
double maxval(const std::vector<double>& vec);
double minval(const std::vector<double>& vec);
void histinfo_fixbar(const std::vector<double>& xdist,
                     std::vector<int>& numinbar,
                     double xmin, double xmax);

}  // namespace coulomb
