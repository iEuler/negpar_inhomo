#include "utils.h"

#include <cstdint>
#include <limits>
#include <omp.h>
#include <random>
#include <stdexcept>

namespace coulomb {

namespace {

struct ThreadRandomState {
  std::mt19937 engine;
  std::uniform_real_distribution<> uniform{0, 1};
  std::normal_distribution<> normal{0, 1};
  const RandomContext* owner = nullptr;
  std::uint32_t seed = 0;
  std::uint32_t thread_id = std::numeric_limits<std::uint32_t>::max();
  std::uint64_t generation = 0;
};

ThreadRandomState& random_state(RandomContext& context) {
  thread_local ThreadRandomState state;
  const auto thread_id = static_cast<std::uint32_t>(omp_get_thread_num());
  if (state.owner != &context || state.seed != context.seed ||
      state.generation != context.generation ||
      state.thread_id != thread_id) {
    std::seed_seq sequence{context.seed, thread_id};
    state.engine.seed(sequence);
    state.uniform.reset();
    state.normal.reset();
    state.owner = &context;
    state.seed = context.seed;
    state.thread_id = thread_id;
    state.generation = context.generation;
  }
  return state;
}

}  // namespace

// unsigned int MYRANDOM_X = 1234567, MYRANDOM_M = 1<<30, MYRANDOM_A = 1664525,
// MYRANDOM_C = 1013904223;
//
// ========================================================================
// random number generator

double myrand(RandomContext& context) {
  double u = -1;
  auto& state = random_state(context);
  while ((u <= 0) || (u >= 1)) u = state.uniform(state.engine);
  return u;
}
// ========================================================================
// random number generator for standard normal distribution

double myrandn(RandomContext& context) {
  auto& state = random_state(context);
  return state.normal(state.engine);
}

void reseed_random(RandomContext& context, std::uint32_t seed) {
  context.seed = seed;
  ++context.generation;
}

std::uint32_t generate_random_seed() {
  std::random_device random_device;
  return random_device();
}

// ========================================================================
// ----------------------- a random permutation ----------------- */

std::vector<int> myrandperm(int Nin, int Nout, RandomContext& context) {
  // extract a nonrepeated length Nout vector p from (1:Nin), Nout=p.size()<=Nin

  if (Nout > Nin)
    throw std::runtime_error("Nout [" + std::to_string(Nout) +
                             "] must not be larger than Nin [" +
                             std::to_string(Nin) + "].");

  std::vector<int> p(Nout);
  std::vector<int> q(Nin);
  for (int j = 0; j < Nin; j++) q[j] = j + 1;

  int Nnow;

  for (int j = 0; j < Nout; j++) {
    Nnow = Nin - j;
    double u = myrand(context);
    int k = (int)(Nnow * u);
    p[j] = q[k];
    q[k] = q[Nnow - 1];
  }
  return p;
}

// ========================================================================
// ----------------------- randomly take integer part ----------------- */

int myfloor(double x, RandomContext& context) {
  int n = (int)(x);
  double y = x - n;

  if (myrand(context) < y) n++;

  return n;
}

// ========================================================================
// ----------------------- maxval ----------------- */
// return the max value of vec with size N
double maxval(const std::vector<double>& vec) {
  int N = static_cast<int>(vec.size());
  double a = vec[0];
  for (int k = 1; k < N; k++) a = (a > vec[k] ? a : vec[k]);
  return a;
}
// ========================================================================
// ----------------------- maxval ----------------- */
// return the max value of vec with size N
double minval(const std::vector<double>& vec) {
  int N = static_cast<int>(vec.size());
  double a = vec[0];
  for (int k = 1; k < N; k++) a = (a < vec[k] ? a : vec[k]);
  return a;
}

// ========================================================================
// ----------------------- histogram info on fixed bars ----------------- */

void histinfo_fixbar(const std::vector<double>& xdist,
                     std::vector<int>& numinbar,
                     double xmin, double xmax) {
  // the histogram info on fixed bars, ranging over [xmin xmax], with number
  // Nbar equivalent to the histc command in Matlab

  int Ndist = static_cast<int>(xdist.size());
  int Nbar = static_cast<int>(numinbar.size());

  double dx = (xmax - xmin) / Nbar;

  for (int kbar = 0; kbar < Nbar; kbar++) numinbar[kbar] = 0;

  double x;
  int kx;

  for (int j = 0; j < Ndist; j++) {
    x = xdist[j];
    kx = (int)((x - xmin) / dx);
    if (kx < 0) kx = 0;
    if (kx >= Nbar) kx = Nbar - 1;
    numinbar[kx]++;
  }
}
}  // namespace coulomb
