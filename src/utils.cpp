#include "utils.h"

#include <random>
#include <stdexcept>

namespace coulomb {
// initialize random number generators
std::random_device RANDOMDEVICE0;
std::mt19937 RANDOM_GEN(RANDOMDEVICE0());

std::uniform_real_distribution<> DIS_UNIFORM(0, 1);
std::normal_distribution<> DIS_NORMAL(0, 1);

// unsigned int MYRANDOM_X = 1234567, MYRANDOM_M = 1<<30, MYRANDOM_A = 1664525,
// MYRANDOM_C = 1013904223;
//
// ========================================================================
// random number generator

double myrand() {
  double u = -1;
  while ((u <= 0) || (u >= 1)) u = DIS_UNIFORM(RANDOM_GEN);
  return u;
}
// ========================================================================
// random number generator for standard normal distribution

double myrandn() { return DIS_NORMAL(RANDOM_GEN); }

// ========================================================================
// ----------------------- a random permutation ----------------- */

std::vector<int> myrandperm(int Nin, int Nout) {
  // extract a nonrepeated length Nout vector p from (1:Nin), Nout=p.size()<=Nin

  if (Nout > Nin)
    throw std::runtime_error("Nout [" + std::to_string(Nout) +
                             "] must not be larger than Nin [" +
                             std::to_string(Nin) + "].");

  std::vector<int> p(Nout);
  vector<int> q(Nin);
  for (int j = 0; j < Nin; j++) q[j] = j + 1;

  int Nnow;

  for (int j = 0; j < Nout; j++) {
    Nnow = Nin - j;
    double u = myrand();
    int k = (int)(Nnow * u);
    p[j] = q[k];
    q[k] = q[Nnow - 1];
  }
  return p;
}

// ========================================================================
// ----------------------- randomly take integer part ----------------- */

int myfloor(double x) {
  int n = (int)(x);
  double y = x - n;

  if (myrand() < y) n++;

  return n;
}

// ========================================================================
// ----------------------- maxval ----------------- */
// return the max value of vec with size N
double maxval(const vector<double>& vec) {
  int N = static_cast<int>(vec.size());
  double a = vec[0];
  for (int k = 1; k < N; k++) a = (a > vec[k] ? a : vec[k]);
  return a;
}
// ========================================================================
// ----------------------- maxval ----------------- */
// return the max value of vec with size N
double minval(const vector<double>& vec) {
  int N = static_cast<int>(vec.size());
  double a = vec[0];
  for (int k = 1; k < N; k++) a = (a < vec[k] ? a : vec[k]);
  return a;
}

// ========================================================================
// ----------------------- histogram info on fixed bars ----------------- */

void histinfo_fixbar(const vector<double>& xdist, vector<int>& numinbar,
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