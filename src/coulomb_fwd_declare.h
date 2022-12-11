#pragma once

#include <fftw3.h>

#include <random>
#include <string>
#include <vector>

using namespace std;

namespace coulomb {
template <class T>
void save_macro(const std::vector<T>& macro, std::string filename);

// class fwd declare
class ParaClass;
class IniValClass;
class NumericGridClass;
class Particle1d3d;
class ParticleGroup;
class NeParticleGroup;
}  // namespace coulomb