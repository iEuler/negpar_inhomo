#pragma once

#include <fftw3.h>

#include <random>

using namespace std;

// Global variables
#define pi 3.1415926535897932
int FLAG_PRECOMPUTE_ALPHA_U = 1;
int K_SAVE_TIME = 0;                  // record the file numbers
bool FLAG_FILENAME_WITH_NUM = false;  // append number to the filename

bool FLAG_SAVEFLUX = false;
int NUM_MOVED = 0;
int NUM_RESAMPLE = 0;

double resample_spatial_ratio = 1.0;

double SYNC_TIME = 0.0;

// global variables on time
std::clock_t t0_all, t1_all, t0_coll, t1_coll, t0_adve, t1_adve, t0_resamp,
    t1_resamp;

template <class T>
void save_macro(const std::vector<T> &macro, std::string filename);

// class fwd declare
class ParaClass;
class IniValClass;
class NumericGridClass;
class Particle1d3d;
class ParticleGroup;
class NeParticleGroup;