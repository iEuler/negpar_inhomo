#pragma once

#include "coulomb_fwd_declare.h"

double myrand();
double myrandn();
std::vector<int> myrandperm(int Nin, int Nout);
int myfloor(double x);
double maxval(const vector<double>& vec);
double minval(const vector<double>& vec);
void histinfo_fixbar(const vector<double>& xdist, vector<int>& numinbar,
                     double xmin, double xmax);