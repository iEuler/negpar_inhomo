#pragma once

#include "ParticleGroup.h"
#include "RandomContext.h"
#include "TensorTypes.h"

namespace coulomb::resampling {

VectorBool3D filterFourierCoeff(VectorComplex3D &Fouriercoeff);

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
Vector3D upperBoundFunc(const Vector3D &fc);

std::vector<double> getValuesByLoc(const std::vector<Vector3D> &fvecs, int kx,
                                   int ky, int kz);

double fvalueFromFFT(const std::vector<double> &Sf,
                     const VectorComplex3D &Fouriercoeff,
                     const std::vector<std::complex<double>> &ifreq1,
                     const std::vector<std::complex<double>> &ifreq2,
                     const std::vector<std::complex<double>> &ifreq3,
                     const VectorBool3D &flag_Fouriercoeff);

void acceptSampled(const std::vector<double> &Sf, NeParticleGroup &S_x_incell,
                   double fval, double &maxf, RandomContext& random);

}  // namespace coulomb::resampling
