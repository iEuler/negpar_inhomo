#pragma once

#include "Classes.h"

namespace coulomb {

NeParticleGroup interp3dRenormalize(NeParticleGroup &S_x);

void interp3dInvertRenormalize(std::vector<Particle1d3d> &Sp,
                               const std::vector<double> &xyz_minmax);
/**
 the frequency grids
*/
int frequency(int kth, size_t Nfreq);

int frequencyInverse(int kfreq, size_t Nfreq);

std::vector<double> interpFrequencies(size_t Nfreq);

std::vector<size_t> augmentedLocation(size_t Nfreq, size_t augFactor);

std::vector<std::complex<double>> interpFrequenciesComplex(size_t Nfreq);

VectorBool3D filterFourierCoeff(VectorComplex3D &Fouriercoeff);

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
Vector3D upperBoundFunc(const Vector3D &fc);

std::vector<double> getValuesByLoc(const std::vector<Vector3D> &fvecs, int kx,
                                   int ky, int kz);

double fvalueApproxFromDeriv(double deltax, double deltay, double deltaz,
                             const std::vector<double> &fDeriv);

double fvalueFromFFT(const std::vector<double> &Sf,
                     const VectorComplex3D &Fouriercoeff,
                     const std::vector<std::complex<double>> &ifreq1,
                     const std::vector<std::complex<double>> &ifreq2,
                     const std::vector<std::complex<double>> &ifreq3,
                     const VectorBool3D &flag_Fouriercoeff);

void acceptSampled(const std::vector<double> &Sf, NeParticleGroup &S_x_incell,
                   double fval, double &maxf,
                   bool sampleFromFullDistribution = false);

void addMaxwellian(const NeParticleGroup &S_x, double Neff,
                   std::vector<Vector3D> &fDerivatives, int Nfreq,
                   int augFactor, double dxSpace);
}  // namespace coulomb