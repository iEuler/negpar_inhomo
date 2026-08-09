#pragma once

// Deterministic Fourier numerics used by full-particle reconstruction.

#include <complex>
#include <vector>

namespace coulomb {

class NeParticleGroup;

std::vector<std::complex<double>> interp3d_fft_approx(NeParticleGroup& groups,
                                                      int nfreq1, int nfreq2,
                                                      int nfreq3);
void filter_Fourier(std::vector<std::complex<double>>& Fouriercoeff,
                    std::vector<int>& flags, int size);
std::vector<double> interp3d_fcoarse(
    const std::vector<std::complex<double>>& Fouriercoeff, int nfreq1,
    int nfreq2, int nfreq3);
std::vector<std::vector<double>> interp3d_fxyz(
    const std::vector<std::complex<double>>& Fouriercoeff, int nfreq1,
    int nfreq2, int nfreq3, int augmentation_factor);
std::vector<double> func_fourierupper3d(int count,
                                        const std::vector<double>& values);
std::vector<double> getKthValues(
    const std::vector<std::vector<double>>& values, int index);
void addMaxwellian(double rhoM, std::vector<double> uM,
                   std::vector<double> TM, double effective_particles,
                   std::vector<std::vector<double>>& derivatives, int frequency,
                   int augmentation_factor);
}  // namespace coulomb
