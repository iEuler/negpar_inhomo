#pragma once

// Full-particle Fourier reconstruction used by resampling orchestration.

#include <complex>
#include <vector>

namespace coulomb {

class NeParticleGroup;
struct RandomContext;

void interp3d_fft_approx_terms(
    std::vector<std::complex<double>>& Fouriercoeff, std::vector<double>& f,
    int nfreq1, int nfreq2, int nfreq3, int augmentation_factor, int orderx,
    int ordery, int orderz);
std::vector<std::complex<double>> interp3d_fft_approx(NeParticleGroup& groups,
                                                      int nfreq1, int nfreq2,
                                                      int nfreq3);
void filter_Fourier(std::vector<std::complex<double>>& Fouriercoeff,
                    std::vector<int>& flags, int size);
void resampleF_acceptsampled(const std::vector<double>& sample,
                             NeParticleGroup& groups, double value,
                             double& max_value, RandomContext& random);
std::vector<double> interp3d_fcoarse(
    const std::vector<std::complex<double>>& Fouriercoeff, int nfreq1,
    int nfreq2, int nfreq3);
std::vector<double> interp3d_fxyz_terms(
    const std::vector<std::complex<double>>& Fouriercoeff, int nfreq1,
    int nfreq2, int nfreq3, int augmentation_factor, int orderx, int ordery,
    int orderz);
std::vector<std::vector<double>> interp3d_fxyz(
    const std::vector<std::complex<double>>& Fouriercoeff, int nfreq1,
    int nfreq2, int nfreq3, int augmentation_factor);
std::vector<double> func_fourierupper3d(int count,
                                        const std::vector<double>& values);
std::vector<double> getKthValues(
    const std::vector<std::vector<double>>& values, int index);
void addMaxwellian_terms(double rhoM, std::vector<double> uM,
                         std::vector<double> TM, double effective_particles,
                         std::vector<double>& values, int frequency,
                         int augmentation_factor, int orderx, int ordery,
                         int orderz);
void addMaxwellian(double rhoM, std::vector<double> uM,
                   std::vector<double> TM, double effective_particles,
                   std::vector<std::vector<double>>& derivatives, int frequency,
                   int augmentation_factor);
NeParticleGroup resample_F_from_MPN(NeParticleGroup& groups, int frequency,
                                    double effective_particles,
                                    double effective_full_particles,
                                    double dx_space, RandomContext& random);

}  // namespace coulomb
