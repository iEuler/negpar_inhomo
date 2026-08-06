#pragma once

#include <complex>
#include <vector>

#include "Classes.h"
#include "Output.h"
#include "_global_variables.h"

namespace coulomb {

std::vector<double> interp3d_xyzminmax(NeParticleGroup& groups);
void interp3d_renormalize(NeParticleGroup& groups,
                          NeParticleGroup& normalized_groups);
void interp3d_rescale(std::vector<Particle1d3d>& particles, int count,
                      const std::vector<double>& xyz_minmax);
std::vector<std::complex<double>> interp_ifreq(int frequency_count);
int freq_sequence(int index, int frequency_count);
int freq_sequence_reverse(int frequency, int frequency_count);
std::vector<double> interp_freq(int frequency_count);
std::vector<int> interp_freq_aug(int frequency_count, int augmentation_factor);
std::vector<std::complex<double>> interp3d_fft(NeParticleGroup& groups,
                                               int nfreq1, int nfreq2,
                                               int nfreq3);
void interp3d_fft_approx_terms(
    std::vector<std::complex<double>>& Fouriercoeff, std::vector<double>& f,
    int nfreq1, int nfreq2, int nfreq3, int augmentation_factor, int orderx,
    int ordery, int orderz);
std::vector<std::complex<double>> interp3d_fft_approx(NeParticleGroup& groups,
                                                      int nfreq1, int nfreq2,
                                                      int nfreq3);
void filter_Fourier(std::vector<std::complex<double>>& Fouriercoeff,
                    std::vector<int>& flags, int size);
double interp3d_fvalue(
    const std::vector<double>& sample, 
    const std::vector<std::complex<double>>& Fouriercoeff,
    const std::vector<std::complex<double>>& ifreq1,
    const std::vector<std::complex<double>>& ifreq2,
    const std::vector<std::complex<double>>& ifreq3,
    std::vector<int>& flags, int nfreq1, int nfreq2, int nfreq3);
double interp3d_fvalue_approx(double deltax, double deltay, double deltaz,
                              const std::vector<double>& derivatives);
void interp3d_acceptsampled(const std::vector<double>& sample,
                            NeParticleGroup& groups, double value,
                            double& max_value, RandomContext& random);
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
NeParticleGroup samplefromfourier3d(NeParticleGroup& groups, int frequency,
                                    RandomContext& random);
void sampleF(NeParticleGroup& groups, double new_effective_particles,
             double old_effective_particles, RandomContext& random);
void sampleF_inhomo(std::vector<NeParticleGroup>& groups,
                    NumericGridClass& grid, ParaClass& parameters,
                    RandomContext& random);
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
bool particleresample_homo(NeParticleGroup& groups,
                           const ParaClass& parameters,
                           SimulationState& state);
void resampleF_homo(NeParticleGroup& groups, double new_effective_particles,
                    double effective_particles, int frequency, double dx_space,
                    RandomContext& random);
void resampleF_inhomo(std::vector<NeParticleGroup>& groups,
                      double new_effective_particles, NumericGridClass& grid,
                      int frequency, SimulationState& state);
void resampleF_keeptotalmass(std::vector<NeParticleGroup>& groups,
                             NumericGridClass& grid, int old_count,
                             RandomContext& random);
void particleresample_inhomo(std::vector<NeParticleGroup>& groups,
                             NumericGridClass& grid, ParaClass& parameters,
                             SimulationState& state);
void sync_coarse(std::vector<NeParticleGroup>& groups, NumericGridClass& grid,
                 ParaClass& parameters, SimulationState& state);

}  // namespace coulomb
