#pragma once

// Deterministic Fourier numerics used by full-particle reconstruction.

#include <complex>
#include <vector>

namespace coulomb {

class NeParticleGroup;

class FullParticleFourier {
public:
  static std::vector<std::complex<double>>
  approximate_transform(NeParticleGroup &groups, int frequency_x,
                        int frequency_y, int frequency_z);
  static void filter(std::vector<std::complex<double>> &coefficients,
                     std::vector<int> &flags, int size);
  static std::vector<double>
  interpolate_coarse(const std::vector<std::complex<double>> &coefficients,
                     int frequency_x, int frequency_y, int frequency_z);
  static std::vector<std::vector<double>>
  interpolate_derivatives(const std::vector<std::complex<double>> &coefficients,
                          int frequency_x, int frequency_y, int frequency_z,
                          int augmentation_factor);
  static std::vector<double> upper_bound(int count,
                                         const std::vector<double> &values);
  static std::vector<double>
  values_at(const std::vector<std::vector<double>> &values, int index);
  static void add_maxwellian(double density, std::vector<double> velocity,
                             std::vector<double> temperature,
                             double effective_particles,
                             std::vector<std::vector<double>> &derivatives,
                             int frequency, int augmentation_factor);
};
} // namespace coulomb
