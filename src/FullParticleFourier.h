#pragma once

// Deterministic Fourier numerics used by full-particle reconstruction.

#include <complex>
#include <vector>

namespace coulomb {

class NeParticleGroup;

class FullParticleFourier {
public:
  std::vector<std::complex<double>>
  approximate_transform(NeParticleGroup &groups, int frequency_x,
                        int frequency_y, int frequency_z);
  void filter(std::vector<std::complex<double>> &coefficients,
              std::vector<int> &flags, int size);
  std::vector<double>
  interpolate_coarse(const std::vector<std::complex<double>> &coefficients,
                     int frequency_x, int frequency_y, int frequency_z);
  std::vector<std::vector<double>>
  interpolate_derivatives(const std::vector<std::complex<double>> &coefficients,
                          int frequency_x, int frequency_y, int frequency_z,
                          int augmentation_factor);
  std::vector<double> upper_bound(int count, const std::vector<double> &values);
  std::vector<double> values_at(const std::vector<std::vector<double>> &values,
                                int index);
  void add_maxwellian(double density, std::vector<double> velocity,
                      std::vector<double> temperature,
                      double effective_particles,
                      std::vector<std::vector<double>> &derivatives,
                      int frequency, int augmentation_factor);
};
} // namespace coulomb
