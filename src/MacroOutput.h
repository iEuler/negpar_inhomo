#pragma once

#include <complex>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

#include "OutputPaths.h"
#include "SimulationState.h"

namespace coulomb {

class NumericGridClass;
class ParticleGroup;
class NeParticleGroup;

class MacroOutput {
public:
  void save_complex(int size, std::complex<double> *values,
                    const std::string &filename, const SimulationState &state);
  void save_2d(int rows, int columns,
               const std::vector<std::vector<double>> &values,
               const std::string &filename, const SimulationState &state);
  void save_rhouT(int size, const std::vector<double> &rho,
                  const std::vector<double> &u1, const std::vector<double> &u2,
                  const std::vector<double> &u3,
                  const std::vector<double> &temperature,
                  const SimulationState &state);
  void save_rhouT(const std::vector<ParticleGroup> &groups,
                  const NumericGridClass &grid, const SimulationState &state);
  void save_rhouT(std::vector<NeParticleGroup> &groups,
                  const NumericGridClass &grid, const SimulationState &state);
  void save_rhouT_F(std::vector<NeParticleGroup> &groups,
                    const NumericGridClass &grid, const SimulationState &state);
  void save_rhouT_maxwellian(const std::vector<NeParticleGroup> &groups,
                             const NumericGridClass &grid,
                             const SimulationState &state);
  void save_E(std::vector<NeParticleGroup> &groups,
              const NumericGridClass &grid, const SimulationState &state);
  void save_E(std::vector<ParticleGroup> &groups, const NumericGridClass &grid,
              const SimulationState &state);
  void save_NpNn(std::vector<NeParticleGroup> &groups,
                 const NumericGridClass &grid, const SimulationState &state);
  void save_d_rhouT(std::vector<NeParticleGroup> &groups,
                    const NumericGridClass &grid, const SimulationState &state);
  void save_dx_rhouT_M(std::vector<NeParticleGroup> &groups,
                       const NumericGridClass &grid,
                       const SimulationState &state);
  void save_m012_PN(std::vector<NeParticleGroup> &groups,
                    const NumericGridClass &grid, const SimulationState &state);
  void save_macro_evolution(std::vector<NeParticleGroup> &groups,
                            const NumericGridClass &grid,
                            const SimulationState &state);

  template <class T>
  void save_macro(const std::vector<T> &macro, const std::string &filename,
                  const SimulationState &state) {
    const std::string suffix = ".txt";
    const std::string filename_with_suffix =
        filename +
        (state.filenameWithNumber
             ? OutputPaths(state).format_index(state.saveIndex)
             : "") +
        suffix;
    const std::string path = OutputPaths(state).resolve(filename_with_suffix);

    std::ofstream file(path);
    if (!file)
      throw std::runtime_error("Unable to open output file: " + path);
    file << std::setprecision(15);
    for (const auto &value : macro)
      file << value << '\n';
  }
};

} // namespace coulomb
