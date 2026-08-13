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
	void saveComplex(int size, std::complex<double>* values,
					 const std::string& filename, const SimulationState& state);
	void save2D(int rows, int columns,
				const std::vector<std::vector<double>>& values,
				const std::string& filename, const SimulationState& state);
	void saveRhouT(int size, const std::vector<double>& rho,
				   const std::vector<double>& u1, const std::vector<double>& u2,
				   const std::vector<double>& u3,
				   const std::vector<double>& temperature,
				   const SimulationState& state);
	void saveRhouT(const std::vector<ParticleGroup>& groups,
				   const NumericGridClass& grid, const SimulationState& state);
	void saveRhouT(std::vector<NeParticleGroup>& groups,
				   const NumericGridClass& grid, const SimulationState& state);
	void saveRhouTF(std::vector<NeParticleGroup>& groups,
					const NumericGridClass& grid, const SimulationState& state);
	void saveRhouTMaxwellian(const std::vector<NeParticleGroup>& groups,
							 const NumericGridClass& grid,
							 const SimulationState& state);
	void saveE(std::vector<NeParticleGroup>& groups,
			   const NumericGridClass& grid, const SimulationState& state);
	void saveE(std::vector<ParticleGroup>& groups, const NumericGridClass& grid,
			   const SimulationState& state);
	void saveNpNn(std::vector<NeParticleGroup>& groups,
				  const NumericGridClass& grid, const SimulationState& state);
	void saveDRhouT(std::vector<NeParticleGroup>& groups,
					const NumericGridClass& grid, const SimulationState& state);
	void saveDxRhouTM(std::vector<NeParticleGroup>& groups,
					  const NumericGridClass& grid,
					  const SimulationState& state);
	void saveM012Pn(std::vector<NeParticleGroup>& groups,
					const NumericGridClass& grid, const SimulationState& state);
	void saveMacroEvolution(std::vector<NeParticleGroup>& groups,
							const NumericGridClass& grid,
							const SimulationState& state);

	template <class T>
	void saveMacro(const std::vector<T>& macro, const std::string& filename,
				   const SimulationState& state) {
		const std::string suffix = ".txt";
		const std::string filenameWithSuffix =
			filename +
			(state.filenameWithNumber
				 ? OutputPaths(state).formatIndex(state.saveIndex)
				 : "") +
			suffix;
		const std::string path = OutputPaths(state).resolve(filenameWithSuffix);

		std::ofstream file(path);
		if (!file)
			throw std::runtime_error("Unable to open output file: " + path);
		file << std::setprecision(15);
		for (const auto& value : macro)
			file << value << '\n';
	}
};

} // namespace coulomb
