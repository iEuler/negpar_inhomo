#include "MacroOutput.h"

#include <fstream>
#include <stdexcept>
#include <vector>

#include "Grid.h"
#include "Moments.h"
#include "ParticleGroup.h"
#include "ParticleOutput.h"

namespace coulomb {
void MacroOutput::saveMacroEvolution(std::vector<NeParticleGroup>& groups,
									 const NumericGridClass& grid,
									 const SimulationState& state) {
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].computeMoments();
	MomentOperations{}.updateMacro(groups, grid);
	MacroOutput::saveRhouT(groups, grid, state);
	MacroOutput::saveRhouTF(groups, grid, state);
	MacroOutput::saveE(groups, grid, state);
	ParticleOutput{}.saveDistribution(groups, grid, state);
	MacroOutput::saveRhouTMaxwellian(groups, grid, state);

	std::ofstream file(OutputPaths(state).resolve("numOfDist.txt"));
	if (!file)
		throw std::runtime_error("Unable to open distribution count output");
	file << state.saveIndex;
}

void MacroOutput::saveComplex(int size, std::complex<double>* values,
							  const std::string& filename,
							  const SimulationState& state) {
	if (size < 0 || (size > 0 && values == nullptr))
		throw std::invalid_argument(
			"MacroOutput::saveComplex received invalid input");

	const auto realPath = OutputPaths(state).resolve(filename + "_r.txt");
	const auto imaginaryPath = OutputPaths(state).resolve(filename + "_i.txt");
	std::ofstream realFile(realPath);
	std::ofstream imaginaryFile(imaginaryPath);
	if (!realFile || !imaginaryFile)
		throw std::runtime_error("Unable to open complex output files");
	for (int index = 0; index < size; ++index) {
		realFile << values[index].real() << '\n';
		imaginaryFile << values[index].imag() << '\n';
	}
}

void MacroOutput::save2D(int rows, int columns,
						 const std::vector<std::vector<double>>& values,
						 const std::string& filename,
						 const SimulationState& state) {
	if (rows < 0 || columns < 0 ||
		static_cast<std::size_t>(rows) != values.size())
		throw std::invalid_argument(
			"MacroOutput::save2d dimensions do not match rows");
	for (const auto& row : values)
		if (static_cast<std::size_t>(columns) != row.size())
			throw std::invalid_argument(
				"MacroOutput::save2d dimensions do not match columns");

	const auto numberedName =
		filename +
		(state.filenameWithNumber
			 ? OutputPaths(state).formatIndex(state.saveIndex)
			 : "") +
		".txt";
	const auto path = OutputPaths(state).resolve(numberedName);
	std::ofstream file(path);
	if (!file)
		throw std::runtime_error("Unable to open output file: " + path);
	for (const auto& row : values) {
		for (const auto value : row)
			file << value << ' ';
		file << '\n';
	}
}

void MacroOutput::saveRhouT(int size, const std::vector<double>& rho,
							const std::vector<double>& u1,
							const std::vector<double>& u2,
							const std::vector<double>& u3,
							const std::vector<double>& temperature,
							const SimulationState& state) {
	if (size < 0 || static_cast<std::size_t>(size) != rho.size() ||
		u1.size() != rho.size() || u2.size() != rho.size() ||
		u3.size() != rho.size() || temperature.size() != rho.size())
		throw std::invalid_argument(
			"MacroOutput::saveRhouT input size mismatch");
	MacroOutput::saveMacro(rho, "rho", state);
	MacroOutput::saveMacro(u1, "u1", state);
	MacroOutput::saveMacro(u2, "u2", state);
	MacroOutput::saveMacro(u3, "u3", state);
	MacroOutput::saveMacro(temperature, "tprt", state);
}

void MacroOutput::saveRhouT(const std::vector<ParticleGroup>& groups,
							const NumericGridClass& grid,
							const SimulationState& state) {
	const int size = grid.nx;
	std::vector<double> rho(size), u1(size), u2(size), u3(size),
		temperature(size);
	MomentOperations{}.computePrimitive(size, groups, grid.neff, rho, u1, u2,
										u3, temperature);
	MacroOutput::saveRhouT(size, rho, u1, u2, u3, temperature, state);
}

void MacroOutput::saveRhouT(std::vector<NeParticleGroup>& groups,
							const NumericGridClass& grid,
							const SimulationState& state) {
	MomentOperations{}.updatePrimitive(groups, grid);
	const int size = grid.nx;
	std::vector<double> rho(size), u1(size), temperature(size);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = groups[cell].rho;
		u1[cell] = groups[cell].u1;
		temperature[cell] = groups[cell].tprt;
	}
	MacroOutput::saveMacro(rho, "rho", state);
	MacroOutput::saveMacro(u1, "u1", state);
	MacroOutput::saveMacro(temperature, "tprt", state);
}

void MacroOutput::saveRhouTF(std::vector<NeParticleGroup>& groups,
							 const NumericGridClass& grid,
							 const SimulationState& state) {
	MomentOperations{}.updateFullPrimitive(groups, grid);
	const int size = grid.nx;
	std::vector<double> rho(size), u1(size), temperature(size);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = groups[cell].rhoF;
		u1[cell] = groups[cell].u1F;
		temperature[cell] = groups[cell].tprtF;
	}
	MacroOutput::saveMacro(rho, "rhoF", state);
	MacroOutput::saveMacro(u1, "u1F", state);
	MacroOutput::saveMacro(temperature, "tprtF", state);
}

void MacroOutput::saveRhouTMaxwellian(
	const std::vector<NeParticleGroup>& groups, const NumericGridClass& grid,
	const SimulationState& state) {
	const int size = grid.nx;
	std::vector<double> rho(size), u1(size), temperature(size);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = groups[cell].rhoM;
		u1[cell] = groups[cell].u1M;
		temperature[cell] = groups[cell].tprtM;
	}
	MacroOutput::saveMacro(rho, "rhoM", state);
	MacroOutput::saveMacro(u1, "u1M", state);
	MacroOutput::saveMacro(temperature, "tprtM", state);
}

void MacroOutput::saveE(std::vector<NeParticleGroup>& groups,
						const NumericGridClass& grid,
						const SimulationState& state) {
	MomentOperations{}.updatePrimitive(groups, grid);
	const int size = grid.nx;
	std::vector<double> electricField(size), electricFieldFull(size);
	for (int cell = 0; cell < size; ++cell) {
		electricField[cell] = groups[cell].elecField;
		electricFieldFull[cell] = groups[cell].elecFieldF;
	}
	MacroOutput::saveMacro(electricField, "elecfield", state);
	MacroOutput::saveMacro(electricFieldFull, "elecfieldF", state);
}

void MacroOutput::saveE(std::vector<ParticleGroup>& groups,
						const NumericGridClass& grid,
						const SimulationState& state) {
	const int size = grid.nx;
	std::vector<double> electricField(size);
	for (int cell = 0; cell < size; ++cell)
		electricField[cell] = groups[cell].elecField;
	MacroOutput::saveMacro(electricField, "E", state);
}

void MacroOutput::saveNpNn(std::vector<NeParticleGroup>& groups,
						   const NumericGridClass& grid,
						   const SimulationState& state) {
	const int size = grid.nx;
	std::vector<int> positive(size), negative(size);
	for (int cell = 0; cell < size; ++cell) {
		positive[cell] = groups[cell].size(ParticleKind::Positive);
		negative[cell] = groups[cell].size(ParticleKind::Negative);
	}
	MacroOutput::saveMacro(positive, "Np", state);
	MacroOutput::saveMacro(negative, "Nn", state);
}

void MacroOutput::saveDRhouT(std::vector<NeParticleGroup>& groups,
							 const NumericGridClass& grid,
							 const SimulationState& state) {
	const int size = grid.nx;
	std::vector<double> densityChange(size), momentumChange(size),
		energyChange(size);
	for (int cell = 0; cell < size; ++cell) {
		densityChange[cell] = groups[cell].drho;
		momentumChange[cell] = groups[cell].dm1;
		energyChange[cell] = groups[cell].dEnergy;
	}
	MacroOutput::saveMacro(densityChange, "drho", state);
	MacroOutput::saveMacro(momentumChange, "dm1", state);
	MacroOutput::saveMacro(energyChange, "denergy", state);
}

void MacroOutput::saveDxRhouTM(std::vector<NeParticleGroup>& groups,
							   const NumericGridClass& grid,
							   const SimulationState& state) {
	const int size = grid.nx;
	std::vector<double> densityGradient(size), velocityGradient(size),
		temperatureGradient(size);
	for (int cell = 0; cell < size; ++cell) {
		densityGradient[cell] = groups[cell].dxRhoM;
		velocityGradient[cell] = groups[cell].dxU1M;
		temperatureGradient[cell] = groups[cell].dxTprtM;
	}
	MacroOutput::saveMacro(densityGradient, "dxrho", state);
	MacroOutput::saveMacro(velocityGradient, "dxu1", state);
	MacroOutput::saveMacro(temperatureGradient, "dxTprt", state);
}

void MacroOutput::saveM012Pn(std::vector<NeParticleGroup>& groups,
							 const NumericGridClass& /*grid*/,
							 const SimulationState& state) {
	const int size = static_cast<int>(groups.size());
	std::vector<double> zeroth(size), first(size), second(size);
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].computeMoments();
		zeroth[cell] =
			groups[cell].positiveMoments.m0 - groups[cell].negativeMoments.m0;
		first[cell] =
			groups[cell].positiveMoments.m11 - groups[cell].negativeMoments.m11;
		second[cell] =
			groups[cell].positiveMoments.m2 - groups[cell].negativeMoments.m2;
	}
	MacroOutput::saveMacro(zeroth, "m0", state);
	MacroOutput::saveMacro(first, "m1", state);
	MacroOutput::saveMacro(second, "m2", state);
}

} // namespace coulomb
