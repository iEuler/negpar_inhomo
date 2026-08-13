#include "ElectricField.h"

#include <complex>

#include "Constants.h"
#include "FFT.h"
#include "Grid.h"
#include "ParticleGroup.h"

namespace coulomb {

std::vector<double>
ElectricFieldSolver::solvePoisson(const std::vector<double>& rho,
								  double lambda) {
	const int gridSize = gridRef.nx;
	const double domainSize = gridRef.xmax - gridRef.xmin;
	std::vector<double> electricField(gridSize);
	Fft1D fftCalculator(gridSize);
	const auto rhoFft = fftCalculator.fft(rho);

	std::vector<double> waveNumbers(gridSize);
	for (int index = 0; index < gridSize / 2 + 1; ++index)
		waveNumbers[index] = index * 2.0 * pi / domainSize;
	for (int index = gridSize / 2 + 1; index < gridSize; ++index)
		waveNumbers[index] = (index - gridSize) * 2.0 * pi / domainSize;

	std::vector<std::complex<double>> electricFft(gridSize);
	electricFft[0] = {0.0, 0.0};
	for (int index = 1; index < gridSize; ++index) {
		electricFft[index] = {
			rhoFft[index].imag() / waveNumbers[index] / gridSize,
			-rhoFft[index].real() / waveNumbers[index] / gridSize};
	}

	const auto electric = fftCalculator.ifft(electricFft);
	for (int index = 0; index < gridSize; ++index)
		electricField[index] = lambda * electric[index];
	return electricField;
}

void ElectricFieldSolver::update(std::vector<ParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].computeMoments();

	std::vector<double> rho(grid.nx);
	for (int cell = 0; cell < grid.nx; ++cell)
		rho[cell] = groups[cell].moments.m0 * grid.neff / grid.dx;

	const auto electricField = solvePoisson(rho, 10.0);
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].elecField = electricField[cell];
}

void ElectricFieldSolver::update(std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].computeMoments();

	std::vector<double> rho(grid.nx);
	for (int cell = 0; cell < grid.nx; ++cell)
		rho[cell] = groups[cell].rhoM + (groups[cell].positiveMoments.m0 -
										 groups[cell].negativeMoments.m0) *
											grid.neff / grid.dx;

	auto electricField = solvePoisson(rho, grid.lambdaPoisson);
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].elecField = electricField[cell];

	for (int cell = 0; cell < grid.nx; ++cell)
		rho[cell] = groups[cell].fullMoments.m0 * grid.neffF / grid.dx;
	electricField = solvePoisson(rho, grid.lambdaPoisson);
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].elecFieldF = electricField[cell];
}

void ElectricFieldSolver::updatePic(std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].computeMoments();

	std::vector<double> rho(grid.nx);
	for (int cell = 0; cell < grid.nx; ++cell)
		rho[cell] = groups[cell].fullMoments.m0 * grid.neffF / grid.dx;
	const auto electricField = solvePoisson(rho, grid.lambdaPoisson);
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].elecFieldF = electricField[cell];
}

void ElectricFieldSolver::updateFromCoarse(
	std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].computeMoments();

	std::vector<double> rho(grid.nx);
	for (int cell = 0; cell < grid.nx; ++cell)
		rho[cell] = groups[cell].fullMoments.m0 * grid.neffF / grid.dx;
	const auto electricField = solvePoisson(rho, grid.lambdaPoisson);
	for (int cell = 0; cell < grid.nx; ++cell)
		groups[cell].elecField = electricField[cell];
}

void ElectricFieldSolver::clear(std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell) {
		groups[cell].computeMoments();
		groups[cell].elecField = 0.0;
		groups[cell].elecFieldF = 0.0;
	}
}

void ElectricFieldSolver::updateFromDensity(
	std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell) {
		groups[cell].computeMoments();
		groups[cell].elecField = groups[cell].rho;
		groups[cell].elecFieldF = groups[cell].rhoF;
	}
}

} // namespace coulomb
