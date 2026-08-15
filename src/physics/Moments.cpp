#include "Moments.h"

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationState.h"

#include "MacroOutput.h"
#include "Numerics.h"
#include "WeightedHdpCoupling.h"

#include <stdexcept>

namespace coulomb {

void MomentOperations::updateMacro(std::vector<NeParticleGroup>& groups,
								   const NumericGridClass& grid) {
	for (auto& group : groups)
		group.computeMoments();
	MomentOperations::updatePrimitive(groups, grid);
	MomentOperations::updateMaxwellianDerivatives(groups, grid);
}

void MomentOperations::computeMacroChange(std::vector<NeParticleGroup>& groups,
										  const NumericGridClass& grid,
										  SimulationState& state,
										  HdpCouplingMode couplingMode) {
	const int size = grid.nx;
	std::vector<double> density(size), velocity(size), temperature(size);
	std::vector<double> densityChange(size), momentumChange(size),
		energyChange(size);
	for (int cell = 0; cell < size; ++cell) {
		density[cell] = groups[cell].rhoM;
		velocity[cell] = groups[cell].u1M;
		temperature[cell] = groups[cell].tprtM;
	}

	Numerics{}.advanceKineticEuler(density, velocity, temperature, size,
								   grid.dx, grid.dt, grid.bdryX, densityChange,
								   momentumChange, energyChange);
	if (state.saveFlux) {
		MacroOutput{}.saveMacro(densityChange, "drho_euler", state);
		MacroOutput{}.saveMacro(momentumChange, "dm1_euler", state);
		MacroOutput{}.saveMacro(energyChange, "denergy_euler", state);
	}
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].drho = densityChange[cell];
		groups[cell].dm1 = momentumChange[cell];
		groups[cell].dEnergy = energyChange[cell];
	}

	std::tie(densityChange, momentumChange, energyChange) =
		MomentOperations::momentChange(groups.data(), grid, couplingMode);
	if (state.saveFlux) {
		MacroOutput{}.saveMacro(densityChange, "drhoG", state);
		MacroOutput{}.saveMacro(momentumChange, "dm1G", state);
		MacroOutput{}.saveMacro(energyChange, "denergyG", state);
	}
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].drhoG = densityChange[cell];
		groups[cell].dm1G = momentumChange[cell];
		groups[cell].dEnergyG = energyChange[cell];
		groups[cell].drho += densityChange[cell];
		groups[cell].dm1 += momentumChange[cell];
		groups[cell].dEnergy += energyChange[cell];
		groups[cell].dm1 -=
			grid.dt * groups[cell].rhoO * groups[cell].elecField;
		groups[cell].dEnergy -= grid.dt * groups[cell].rhoO * groups[cell].u1O *
								groups[cell].elecField;
	}
}

void MomentOperations::primitiveToConserved(const double& rho, const double* u,
											const double& temperature,
											double* momentum, double& energy) {
	if (u == nullptr || momentum == nullptr)
		throw std::invalid_argument("MomentOperations::primitiveToConserved "
									"velocity or momentum output is null");
	for (int component = 0; component < 3; ++component)
		momentum[component] = rho * u[component];

	double velocitySquared = 0.0;
	for (int component = 0; component < 3; ++component)
		velocitySquared += u[component] * u[component];
	energy = 0.5 * rho * (velocitySquared + 3.0 * temperature);
}

void MomentOperations::conservedToPrimitive(const double& rho,
											const double* momentum,
											const double& energy, double* u,
											double& temperature) {
	if (momentum == nullptr || u == nullptr)
		throw std::invalid_argument("MomentOperations::conservedToPrimitive "
									"momentum or velocity output is null");
	for (int component = 0; component < 3; ++component)
		u[component] = momentum[component] / rho;

	double velocitySquared = 0.0;
	for (int component = 0; component < 3; ++component)
		velocitySquared += u[component] * u[component];
	temperature = (energy * 2.0 / rho - velocitySquared) / 3.0;
}

void MomentOperations::primitiveToConserved(
	int gridSize, const std::vector<double>& rho, const std::vector<double>& u1,
	const std::vector<double>& temperature, std::vector<double>& m1,
	std::vector<double>& energy) {
	if (gridSize < 0 || static_cast<std::size_t>(gridSize) != rho.size() ||
		u1.size() != rho.size() || temperature.size() != rho.size() ||
		m1.size() != rho.size() || energy.size() != rho.size())
		throw std::invalid_argument(
			"MomentOperations::primitiveToConserved input size mismatch");
	double velocity[3] = {0.0, 0.0, 0.0};
	double momentum[3] = {0.0, 0.0, 0.0};
	double energyAtCell = 0.0;
	for (int cell = 0; cell < gridSize; ++cell) {
		velocity[0] = u1[cell];
		MomentOperations::primitiveToConserved(
			rho[cell], velocity, temperature[cell], momentum, energyAtCell);
		m1[cell] = momentum[0];
		energy[cell] = energyAtCell;
	}
}

void MomentOperations::conservedToPrimitive(int gridSize,
											const std::vector<double>& rho,
											const std::vector<double>& m1,
											const std::vector<double>& energy,
											std::vector<double>& u1,
											std::vector<double>& temperature) {
	if (gridSize < 0 || static_cast<std::size_t>(gridSize) != rho.size() ||
		m1.size() != rho.size() || energy.size() != rho.size() ||
		u1.size() != rho.size() || temperature.size() != rho.size())
		throw std::invalid_argument(
			"MomentOperations::conservedToPrimitive input size mismatch");
	double velocity[3] = {0.0, 0.0, 0.0};
	double momentum[3] = {0.0, 0.0, 0.0};
	double temperatureAtCell = 0.0;
	for (int cell = 0; cell < gridSize; ++cell) {
		momentum[0] = m1[cell];
		MomentOperations::conservedToPrimitive(
			rho[cell], momentum, energy[cell], velocity, temperatureAtCell);
		u1[cell] = velocity[0];
		temperature[cell] = temperatureAtCell;
	}
}

void MomentOperations::particleToConserved(const ParticleGroup& group,
										   double effectiveParticles,
										   double& rho, double* momentum,
										   double& energy) {
	if (momentum == nullptr)
		throw std::invalid_argument(
			"MomentOperations::particleToConserved momentum output is null");
	rho = group.moments.m0 * effectiveParticles;
	momentum[0] = group.moments.m11 * effectiveParticles;
	momentum[1] = group.moments.m12 * effectiveParticles;
	momentum[2] = group.moments.m13 * effectiveParticles;
	energy = 0.5 * group.moments.m2 * effectiveParticles;
}

void MomentOperations::computePrimitive(
	int gridSize, const std::vector<ParticleGroup>& groups,
	double effectiveParticles, std::vector<double>& rho,
	std::vector<double>& u1, std::vector<double>& u2, std::vector<double>& u3,
	std::vector<double>& temperature) {
	if (gridSize < 0 || static_cast<std::size_t>(gridSize) != groups.size() ||
		rho.size() != groups.size() || u1.size() != groups.size() ||
		u2.size() != groups.size() || u3.size() != groups.size() ||
		temperature.size() != groups.size())
		throw std::invalid_argument(
			"MomentOperations::computePrimitive input size mismatch");
	for (int cell = 0; cell < gridSize; ++cell) {
		double density = 0.0;
		double momentum[3] = {0.0, 0.0, 0.0};
		double energy = 0.0;
		double velocity[3] = {0.0, 0.0, 0.0};
		double cellTemperature = 0.0;
		MomentOperations::particleToConserved(groups[cell], effectiveParticles,
											  density, momentum, energy);
		MomentOperations::conservedToPrimitive(density, momentum, energy,
											   velocity, cellTemperature);
		rho[cell] = density;
		u1[cell] = velocity[0];
		u2[cell] = velocity[1];
		u3[cell] = velocity[2];
		temperature[cell] = cellTemperature;
	}
}

void MomentOperations::updatePrimitive(std::vector<NeParticleGroup>& groups,
									   const NumericGridClass& grid) {
	const int size = grid.nx;
	const double effectiveParticles = grid.neff;
	const double dx = grid.dx;
	std::vector<double> rhoM(size), u1M(size), temperatureM(size);
	std::vector<double> rhoPn(size), momentumPn(size), energyPn(size);
	std::vector<double> rho(size), momentum(size), energy(size);
	std::vector<double> u1(size), temperature(size);

	for (int cell = 0; cell < size; ++cell) {
		rhoM[cell] = groups[cell].rhoM;
		u1M[cell] = groups[cell].u1M;
		temperatureM[cell] = groups[cell].tprtM;
		rhoPn[cell] =
			effectiveParticles / dx *
			(groups[cell].positiveMoments.m0 - groups[cell].negativeMoments.m0);
		momentumPn[cell] = effectiveParticles / dx *
						   (groups[cell].positiveMoments.m11 -
							groups[cell].negativeMoments.m11);
		energyPn[cell] =
			0.5 * effectiveParticles / dx *
			(groups[cell].positiveMoments.m2 - groups[cell].negativeMoments.m2);
		rho[cell] = rhoM[cell];
	}
	MomentOperations::primitiveToConserved(size, rhoM, u1M, temperatureM,
										   momentum, energy);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] += rhoPn[cell];
		momentum[cell] += momentumPn[cell];
		energy[cell] += energyPn[cell];
	}
	MomentOperations::conservedToPrimitive(size, rho, momentum, energy, u1,
										   temperature);
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].rho = rho[cell];
		groups[cell].u1 = u1[cell];
		groups[cell].tprt = temperature[cell];
	}
}

void MomentOperations::updateFullPrimitive(std::vector<NeParticleGroup>& groups,
										   const NumericGridClass& grid) {
	const int size = grid.nx;
	const double effectiveParticles = grid.neffF;
	std::vector<double> rho(size), momentum(size), energy(size);
	std::vector<double> u1(size), temperature(size);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = effectiveParticles / grid.dx * groups[cell].fullMoments.m0;
		momentum[cell] =
			effectiveParticles / grid.dx * groups[cell].fullMoments.m11;
		energy[cell] =
			0.5 * effectiveParticles / grid.dx * groups[cell].fullMoments.m2;
	}
	MomentOperations::conservedToPrimitive(size, rho, momentum, energy, u1,
										   temperature);
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].rhoF = rho[cell];
		groups[cell].u1F = u1[cell];
		groups[cell].tprtF = temperature[cell];
	}
}

void MomentOperations::updateMaxwellianDerivatives(
	std::vector<NeParticleGroup>& groups, const NumericGridClass& grid) {
	const int size = grid.nx;
	std::vector<double> rho(size), u1(size), temperature(size);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = groups[cell].rhoM;
		u1[cell] = groups[cell].u1M;
		temperature[cell] = groups[cell].tprtM;
	}
	auto dxRho = Numerics{}.centralDifference(rho, size, grid.bdryX);
	auto dxU1 = Numerics{}.centralDifference(u1, size, grid.bdryX);
	auto dxTemperature =
		Numerics{}.centralDifference(temperature, size, grid.bdryX);
	for (int cell = 0; cell < size; ++cell) {
		dxRho[cell] /= grid.dx;
		dxU1[cell] /= grid.dx;
		dxTemperature[cell] /= grid.dx;
		groups[cell].dxRhoM = dxRho[cell];
		groups[cell].dxU1M = dxU1[cell];
		groups[cell].dxTprtM = dxTemperature[cell];
	}
}

void MomentOperations::updateMaxwellian(std::vector<NeParticleGroup>& groups,
										const NumericGridClass& grid) {
	const int size = grid.nx;
	std::vector<double> rho(size), velocity(size), temperature(size);
	std::vector<double> momentum(size), energy(size);

	for (int cell = 0; cell < size; ++cell) {
		rho[cell] = groups[cell].rhoM;
		velocity[cell] = groups[cell].u1M;
		temperature[cell] = groups[cell].tprtM;
	}

	MomentOperations::primitiveToConserved(size, rho, velocity, temperature,
										   momentum, energy);
	for (int cell = 0; cell < size; ++cell) {
		rho[cell] -= groups[cell].drho;
		momentum[cell] -= groups[cell].dm1;
		energy[cell] -= groups[cell].dEnergy;
	}

	MomentOperations::conservedToPrimitive(size, rho, momentum, energy,
										   velocity, temperature);
	for (int cell = 0; cell < size; ++cell) {
		if (!(temperature[cell] >= 0.0))
			throw std::runtime_error(
				"Maxwellian update produced invalid temperature");
		groups[cell].rhoM = rho[cell];
		groups[cell].u1M = velocity[cell];
		groups[cell].tprtM = temperature[cell];
	}
}

void MomentOperations::computeKineticMacroChange(
	std::vector<NeParticleGroup>& groups, const NumericGridClass& grid) {
	const int size = grid.nx;
	std::vector<double> density(size), velocity(size), temperature(size);
	std::vector<double> densityChange(size), momentumChange(size),
		energyChange(size);
	for (int cell = 0; cell < size; ++cell) {
		density[cell] = groups[cell].rhoM;
		velocity[cell] = groups[cell].u1M;
		temperature[cell] = groups[cell].tprtM;
	}

	Numerics{}.advanceKineticEuler(density, velocity, temperature, size,
								   grid.dx, grid.dt, grid.bdryX, densityChange,
								   momentumChange, energyChange);
	for (int cell = 0; cell < size; ++cell) {
		groups[cell].drho = densityChange[cell];
		groups[cell].dm1 = momentumChange[cell];
		groups[cell].dEnergy = energyChange[cell];
		groups[cell].drhoG = 0.0;
		groups[cell].dm1G = 0.0;
		groups[cell].dEnergyG = 0.0;
	}
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
MomentOperations::momentChange(const NeParticleGroup* groups,
							   const NumericGridClass& grid,
							   HdpCouplingMode couplingMode) {
	const int size = grid.nx;
	if (size < 0)
		throw std::invalid_argument(
			"MomentOperations::momentChange grid size is negative");
	if (size > 0 && groups == nullptr)
		throw std::invalid_argument(
			"MomentOperations::momentChange groups pointer is null");
	const double effectiveParticles = grid.neff;
	const double dx = grid.dx;
	std::vector<double> rhoFlux(size), momentumFlux(size), energyFlux(size);
	std::vector<double> rho(size), momentum(size);
	for (int cell = 0; cell < size; ++cell) {
		if (couplingMode == HdpCouplingMode::VarianceWeighted) {
			rhoFlux[cell] =
				WeightedHdpCoupling::densityFlux(groups[cell], grid);
			momentumFlux[cell] =
				WeightedHdpCoupling::momentumFlux(groups[cell], grid);
			energyFlux[cell] =
				WeightedHdpCoupling::energyFlux(groups[cell], grid);
			rho[cell] = WeightedHdpCoupling::density(groups[cell], grid);
			momentum[cell] =
				WeightedHdpCoupling::momentumX(groups[cell], grid);
		} else {
			rhoFlux[cell] = effectiveParticles / dx *
							(groups[cell].positiveMoments.m11 -
							 groups[cell].negativeMoments.m11);
			momentumFlux[cell] = effectiveParticles / dx *
								 (groups[cell].positiveMoments.m21 -
								  groups[cell].negativeMoments.m21);
			energyFlux[cell] = 0.5 * effectiveParticles / dx *
							   (groups[cell].positiveMoments.m31 -
								groups[cell].negativeMoments.m31);
			rho[cell] =
				effectiveParticles / dx *
				(groups[cell].positiveMoments.m0 -
				 groups[cell].negativeMoments.m0);
			momentum[cell] = effectiveParticles / dx *
							 (groups[cell].positiveMoments.m11 -
							  groups[cell].negativeMoments.m11);
		}
	}
	auto drho = Numerics{}.centralDifference(rhoFlux, size, grid.bdryX);
	auto dm1 = Numerics{}.centralDifference(momentumFlux, size, grid.bdryX);
	auto dEnergy = Numerics{}.centralDifference(energyFlux, size, grid.bdryX);
	const double timeSpaceRatio = grid.dt / grid.dx;
	for (int cell = 0; cell < size; ++cell) {
		drho[cell] *= timeSpaceRatio;
		dm1[cell] *= timeSpaceRatio;
		dEnergy[cell] *= timeSpaceRatio;
		dm1[cell] -= grid.dt * rho[cell] * groups[cell].elecField;
		dEnergy[cell] -= grid.dt * momentum[cell] * groups[cell].elecField;
	}
	return {std::move(drho), std::move(dm1), std::move(dEnergy)};
}

} // namespace coulomb
