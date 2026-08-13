
#include <ctime>
#include <iostream>
#include <vector>

#include "Diagnostics.h"
#include "ElectricField.h"
#include "Initialization.h"
#include "MacroOutput.h"
#include "Moments.h"
#include "RunMetadataOutput.h"
#include "RunOptions.h"
#include "Simulation.h"
#include "SimulationSteps.h"

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationConfig.h"
#include "SimulationTypes.h"

// #define I complex<double>(0.,1.)

// g++ *.cpp -o out -lfftw3 -std=c++17 -fopenmp

// ===============================================================
// Main program

using namespace coulomb;
using std::cout;
using std::endl;

namespace {

struct SimulationHistory {
	std::vector<double> electricEnergy;
	std::vector<double> electricEnergyFull;
	std::vector<double> totalEnergy;
	std::vector<double> totalEnergyFull;
	std::vector<double> fullEffectiveParticleCount;
	std::vector<double> distributionTimes;
	std::vector<double> advectionCpuTime;
	std::vector<double> collisionCpuTime;
	std::vector<double> totalCpuTime;
	std::vector<double> resamplingCpuTime;
	std::vector<int> positiveParticleCount;
	std::vector<int> negativeParticleCount;
	std::vector<int> fullParticleCount;
	std::vector<int> resamplingCount;

	void recordState(const std::vector<NeParticleGroup>& groups,
					 const NumericGridClass& grid, SimulationState& state) {
		Diagnostics diagnostics(grid);
		electricEnergy.push_back(diagnostics.electricEnergy(groups));
		electricEnergyFull.push_back(diagnostics.fullElectricEnergy(groups));
		totalEnergy.push_back(diagnostics.totalEnergy(groups));
		totalEnergyFull.push_back(diagnostics.fullTotalEnergy(groups));
		positiveParticleCount.push_back(
			diagnostics.particleCount(groups, grid.nx, ParticleKind::Positive));
		negativeParticleCount.push_back(
			diagnostics.particleCount(groups, grid.nx, ParticleKind::Negative));
		fullParticleCount.push_back(
			diagnostics.particleCount(groups, grid.nx, ParticleKind::Full));
		fullEffectiveParticleCount.push_back(grid.neffF);
		resamplingCount.push_back(state.resampleCount);
		state.resampleCount = 0;
	}

	void recordTiming(const SimulationState& state) {
		totalCpuTime.push_back(static_cast<float>(state.t1All - state.t0All) /
							   CLOCKS_PER_SEC);
		advectionCpuTime.push_back(
			static_cast<float>(state.t1Advection - state.t0Advection) /
			CLOCKS_PER_SEC);
		collisionCpuTime.push_back(
			static_cast<float>(state.t1Collision - state.t0Collision) /
			CLOCKS_PER_SEC);
		resamplingCpuTime.push_back(
			static_cast<float>(state.t1Resampling - state.t0Resampling) /
			CLOCKS_PER_SEC);
	}

	void savePartial(const SimulationState& state) const {
		MacroOutput{}.saveMacro(electricEnergy, "elec_energy", state);
		MacroOutput{}.saveMacro(electricEnergyFull, "elec_energy_F", state);
		MacroOutput{}.saveMacro(positiveParticleCount, "Np_rec", state);
		MacroOutput{}.saveMacro(negativeParticleCount, "Nn_rec", state);
		MacroOutput{}.saveMacro(fullEffectiveParticleCount, "Neff_F_rec",
								state);
		MacroOutput{}.saveMacro(resamplingCount, "num_resample", state);
	}

	void saveAll(const SimulationState& state) const {
		MacroOutput{}.saveMacro(electricEnergy, "elec_energy", state);
		MacroOutput{}.saveMacro(electricEnergyFull, "elec_energy_F", state);
		MacroOutput{}.saveMacro(totalEnergy, "totalEnergy", state);
		MacroOutput{}.saveMacro(totalEnergyFull, "total_energy_F", state);
		MacroOutput{}.saveMacro(totalCpuTime, "cputime_all", state);
		MacroOutput{}.saveMacro(advectionCpuTime, "cputime_adve", state);
		MacroOutput{}.saveMacro(collisionCpuTime, "cputime_coll", state);
		MacroOutput{}.saveMacro(resamplingCpuTime, "cputime_resamp", state);
		MacroOutput{}.saveMacro(positiveParticleCount, "Np_rec", state);
		MacroOutput{}.saveMacro(negativeParticleCount, "Nn_rec", state);
		MacroOutput{}.saveMacro(fullParticleCount, "Nf_rec", state);
		MacroOutput{}.saveMacro(fullEffectiveParticleCount, "Neff_F_rec",
								state);
		MacroOutput{}.saveMacro(distributionTimes, "time_dist", state);
		MacroOutput{}.saveMacro(resamplingCount, "num_resample", state);
	}
};

class SimulationRunner {
  public:
	SimulationRunner(const RunOptions& options, SimulationState& state)
		: options(options), state(state), parameters(),
		  grid(100, parameters.method), groups(grid.nx) {}

	int run() {
		initialize();
		for (int step = 0; step < grid.nt; ++step) {
			std::cout << "step " << step << '/' << grid.nt << endl;
			saveDistributionIfDue(step);
			history.recordState(groups, grid, state);
			advanceOneStep();
			saveCheckpointIfDue(step);
		}
		finalize();
		cout << "Finished" << endl;
		return 0;
	}

  private:
	void initialize() {
		if (options.steps)
			grid.nt = *options.steps;
		parameters.dt = grid.dt;
		grid.lambdaPoisson = parameters.lambdaPoisson;

		RunMetadataOutput{}.saveGrid(grid, state);
		RunMetadataOutput{}.saveParameters(parameters, grid, state);
		Initialization{}.initialize(grid, groups, state);
		parameters.lambdaPoisson = grid.lambdaPoisson;

		cout << "method = " << SimulationTypes{}.methodName(parameters.method)
			 << endl;
		state.filenameWithNumber = false;
		MomentOperations{}.updateMacro(groups, grid);
		ElectricFieldSolver(grid).update(groups);
		state.syncTime = 0;
	}

	void saveDistributionIfDue(int step) {
		if (step < nextDistributionStep)
			return;

		state.filenameWithNumber = true;
		MacroOutput{}.saveMacroEvolution(groups, grid, state);
		++state.saveIndex;
		nextDistributionStep += 40;
		state.filenameWithNumber = false;
		history.distributionTimes.push_back(step * grid.dt);
	}

	void advanceOneStep() {
		state.t0All = clock();
		if (parameters.method == SimulationMethod::HDP)
			SimulationSteps(grid, parameters, state).advanceHdp(groups);
		else
			SimulationSteps(grid, parameters, state).advancePic(groups);
		state.syncTime += grid.dt;
		state.t1All = clock();
		history.recordTiming(state);
	}

	void saveCheckpointIfDue(int step) {
		if (step < nextCheckpointStep)
			return;

		state.filenameWithNumber = false;
		history.saveAll(state);
		state.filenameWithNumber = true;
		nextCheckpointStep += 50;
	}

	void finalize() {
		state.filenameWithNumber = false;
		history.savePartial(state);
		MacroOutput{}.saveMacroEvolution(groups, grid, state);
		state.filenameWithNumber = false;
		history.saveAll(state);
	}

	const RunOptions& options;
	SimulationState& state;
	ParaClass parameters;
	NumericGridClass grid;
	std::vector<NeParticleGroup> groups;
	SimulationHistory history;
	int nextCheckpointStep{};
	int nextDistributionStep{};
};

} // namespace

coulomb::Simulation::Simulation(coulomb::RunOptions options,
								coulomb::SimulationState& state)
	: optionsValue(options), stateRef(state) {}

int coulomb::Simulation::run() {
	return SimulationRunner(optionsValue, stateRef).run();
}
