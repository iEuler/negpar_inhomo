
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
	std::vector<double> electric_energy;
	std::vector<double> electric_energy_full;
	std::vector<double> total_energy;
	std::vector<double> total_energy_full;
	std::vector<double> full_effective_particle_count;
	std::vector<double> distribution_times;
	std::vector<double> advection_cpu_time;
	std::vector<double> collision_cpu_time;
	std::vector<double> total_cpu_time;
	std::vector<double> resampling_cpu_time;
	std::vector<int> positive_particle_count;
	std::vector<int> negative_particle_count;
	std::vector<int> full_particle_count;
	std::vector<int> resampling_count;

	void record_state(const std::vector<NeParticleGroup>& groups,
					  const NumericGridClass& grid, SimulationState& state) {
		Diagnostics diagnostics(grid);
		electric_energy.push_back(diagnostics.electric_energy(groups));
		electric_energy_full.push_back(
			diagnostics.full_electric_energy(groups));
		total_energy.push_back(diagnostics.total_energy(groups));
		total_energy_full.push_back(diagnostics.full_total_energy(groups));
		positive_particle_count.push_back(diagnostics.particle_count(
			groups, grid.Nx, ParticleKind::Positive));
		negative_particle_count.push_back(diagnostics.particle_count(
			groups, grid.Nx, ParticleKind::Negative));
		full_particle_count.push_back(
			diagnostics.particle_count(groups, grid.Nx, ParticleKind::Full));
		full_effective_particle_count.push_back(grid.Neff_F);
		resampling_count.push_back(state.resampleCount);
		state.resampleCount = 0;
	}

	void record_timing(const SimulationState& state) {
		total_cpu_time.push_back(static_cast<float>(state.t1All - state.t0All) /
								 CLOCKS_PER_SEC);
		advection_cpu_time.push_back(
			static_cast<float>(state.t1Advection - state.t0Advection) /
			CLOCKS_PER_SEC);
		collision_cpu_time.push_back(
			static_cast<float>(state.t1Collision - state.t0Collision) /
			CLOCKS_PER_SEC);
		resampling_cpu_time.push_back(
			static_cast<float>(state.t1Resampling - state.t0Resampling) /
			CLOCKS_PER_SEC);
	}

	void save_partial(const SimulationState& state) const {
		MacroOutput{}.save_macro(electric_energy, "elec_energy", state);
		MacroOutput{}.save_macro(electric_energy_full, "elec_energy_F", state);
		MacroOutput{}.save_macro(positive_particle_count, "Np_rec", state);
		MacroOutput{}.save_macro(negative_particle_count, "Nn_rec", state);
		MacroOutput{}.save_macro(full_effective_particle_count, "Neff_F_rec",
								 state);
		MacroOutput{}.save_macro(resampling_count, "num_resample", state);
	}

	void save_all(const SimulationState& state) const {
		MacroOutput{}.save_macro(electric_energy, "elec_energy", state);
		MacroOutput{}.save_macro(electric_energy_full, "elec_energy_F", state);
		MacroOutput{}.save_macro(total_energy, "total_energy", state);
		MacroOutput{}.save_macro(total_energy_full, "total_energy_F", state);
		MacroOutput{}.save_macro(total_cpu_time, "cputime_all", state);
		MacroOutput{}.save_macro(advection_cpu_time, "cputime_adve", state);
		MacroOutput{}.save_macro(collision_cpu_time, "cputime_coll", state);
		MacroOutput{}.save_macro(resampling_cpu_time, "cputime_resamp", state);
		MacroOutput{}.save_macro(positive_particle_count, "Np_rec", state);
		MacroOutput{}.save_macro(negative_particle_count, "Nn_rec", state);
		MacroOutput{}.save_macro(full_particle_count, "Nf_rec", state);
		MacroOutput{}.save_macro(full_effective_particle_count, "Neff_F_rec",
								 state);
		MacroOutput{}.save_macro(distribution_times, "time_dist", state);
		MacroOutput{}.save_macro(resampling_count, "num_resample", state);
	}
};

class SimulationRunner {
  public:
	SimulationRunner(const RunOptions& options, SimulationState& state)
		: options_(options), state_(state), parameters_(),
		  grid_(100, parameters_.method), groups_(grid_.Nx) {}

	int run() {
		initialize();
		for (int step = 0; step < grid_.Nt; ++step) {
			std::cout << "step " << step << '/' << grid_.Nt << endl;
			save_distribution_if_due(step);
			history_.record_state(groups_, grid_, state_);
			advance_one_step();
			save_checkpoint_if_due(step);
		}
		finalize();
		cout << "Finished" << endl;
		return 0;
	}

  private:
	void initialize() {
		if (options_.steps)
			grid_.Nt = *options_.steps;
		parameters_.dt = grid_.dt;
		grid_.lambda_Poisson = parameters_.lambda_Poisson;

		RunMetadataOutput{}.save_grid(grid_, state_);
		RunMetadataOutput{}.save_parameters(parameters_, grid_, state_);
		Initialization{}.initialize(grid_, groups_, state_);
		parameters_.lambda_Poisson = grid_.lambda_Poisson;

		cout << "method = " << SimulationTypes{}.method_name(parameters_.method)
			 << endl;
		state_.filenameWithNumber = false;
		MomentOperations{}.update_macro(groups_, grid_);
		ElectricFieldSolver(grid_).update(groups_);
		state_.syncTime = 0;
	}

	void save_distribution_if_due(int step) {
		if (step < next_distribution_step_)
			return;

		state_.filenameWithNumber = true;
		MacroOutput{}.save_macro_evolution(groups_, grid_, state_);
		++state_.saveIndex;
		next_distribution_step_ += 40;
		state_.filenameWithNumber = false;
		history_.distribution_times.push_back(step * grid_.dt);
	}

	void advance_one_step() {
		state_.t0All = clock();
		if (parameters_.method == SimulationMethod::HDP)
			SimulationSteps(grid_, parameters_, state_).advance_hdp(groups_);
		else
			SimulationSteps(grid_, parameters_, state_).advance_pic(groups_);
		state_.syncTime += grid_.dt;
		state_.t1All = clock();
		history_.record_timing(state_);
	}

	void save_checkpoint_if_due(int step) {
		if (step < next_checkpoint_step_)
			return;

		state_.filenameWithNumber = false;
		history_.save_all(state_);
		state_.filenameWithNumber = true;
		next_checkpoint_step_ += 50;
	}

	void finalize() {
		state_.filenameWithNumber = false;
		history_.save_partial(state_);
		MacroOutput{}.save_macro_evolution(groups_, grid_, state_);
		state_.filenameWithNumber = false;
		history_.save_all(state_);
	}

	const RunOptions& options_;
	SimulationState& state_;
	ParaClass parameters_;
	NumericGridClass grid_;
	std::vector<NeParticleGroup> groups_;
	SimulationHistory history_;
	int next_checkpoint_step_{};
	int next_distribution_step_{};
};

} // namespace

coulomb::Simulation::Simulation(coulomb::RunOptions options,
								coulomb::SimulationState& state)
	: options_(options), state_(state) {}

int coulomb::Simulation::run() {
	return SimulationRunner(options_, state_).run();
}
