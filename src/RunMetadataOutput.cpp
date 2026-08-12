#include "RunMetadataOutput.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "Grid.h"
#include "MacroOutput.h"
#include "OutputPaths.h"
#include "SimulationConfig.h"

namespace coulomb {
void RunMetadataOutput::save_grid(const NumericGridClass& grid,
								  const SimulationState& state) {
	MacroOutput{}.save_macro(grid.x, "x", state);
	MacroOutput{}.save_macro(grid.vx, "v", state);
}

void RunMetadataOutput::save_parameters(const ParaClass& parameters,
										const NumericGridClass& grid,
										const SimulationState& state) {
	std::ofstream file(OutputPaths(state).resolve("parameter.txt"));
	if (!file)
		throw std::runtime_error("Unable to open parameter output file");
	file << std::setprecision(15) << "coeff_binarycoll "
		 << parameters.coeff_binarycoll << '\n'
		 << "resample_ratio " << parameters.resample_ratio << '\n'
		 << "Npickup_neg " << parameters.Npickup_neg << '\n'
		 << "Nfreq " << parameters.Nfreq << '\n'
		 << "xmax " << grid.xmax << '\n'
		 << "xmin " << grid.xmin << '\n'
		 << "vmax " << grid.vmax << '\n'
		 << "vmin " << grid.vmin << '\n'
		 << "tmax " << grid.tmax << '\n'
		 << "Nx " << grid.Nx << '\n'
		 << "Nt " << grid.Nt << '\n'
		 << "Nv " << grid.Nv << '\n'
		 << "dx " << grid.dx << '\n'
		 << "dv " << grid.dv << '\n'
		 << "dt " << grid.dt << '\n'
		 << "Neff " << grid.Neff << '\n'
		 << "Neff_F " << grid.Neff_F << '\n'
		 << "collision_kernel "
		 << SimulationTypes{}.collision_name(parameters.collisionType) << '\n'
		 << "lambda_Poisson " << parameters.lambda_Poisson << '\n'
		 << "resample_spatial_ratio " << parameters.resample_spatial_ratio
		 << '\n'
		 << "sync_time_interval " << parameters.sync_time_interval << '\n'
		 << "resample_sync_ratio " << parameters.resample_sync_ratio << '\n';

	std::ofstream second_file(OutputPaths(state).resolve("parameter2.txt"));
	if (!second_file)
		throw std::runtime_error(
			"Unable to open secondary parameter output file");
	second_file << std::setprecision(15) << "method_binarycoll "
				<< SimulationTypes{}.binary_collision_name(
					   parameters.method_binarycoll)
				<< '\n'
				<< "bdry_x " << SimulationTypes{}.encode_boundary(grid.bdry_x)
				<< '\n'
				<< "bdry_v " << SimulationTypes{}.encode_boundary(grid.bdry_v)
				<< '\n'
				<< "method " << SimulationTypes{}.method_name(parameters.method)
				<< '\n';
}

void RunMetadataOutput::save_initial_conditions(IniValClass& initial_data,
												const SimulationState& state) {
	std::ofstream file(OutputPaths(state).resolve("parameter.txt"),
					   std::ios::app);
	if (!file)
		throw std::runtime_error("Unable to append parameter output");
	file << std::setprecision(15);
	if (initial_data.probname == "LandauDamping") {
		file << "LD_alpha " << initial_data.LD_alpha << '\n'
			 << "totalmass " << initial_data.totalmass << '\n';
	} else if (initial_data.probname == "BumpOnTail") {
		file << "BOT_beta " << initial_data.BOT_beta << '\n'
			 << "BOT_rho0 " << initial_data.BOT_rho0 << '\n'
			 << "BOT_Tprt " << initial_data.BOT_Tprt << '\n'
			 << "BOT_dTprt " << initial_data.BOT_dTprt << '\n'
			 << "BOT_Tx " << initial_data.BOT_Tx << '\n'
			 << "BOT_ub " << initial_data.BOT_ub << '\n';
	}
	std::ofstream second_file(OutputPaths(state).resolve("parameter2.txt"),
							  std::ios::app);
	if (!second_file)
		throw std::runtime_error("Unable to append secondary parameter output");
	second_file << "problem_name " << initial_data.probname << '\n';
}

} // namespace coulomb
