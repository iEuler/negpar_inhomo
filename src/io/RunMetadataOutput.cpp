#include "RunMetadataOutput.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "Grid.h"
#include "MacroOutput.h"
#include "OutputPaths.h"
#include "SimulationConfig.h"

namespace coulomb {
void RunMetadataOutput::saveGrid(const NumericGridClass& grid,
								 const SimulationState& state) {
	MacroOutput{}.saveMacro(grid.x, "x", state);
	MacroOutput{}.saveMacro(grid.vx, "v", state);
}

void RunMetadataOutput::saveParameters(const ParaClass& parameters,
									   const NumericGridClass& grid,
									   const SimulationState& state) {
	std::ofstream file(OutputPaths(state).resolve("parameter.txt"));
	if (!file)
		throw std::runtime_error("Unable to open parameter output file");
	file << std::setprecision(15) << "coeffBinarycoll "
		 << parameters.coeffBinaryColl << '\n'
		 << "resampleRatio " << parameters.resampleRatio << '\n'
		 << "npickupNeg " << parameters.nPickupNeg << '\n'
		 << "nfreq " << parameters.nfreq << '\n'
		 << "xmax " << grid.xmax << '\n'
		 << "xmin " << grid.xmin << '\n'
		 << "vmax " << grid.vmax << '\n'
		 << "vmin " << grid.vmin << '\n'
		 << "tmax " << grid.tmax << '\n'
		 << "nx " << grid.nx << '\n'
		 << "nt " << grid.nt << '\n'
		 << "nv " << grid.nv << '\n'
		 << "dx " << grid.dx << '\n'
		 << "dv " << grid.dv << '\n'
		 << "dt " << grid.dt << '\n'
		 << "neff " << grid.neff << '\n'
		 << "neffF " << grid.neffF << '\n'
		 << "collision_kernel "
		 << SimulationTypes{}.collisionName(parameters.collisionType) << '\n'
		 << "lambdaPoisson " << parameters.lambdaPoisson << '\n'
		 << "resampleSpatialRatio " << parameters.resampleSpatialRatio << '\n'
		 << "syncTimeInterval " << parameters.syncTimeInterval << '\n'
		 << "resampleSyncRatio " << parameters.resampleSyncRatio << '\n'
		 << "hdpCoupling "
		 << SimulationTypes{}.hdpCouplingName(parameters.hdpCouplingMode) << '\n'
		 << "effectiveWeightPolicy "
		 << SimulationTypes{}.effectiveWeightPolicyName(
				parameters.effectiveWeightPolicy)
		 << '\n'
		 << "collisionCoupling "
		 << SimulationTypes{}.collisionCouplingName(
				parameters.collisionCoupling)
		 << '\n'
		 << "projectionMode "
		 << SimulationTypes{}.projectionModeName(parameters.projectionMode)
		 << '\n'
		 << "deltaMMode "
		 << SimulationTypes{}.deltaMModeName(parameters.deltaMMode) << '\n'
		 << "weightedFourierResampling "
		 << (parameters.weightedFourierResampling ? "true" : "false") << '\n'
		 << "partialResampling "
		 << (parameters.partialResampling ? "true" : "false") << '\n'
		 << "partialResamplingCutoff "
		 << parameters.partialResamplingCutoff << '\n'
		 << "conserveWeightedMoments "
		 << (parameters.conserveWeightedMoments ? "true" : "false") << '\n'
		 << "signedWeightMin " << parameters.signedWeightMin << '\n'
		 << "signedWeightMax " << parameters.signedWeightMax << '\n'
		 << "fullWeightMin " << parameters.fullWeightMin << '\n'
		 << "fullWeightMax " << parameters.fullWeightMax << '\n'
		 << "cpuCostConstant " << parameters.cpuCostConstant << '\n'
		 << "cpuCostCollisionCoefficient "
		 << parameters.cpuCostCollisionCoefficient << '\n'
		 << "coulombBgkHybrid "
		 << (parameters.coulombBgkHybrid ? "true" : "false") << '\n'
		 << "bgkStrength " << parameters.bgkStrength << '\n';

	std::ofstream secondFile(OutputPaths(state).resolve("parameter2.txt"));
	if (!secondFile)
		throw std::runtime_error(
			"Unable to open secondary parameter output file");
	secondFile << std::setprecision(15) << "methodBinarycoll "
			   << SimulationTypes{}.binaryCollisionName(
					  parameters.methodBinaryColl)
			   << '\n'
			   << "bdryX " << SimulationTypes{}.encodeBoundary(grid.bdryX)
			   << '\n'
			   << "bdryV " << SimulationTypes{}.encodeBoundary(grid.bdryV)
			   << '\n'
			   << "method " << SimulationTypes{}.methodName(parameters.method)
			   << '\n';
}

void RunMetadataOutput::saveInitialConditions(IniValClass& initialData,
											  const SimulationState& state) {
	std::ofstream file(OutputPaths(state).resolve("parameter.txt"),
					   std::ios::app);
	if (!file)
		throw std::runtime_error("Unable to append parameter output");
	file << std::setprecision(15);
	if (initialData.problemName == "LandauDamping") {
		file << "ldAlpha " << initialData.ldAlpha << '\n'
			 << "totalMass " << initialData.totalMass << '\n';
	} else if (initialData.problemName == "BumpOnTail") {
		file << "botBeta " << initialData.botBeta << '\n'
			 << "botRho0 " << initialData.botRho0 << '\n'
			 << "botTprt " << initialData.botTprt << '\n'
			 << "botDTprt " << initialData.botDTprt << '\n'
			 << "botTx " << initialData.botTx << '\n'
			 << "botUb " << initialData.botUb << '\n';
	}
	std::ofstream secondFile(OutputPaths(state).resolve("parameter2.txt"),
							 std::ios::app);
	if (!secondFile)
		throw std::runtime_error("Unable to append secondary parameter output");
	secondFile << "problem_name " << initialData.problemName << '\n';
}

} // namespace coulomb
