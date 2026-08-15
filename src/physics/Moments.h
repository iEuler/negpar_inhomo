#pragma once

#include <tuple>
#include <vector>

#include "SimulationTypes.h"

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;
class ParticleGroup;
struct SimulationState;

class MomentOperations {
  public:
	void primitiveToConserved(const double& rho, const double* velocity,
							  const double& temperature, double* momentum,
							  double& energy);
	void conservedToPrimitive(const double& rho, const double* momentum,
							  const double& energy, double* velocity,
							  double& temperature);

	void primitiveToConserved(int gridSize, const std::vector<double>& rho,
							  const std::vector<double>& velocity,
							  const std::vector<double>& temperature,
							  std::vector<double>& momentum,
							  std::vector<double>& energy);
	void conservedToPrimitive(int gridSize, const std::vector<double>& rho,
							  const std::vector<double>& momentum,
							  const std::vector<double>& energy,
							  std::vector<double>& velocity,
							  std::vector<double>& temperature);

	void particleToConserved(const ParticleGroup& group,
							 double effectiveParticles, double& rho,
							 double* momentum, double& energy);
	void computePrimitive(int gridSize,
						  const std::vector<ParticleGroup>& groups,
						  double effectiveParticles, std::vector<double>& rho,
						  std::vector<double>& u1, std::vector<double>& u2,
						  std::vector<double>& u3,
						  std::vector<double>& temperature);

	void updatePrimitive(std::vector<NeParticleGroup>& groups,
						 const NumericGridClass& grid);
	void updateFullPrimitive(std::vector<NeParticleGroup>& groups,
							 const NumericGridClass& grid);
	void updateMaxwellianDerivatives(std::vector<NeParticleGroup>& groups,
									 const NumericGridClass& grid);
	void updateMacro(std::vector<NeParticleGroup>& groups,
					 const NumericGridClass& grid);
	void updateMaxwellian(std::vector<NeParticleGroup>& groups,
						  const NumericGridClass& grid);
	void computeMacroChange(std::vector<NeParticleGroup>& groups,
							const NumericGridClass& grid,
							SimulationState& state,
							HdpCouplingMode couplingMode =
								HdpCouplingMode::Decoupled);
	void computeKineticMacroChange(std::vector<NeParticleGroup>& groups,
								   const NumericGridClass& grid);

	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
	momentChange(const NeParticleGroup* groups, const NumericGridClass& grid,
				 HdpCouplingMode couplingMode =
					 HdpCouplingMode::Decoupled);
};

} // namespace coulomb
