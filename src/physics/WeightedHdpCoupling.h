#pragma once

namespace coulomb {

class NeParticleGroup;
class NumericGridClass;

class WeightedHdpCoupling {
  public:
	static double blendWeight(const NeParticleGroup& group,
							  const NumericGridClass& grid);
	static double density(const NeParticleGroup& group,
						  const NumericGridClass& grid);
	static double momentumX(const NeParticleGroup& group,
							const NumericGridClass& grid);
	static double energy(const NeParticleGroup& group,
						 const NumericGridClass& grid);
	static double densityFlux(const NeParticleGroup& group,
							  const NumericGridClass& grid);
	static double momentumFlux(const NeParticleGroup& group,
							   const NumericGridClass& grid);
	static double energyFlux(const NeParticleGroup& group,
							 const NumericGridClass& grid);
};

} // namespace coulomb
