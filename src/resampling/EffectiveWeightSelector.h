#pragma once

#include <vector>

namespace coulomb {

struct EffectiveWeightCell {
	double omega{};
	double macroMass{};
	int signedParticleCount{};
	int fullParticleCount{};
};

struct EffectiveWeightSelection {
	double signedWeight{};
	double fullWeight{};
};

class EffectiveWeightSelector {
  public:
	EffectiveWeightSelection select(
		double currentSignedWeight, double currentFullWeight,
		double signedWeightMin, double signedWeightMax, double fullWeightMin,
		double fullWeightMax, double cpuCostConstant,
		double cpuCostCollisionCoefficient, double dt,
		double collisionCoefficient,
		const std::vector<EffectiveWeightCell>& cells) const;

  private:
	static double interpolateQuadratic(const double* x, const double* y,
									   double xValue);
};

} // namespace coulomb
