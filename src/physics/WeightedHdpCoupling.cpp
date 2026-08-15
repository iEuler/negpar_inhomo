#include "WeightedHdpCoupling.h"

#include <algorithm>
#include <cmath>

#include "Grid.h"
#include "ParticleGroup.h"

namespace coulomb {
namespace {

double signedDensity(const NeParticleGroup& group,
					 const NumericGridClass& grid) {
	return group.rhoM + grid.neff / grid.dx *
							  (group.positiveMoments.m0 -
							   group.negativeMoments.m0);
}

double fullDensity(const NeParticleGroup& group,
				   const NumericGridClass& grid) {
	return grid.neffF / grid.dx * group.fullMoments.m0;
}

double signedMomentumX(const NeParticleGroup& group,
					   const NumericGridClass& grid) {
	return group.rhoM * group.u1M + grid.neff / grid.dx *
										   (group.positiveMoments.m11 -
											group.negativeMoments.m11);
}

double fullMomentumX(const NeParticleGroup& group,
					 const NumericGridClass& grid) {
	return grid.neffF / grid.dx * group.fullMoments.m11;
}

double signedEnergy(const NeParticleGroup& group,
					const NumericGridClass& grid) {
	const double maxwellianEnergy =
		0.5 * group.rhoM *
		(group.u1M * group.u1M + group.u2M * group.u2M +
		 group.u3M * group.u3M + 3.0 * group.tprtM);
	return maxwellianEnergy + 0.5 * grid.neff / grid.dx *
										   (group.positiveMoments.m2 -
											group.negativeMoments.m2);
}

double fullEnergy(const NeParticleGroup& group,
				  const NumericGridClass& grid) {
	return 0.5 * grid.neffF / grid.dx * group.fullMoments.m2;
}

template <typename SignedValue, typename FullValue>
double blend(const NeParticleGroup& group, const NumericGridClass& grid,
			 SignedValue signedValue, FullValue fullValue) {
	const double weight = WeightedHdpCoupling::blendWeight(group, grid);
	return (1.0 - weight) * signedValue(group, grid) +
		   weight * fullValue(group, grid);
}

double signedDensityFlux(const NeParticleGroup& group,
						 const NumericGridClass& grid) {
	return grid.neff / grid.dx *
		   (group.positiveMoments.m11 - group.negativeMoments.m11);
}

double fullDensityFlux(const NeParticleGroup& group,
					   const NumericGridClass& grid) {
	return grid.neffF / grid.dx * group.fullMoments.m11;
}

double signedMomentumFlux(const NeParticleGroup& group,
						  const NumericGridClass& grid) {
	return grid.neff / grid.dx *
		   (group.positiveMoments.m21 - group.negativeMoments.m21);
}

double fullMomentumFlux(const NeParticleGroup& group,
						const NumericGridClass& grid) {
	return grid.neffF / grid.dx * group.fullMoments.m21;
}

double signedEnergyFlux(const NeParticleGroup& group,
						const NumericGridClass& grid) {
	return 0.5 * grid.neff / grid.dx *
		   (group.positiveMoments.m31 - group.negativeMoments.m31);
}

double fullEnergyFlux(const NeParticleGroup& group,
					  const NumericGridClass& grid) {
	return 0.5 * grid.neffF / grid.dx * group.fullMoments.m31;
}

} // namespace

double WeightedHdpCoupling::blendWeight(const NeParticleGroup& group,
										const NumericGridClass& grid) {
	const double signedVariance =
		grid.neff * grid.neff *
		static_cast<double>(group.size(ParticleKind::Positive) +
							group.size(ParticleKind::Negative));
	const double fullVariance =
		grid.neffF * grid.neffF *
		static_cast<double>(group.size(ParticleKind::Full));
	const double denominator = signedVariance + fullVariance;
	if (!(denominator > 0.0) || !std::isfinite(denominator))
		return 0.0;
	return std::clamp(signedVariance / denominator, 0.0, 1.0);
}

double WeightedHdpCoupling::density(const NeParticleGroup& group,
									const NumericGridClass& grid) {
	return blend(group, grid, signedDensity, fullDensity);
}

double WeightedHdpCoupling::momentumX(const NeParticleGroup& group,
									  const NumericGridClass& grid) {
	return blend(group, grid, signedMomentumX, fullMomentumX);
}

double WeightedHdpCoupling::energy(const NeParticleGroup& group,
								   const NumericGridClass& grid) {
	return blend(group, grid, signedEnergy, fullEnergy);
}

double WeightedHdpCoupling::densityFlux(const NeParticleGroup& group,
										const NumericGridClass& grid) {
	return blend(group, grid, signedDensityFlux, fullDensityFlux);
}

double WeightedHdpCoupling::momentumFlux(const NeParticleGroup& group,
										 const NumericGridClass& grid) {
	return blend(group, grid, signedMomentumFlux, fullMomentumFlux);
}

double WeightedHdpCoupling::energyFlux(const NeParticleGroup& group,
									   const NumericGridClass& grid) {
	return blend(group, grid, signedEnergyFlux, fullEnergyFlux);
}

} // namespace coulomb
