#include "ParticleGroupOperations.h"

#include "ParticleGroup.h"
#include "RandomSampling.h"

namespace coulomb {

void ParticleGroupOperations::mergeSigned(NeParticleGroup& sX,
										  const NeParticleGroup& sXNew) {
	auto& sp = sXNew.list(ParticleKind::Positive);
	auto& sn = sXNew.list(ParticleKind::Negative);
	for (const auto& sone : sp)
		sX.pushBack(sone, ParticleKind::Positive);
	for (const auto& sone : sn)
		sX.pushBack(sone, ParticleKind::Negative);
}

void ParticleGroupOperations::mergeFull(NeParticleGroup& sX,
										const NeParticleGroup& sXNew) {
	auto& sf = sXNew.list(ParticleKind::Full);
	for (const auto& sone : sf)
		sX.pushBack(sone, ParticleKind::Full);
}

void ParticleGroupOperations::merge(NeParticleGroup& sX,
									const NeParticleGroup& sXNew,
									const std::vector<ParticleKind>& parTypes) {
	for (const auto kind : parTypes) {
		auto& sp = sXNew.list(kind);
		for (const auto& sone : sp)
			sX.pushBack(sone, kind);
	}
}

void ParticleGroupOperations::assignPositions(NeParticleGroup& sNew,
											  double xmin, double xmax,
											  RandomContext& random) {
	double x1 = xmin, x2 = xmax;
	auto& sp = sNew.list(ParticleKind::Positive);
	auto& sn = sNew.list(ParticleKind::Negative);
	auto& sf = sNew.list(ParticleKind::Full);
	for (int kp = 0; kp < sNew.size(ParticleKind::Positive); kp++)
		sp[kp].setPosition(RandomSampling(random).uniform() * (x2 - x1) + x1);
	for (int kp = 0; kp < sNew.size(ParticleKind::Negative); kp++)
		sn[kp].setPosition(RandomSampling(random).uniform() * (x2 - x1) + x1);
	for (int kp = 0; kp < sNew.size(ParticleKind::Full); kp++)
		sf[kp].setPosition(RandomSampling(random).uniform() * (x2 - x1) + x1);
}

} // namespace coulomb
