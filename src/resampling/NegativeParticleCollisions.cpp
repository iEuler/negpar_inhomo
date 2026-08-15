#include "NegativeParticleCollisions.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <omp.h>

#include "Collisions.h"
#include "Grid.h"
#include "NegativeParticleSampling.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {

void NegativeParticleCollisions::collideWithFull(NeParticleGroup& groups) {
	const auto& parameters = parametersRef;
	auto& random = randomContext;
	const int fullCount = groups.size(ParticleKind::Full);
	const int positiveCount = groups.size(ParticleKind::Positive);
	const int negativeCount = groups.size(ParticleKind::Negative);
	if (fullCount < positiveCount + negativeCount) {
		std::cout << "Too few F particles." << std::endl;
		std::cout << "(" << positiveCount << ", " << negativeCount << ", "
				  << fullCount << ") " << std::endl;
	}

	auto& positive = groups.list(ParticleKind::Positive);
	auto& negative = groups.list(ParticleKind::Negative);
	auto& full = groups.list(ParticleKind::Full);
	CollisionOperator collision(parameters, random);
	if (parameters.collisionCoupling == CollisionCoupling::Linearized) {
		const double sqrtTemperature = std::sqrt(groups.tprtM);
		std::array<double, 3> maxwellianVelocity{};
		auto sampleMaxwellianVelocity = [&]() {
			maxwellianVelocity[0] =
				groups.u1M + sqrtTemperature * RandomSampling(random).normal();
			maxwellianVelocity[1] =
				groups.u2M + sqrtTemperature * RandomSampling(random).normal();
			maxwellianVelocity[2] =
				groups.u3M + sqrtTemperature * RandomSampling(random).normal();
			return std::vector<double>(maxwellianVelocity.begin(),
									   maxwellianVelocity.end());
		};
		for (auto& particle : positive) {
			const auto velocities =
				collision.collidePair(particle.velocity(), sampleMaxwellianVelocity());
			particle.setVelocity(velocities.first);
		}
		for (auto& particle : negative) {
			const auto velocities =
				collision.collidePair(particle.velocity(), sampleMaxwellianVelocity());
			particle.setVelocity(velocities.first);
		}
		return;
	}

	const auto permutation = RandomSampling(random).permutation(
		fullCount, positiveCount + negativeCount);
	for (int index = 0; index < positiveCount; ++index) {
		const int fullIndex = permutation[index] - 1;
		const auto velocities =
			collision.collidePair(positive[index].velocity(), full[fullIndex].velocity());
		positive[index].setVelocity(velocities.first);
	}
	for (int index = 0; index < negativeCount; ++index) {
		const int fullIndex = permutation[index + positiveCount] - 1;
		const auto velocities =
			collision.collidePair(negative[index].velocity(), full[fullIndex].velocity());
		negative[index].setVelocity(velocities.first);
	}
}

void NegativeParticleCollisions::collideHomogeneous(NeParticleGroup& sX) {
	const auto& para = parametersRef;
	const double neff = gridRef.neff;
	auto& random = randomContext;
	NeParticleGroup sXNew;

	if (para.deltaMMode == DeltaMMode::Enabled)
		NegativeParticleSampling{}.sampleDelta(sX, sXNew, para, neff, random);
	ParticleGroupOperations{}.assignPositions(sXNew, sX.getXMin(), sX.getXMax(),
											  random);
	collideWithFull(sX);
	ParticleGroupOperations{}.mergeSigned(sX, sXNew);

	auto& sf = sX.list(ParticleKind::Full);
	CollisionOperator(para, random)
		.collideHomogeneous(sf, sX.size(ParticleKind::Full));
}

void NegativeParticleCollisions::collide(std::vector<NeParticleGroup>& sX) {
	const auto& grid = gridRef;
	const auto& para = parametersRef;
	NegativeParticleSampling{}.updateBounds(sX, grid, para);
	for (int kx = 0; kx < grid.nx; kx++)
		collideHomogeneous(sX[kx]);
}

void NegativeParticleCollisions::collideParallel(
	std::vector<NeParticleGroup>& sX) {
	const auto& grid = gridRef;
	const auto& para = parametersRef;
	NegativeParticleSampling{}.updateBounds(sX, grid, para);
#pragma omp parallel if (para.flagUseOpenMp)
	{
#pragma omp for
		for (int kx = 0; kx < grid.nx; kx++)
			collideHomogeneous(sX[kx]);
	}
}

void NegativeParticleCollisions::collideBgkHomogeneous(
	NeParticleGroup& groups) {
	auto& parameters = parametersRef;
	auto& random = randomContext;
	const int positiveCount = groups.size(ParticleKind::Positive);
	const int negativeCount = groups.size(ParticleKind::Negative);
	const int fullCount = groups.size(ParticleKind::Full);
	const double collisionScale =
		parameters.bgkStrength > 0.0 ? parameters.bgkStrength : 1.0;
	const double changeRate = std::clamp(
		parameters.dt * parameters.coeffBinaryColl * collisionScale, 0.0, 1.0);
	const int positiveRemove = RandomSampling(random).stochasticFloor(
		positiveCount * changeRate);
	const int negativeRemove = RandomSampling(random).stochasticFloor(
		negativeCount * changeRate);

	for (int index = 0; index < positiveRemove; ++index) {
		const int removeIndex =
			static_cast<int>(RandomSampling(random).uniform() *
							 groups.size(ParticleKind::Positive));
		groups.erase(removeIndex, ParticleKind::Positive);
	}
	for (int index = 0; index < negativeRemove; ++index) {
		const int removeIndex =
			static_cast<int>(RandomSampling(random).uniform() *
							 groups.size(ParticleKind::Negative));
		groups.erase(removeIndex, ParticleKind::Negative);
	}

	const double sqrtTemperature = std::sqrt(groups.tprtM);
	std::array<double, 3> velocity{};
	auto& full = groups.list(ParticleKind::Full);
	for (int index = 0; index < fullCount; ++index) {
		if (RandomSampling(random).uniform() < changeRate) {
			velocity[0] =
				groups.u1M + sqrtTemperature * RandomSampling(random).normal();
			velocity[1] =
				groups.u2M + sqrtTemperature * RandomSampling(random).normal();
			velocity[2] =
				groups.u3M + sqrtTemperature * RandomSampling(random).normal();
			full[index].setVelocity(velocity);
		}
	}
}

void NegativeParticleCollisions::collideBgk(
	std::vector<NeParticleGroup>& groups) {
	const auto& grid = gridRef;
	for (int cell = 0; cell < grid.nx; ++cell)
		collideBgkHomogeneous(groups[cell]);
}

} // namespace coulomb
