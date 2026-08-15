#include "ParticleResampling.h"

#include "Grid.h"
#include "ParticleGroup.h"
#include "SimulationConfig.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "Diagnostics.h"
#include "FullParticleSampling.h"
#include "MacroOutput.h"
#include "ParticleConservation.h"
#include "ParticlePartition.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "Resampler.h"
#include "EffectiveWeightSelector.h"

namespace coulomb {
using std::cout;
using std::endl;
using std::vector;

bool ParticleResampling::resampleHomogeneous(
	NeParticleGroup& sX, double previousSignedWeight,
	double previousFullWeight) {
	const auto& para = parametersRef;
	auto& state = stateRef;
	// // cout << " resample 0" << endl;

	state.resampleCount++;

	int npOld = sX.size(ParticleKind::Positive);
	int nnOld = sX.size(ParticleKind::Negative);

	// int Nmax = 2*max(S_x . size('p'), S_x . size('n'));

	// cout << " resample 1" << endl;

	// resample particles
	resampling::FourierResamplerConfig config;
	config.frequencyCount = static_cast<std::size_t>(para.nfreq);
	config.useApproximation = true;
	config.fullParticleWeight = gridRef.neffF;
	config.weightedCoupling = para.weightedFourierResampling;
	const auto resampleSigned = [&](const NeParticleGroup& source) {
		if (source.size(ParticleKind::Positive) == 0 &&
			source.size(ParticleKind::Negative) == 0) {
			NeParticleGroup empty = source;
			empty.clear(ParticleKind::Positive);
			empty.clear(ParticleKind::Negative);
			empty.clear(ParticleKind::Full);
			return empty;
		}
		return resampling::FourierResampler(source, config).resample(
			state.random);
	};
	const auto conserveSignedMoments = [&](const NeParticleGroup& source,
										 NeParticleGroup& reconstructed) {
		if (!para.conserveWeightedMoments ||
			reconstructed.size(ParticleKind::Positive) == 0 ||
			reconstructed.size(ParticleKind::Negative) == 0)
			return;
		auto target = source;
		target.computeMoments();
		ParticleConservation{}.enforce(
			gridRef.neff *
				(target.positiveMoments.m0 - target.negativeMoments.m0),
			gridRef.neff *
				(target.positiveMoments.m11 - target.negativeMoments.m11),
			gridRef.neff *
				(target.positiveMoments.m12 - target.negativeMoments.m12),
			gridRef.neff *
				(target.positiveMoments.m13 - target.negativeMoments.m13),
			gridRef.neff *
				(target.positiveMoments.m21 - target.negativeMoments.m21),
			gridRef.neff *
				(target.positiveMoments.m22 - target.negativeMoments.m22),
			gridRef.neff *
				(target.positiveMoments.m23 - target.negativeMoments.m23),
			reconstructed, gridRef.neff, true, state.random);
	};
	const auto adjustTail = [&](NeParticleGroup& tail, ParticleKind kind,
								double oldWeight, double newWeight) {
		if (!(oldWeight > 0.0) || !(newWeight > 0.0) ||
			std::abs(oldWeight - newWeight) < 1e-14)
			return;
		if (newWeight >= oldWeight) {
			const double removeRate = 1.0 - oldWeight / newWeight;
			const int removeCount = std::min(
				tail.size(kind),
				RandomSampling(state.random).stochasticFloor(
					removeRate * tail.size(kind)));
			for (int index = 0; index < removeCount; ++index) {
				if (tail.size(kind) == 0)
					break;
				const int removeIndex = static_cast<int>(
					RandomSampling(state.random).uniform() * tail.size(kind));
				tail.erase(removeIndex, kind);
			}
			return;
		}
		const int duplicateCount = RandomSampling(state.random).stochasticFloor(
			(oldWeight / newWeight - 1.0) * tail.size(kind));
		const auto original = tail.list(kind);
		for (int index = 0; index < duplicateCount; ++index) {
			if (original.empty())
				break;
			const int sourceIndex = static_cast<int>(
				RandomSampling(state.random).uniform() * original.size());
			tail.pushBack(original[sourceIndex], kind);
		}
	};

	if (para.partialResampling) {
		const auto partition =
			ParticlePartitioning::split(sX, para.partialResamplingCutoff);
		auto coreNew = resampleSigned(partition.core);
		conserveSignedMoments(partition.core, coreNew);
		ParticleGroupOperations{}.assignPositions(
			coreNew, partition.core.getXMin(), partition.core.getXMax(),
			state.random);
		auto merged = partition.tail;
		if (previousSignedWeight > 0.0)
			for (const auto kind :
				 {ParticleKind::Positive, ParticleKind::Negative})
				adjustTail(merged, kind, previousSignedWeight, gridRef.neff);
		if (previousFullWeight > 0.0)
			adjustTail(merged, ParticleKind::Full, previousFullWeight,
					   gridRef.neffF);
		ParticleGroupOperations{}.mergeSigned(merged, coreNew);
		ParticleGroupOperations{}.mergeFull(merged, partition.core);
		const int npNew = merged.size(ParticleKind::Positive);
		const int nnNew = merged.size(ParticleKind::Negative);
		if (npNew < npOld && nnNew < nnOld) {
			sX = std::move(merged);
			sX.isResampled = true;
			return true;
		}
		cout << "New partial sampled particles rejected." << endl;
		return false;
	}

	auto sXNew = resampleSigned(sX);
	conserveSignedMoments(sX, sXNew);

	int npNew = sXNew.size(ParticleKind::Positive);
	int nnNew = sXNew.size(ParticleKind::Negative);

	// cout << " Resample finished." << endl;
	// cout << "After resampling N = (" << ptr_S_x_new .size('p') << ", " <<
	// ptr_S_x_new .size('n') << ");" << endl;
	ParticleGroupOperations{}.assignPositions(sXNew, sX.getXMin(), sX.getXMax(),
											  state.random);

	sX.isResampled = true;

	// replace old particles by new sampled particles

	// save_particles(S_x, ptr_S_x_new);

	// Replace the original particles by new sampled particles
	if ((npNew < npOld) && (nnNew < nnOld)) {
		// cout << "Replace by new sampled particles" << endl;
		sX.clear(ParticleKind::Positive);
		sX.clear(ParticleKind::Negative);
		ParticleGroupOperations{}.mergeSigned(sX, sXNew);
		return true;
	} else {
		cout << "New sampled particles rejected." << endl;
		return false;
	}
}

// void ParticleResampling::resample(NeParticleGroup * S_x, NumericGridClass &
// grid, ParaClass & para, MultlLevelGroup * MLsol) {
void ParticleResampling::resample(std::vector<NeParticleGroup>& sX) {
	auto& grid = gridRef;
	auto& para = parametersRef;
	auto& state = stateRef;
	bool needGlobalResample = false;

	bool flagResampleSuccess = true;

	for (int kx = 0; kx < grid.nx; kx++) {
		if ((sX[kx].size(ParticleKind::Positive) +
			 sX[kx].size(ParticleKind::Negative)) >=
			sX[kx].size(ParticleKind::Full))
			needGlobalResample = true;
	}
	if (needGlobalResample) {
		double resampleSpatialRatio = 0;
		for (int kx = 0; kx < grid.nx; kx++) {
			resampleSpatialRatio += (sX[kx].size(ParticleKind::Positive) +
									 sX[kx].size(ParticleKind::Negative)) /
									sX[kx].size(ParticleKind::Full);
		}
		resampleSpatialRatio /= grid.nx;

		resampleSpatialRatio = para.resampleSpatialRatio;

		// for (int kx = 0; kx < grid.nx; kx ++) {
		int kx = 0;
		while ((flagResampleSuccess) && (kx < grid.nx)) {
			if ((sX[kx].size(ParticleKind::Positive) +
				 sX[kx].size(ParticleKind::Negative)) >=
				resampleSpatialRatio * sX[kx].size(ParticleKind::Full)) {
				cout << "Particles resampling: ( "
					 << sX[kx].size(ParticleKind::Positive) << ", "
					 << sX[kx].size(ParticleKind::Negative) << ", "
					 << sX[kx].size(ParticleKind::Full) << ") " << endl;

				flagResampleSuccess = resampleHomogeneous(sX[kx]);

				cout << "After resampling: ( "
					 << sX[kx].size(ParticleKind::Positive) << ", "
					 << sX[kx].size(ParticleKind::Negative) << ", "
					 << sX[kx].size(ParticleKind::Full) << ") " << endl;
			}
			kx++;
		}
	}

	if (!flagResampleSuccess) {
		resampleFull(sX, grid.neffF / 2, para.nfreq);

		int nx = grid.nx;
		vector<double> rho(nx), rhoF(nx);

		for (int kx = 0; kx < nx; kx++)
			sX[kx].computeMoments();

		for (int kx = 0; kx < nx; kx++)
			rho[kx] = sX[kx].rhoM +
					  (sX[kx].positiveMoments.m0 - sX[kx].negativeMoments.m0) *
						  grid.neff / grid.dx;

		for (int kx = 0; kx < nx; kx++)
			rhoF[kx] = sX[kx].fullMoments.m0 * grid.neffF / grid.dx;

		MacroOutput{}.saveMacro<double>(rho, "rho_test", state);
		MacroOutput{}.saveMacro<double>(rhoF, "rhoF_test", state);
	}
}

void ParticleResampling::resampleFullHomogeneous(NeParticleGroup& sX,
												 double neffFNew, double neff,
												 int nfreq, double dxSpace) {
	auto& random = stateRef.random;
	// resample particles
	auto sXNew = FullParticleSampling{}.resample(sX, nfreq, neff, neffFNew,
												 dxSpace, random);

	ParticleGroupOperations{}.assignPositions(sXNew, sX.getXMin(), sX.getXMax(),
											  random);

	// replace old particles by new sampled particles

	sX.clear(ParticleKind::Full);
	ParticleGroupOperations{}.mergeFull(sX, sXNew);
}

void ParticleResampling::resampleFull(std::vector<NeParticleGroup>& sX,
									  double neffFNew, int nfreq) {
	auto& grid = gridRef;
	auto& state = stateRef;
	for (int kx = 0; kx < grid.nx; kx++) {
		resampleFullHomogeneous(sX[kx], neffFNew, grid.neff, nfreq, grid.dx);
	}

	grid.neffF = neffFNew;
	state.syncTime = 0;

	cout << "F particle resampled." << endl;
}

void ParticleResampling::resampleFullPreservingMass(
	std::vector<NeParticleGroup>& sX, int nfOld) {
	auto& grid = gridRef;
	auto& random = stateRef.random;
	int nfNew =
		Diagnostics(grid).particleCount(sX, grid.nx, ParticleKind::Full);
	if (nfNew > nfOld) {
		double neffFNew = grid.neffF;
		double totalMass = 0;

		for (int kx = 0; kx < grid.nx; kx++) {
			double mass = sX[kx].rhoM * grid.dx;
			totalMass += mass;
			int nk = (int)(mass / neffFNew);

			int nkRemove = sX[kx].size(ParticleKind::Full) - nk;

			for (int kp = 0; kp < nkRemove; kp++) {
				int kRemove = (int)(RandomSampling(random).uniform() *
									sX[kx].size(ParticleKind::Full));
				sX[kx].erase(kRemove, ParticleKind::Full);
			}
		}

		grid.neffF = totalMass / Diagnostics(grid).particleCount(
									 sX, grid.nx, ParticleKind::Full);
	}
}

void ParticleResampling::synchronizeCoarse(std::vector<NeParticleGroup>& sX) {
	auto& grid = gridRef;
	auto& para = parametersRef;
	auto& state = stateRef;
	if (para.collisionType == CollisionType::Coulomb) {
		if (state.syncTime > para.syncTimeInterval) {
			cout << "Start resample F" << endl;

			// cout << "First resample P and N" << endl;
			state.syncTime = 0;

			if (para.effectiveWeightPolicy ==
				EffectiveWeightPolicy::QuadraticAdaptive) {
				const auto originalGroups = sX;
				const double oldSignedWeight = grid.neff;
				const double oldFullWeight = grid.neffF;
				const int oldResampleCount = state.resampleCount;
				std::vector<EffectiveWeightCell> cells;
				cells.reserve(grid.nx);
				for (const auto& group : sX) {
					const double signedVariance =
						oldSignedWeight * oldSignedWeight *
						static_cast<double>(
							group.size(ParticleKind::Positive) +
							group.size(ParticleKind::Negative));
					const double fullVariance =
						oldFullWeight * oldFullWeight *
						static_cast<double>(group.size(ParticleKind::Full));
					const double denominator = signedVariance + fullVariance;
					cells.push_back(
						{denominator > 0.0 ? signedVariance / denominator : 0.0,
						 group.rhoM * grid.dx,
						 group.size(ParticleKind::Positive) +
							 group.size(ParticleKind::Negative),
						 group.size(ParticleKind::Full)});
				}
				const auto selection = EffectiveWeightSelector{}.select(
					oldSignedWeight, oldFullWeight, para.signedWeightMin,
					para.signedWeightMax, para.fullWeightMin,
					para.fullWeightMax, para.cpuCostConstant,
					para.cpuCostCollisionCoefficient, para.dt,
					para.coeffBinaryColl, cells);
				grid.neff = selection.signedWeight;
				grid.neffF = selection.fullWeight;

				bool adaptiveSuccess = true;
				for (auto& group : sX)
					adaptiveSuccess =
						ParticleResampling(grid, para, state)
							.resampleHomogeneous(group, oldSignedWeight,
												 oldFullWeight) &&
						adaptiveSuccess;
				if (adaptiveSuccess &&
					std::abs(selection.fullWeight - oldFullWeight) > 1e-14)
					resampleFull(sX, selection.fullWeight, para.nfreq);
				if (!adaptiveSuccess) {
					sX = originalGroups;
					grid.neff = oldSignedWeight;
					grid.neffF = oldFullWeight;
					state.resampleCount = oldResampleCount;
					state.syncTime = para.syncTimeInterval;
				}
				return;
			}

			resample(sX);

			cout << "P and N resampled" << endl;

			/*
			double totalMass = 0;
			for (int kx = 0; kx < grid.nx; kx ++)  totalMass += (S_x+kx) . rhoM;
			totalMass *= grid.dx;

			int Nd = 0;
			for (int kx = 0; kx < grid.nx; kx ++)  Nd += ( (S_x+kx) . size('p')
			+ (S_x+kx) . size('n') ); int Nf = (int)(1.2*Nd);

			double Neff_F_new = totalMass/Nf;
			if (Neff_F_new < grid.neffF) Neff_F_new = grid.neffF;
			*/

			double neffFNew = 100;
			for (int kx = 0; kx < grid.nx; kx++) {
				int nOne = (sX[kx].size(ParticleKind::Positive) +
							sX[kx].size(ParticleKind::Negative));
				double neffFOne = (sX[kx].rhoM) * grid.dx / nOne / 1.1;
				if (neffFNew > neffFOne)
					neffFNew = neffFOne;
			}

			if (neffFNew < grid.neffF)
				neffFNew = grid.neffF;

			// cout << "s resample F" << endl;

			int nfOld = Diagnostics(grid).particleCount(sX, grid.nx,
														ParticleKind::Full);

			resampleFull(sX, neffFNew, para.nfreq);
			cout << "F resampled" << endl;
			resampleFullPreservingMass(sX, nfOld);
		}
	}
}

} // namespace coulomb
