#include "FullParticleSampling.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "FullParticleFourier.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"

namespace coulomb {
using std::abs;
using std::sqrt;
using std::vector;

namespace {

void resampleFAcceptSampled(const std::vector<double>& sf,
							NeParticleGroup& ptrSXInCell, double fval,
							double& maxF, RandomContext& random) {
	if (abs(fval) > maxF) {
		// keep sampled particles with rate maxF/maxF_new
		double keepRate = maxF / (1.5 * abs(fval));
		maxF = 1.5 * abs(fval);

		int npRemove = RandomSampling(random).stochasticFloor(
			(1 - keepRate) * ptrSXInCell.size(ParticleKind::Full));

		for (int kp = 0; kp < npRemove; kp++) {
			int kRemove = (int)(RandomSampling(random).uniform() *
								ptrSXInCell.size(ParticleKind::Full));
			ptrSXInCell.erase(kRemove, ParticleKind::Full);
		}
	}

	// Full particles represent a nonnegative density. Negative spectral
	// interpolation values must not be turned into extra positive mass.
	if (maxF > 0.0 && RandomSampling(random).uniform() < std::max(0.0, fval / maxF)) {
		std::vector<double> wrapped = sf;
		for (auto& component : wrapped) {
			component = std::fmod(component, 2.0 * pi);
			if (component < 0.0) component += 2.0 * pi;
		}
		ptrSXInCell.pushBack(Particle1D3D(wrapped), ParticleKind::Full);
	}
}

} // namespace

NeParticleGroup FullParticleSampling::resample(NeParticleGroup& sX, int nfreq,
											   double neff, double neffF,
											   double dxSpace,
											   RandomContext& random) {
	NeParticleGroup sXNew;
	/* Normalize particle velocity to [0 2*pi] */
	sX.setXyzRange();
	if (!(sX.tprtM > 0.0) || !std::isfinite(sX.tprtM) ||
		!(neffF > 0.0) || !std::isfinite(neffF))
		throw std::invalid_argument("Full reconstruction requires positive finite temperature and weight");
	// The Maxwellian remains broad even when the signed population is empty
	// or concentrated near one velocity. Signed-only bounds truncate its mass.
	const double thermalRadius = 6.0 * std::sqrt(sX.tprtM);
	const double mean[3] = {sX.u1M, sX.u2M, sX.u3M};
	for (int component = 0; component < 3; ++component) {
		sX.xyzMinMax[2 * component] = std::min(sX.xyzMinMax[2 * component], mean[component] - thermalRadius);
		sX.xyzMinMax[2 * component + 1] = std::max(sX.xyzMinMax[2 * component + 1], mean[component] + thermalRadius);
	}

	auto sXRenormalized = resampling::ResamplingVelocity{}.normalizeSigned(sX);

	const auto ifreq =
		resampling::ResamplingNumerics{}.imaginaryFrequencies(nfreq);
	vector<double> interpX(nfreq);
	for (int kx = 0; kx < nfreq; kx++)
		interpX[kx] = kx * 2 * pi / nfreq;

	vector<int> flagFouriercoeff(nfreq * nfreq * nfreq);

	auto fourierCoeff = FullParticleFourier{}.approximateTransform(
		sXRenormalized, nfreq, nfreq, nfreq);
	FullParticleFourier{}.filter(fourierCoeff, flagFouriercoeff,
								 nfreq * nfreq * nfreq);

	const auto fcoarse = FullParticleFourier{}.interpolateCoarse(
		fourierCoeff, nfreq, nfreq, nfreq);

	int augFactor = 2;
	auto fDerivatives = FullParticleFourier{}.interpolateDerivatives(
		fourierCoeff, nfreq, nfreq, nfreq, augFactor);

	vector<double> uM(3);
	vector<double> tm(3);
	double rhoM = sX.rhoM * dxSpace;
	uM[0] = sXRenormalized.u1M;
	uM[1] = sXRenormalized.u2M;
	uM[2] = sXRenormalized.u3M;
	tm[0] = sXRenormalized.t1M;
	tm[1] = sXRenormalized.t2M;
	tm[2] = sXRenormalized.t3M;

	FullParticleFourier{}.addMaxwellian(rhoM, uM, tm, neff, fDerivatives, nfreq,
										augFactor);
	const auto f = fDerivatives[0];

	const auto fUp = FullParticleFourier{}.upperBound(augFactor * nfreq, f);

	double dxaug = 2.0 * pi / nfreq / augFactor;
	vector<double> interpXaug(nfreq * augFactor);
	for (int kx = 0; kx < nfreq * augFactor; kx++)
		interpXaug[kx] = kx * 2 * pi / nfreq / augFactor;

	for (int kx = 0; kx < augFactor * nfreq; kx++) {
		for (int ky = 0; ky < augFactor * nfreq; ky++) {
			for (int kz = 0; kz < augFactor * nfreq; kz++) {
				int kk = kz + augFactor * nfreq * (ky + augFactor * nfreq * kx);

				double xc = interpXaug[kx];
				double yc = interpXaug[ky];
				double zc = interpXaug[kz];
				double fcc = fUp[kk];

				double maxF = 1.5 * abs(fcc);
				int nIncell = RandomSampling(random).stochasticFloor(
					maxF * dxaug * dxaug * dxaug / neffF);

				int kVirtual = 0;
				NeParticleGroup sXInCell;

				while (kVirtual < nIncell) {
					double deltax =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					double deltay =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					double deltaz =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					std::vector<double> sf{xc + deltax, yc + deltay,
										   zc + deltaz};

					const auto fDeriv =
						FullParticleFourier{}.valuesAt(fDerivatives, kk);
					double fval = resampling::ResamplingNumerics{}
									  .evaluateQuadraticTaylor(deltax, deltay,
															   deltaz, fDeriv);

					resampleFAcceptSampled(sf, sXInCell, fval, maxF, random);

					nIncell = RandomSampling(random).stochasticFloor(
						maxF / (neffF / (dxaug * dxaug * dxaug)));
					kVirtual++;
				}

				ParticleGroupOperations{}.mergeFull(sXNew, sXInCell);
			}
		}
	}

	auto& spSampled = sXNew.list(ParticleKind::Full);
	const auto& xyzMinMax = sX.xyzMinMax;
	resampling::ResamplingVelocity{}.restore(spSampled, xyzMinMax);

	std::cout << "# resampled F = " << sXNew.size(ParticleKind::Full)
			  << std::endl;

	return sXNew;
}

} // namespace coulomb
