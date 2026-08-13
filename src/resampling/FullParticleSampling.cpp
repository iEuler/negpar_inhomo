#include "FullParticleSampling.h"

#include <cmath>
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

	// accept this particle with rate abs(fval/maxF)
	if (RandomSampling(random).uniform() < (abs(fval / maxF))) {
		double sumSfPiSq = 0.;
		for (int kv = 0; kv < 3; kv++)
			sumSfPiSq += (sf[kv] - pi) * (sf[kv] - pi);
		if (sqrt(sumSfPiSq) < pi) {
			Particle1D3D sOne({sf[0], sf[1], sf[2]});
			ptrSXInCell.pushBack(sOne, ParticleKind::Full);
		}
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
