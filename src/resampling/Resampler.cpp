#include "Resampler.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Constants.h"
#include "FFT.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "ResamplerHelper.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"

namespace coulomb::resampling {

using std::abs;
using std::complex;
using std::exp;
using std::max;
using std::to_string;

FourierResampler::FourierResampler(const NeParticleGroup& particles,
								   FourierResamplerConfig config)
	: particlesValue(particles), neff(config.effectiveParticleWeight),
	  nfreq(config.frequencyCount), useApproximation(config.useApproximation),
	  maxSamplingAttempts(config.maxSamplingAttempts) {
	if (!(neff > 0.0))
		throw std::invalid_argument(
			"Resampler effective particle weight must be positive");
	if (nfreq < 2 || nfreq % 2 != 0)
		throw std::invalid_argument(
			"Resampler frequency count must be an even value of at least 2");
	if (nfreq >
		static_cast<size_t>(std::numeric_limits<int>::max()) / augFactor)
		throw std::invalid_argument("Resampler frequency grid is too large");
	if (maxSamplingAttempts == 0)
		throw std::invalid_argument(
			"Resampler sampling-attempt budget must be positive");
}

Vector3D
FourierResampler::funcOnAugGrid(const VectorComplex3D& fourierCoeff) const {
	auto fft3D = Fft3D(nfreq * augFactor, nfreq * augFactor, nfreq * augFactor);
	return fft3D.ifft(fourierCoeff);
}

Vector3D
FourierResampler::derivativesFromFftOneTerm(const VectorComplex3D& fourierCoeff,
											int orderx, int ordery,
											int orderz) const {
	const auto augmentedSize = augFactor * nfreq;
	auto fSaug = VectorComplex3D(
		augmentedSize,
		std::vector(augmentedSize,
					std::vector<std::complex<double>>(augmentedSize)));

	// 1i *freq
	const auto freq1 = resampling::ResamplingNumerics{}.frequencies(nfreq);
	const auto freq2 = resampling::ResamplingNumerics{}.frequencies(nfreq);
	const auto freq3 = resampling::ResamplingNumerics{}.frequencies(nfreq);

	const auto loc1 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);
	const auto loc2 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);
	const auto loc3 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);

	for (int kk1 = 0; kk1 < nfreq; kk1++) {
		for (int kk2 = 0; kk2 < nfreq; kk2++) {
			for (int kk3 = 0; kk3 < nfreq; kk3++) {
				size_t kk1Aug = loc1[kk1];
				size_t kk2Aug = loc2[kk2];
				size_t kk3Aug = loc3[kk3];

				double freq = 1.0;
				for (int kx = 0; kx < orderx; kx++)
					freq *= freq1[kk1];
				for (int kx = 0; kx < ordery; kx++)
					freq *= freq2[kk2];
				for (int kx = 0; kx < orderz; kx++)
					freq *= freq3[kk3];

				fSaug[kk1Aug][kk2Aug][kk3Aug] =
					freq * fourierCoeff[kk1][kk2][kk3];
			}
		}
	}

	return funcOnAugGrid(fSaug);
}

std::vector<Vector3D> FourierResampler::derivativesFromFft(
	const VectorComplex3D& fourierCoeff) const {
	const auto orders = std::vector<std::vector<int>>{
		{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
		{0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

	std::vector<Vector3D> fDerivatives;
	for (const auto& order : orders) {
		fDerivatives.push_back(derivativesFromFftOneTerm(fourierCoeff, order[0],
														 order[1], order[2]));
	}
	return fDerivatives;
}

VectorComplex3D FourierResampler::fft3DApproxOneterm(const Vector3D& f,
													 int orderx, int ordery,
													 int orderz) const {
	const auto n = augFactor * nfreq;
	auto fourierCoeff =
		std::vector(n, std::vector(n, std::vector<std::complex<double>>(n)));

	auto fft3D = Fft3D(n, n, n);

	const auto fSaug = fft3D.fft(f);

	// 1i *freq
	const auto freq1 = resampling::ResamplingNumerics{}.frequencies(nfreq);
	const auto freq2 = resampling::ResamplingNumerics{}.frequencies(nfreq);
	const auto freq3 = resampling::ResamplingNumerics{}.frequencies(nfreq);

	const auto loc1 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);
	const auto loc2 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);
	const auto loc3 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq, augFactor);

	for (size_t kk1 = 0; kk1 < nfreq; kk1++) {
		for (size_t kk2 = 0; kk2 < nfreq; kk2++) {
			for (size_t kk3 = 0; kk3 < nfreq; kk3++) {
				size_t kk1Aug = loc1[kk1];
				size_t kk2Aug = loc2[kk2];
				size_t kk3Aug = loc3[kk3];

				double freq = neff;
				for (int kx = 0; kx < orderx; kx++)
					freq *= freq1[kk1];
				for (int kx = 0; kx < ordery; kx++)
					freq *= freq2[kk2];
				for (int kx = 0; kx < orderz; kx++)
					freq *= freq3[kk3];

				if ((orderx + ordery + orderz) == 2) {
					if ((orderx == 2) || (ordery == 2) || (orderz == 2))
						freq *= -.5;
					else
						freq *= -1;
				}

				const auto augOne = fSaug[kk1Aug][kk2Aug][kk3Aug];
				if ((orderx + ordery + orderz) == 1) {
					fourierCoeff[kk1][kk2][kk3] =
						freq * complex<double>(augOne.imag(), -augOne.real());
				} else {
					fourierCoeff[kk1][kk2][kk3] = freq * augOne;
				}
			}
		}
	}

	return fourierCoeff;
}

NeParticleGroup FourierResampler::resample(RandomContext& random) const {
	NeParticleGroup sXNew;
	auto sX = particlesValue;

	const auto parTypes = std::vector<ParticleKind>{ParticleKind::Positive,
													ParticleKind::Negative};

	/* Normalize particle velocity to [0 2*pi] */
	sX.setXyzRange({ParticleKind::Positive, ParticleKind::Negative});
	auto sXRenormalized = resampling::ResamplingVelocity{}.normalizeSigned(sX);

	/* Prepare the grids in physical space and frequence space */
	// double dx = 2.0*pi/nfreq;

	const auto ifreq =
		resampling::ResamplingNumerics{}.imaginaryFrequencies(nfreq);
	std::vector<double> interpX(nfreq);
	for (int kx = 0; kx < nfreq; kx++)
		interpX[kx] = kx * 2 * pi / nfreq;

	/* Compute the Fourier coefficient */
	VectorComplex3D fourierCoeff =
		useApproximation ? fft3DApprox(sXRenormalized) : fft3D(sXRenormalized);

	// Apply the filter on Fourier coefficients
	auto flagFouriercoeff =
		ResamplerHelper(random).filterFourierCoefficients(fourierCoeff);

	// cout << " F coeff computed " << endl;

	/* Compute a coarse interpolation in physical space */
	//  const auto fcoarse =
	//  FullParticleFourier{}.interpolateCoarse(Fouriercoeff, nfreq, nfreq,
	//  nfreq);

	auto fDerivatives = derivativesFromFft(fourierCoeff);

	/* evaluate the upperbound of f */
	ResamplerHelper helper(random);
	const auto fUp = helper.upperBound(fDerivatives[0]);

	/* refined x grid */
	double dxaug = 2.0 * pi / nfreq / augFactor;
	std::vector<double> interpXaug(nfreq * augFactor);
	for (int kx = 0; kx < nfreq * augFactor; kx++)
		interpXaug[kx] = kx * 2 * pi / nfreq / augFactor;

	const auto& f = fDerivatives[0];
	/* create a NeParticleGroup to host the P and N particles in current cell */

	/* Start sampling */

	size_t samplingAttempts = 0;

	for (int kx = 0; kx < augFactor * nfreq; kx++) {
		for (int ky = 0; ky < augFactor * nfreq; ky++) {
			for (int kz = 0; kz < augFactor * nfreq; kz++) {
				double xc = interpXaug[kx];
				double yc = interpXaug[ky];
				double zc = interpXaug[kz];

				double fcc = fUp[kx][ky][kz];

				if (fcc < std::abs(f[kx][ky][kz]))
					throw std::exception("small bound!");

				double maxF = 1.5 * fcc;
				int nIncell = RandomSampling(random).stochasticFloor(
					maxF * dxaug * dxaug * dxaug / neff);

				int kVirtual = 0;
				NeParticleGroup sXInCell;
				const auto fDeriv = helper.valuesAt(fDerivatives, kx, ky, kz);

				while (kVirtual < nIncell) {
					if (++samplingAttempts > maxSamplingAttempts)
						throw std::runtime_error(
							"Resampler exceeded its sampling-attempt budget at "
							"cell (" +
							std::to_string(kx) + ", " + std::to_string(ky) +
							", " + std::to_string(kz) + ") with target " +
							std::to_string(nIncell) + " and envelope " +
							std::to_string(maxF));
					// create a particle in the cell
					// Sample offsets from the explicit RandomContext below.
					double deltax =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					double deltay =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					double deltaz =
						RandomSampling(random).uniform() * dxaug - 0.5 * dxaug;
					std::vector<double> sf{xc + deltax, yc + deltay,
										   zc + deltaz};

					// compute f at this point
					double fval = useApproximation
									  ? resampling::ResamplingNumerics{}
											.evaluateQuadraticTaylor(
												deltax, deltay, deltaz, fDeriv)
									  : helper.valueFromFft(sf, fourierCoeff,
															ifreq, ifreq, ifreq,
															flagFouriercoeff);

					// reset current cell if fval>maxF, otherwise continue
					// sampling in current cell
					helper.acceptSample(sf, sXInCell, fval, maxF);

					// reset N_incell if maxF is changed
					nIncell = RandomSampling(random).stochasticFloor(
						maxF * dxaug * dxaug * dxaug / neff);
					kVirtual++;
				}

				::coulomb::ParticleGroupOperations{}.merge(sXNew, sXInCell,
														   parTypes);
			}
		}
	}

	// cout << "Resampled." << endl;

	// rescale to the original coordinates

	const auto& xyzMinMax = sX.xyzMinMax;
	for (const auto parType : parTypes) {
		auto& spSampled = sXNew.list(parType);
		resampling::ResamplingVelocity{}.restore(spSampled, xyzMinMax);
	}

	// cout << "Rescaled." << endl;

	return sXNew;
	// return std::make_shared<NeParticleGroup>(S_x_new);
}

VectorComplex3D FourierResampler::fft3D(NeParticleGroup& sX) const {
	auto fourierCoeff = std::vector(
		nfreq, std::vector(nfreq, std::vector<std::complex<double>>(nfreq)));

	int np = sX.size(ParticleKind::Positive);
	int nn = sX.size(ParticleKind::Negative);

	auto& sp = sX.list(ParticleKind::Positive);
	auto& sn = sX.list(ParticleKind::Negative);

	// double Neff_temp = 1./Np;

	const auto ifreq1 =
		resampling::ResamplingNumerics{}.imaginaryFrequencies(nfreq);
	const auto ifreq2 =
		resampling::ResamplingNumerics{}.imaginaryFrequencies(nfreq);
	const auto ifreq3 =
		resampling::ResamplingNumerics{}.imaginaryFrequencies(nfreq);

	double cubic2Pi = 8.0 * pi * pi * pi;

	double coeffFft = 1. / cubic2Pi * nfreq * nfreq * nfreq;
	double maxFS = 0.0;

	// the (i,j,k)-th element of the array with size (nx,Ny,Nz), you would use
	// the expression an_array[k + Nz * (j + Ny * i)].

	for (int kk1 = 0; kk1 < nfreq; kk1++) {
		for (int kk2 = 0; kk2 < nfreq; kk2++) {
			for (int kk3 = 0; kk3 < nfreq; kk3++) {
				auto coeff = complex<double>(0., 0.);
				for (int kp = 0; kp < np; kp++) {
					auto& vp = sp[kp].velocity();
					complex<double> expterm = vp[0] * ifreq1[kk1] +
											  vp[1] * ifreq2[kk2] +
											  vp[2] * ifreq3[kk3];
					coeff += exp(-expterm);
				}
				for (int kn = 0; kn < nn; kn++) {
					auto& vp = sn[kn].velocity();
					complex<double> expterm = vp[0] * ifreq1[kk1] +
											  vp[1] * ifreq2[kk2] +
											  vp[2] * ifreq3[kk3];
					coeff -= exp(-expterm);
				}
				// Fouriercoeff[kk] *= neff * coeff_fft;
				coeff *= neff * coeffFft / (nfreq * nfreq * nfreq);
				maxFS = max(maxFS, abs(coeff));
				fourierCoeff[kk1][kk2][kk3] = coeff;
			}
		}
	}

	return fourierCoeff;
}

VectorComplex3D FourierResampler::fft3DApprox(NeParticleGroup& sX) const {
	auto fourierCoeff =
		std::vector(nfreq, std::vector(nfreq, std::vector<std::complex<double>>(
												  nfreq, {0.0, 0.0})));

	int np = sX.size(ParticleKind::Positive);
	int nn = sX.size(ParticleKind::Negative);

	auto& sp = sX.list(ParticleKind::Positive);
	auto& sn = sX.list(ParticleKind::Negative);

	double cubic2Pi = 8.0 * pi * pi * pi;

	double coeffFft = 1. / cubic2Pi;

	// for the (i,j,k)-th element of the array with size (nx,Ny,Nz), use the
	// expression an_array[k + Nz * (j + Ny * i)].

	// create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
	double dx = 2.0 * pi / augFactor / nfreq;
	double dy = 2.0 * pi / augFactor / nfreq;
	double dz = 2.0 * pi / augFactor / nfreq;

	const auto n = augFactor * nfreq;
	Vector3D f = std::vector(n, std::vector(n, std::vector<double>(n, 0.0)));
	auto fx = f;
	auto fy = f;
	auto fz = f;
	auto fxx = f;
	auto fyy = f;
	auto fzz = f;
	auto fxy = f;
	auto fxz = f;
	auto fyz = f;

	auto sizeF = n * n * n;

	for (int kp = 0; kp < np; kp++) {
		double x0 = sp[kp].velocity(0);
		double y0 = sp[kp].velocity(1);
		double z0 = sp[kp].velocity(2);
		int xloc = (int)(floor(x0 / dx + 0.5));
		int yloc = (int)(floor(y0 / dy + 0.5));
		int zloc = (int)(floor(z0 / dz + 0.5));
		if (xloc >= n)
			xloc--;
		if (yloc >= n)
			yloc--;
		if (zloc >= n)
			zloc--;
		double xdelta = x0 - xloc * dx;
		double ydelta = y0 - yloc * dy;
		double zdelta = z0 - zloc * dz;

		size_t loc = zloc + n * (yloc + n * xloc);
		if ((loc >= sizeF) || (loc < 0)) {
			std::string errMsg =
				"error: in approximation. Particle moved out of range. x =  (" +
				to_string(x0) + ", " + to_string(y0) + ", " + to_string(z0) +
				"), dx = " + to_string(dx);
			throw std::exception(errMsg.c_str());
		}

		f[xloc][yloc][zloc]++;
		fx[xloc][yloc][zloc] += xdelta;
		fy[xloc][yloc][zloc] += ydelta;
		fz[xloc][yloc][zloc] += zdelta;
		fxx[xloc][yloc][zloc] += xdelta * xdelta;
		fyy[xloc][yloc][zloc] += ydelta * ydelta;
		fzz[xloc][yloc][zloc] += zdelta * zdelta;
		fxy[xloc][yloc][zloc] += xdelta * ydelta;
		fyz[xloc][yloc][zloc] += ydelta * zdelta;
		fxz[xloc][yloc][zloc] += zdelta * xdelta;
	}

	// cout << "Approx 2" << endl;

	for (int kp = 0; kp < nn; kp++) {
		double x0 = sn[kp].velocity(0);
		double y0 = sn[kp].velocity(1);
		double z0 = sn[kp].velocity(2);
		int xloc = (int)(floor(x0 / dx + 0.5));
		int yloc = (int)(floor(y0 / dy + 0.5));
		int zloc = (int)(floor(z0 / dz + 0.5));
		if (xloc >= n)
			xloc--;
		if (yloc >= n)
			yloc--;
		if (zloc >= n)
			zloc--;
		double xdelta = x0 - xloc * dx;
		double ydelta = y0 - yloc * dy;
		double zdelta = z0 - zloc * dz;

		size_t loc = zloc + n * (yloc + n * xloc);
		if ((loc >= sizeF) || (loc < 0)) {
			std::string errMsg =
				"error: in approximation. Particle moved out of range. x =  (" +
				to_string(x0) + ", " + to_string(y0) + ", " + to_string(z0) +
				"), dx = " + to_string(dx);
			throw std::exception(errMsg.c_str());
		}

		f[xloc][yloc][zloc]--;
		fx[xloc][yloc][zloc] -= xdelta;
		fy[xloc][yloc][zloc] -= ydelta;
		fz[xloc][yloc][zloc] -= zdelta;
		fxx[xloc][yloc][zloc] -= xdelta * xdelta;
		fyy[xloc][yloc][zloc] -= ydelta * ydelta;
		fzz[xloc][yloc][zloc] -= zdelta * zdelta;
		fxy[xloc][yloc][zloc] -= xdelta * ydelta;
		fyz[xloc][yloc][zloc] -= ydelta * zdelta;
		fxz[xloc][yloc][zloc] -= zdelta * xdelta;
	}

	const auto orders = std::vector<std::vector<int>>{
		{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
		{0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

	const auto fDerivs = std::vector{
		std::move(f),	std::move(fx),	std::move(fy),	std::move(fz),
		std::move(fxx), std::move(fyy), std::move(fzz), std::move(fxy),
		std::move(fxz), std::move(fyz),
	};
	for (size_t korder = 0; korder < orders.size(); ++korder) {
		const auto order = orders[korder];
		const auto fourierCoeffOneTerm =
			fft3DApproxOneterm(fDerivs[korder], order[0], order[1], order[2]);
		for (size_t kx = 0; kx < nfreq; kx++) {
			for (size_t ky = 0; ky < nfreq; ky++) {
				for (size_t kz = 0; kz < nfreq; kz++) {
					fourierCoeff[kx][ky][kz] += fourierCoeffOneTerm[kx][ky][kz];
				}
			}
		}
	}

	for (size_t kx = 0; kx < nfreq; kx++) {
		for (size_t ky = 0; ky < nfreq; ky++) {
			for (size_t kz = 0; kz < nfreq; kz++) {
				fourierCoeff[kx][ky][kz] *= coeffFft;
			}
		}
	}

	return fourierCoeff;
}

} // namespace coulomb::resampling
