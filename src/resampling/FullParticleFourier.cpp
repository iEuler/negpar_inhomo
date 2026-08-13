#include "FullParticleFourier.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>

#include "Constants.h"
#include "FFT.h"
#include "ParticleGroup.h"
#include "ResamplingNumerics.h"

namespace coulomb {
using std::abs;
using std::complex;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::min;
using std::pow;
using std::sqrt;
using std::vector;

namespace {

/* ======================================================== *\
		Use Fourier transform for 3D interpolation
\* ======================================================== */

// xyzMinMax = [xmin, xmax, ymin, ymax, zmin, zmax];

void interp3dFftApproxTerms(std::vector<std::complex<double>>& fourierCoeff,
							vector<double>& f, int nfreq1, int nfreq2,
							int nfreq3, int augFactor, int orderx, int ordery,
							int orderz) {
	int sizeF = augFactor * augFactor * augFactor * nfreq1 * nfreq2 * nfreq3;

	// 1i *freq
	const auto freq1 = resampling::ResamplingNumerics{}.frequencies(nfreq1);
	const auto freq2 = resampling::ResamplingNumerics{}.frequencies(nfreq2);
	const auto freq3 = resampling::ResamplingNumerics{}.frequencies(nfreq3);

	const auto loc1 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq1, augFactor);
	const auto loc2 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq2, augFactor);
	const auto loc3 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq3, augFactor);

	FftwBuffer fin(static_cast<fftw_complex*>(
		fftw_malloc(static_cast<std::size_t>(sizeF) * sizeof(fftw_complex))));
	FftwBuffer fSaug(static_cast<fftw_complex*>(
		fftw_malloc(static_cast<std::size_t>(sizeF) * sizeof(fftw_complex))));
	if (!fin || !fSaug)
		throw std::runtime_error(
			"Unable to allocate interpolated FFTW buffers");

	FftwPlan plan3dFt(fftw_plan_dft_3d(
		augFactor * nfreq1, augFactor * nfreq2, augFactor * nfreq3, fin.get(),
		fSaug.get(), FFTW_FORWARD, FFTW_ESTIMATE));
	if (!plan3dFt)
		throw std::runtime_error("Unable to create interpolated FFTW plan");

	for (int kfc = 0; kfc < sizeF; kfc++) {
		fin[kfc][0] = f[kfc];
		fin[kfc][1] = 0.;
	}

	fftw_execute(plan3dFt.get());

	for (int kk1 = 0; kk1 < nfreq1; kk1++) {
		for (int kk2 = 0; kk2 < nfreq2; kk2++) {
			for (int kk3 = 0; kk3 < nfreq3; kk3++) {
				int kk = kk3 + nfreq3 * (kk2 + nfreq2 * kk1);

				const std::size_t kk1Aug = loc1[kk1];
				const std::size_t kk2Aug = loc2[kk2];
				const std::size_t kk3Aug = loc3[kk3];
				const std::size_t kkaug =
					kk3Aug + static_cast<std::size_t>(augFactor) * nfreq3 *
								 (kk2Aug + static_cast<std::size_t>(augFactor) *
											   nfreq2 * kk1Aug);

				double freq = 1.0;
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

				if ((orderx + ordery + orderz) == 1) {
					fourierCoeff[kk] +=
						freq *
						complex<double>(fSaug[kkaug][1], -fSaug[kkaug][0]);
				} else {
					fourierCoeff[kk] += freq * complex<double>(fSaug[kkaug][0],
															   fSaug[kkaug][1]);
				}
			}
		}
	}
}

} // namespace

std::vector<std::complex<double>>
FullParticleFourier::approximateTransform(NeParticleGroup& sX, int nfreq1,
										  int nfreq2, int nfreq3) {
	std::vector<std::complex<double>> fourierCoeff(nfreq1 * nfreq2 * nfreq3,
												   {0., 0.});
	int augFactor = 2;

	int np = sX.size(ParticleKind::Positive);
	int nn = sX.size(ParticleKind::Negative);

	auto& sp = sX.list(ParticleKind::Positive);
	auto& sn = sX.list(ParticleKind::Negative);

	double cubic2Pi = 8.0 * pi * pi * pi;

	double coeffFft = 1. / cubic2Pi;

	// for the (i,j,k)-th element of the array with size (nx,Ny,Nz), use the
	// expression an_array[k + Nz * (j + Ny * i)].

	// create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
	int sizeF = augFactor * augFactor * augFactor * nfreq1 * nfreq2 * nfreq3;
	double dx = 2.0 * pi / augFactor / nfreq1;
	double dy = 2.0 * pi / augFactor / nfreq2;
	double dz = 2.0 * pi / augFactor / nfreq3;

	vector<double> f(sizeF);
	vector<double> fx(sizeF);
	vector<double> fy(sizeF);
	vector<double> fz(sizeF);
	vector<double> fxx(sizeF);
	vector<double> fyy(sizeF);
	vector<double> fzz(sizeF);
	vector<double> fxy(sizeF);
	vector<double> fxz(sizeF);
	vector<double> fyz(sizeF);

	// cout << "Approx 1" << endl;

	for (int kk = 0; kk < sizeF; kk++) {
		f[kk] = 0.;
		fx[kk] = 0.;
		fy[kk] = 0.;
		fz[kk] = 0.;
		fxx[kk] = 0.;
		fyy[kk] = 0.;
		fzz[kk] = 0.;
		fxy[kk] = 0.;
		fxz[kk] = 0.;
		fyz[kk] = 0.;
	}

	for (int kp = 0; kp < np; kp++) {
		double x0 = sp[kp].velocity(0);
		double y0 = sp[kp].velocity(1);
		double z0 = sp[kp].velocity(2);
		int xloc = (int)(floor(x0 / dx + 0.5));
		int yloc = (int)(floor(y0 / dy + 0.5));
		int zloc = (int)(floor(z0 / dz + 0.5));
		if (xloc >= augFactor * nfreq1)
			xloc--;
		if (yloc >= augFactor * nfreq2)
			yloc--;
		if (zloc >= augFactor * nfreq3)
			zloc--;
		double xdelta = x0 - xloc * dx;
		double ydelta = y0 - yloc * dy;
		double zdelta = z0 - zloc * dz;

		int loc =
			zloc + augFactor * nfreq3 * (yloc + augFactor * nfreq2 * xloc);

		if ((loc >= sizeF) || (loc < 0))
			cout << x0 << ' ' << y0 << ' ' << z0 << ' ' << dx << ' ' << loc
				 << endl;

		f[loc]++;
		fx[loc] += xdelta;
		fy[loc] += ydelta;
		fz[loc] += zdelta;
		fxx[loc] += xdelta * xdelta;
		fyy[loc] += ydelta * ydelta;
		fzz[loc] += zdelta * zdelta;
		fxy[loc] += xdelta * ydelta;
		fyz[loc] += ydelta * zdelta;
		fxz[loc] += zdelta * xdelta;
	}

	// cout << "Approx 2" << endl;

	for (int kp = 0; kp < nn; kp++) {
		double x0 = sn[kp].velocity(0);
		double y0 = sn[kp].velocity(1);
		double z0 = sn[kp].velocity(2);
		// int xloc = floor(x0/dx);
		// int yloc = floor(y0/dy);
		// int zloc = floor(z0/dz);
		int xloc = (int)(floor(x0 / dx + 0.5));
		int yloc = (int)(floor(y0 / dy + 0.5));
		int zloc = (int)(floor(z0 / dz + 0.5));
		if (xloc >= augFactor * nfreq1)
			xloc--;
		if (yloc >= augFactor * nfreq2)
			yloc--;
		if (zloc >= augFactor * nfreq3)
			zloc--;
		double xdelta = x0 - xloc * dx;
		double ydelta = y0 - yloc * dy;
		double zdelta = z0 - zloc * dz;

		int loc =
			zloc + augFactor * nfreq3 * (yloc + augFactor * nfreq2 * xloc);

		if ((loc >= sizeF) || (loc < 0)) {
			cout
				<< "error: in approximation. Particle moved out of range. kx = "
				<< xloc << ' ' << yloc << ' ' << zloc << ' ' << loc << endl;
			exit(0);
		}

		f[loc]--;
		fx[loc] -= xdelta;
		fy[loc] -= ydelta;
		fz[loc] -= zdelta;
		fxx[loc] -= xdelta * xdelta;
		fyy[loc] -= ydelta * ydelta;
		fzz[loc] -= zdelta * zdelta;
		fxy[loc] -= xdelta * ydelta;
		fyz[loc] -= ydelta * zdelta;
		fxz[loc] -= zdelta * xdelta;
	}
	// cout << "Approx 3" << endl;

	interp3dFftApproxTerms(fourierCoeff, f, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 0, 0);

	interp3dFftApproxTerms(fourierCoeff, fx, nfreq1, nfreq2, nfreq3, augFactor,
						   1, 0, 0);
	interp3dFftApproxTerms(fourierCoeff, fy, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 1, 0);
	interp3dFftApproxTerms(fourierCoeff, fz, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 0, 1);
	interp3dFftApproxTerms(fourierCoeff, fxx, nfreq1, nfreq2, nfreq3, augFactor,
						   2, 0, 0);
	interp3dFftApproxTerms(fourierCoeff, fyy, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 2, 0);
	interp3dFftApproxTerms(fourierCoeff, fzz, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 0, 2);
	interp3dFftApproxTerms(fourierCoeff, fxy, nfreq1, nfreq2, nfreq3, augFactor,
						   1, 1, 0);
	interp3dFftApproxTerms(fourierCoeff, fxz, nfreq1, nfreq2, nfreq3, augFactor,
						   1, 0, 1);
	interp3dFftApproxTerms(fourierCoeff, fyz, nfreq1, nfreq2, nfreq3, augFactor,
						   0, 1, 1);

	for (int kk = 0; kk < nfreq1 * nfreq2 * nfreq3; kk++)
		fourierCoeff[kk] *= coeffFft;

	return fourierCoeff;
	// cout << "Approx finished." << endl;
}

void FullParticleFourier::filter(
	std::vector<std::complex<double>>& /*Fouriercoeff*/,
	vector<int>& flagFouriercoeff, int sizeFc) {
	// double thres = 10.0;
	for (int k = 0; k < sizeFc; k++) {
		flagFouriercoeff[k] = 1;
		/*
		double abs_FC = abs(Fouriercoeff[k]);
		if (abs_FC < thres) {
		  Fouriercoeff[k] *= 0.;
		  flag_Fouriercoeff[k] = 0;
		}
		*/
	}
}

/******************************************************************/
/* ---------- Use Fourier transform for 3D interpolation -------- */
/******************************************************************/

/*
  Find the coarse approximation with the given Fourier coefficients
  Need to include 'fftw3.f'
*/
vector<double> FullParticleFourier::interpolateCoarse(
	const std::vector<std::complex<double>>& fourierCoeff, int nfreq1,
	int nfreq2, int nfreq3) {
	// double Lcubic = (double) (Nfreq1*Nfreq2*Nfreq3);

	vector<double> fcoarse(nfreq1 * nfreq2 * nfreq3);

	const auto elementCount =
		static_cast<std::size_t>(nfreq1) * nfreq2 * nfreq3;
	FftwBuffer fc(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	FftwBuffer fcoarseC(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	if (!fc || !fcoarseC)
		throw std::runtime_error("Unable to allocate coarse FFTW buffers");

	for (int kfc = 0; kfc < nfreq1 * nfreq2 * nfreq3; kfc++) {
		fc[kfc][0] = fourierCoeff[kfc].real();
		fc[kfc][1] = fourierCoeff[kfc].imag();
	}

	// use fftw to obtain an estimation of f_p - f_n
	// // cout << " resample 1.3" << endl;
	FftwPlan plan3dIft(fftw_plan_dft_3d(nfreq1, nfreq2, nfreq3, fc.get(),
										fcoarseC.get(), FFTW_BACKWARD,
										FFTW_ESTIMATE));
	if (!plan3dIft)
		throw std::runtime_error("Unable to create coarse FFTW plan");

	fftw_execute(plan3dIft.get());

	for (int kfc = 0; kfc < nfreq1 * nfreq2 * nfreq3; kfc++) {
		// fcoarse[kfc] = fcoarse_c[kfc][0] / Lcubic;
		fcoarse[kfc] = fcoarseC[kfc][0];
	}
	return fcoarse;
}

/** The derivatives of f on grids
 */
// CONTINUE HERE

namespace {

vector<double>
interp3dFxyzTerms(const std::vector<std::complex<double>>& fourierCoeff,
				  int nfreq1, int nfreq2, int nfreq3, int augFactor, int orderx,
				  int ordery, int orderz) {
	int sizeF = augFactor * augFactor * augFactor * nfreq1 * nfreq2 * nfreq3;
	std::vector<complex<double>> fSaug(sizeF, {0., 0.});

	// 1i *freq
	const auto freq1 = resampling::ResamplingNumerics{}.frequencies(nfreq1);
	const auto freq2 = resampling::ResamplingNumerics{}.frequencies(nfreq2);
	const auto freq3 = resampling::ResamplingNumerics{}.frequencies(nfreq3);

	const auto loc1 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq1, augFactor);
	const auto loc2 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq2, augFactor);
	const auto loc3 =
		resampling::ResamplingNumerics{}.augmentedLocations(nfreq3, augFactor);

	for (int kk1 = 0; kk1 < nfreq1; kk1++) {
		for (int kk2 = 0; kk2 < nfreq2; kk2++) {
			for (int kk3 = 0; kk3 < nfreq3; kk3++) {
				int kk = kk3 + nfreq3 * (kk2 + nfreq2 * kk1);

				const std::size_t kk1Aug = loc1[kk1];
				const std::size_t kk2Aug = loc2[kk2];
				const std::size_t kk3Aug = loc3[kk3];
				const std::size_t kkaug =
					kk3Aug + static_cast<std::size_t>(augFactor) * nfreq3 *
								 (kk2Aug + static_cast<std::size_t>(augFactor) *
											   nfreq2 * kk1Aug);

				double freq = 1.0;
				for (int kx = 0; kx < orderx; kx++)
					freq *= freq1[kk1];
				for (int kx = 0; kx < ordery; kx++)
					freq *= freq2[kk2];
				for (int kx = 0; kx < orderz; kx++)
					freq *= freq3[kk3];

				fSaug[kkaug] = freq * fourierCoeff[kk];
			}
		}
	}

	return FullParticleFourier{}.interpolateCoarse(
		fSaug, augFactor * nfreq1, augFactor * nfreq2, augFactor * nfreq3);
}

} // namespace

std::vector<std::vector<double>> FullParticleFourier::interpolateDerivatives(
	const std::vector<std::complex<double>>& fourierCoeff, int nfreq1,
	int nfreq2, int nfreq3, int augFactor) {
	const auto f = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									 augFactor, 0, 0, 0);
	const auto fx = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									  augFactor, 1, 0, 0);
	const auto fy = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									  augFactor, 0, 1, 0);
	const auto fz = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									  augFactor, 0, 0, 1);
	const auto fxx = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 2, 0, 0);
	const auto fyy = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 0, 2, 0);
	const auto fzz = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 0, 0, 2);
	const auto fxy = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 1, 1, 0);
	const auto fxz = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 1, 0, 1);
	const auto fyz = interp3dFxyzTerms(fourierCoeff, nfreq1, nfreq2, nfreq3,
									   augFactor, 0, 1, 1);
	return {f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz};
}

// void interp3d_fft_eachlevel(NeParticleGroup * S_x, MultlLevelGroup * MLsol,
// int Nlevel); void interp3d_fft_ml(complex<double> *Fouriercoeff, int *
// flag_Fouriercoeff, MultlLevelGroup * MLsol, int Nlevel);
/*
  Sampled particles are stored in particles_Sp_sampled and particles_Sn_sampled
  with size particles%Np_sampled and particles%Nn_sampled
*/

std::vector<double>
FullParticleFourier::valuesAt(const std::vector<std::vector<double>>& fvecs,
							  int k) {
	std::vector<double> result;
	result.reserve(fvecs.size());

	for (const auto& fvec : fvecs) {
		result.push_back(fvec[k]);
	}
	return result;
}

vector<double> FullParticleFourier::upperBound(int n,
											   const vector<double>& fc) {
	vector<double> fUp(n * n * n);
	// Find f_up>abs(fc)
	double fAll[8];
	int kk[8];

	for (int kx = 0; kx < n; kx++) {
		int xr = kx + 1;
		if (kx == n - 1)
			xr = 0;

		for (int ky = 0; ky < n; ky++) {
			int yr = ky + 1;
			if (ky == n - 1)
				yr = 0;

			for (int kz = 0; kz < n; kz++) {
				int zr = kz + 1;
				if (kz == n - 1)
					zr = 0;

				kk[0] = kz + n * (ky + n * kx);
				kk[1] = zr + n * (ky + n * kx);
				kk[2] = kz + n * (yr + n * kx);
				kk[3] = zr + n * (yr + n * kx);
				kk[4] = kz + n * (ky + n * xr);
				kk[5] = zr + n * (ky + n * xr);
				kk[6] = kz + n * (yr + n * xr);
				kk[7] = zr + n * (yr + n * xr);

				double maxFAll = 0.;
				for (int k = 0; k < 8; k++) {
					fAll[k] = abs(fc[kk[k]]);
					maxFAll = max(maxFAll, fAll[k]);
				}

				fUp[kk[0]] = maxFAll;
			}
		}
	}
	return fUp;
}

namespace {

void addMaxwellianTerms(double rhoM, vector<double> uM, vector<double> tm,
						double neff, vector<double>& f, int nfreq,
						int augFactor, int orderx, int ordery, int orderz) {
	double mccCoe = rhoM / sqrt(8.0 * pi * pi * pi * tm[0] * tm[1] * tm[2]);

	vector<double> interpXaug(nfreq * augFactor);
	vector<double> expX(nfreq * augFactor);
	vector<double> expY(nfreq * augFactor);
	vector<double> expZ(nfreq * augFactor);
	for (int kx = 0; kx < nfreq * augFactor; kx++) {
		double xk = kx * 2 * pi / nfreq / augFactor;
		interpXaug[kx] = xk;
		expX[kx] = exp(-(xk - uM[0]) * (xk - uM[0]) / 2 / tm[0]);
		expY[kx] = exp(-(xk - uM[1]) * (xk - uM[1]) / 2 / tm[1]);
		expZ[kx] = exp(-(xk - uM[2]) * (xk - uM[2]) / 2 / tm[2]);
	}

	for (int kx = 0; kx < augFactor * nfreq; kx++) {
		for (int ky = 0; ky < augFactor * nfreq; ky++) {
			for (int kz = 0; kz < augFactor * nfreq; kz++) {
				int kk = kz + augFactor * nfreq * (ky + augFactor * nfreq * kx);

				double mcc = mccCoe * expX[kx] * expY[ky] * expZ[kz];
				double xc = interpXaug[kx];
				double yc = interpXaug[ky];
				double zc = interpXaug[kz];

				if (orderx == 1)
					mcc *= -(xc - uM[0]) / tm[0];
				if (orderx == 2)
					mcc *=
						((xc - uM[0]) * (xc - uM[0]) - tm[0]) / tm[0] / tm[0];
				if (ordery == 1)
					mcc *= -(yc - uM[1]) / tm[1];
				if (ordery == 2)
					mcc *=
						((yc - uM[1]) * (yc - uM[1]) - tm[1]) / tm[1] / tm[1];
				if (orderz == 1)
					mcc *= -(zc - uM[2]) / tm[2];
				if (orderz == 2)
					mcc *=
						((zc - uM[2]) * (zc - uM[2]) - tm[2]) / tm[2] / tm[2];

				f[kk] = neff * f[kk] + mcc;
			}
		}
	}
}

} // namespace

void FullParticleFourier::addMaxwellian(
	double rhoM, vector<double> uM, vector<double> tm, double neff,
	std::vector<vector<double>>& fDerivatives, int nfreq, int augFactor) {
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[0], nfreq, augFactor, 0,
					   0, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[1], nfreq, augFactor, 1,
					   0, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[2], nfreq, augFactor, 0,
					   1, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[3], nfreq, augFactor, 0,
					   0, 1);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[4], nfreq, augFactor, 2,
					   0, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[5], nfreq, augFactor, 0,
					   2, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[6], nfreq, augFactor, 0,
					   0, 2);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[7], nfreq, augFactor, 1,
					   1, 0);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[8], nfreq, augFactor, 1,
					   0, 1);
	addMaxwellianTerms(rhoM, uM, tm, neff, fDerivatives[9], nfreq, augFactor, 0,
					   1, 1);
}

} // namespace coulomb
