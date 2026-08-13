#include "ResamplerHelper.h"

#include <cmath>
#include <complex>
#include <vector>

#include "Constants.h"
#include "RandomSampling.h"

namespace coulomb::resampling {

using std::abs;
using std::complex;
using std::exp;
using std::sqrt;
using std::vector;

VectorBool3D
ResamplerHelper::filterFourierCoefficients(VectorComplex3D& fourierCoeff) {
	// double thres = 10.0;
	const auto n1 = fourierCoeff.size();
	const auto n2 = fourierCoeff.front().size();
	const auto n3 = fourierCoeff.front().front().size();
	auto flagFouriercoeff =
		std::vector(n1, std::vector(n2, std::vector<bool>(n3, true)));

	// for (size_t kk1 = 0; kk1 < n1; kk1++) {
	//   for (size_t kk2 = 0; kk2 < n2; kk2++) {
	//     for (size_t kk3 = 0; kk3 < n3; kk3++) {
	//       const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
	//       flag_Fouriercoeff[kk1][kk2][kk3] = true;
	//     }
	//   }
	// }
	// for (int k = 0; k < size_FC; k++) {
	//   /*
	//   double abs_FC = abs(Fouriercoeff[k]);
	//   if (abs_FC < thres) {
	//     Fouriercoeff[k] *= 0.;
	//     flag_Fouriercoeff[k] = 0;
	//   }
	//   */
	// }
	return flagFouriercoeff;
}

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
Vector3D ResamplerHelper::upperBound(const Vector3D& fc) {
	const auto n = fc.size();
	auto fUp = std::vector(n, std::vector(n, std::vector<double>(n, 0.0)));

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

				double maxFAll = std::abs(fc[kx][ky][kz]);
				maxFAll = std::max(maxFAll, std::abs(fc[kx][ky][zr]));
				maxFAll = std::max(maxFAll, std::abs(fc[kx][yr][kz]));
				maxFAll = std::max(maxFAll, std::abs(fc[xr][ky][kz]));
				maxFAll = std::max(maxFAll, std::abs(fc[kx][yr][zr]));
				maxFAll = std::max(maxFAll, std::abs(fc[xr][yr][kz]));
				maxFAll = std::max(maxFAll, std::abs(fc[xr][ky][zr]));
				maxFAll = std::max(maxFAll, std::abs(fc[xr][yr][zr]));

				fUp[kx][ky][kz] = maxFAll;
			}
		}
	}
	return fUp;
}

std::vector<double>
ResamplerHelper::valuesAt(const std::vector<Vector3D>& fvecs, int kx, int ky,
						  int kz) {
	std::vector<double> result;
	result.reserve(fvecs.size());

	for (const auto& fvec : fvecs) {
		result.push_back(fvec[kx][ky][kz]);
	}
	return result;
}

double
ResamplerHelper::valueFromFft(const std::vector<double>& sf,
							  const VectorComplex3D& fourierCoeff,
							  const std::vector<std::complex<double>>& ifreq1,
							  const std::vector<std::complex<double>>& ifreq2,
							  const std::vector<std::complex<double>>& ifreq3,
							  const VectorBool3D& flagFouriercoeff) {
	std::complex<double> fvalC(0., 0.);

	for (int kk1 = 0; kk1 < ifreq1.size(); kk1++) {
		for (int kk2 = 0; kk2 < ifreq2.size(); kk2++) {
			for (int kk3 = 0; kk3 < ifreq3.size(); kk3++) {
				if (flagFouriercoeff[kk1][kk2][kk3]) {
					fvalC += fourierCoeff[kk1][kk2][kk3] *
							 exp(ifreq1[kk1] * sf[0] + ifreq2[kk2] * sf[1] +
								 ifreq3[kk3] * sf[2]);
				}
			}
		}
	}

	// return real(fval_c)/(Nfreq1*Nfreq2*Nfreq3);
	return fvalC.real();
}

void ResamplerHelper::acceptSample(const std::vector<double>& sf,
								   NeParticleGroup& sXInCell, double fval,
								   double& maxF) {
	RandomSampling randomSampling(randomContext);
	if (abs(fval) > maxF) {
		// keep sampled particles with rate maxF/maxF_new

		double keepRate = maxF / (1.5 * abs(fval));

		maxF = 1.5 * abs(fval);

		const auto removeParticles = [&](ParticleKind kind) {
			const int count = randomSampling.stochasticFloor(
				(1 - keepRate) * sXInCell.size(kind));
			for (int particle = 0; particle < count; ++particle) {
				const int index = static_cast<int>(randomSampling.uniform() *
												   sXInCell.size(kind));
				sXInCell.erase(index, kind);
			}
		};
		removeParticles(ParticleKind::Positive);
		removeParticles(ParticleKind::Negative);
	}

	// accept this particle with rate abs(fval/maxF)
	if (randomSampling.uniform() < (abs(fval / maxF))) {
		double sumSfPiSq = 0.;
		for (int kv = 0; kv < 3; kv++)
			sumSfPiSq += (sf[kv] - pi) * (sf[kv] - pi);
		if (sqrt(sumSfPiSq) < pi) {
			Particle1D3D sOne({sf[0], sf[1], sf[2]});
			const auto kind =
				fval > 0 ? ParticleKind::Positive : ParticleKind::Negative;
			sXInCell.pushBack(sOne, kind);
		}
	}
}

} // namespace coulomb::resampling
