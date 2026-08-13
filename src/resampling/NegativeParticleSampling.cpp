#include "NegativeParticleSampling.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include "Constants.h"
#include "Grid.h"
#include "Numerics.h"
#include "ParticleGroup.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {
using std::abs;
using std::cos;
using std::exp;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::sin;
using std::sqrt;
using std::vector;

/* ======================================================================== *\
		Sample from the source term: the change in Maxwellian due to collisions
		with P and N particles
\* ======================================================================== */

/**
  Evaluate M(v)
  evaluate M(v0) where the moments of Maxwellian is given in moment
*/

double
NegativeParticleSampling::evaluateMaxwellian(const std::vector<double>& v0,
											 const NeParticleGroup& sX) {
	double rho = sX.rhoM;
	double u1 = sX.u1M;
	double u2 = sX.u2M;
	double u3 = sX.u3M;
	double tprt = sX.tprtM;

	double usq = (v0[0] - u1) * (v0[0] - u1) + (v0[1] - u2) * (v0[1] - u2) +
				 (v0[2] - u3) * (v0[2] - u3);
	return rho / pow(sqrt(2.0 * pi * tprt), 3) * exp(-usq / 2.0 / tprt);
}

// ========================================================================

/**
  Evaluate h(v;v1)
  evaluate h(v0;v1) = h(v0), with a source particle at v1
  mode = 0, calculate
  mode = 1, calculate c0
  mode = 2, calculate c1
  mode = 3, calculate c2
*/

double NegativeParticleSampling::evaluateSource(const std::vector<double>& v0,
												const std::vector<double>& v1,
												const NeParticleGroup& sX,
												const ParaClass& para,
												int mode) {
	// double rho = S_x.rhoM;
	double u1 = sX.u1M;
	double u2 = sX.u2M;
	double u3 = sX.u3M;
	double tprt = sX.tprtM;

	double h = 0;

	if (para.methodBinaryColl == BinaryCollisionMethod::TA) {
		double u[3], sqrtU = 0;
		for (int k = 0; k < 3; k++) {
			u[k] = v0[k] - v1[k];
			sqrtU += u[k] * u[k];
		}

		sqrtU = sqrt(sqrtU);
		double sigma2Delta = para.coeffBinaryColl * para.dt / pow(sqrtU, 3);

		double r = 0.5 / sigma2Delta;
		double logr = log(r);
		double maxdelta;
		if (r > 1.0) {
			maxdelta = exp(-0.0026 * logr * logr - 0.4150 * logr + 0.9154);
		} else {
			maxdelta = exp(-0.0006 * logr * logr - 0.2242 * logr + 0.4711);
		}

		double v0Perp[3], v0PerpSq = 0, v0PerpDotU = 0;
		v0Perp[0] = v0[0] - u1;
		v0Perp[1] = v0[1] - u2;
		v0Perp[2] = v0[2] - u3;
		for (int k = 0; k < 3; k++)
			v0PerpDotU += v0Perp[k] * u[k];
		for (int k = 0; k < 3; k++)
			v0Perp[k] -= v0PerpDotU * u[k] / (sqrtU * sqrtU);
		for (int k = 0; k < 3; k++)
			v0PerpSq += v0Perp[k] * v0Perp[k];

		double vperp2 = v0PerpSq / (2.0 * tprt);

		int ndelta = 16;
		vector<double> deltaAll(ndelta);
		for (int j = 0; j < ndelta; j++)
			deltaAll[j] = ((double)j) / ndelta * maxdelta;

		vector<double> coeffSum(ndelta);
		for (int j = 0; j < ndelta; j++) {
			double zeta2p11p5 = pow(sqrt(1.0 + deltaAll[j] * deltaAll[j]), 3);
			coeffSum[j] = sqrt(r / pi) * pow(sqrt(zeta2p11p5), 3) *
						  exp(-r * deltaAll[j] * deltaAll[j] * zeta2p11p5);
		}

		double dxDelta = deltaAll[1] - deltaAll[0];

		coeffSum[0] = coeffSum[0] / 2.0;
		coeffSum[ndelta - 1] = coeffSum[ndelta - 1] / 2.0;

		vector<double> eps2(ndelta);
		for (int j = 0; j < ndelta; j++)
			eps2[j] = sqrtU * sqrtU * deltaAll[j] * deltaAll[j] / 2.0 / tprt;

		double m = NegativeParticleSampling{}.evaluateMaxwellian(v0, sX);
		double hM = 0;

		for (int j = 0; j < ndelta; j++)
			hM += exp(-eps2[j]) *
				  (1.0 + (eps2[j] * vperp2) +
				   .25 * (eps2[j] * vperp2 * eps2[j] * vperp2)) *
				  coeffSum[j];

		hM = hM * dxDelta * 2.0;

		h = m * hM;

		double hh = 0;
		if (mode == 1) { // this h gives c0
			for (int j = 0; j < ndelta; j++)
				hh += exp(-eps2[j]) * coeffSum[j];
			h = sqrtU * sqrtU * (hh * dxDelta * 2.0 - 1.0);
		} else if (mode == 2) { // this h gives c1
			for (int j = 0; j < ndelta; j++)
				hh += exp(-eps2[j]) * eps2[j] * coeffSum[j];
			h = sqrtU * sqrtU * hh * dxDelta * 2.0;
		} else if (mode == 3) { // this h gives c2
			for (int j = 0; j < ndelta; j++)
				hh += exp(-eps2[j]) * eps2[j] * eps2[j] * coeffSum[j];
			h = sqrtU * sqrtU * hh * dxDelta * 2.0;
		}

	} else {
		h = 1.0;
	}

	return h;
}

/**
  Find the lower/upper bound of delta m (v;v1)
  for the delta source at v1, find lower/upper bound of delta m(v;v1) in the
  following form \delta m_n(v;v1) > - alphaNeg * M(v) |v-v1|^2 delta m_p(v;v1)
  < alphaPos * rho_m \detla m = \delta m_p - \delta m_n, with \delta m_p(m) = 0
  if |\delta m_p(m)|<(>) alphaNeg * M(v)
  |\delta m_n(v;v1)| < alphaNeg * M(v)
  |v-v1|^2 delta m_p(v;v1) < alphaPos * rho_m
*/

void NegativeParticleSampling::updateBounds(NeParticleGroup& sX,
											const ParaClass& para) {
	double tprt = sX.tprtM;
	std::vector<double> v1{sqrt(tprt), 0, 0};

	int nepsIn = 20;
	vector<double> epsAll(nepsIn + 1);
	for (int k = 0; k <= nepsIn; k++)
		epsAll[k] = 0.1 / (1 << k);

	double vrange = 3.0 * sqrt(tprt);

	int nepsOut = 40;
	int lengthVAll = 2 * nepsIn + 2 * nepsOut;

	vector<double> vAll1(lengthVAll);

	double dv1 = (v1[0] - epsAll[0] + vrange) / nepsOut;
	double dv2 = (vrange - v1[0] - epsAll[0]) / nepsOut;

	for (int k = 0; k < nepsOut; k++)
		vAll1[k] = -vrange + (k + 1) * dv1;
	for (int k = nepsOut; k < (nepsOut + nepsIn); k++)
		vAll1[k] = v1[0] - epsAll[k - nepsOut + 1];
	for (int k = nepsOut + nepsIn; k < (nepsOut + 2 * nepsIn); k++)
		vAll1[k] = v1[0] + epsAll[k - nepsOut - nepsIn + 1];
	for (int k = nepsOut + 2 * nepsIn; k < lengthVAll; k++)
		vAll1[k] = vrange - (k - nepsOut - 2 * nepsIn + 1) * dv2;

	vector<double> hh(lengthVAll);
	vector<double> mm(lengthVAll);

	// Look for lower bound

	std::vector<double> v0{0, 0, 0};

	for (int kv = 0; kv < lengthVAll; kv++) {
		v0[0] = vAll1[kv];
		hh[kv] = NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para);
		mm[kv] = NegativeParticleSampling{}.evaluateMaxwellian(v0, sX);
	}

	// MacroOutput{}.saveMacro<double>(hh, "hh");

	// delta m = h - m
	vector<double> hhMM(lengthVAll);
	for (int kv = 0; kv < lengthVAll; kv++)
		hhMM[kv] = hh[kv] / mm[kv] - 1;
	double alphaNeg = -*std::min_element(hhMM.begin(), hhMM.end());

	// Look for upper bound

	vector<double> hh0(lengthVAll);
	vector<double> hh1(lengthVAll);
	vector<double> hh2(lengthVAll);
	for (int kv = 0; kv < lengthVAll; kv++) {
		v0[0] = vAll1[kv];

		hh0[kv] =
			NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para, 1);
		hh1[kv] =
			NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para, 2);
		hh2[kv] =
			NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para, 3);
	}

	double beta = 3.0;

	double alphaPos = *std::max_element(hh0.begin(), hh0.end()) +
					  *std::max_element(hh1.begin(), hh1.end()) * beta * beta +
					  *std::max_element(hh2.begin(), hh2.end()) * pow(beta, 4);

	sX.alphaNeg = alphaNeg;
	sX.alphaPos = alphaPos;
	sX.rMax = 6.0 * sqrt(2 * tprt);

	// cout << "alpha = ( " << alphaNeg << ", " << alphaPos << ", " << S_x .
	// rMax << " )" << endl;
}

void NegativeParticleSampling::updateBounds(std::vector<NeParticleGroup>& sX,
											const NumericGridClass& grid,
											const ParaClass& para) {
	double minTprt = sX.front().tprtM;
	int kxMinTprt = 0;
	for (int kx = 1; kx < grid.nx; kx++) {
		if (minTprt > sX[kx].tprtM) {
			minTprt = sX[kx].tprtM;
			kxMinTprt = kx;
		}
	}

	NegativeParticleSampling::updateBounds(sX[kxMinTprt], para);

	double alphaNeg = sX[kxMinTprt].alphaNeg;
	double alphaPos = sX[kxMinTprt].alphaPos;

	for (int kx = 0; kx < grid.nx; kx++) {
		sX[kx].alphaNeg = alphaNeg;
		sX[kx].alphaPos = alphaPos;
		sX[kx].rMax = 6.0 * sqrt(2 * (sX[kx].tprtM));
	}
}

// ========================================================================

struct NegativeParticleSample {
	std::vector<double> velocity;
	ParticleKind kind;
};

/** Sample one accepted particle from the negative part of Delta M. */
std::optional<NegativeParticleSample> samplefromhNeg(NeParticleGroup& sX,
													 const ParaClass& para,
													 double neff,
													 RandomContext& random) {
	std::vector<double> v0{0.0, 0.0, 0.0};

	double alphaNeg = sX.alphaNeg;

	double rhoF = sX.rho;
	double rhoP = sX.positiveMoments.m0 * neff;
	double rhon = sX.negativeMoments.m0 * neff;

	int np, nn;
	np = sX.size(ParticleKind::Positive);
	nn = sX.size(ParticleKind::Negative);

	int npickup = para.nPickupNeg;

	int nNp = min(npickup, np);
	int nNn = min(npickup, nn);

	const auto idp = RandomSampling(random).permutation(np, nNp);
	const auto idn = RandomSampling(random).permutation(nn, nNn);

	v0[0] = RandomSampling(random).normal() * sqrt(sX.tprtM) + sX.u1M;
	v0[1] = RandomSampling(random).normal() * sqrt(sX.tprtM) + sX.u2M;
	v0[2] = RandomSampling(random).normal() * sqrt(sX.tprtM) + sX.u3M;

	double m0 = NegativeParticleSampling{}.evaluateMaxwellian(v0, sX);

	double hp = 0, hn = 0;

	auto& sp = sX.list(ParticleKind::Positive);
	auto& sn = sX.list(ParticleKind::Negative);

	for (int kp = 0; kp < nNp; kp++) {
		auto& v1 = sp[idp[kp] - 1].velocity();
		double h0 =
			NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para) - m0;
		if (h0 < (alphaNeg * m0))
			hp += h0;
	}
	for (int kn = 0; kn < nNn; kn++) {
		auto& v1 = sn[idn[kn] - 1].velocity();
		double h0 =
			NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para) - m0;
		if (h0 < (alphaNeg * m0))
			hn += h0;
	}
	double h = hp * np / (nNp + 1.0e-15) - hn * nn / (nNn + 1.0e-15);
	h = h * neff / rhoF;
	double hbar = max(rhoP, rhon) / rhoF * m0 * alphaNeg;
	double r0 = RandomSampling(random).uniform();
	if (r0 < (abs(h) / hbar)) {
		return NegativeParticleSample{v0, h > 0 ? ParticleKind::Positive
												: ParticleKind::Negative};
	}
	return std::nullopt;
}

// ========================================================================

/**
  the number of virtual particles sampled from h_+
  determine the number of virtual particles sampled from \Delta m_+
  Max_p(j) is the upper bound for h due to source at Sp(idp(j))
*/

int NegativeParticleSampling::estimateVirtualCount(const NeParticleGroup& sX,
												   double neff,
												   RandomContext& random) {
	int np = sX.size(ParticleKind::Positive);
	int nn = sX.size(ParticleKind::Negative);

	double rhoM = sX.rhoM;
	double tprtM = sX.tprtM;

	double rho = rhoM + neff * (np - nn);
	return RandomSampling(random).stochasticFloor(
		4.0 * pi * sX.rMax * sX.alphaPos * rhoM /
		pow(sqrt(2.0 * pi * tprtM), 3) / rho * (np + nn));
}

// ========================================================================

/**
  Sample particles from \Delta M, in homogeneous case
*/

void NegativeParticleSampling::sampleDelta(NeParticleGroup& sX,
										   NeParticleGroup& sXNew,
										   const ParaClass& para, double neff,
										   RandomContext& random) {
	// Sample particles from Delta m
	double alphaNeg = sX.alphaNeg;
	double alphaPos = sX.alphaPos;
	double rMax = sX.rMax;

	int np = sX.size(ParticleKind::Positive);
	int nn = sX.size(ParticleKind::Negative);

	auto& sp = sX.list(ParticleKind::Positive);
	auto& sn = sX.list(ParticleKind::Negative);

	double rhoF = sX.rho;
	double rhoM = sX.rhoM;
	double tprtM = sX.tprtM;
	double maxm = rhoM / pow(sqrt(2.0 * pi * tprtM), 3);

	// Sample from negative part

	double nnegF = max(np, nn) * alphaNeg * rhoM / rhoF;
	int nneg = RandomSampling(random).stochasticFloor(
		nnegF); // Number of virtual particles

	// Particle1D3D * Sp_new = S_x_new . list('p');
	// Particle1D3D * Sn_new = S_x_new . list('n');

	Particle1D3D sOne;

	for (int kneg = 0; kneg < nneg; kneg++) {
		const auto sample = samplefromhNeg(sX, para, neff, random);
		if (sample) {
			sOne.setVelocity(sample->velocity);
			sXNew.pushBack(sOne, sample->kind);
		}
	}

	// cout << "negpart " << COUNT_MYRAND << endl;

	int npos = NegativeParticleSampling::estimateVirtualCount(sX, neff, random);

	// cout << "Npos = " << Npos << endl;

	double rateP = ((double)np) / (np + nn);

	// int kk_test = 418;
	// cout << Npos << ' ' << rate_P << ' ' << COUNT_MYRAND << endl;
	for (int kpos = 0; kpos < npos; kpos++) {
		double rrr = RandomSampling(random).uniform();

		/*
		if (FLAG_CHECK == 1) {
		  if (kpos == 52) {
			cout << COUNT_MYRAND << endl;
			cout << rrr << ' ' << rate_P << endl;
			// std::exit(0);
		  }
		}
					*/

		if (rrr < rateP) {
			// Sample positve particles
			// choose the source particle
			int kp = (int)(np * RandomSampling(random).uniform());
			auto& v1 = sp[kp].velocity();
			// cout << v1[0] << ' '<< v1[1] << ' '<< v1[2] << endl;
			// sample a particle from the nearby
			double r1 = RandomSampling(random).uniform() * rMax;
			double costheta = 2.0 * RandomSampling(random).uniform() - 1.0;
			double sintheta = sqrt(1.0 - costheta * costheta);
			double phi = RandomSampling(random).uniform() * pi * 2.0;
			std::vector<double> v0{v1[0] + r1 * sintheta * cos(phi),
								   v1[1] + r1 * sintheta * sin(phi),
								   v1[2] + r1 * costheta};

			double m0 = NegativeParticleSampling{}.evaluateMaxwellian(v0, sX);

			// if (kpos == kk_test) cout << "test " << COUNT_MYRAND <<' '<<
			// v0[0] << '
			// '<< v0[1] << ' '<< v0[2] << ' ' << M0 << endl;

			if (RandomSampling(random).uniform() < (m0 / maxm)) {
				double h0 =
					NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para);
				double hbar0 = h0 - m0 - alphaNeg * m0;
				if (hbar0 > 0) {
					// check v0 is in the pos zone
					double r2h0 = r1 * r1 * hbar0;
					double rr = RandomSampling(random).uniform();
					if (rr < (r2h0 / (alphaPos * m0))) {
						// accept the virtual particle v0 with suitable rate
						sOne.setVelocity(v0);
						sXNew.pushBack(sOne, ParticleKind::Positive);
						// cout << "pos " <<  COUNT_MYRAND << ' ' << kpos <<
						// endl; cout << v0[0] << ' '<< v0[1] << ' '<< v0[2] <<
						// endl;
					}
				}
			}
		} else {
			// Sample negative particles
			// choose the source particle
			int kn = (int)(nn * RandomSampling(random).uniform());
			auto& v1 = sn[kn].velocity();
			// sample a particle from the nearby
			double r1 = RandomSampling(random).uniform() * rMax;
			double costheta = 2.0 * RandomSampling(random).uniform() - 1.0;
			double sintheta = sqrt(1.0 - costheta * costheta);
			double phi = RandomSampling(random).uniform() * pi * 2.0;
			std::vector<double> v0{v1[0] + r1 * sintheta * cos(phi),
								   v1[1] + r1 * sintheta * sin(phi),
								   v1[2] + r1 * costheta};

			double m0 = NegativeParticleSampling{}.evaluateMaxwellian(v0, sX);

			rrr = RandomSampling(random).uniform();
			if (rrr < (m0 / maxm)) {
				double h0 =
					NegativeParticleSampling{}.evaluateSource(v0, v1, sX, para);
				double hbar0 = h0 - m0 - alphaNeg * m0;
				if (hbar0 > 0) {
					// check v0 is in the pos zone
					double r2h0 = r1 * r1 * hbar0;
					double rr = RandomSampling(random).uniform();
					if (rr < (r2h0 / (alphaPos * m0))) {
						// accept the virtual particle v0 with suitable rate
						sOne.setVelocity(v0);
						sXNew.pushBack(sOne, ParticleKind::Negative);
					}
				}
			}
		}
	}
}

} // namespace coulomb
