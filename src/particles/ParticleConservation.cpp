#include "ParticleConservation.h"

#include <array>
#include <cmath>
#include <iostream>

#include "ParticleGroup.h"
#include "RandomSampling.h"

namespace coulomb {

void ParticleConservation::enforce(double m0, double m11, double m12,
								   double m13, double m21, double m22,
								   double m23, NeParticleGroup& sNew,
								   double neff, bool flagConserveEnergyvector,
								   RandomContext& random) {
	// enforce m0
	double m0Need = m0;
	sNew.computeMoments();
	double m0Actual =
		neff * (sNew.positiveMoments.m0 - sNew.negativeMoments.m0);
	// cout << "before cons = " <<  S_new . positiveMoments.m0 - S_new .
	// negativeMoments.m0;
	int nRemove;
	if (m0Actual < m0Need) {
		nRemove =
			RandomSampling(random).stochasticFloor((m0Need - m0Actual) / neff);
		for (int kp = 0; kp < nRemove; kp++) {
			int kRemove = (int)(RandomSampling(random).uniform() *
								sNew.size(ParticleKind::Negative));
			sNew.erase(kRemove, ParticleKind::Negative);
		}
	} else {
		nRemove =
			RandomSampling(random).stochasticFloor((m0Actual - m0Need) / neff);
		for (int kp = 0; kp < nRemove; kp++) {
			int kRemove = (int)(RandomSampling(random).uniform() *
								sNew.size(ParticleKind::Positive));
			sNew.erase(kRemove, ParticleKind::Positive);
		}
	}

	int np = sNew.size(ParticleKind::Positive);
	int nn = sNew.size(ParticleKind::Negative);

	// enforce m11, m12, m13

	sNew.computeMoments();
	double m1Actual[3], m1Need[3] = {m11, m12, m13};
	m1Actual[0] = neff * (sNew.positiveMoments.m11 - sNew.negativeMoments.m11);
	m1Actual[1] = neff * (sNew.positiveMoments.m12 - sNew.negativeMoments.m12);
	m1Actual[2] = neff * (sNew.positiveMoments.m13 - sNew.negativeMoments.m13);

	std::array<double, 3> v0{};
	double m1Mod[3];
	for (int kv = 0; kv < 3; kv++)
		m1Mod[kv] = -m1Actual[kv] + m1Need[kv];

	if (np > nn) {
		auto& sp = sNew.list(ParticleKind::Positive);
		for (int kp = 0; kp < np; kp++) {
			auto& vkp = sp[kp].velocity();
			for (int kv = 0; kv < 3; kv++)
				v0[kv] = vkp[kv] + m1Mod[kv] / neff / np;
			sp[kp].setVelocity(v0);
		}
	} else {
		auto& sn = sNew.list(ParticleKind::Negative);
		for (int kn = 0; kn < nn; kn++) {
			auto& vkn = sn[kn].velocity();
			for (int kv = 0; kv < 3; kv++)
				v0[kv] = vkn[kv] - m1Mod[kv] / neff / nn;
			sn[kn].setVelocity(v0);
		}
	}

	// enforce m21, m22, m23
	// mu2P*Tp - mu2N*Tn + Np*cp^2 - Nn*cn^2 = Ep_need - En_need

	sNew.computeMoments();

	double cp[3], cn[3], tp[3], tn[3], rhs[3];
	double m2pActual[3], m2nActual[3];
	double m2Need[3] = {m21 / neff, m22 / neff, m23 / neff};

	cp[0] = sNew.positiveMoments.m11;
	cp[1] = sNew.positiveMoments.m12;
	cp[2] = sNew.positiveMoments.m13;
	cn[0] = sNew.negativeMoments.m11;
	cn[1] = sNew.negativeMoments.m12;
	cn[2] = sNew.negativeMoments.m13;
	for (int kv = 0; kv < 3; kv++) {
		cp[kv] /= np;
		cn[kv] /= nn;
	}

	m2pActual[0] = sNew.positiveMoments.m21;
	m2pActual[1] = sNew.positiveMoments.m22;
	m2pActual[2] = sNew.positiveMoments.m23;
	m2nActual[0] = sNew.negativeMoments.m21;
	m2nActual[1] = sNew.negativeMoments.m22;
	m2nActual[2] = sNew.negativeMoments.m23;

	double mu2P[3], mu2N[3];
	for (int kv = 0; kv < 3; kv++) {
		tp[kv] = m2pActual[kv] - np * cp[kv] * cp[kv];
		tn[kv] = m2nActual[kv] - nn * cn[kv] * cn[kv];
		rhs[kv] = m2Need[kv] - np * cp[kv] * cp[kv] + nn * cn[kv] * cn[kv];
	}

	if (flagConserveEnergyvector) {
		for (int kv = 0; kv < 3; kv++) {
			if (np > nn) {
				mu2N[kv] = 1.0;
				mu2P[kv] = 1.0 / tp[kv] * (mu2N[kv] * tn[kv] + rhs[kv]);
				if (mu2P[kv] < 0) {
					mu2P[kv] = 1.0;
					std::cout << "ERROR NOT CONSERVATIVE" << std::endl;
				}
			} else {
				mu2P[kv] = 1.0;
				mu2N[kv] = 1.0 / tn[kv] * (mu2P[kv] * tp[kv] - rhs[kv]);
				if (mu2N[kv] < 0) {
					mu2N[kv] = 1.0;
					std::cout << "ERROR NOT CONSERVATIVE" << std::endl;
				}
			}
		}
	} else {
		double sumRhs = 0., sumTp = 0., sumTn = 0.;
		double mu2NAll, mu2PAll;
		for (int kv = 0; kv < 3; kv++) {
			sumRhs += rhs[kv];
			sumTp += tp[kv];
			sumTn += tn[kv];
		}
		if (np > nn) {
			mu2NAll = 1.0;
			mu2PAll = 1.0 / sumTp * (mu2NAll * sumTn + sumRhs);
			if (mu2PAll < 0)
				mu2PAll = 1.0;
		} else {
			mu2PAll = 1.0;
			mu2NAll = 1.0 / sumTn * (mu2PAll * sumTp - sumRhs);
			if (mu2NAll < 0)
				mu2NAll = 1.0;
		}
		for (int kv = 0; kv < 3; kv++) {
			mu2N[kv] = mu2NAll;
			mu2P[kv] = mu2PAll;
		}
	}

	if (np > nn) {
		auto& sp = sNew.list(ParticleKind::Positive);
		for (int kp = 0; kp < np; kp++) {
			auto& vkp = sp[kp].velocity();
			for (int kv = 0; kv < 3; kv++)
				v0[kv] = std::sqrt(mu2P[kv]) * (vkp[kv] - cp[kv]) + cp[kv];
			sp[kp].setVelocity(v0);
		}
	} else {
		auto& sn = sNew.list(ParticleKind::Negative);
		for (int kn = 0; kn < nn; kn++) {
			auto& vkn = sn[kn].velocity();
			for (int kv = 0; kv < 3; kv++)
				v0[kv] = std::sqrt(mu2N[kv]) * (vkn[kv] - cn[kv]) + cn[kv];
			sn[kn].setVelocity(v0);
		}
	}
	sNew.computeMoments();
}

void ParticleConservation::enforceZero(NeParticleGroup& sNew, double neff,
									   RandomContext& random) {
	ParticleConservation::enforce(0., 0., 0., 0., 0., 0., 0., sNew, neff, false,
								  random);
}

} // namespace coulomb
