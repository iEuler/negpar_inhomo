#include "ParticleInitialization.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Constants.h"
#include "ParticleConservation.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

namespace coulomb {
using std::abs;
using std::cout;
using std::endl;
using std::exp;
using std::pow;
using std::sqrt;
using std::string;
using std::vector;
// Generate a P, N, F particle list with designated distribution

void ParticleInitialization::initialize(NeParticleGroup& sX,
										const IniValClass& initialData,
										double neff, double neffF, double dx,
										RandomContext& random) {
	// initialize_Negpar_size(int &Np, int &Nn, int &Nf);

	string problemName = initialData.problemName;

	double rhoF, rhoP, tprt, maxFOverM;
	rhoF = initialData.tsiRhof;
	rhoP = initialData.tsiRhop;
	tprt = initialData.tsiTprt;
	maxFOverM = initialData.tsiMaxFOverM;

	/*
	int Np = S_x->size('p');
	int Nn = S_x->size('n');
	int Nf = S_x->size('f');
	*/

	double x1, x2;
	x1 = sX.getXMin();
	x2 = sX.getXMax();

	if ((problemName == "LandauDamping") || (problemName == "Efficiency")) {
		// decide the size
		// int Np = 0, Nn = 0;
		int nf = RandomSampling(random).stochasticFloor(initialData.rho * dx /
														neffF);

		// Particle1D3D * Sp = S_x->list('p');
		// Particle1D3D * Sn = S_x->list('n');
		// Particle1D3D * Sf = S_x->list('f');

		// update maxwellian part

		sX.rhoM = initialData.rho;
		sX.u1M = initialData.velocity[0];
		sX.u2M = initialData.velocity[1];
		sX.u3M = initialData.velocity[2];
		sX.tprtM = initialData.tprt;

		// create F particles

		std::vector<double> vf(3);

		double sqrtT = sqrt(initialData.tprt);
		for (int kf = 0; kf < nf; kf++) {
			for (int k = 0; k < 3; k++)
				vf[k] = initialData.velocity[k] +
						sqrtT * RandomSampling(random).normal();

			Particle1D3D sOne(RandomSampling(random).uniform() * (x2 - x1) + x1,
							  vf);
			sX.pushBack(sOne, ParticleKind::Full);
		}

		// create P and N particles

	} else if (problemName == "Delta") {
		// decide the size
		// int Np = 0, Nn = 0;
		int nf = RandomSampling(random).stochasticFloor(initialData.rho * dx /
														neffF);

		// Particle1D3D * Sp = S_x->list('p');
		// Particle1D3D * Sn = S_x->list('n');
		// Particle1D3D * Sf = S_x->list('f');

		// update maxwellian part

		sX.rhoM = initialData.rho;
		sX.u1M = initialData.velocity[0];
		sX.u2M = initialData.velocity[1];
		sX.u3M = initialData.velocity[2];
		sX.tprtM = initialData.tprt;

		// create F particles

		std::vector<double> vf(3);

		for (int kf = 0; kf < nf; kf++) {
			for (int k = 0; k < 3; k++)
				vf[k] = initialData.velocity[k];
			// for (int k=0;k<3;k++) vf[k] = initialData.velocity[k] +
			// A deterministic thermal draw would use
			// RandomSampling(random).normal().

			Particle1D3D sOne(RandomSampling(random).uniform() * (x2 - x1) + x1,
							  vf);
			sX.pushBack(sOne, ParticleKind::Full);
		}

	} else if (problemName == "TwoStreamInstab") {
		// decide the size
		int np, nn, nf;
		np = RandomSampling(random).stochasticFloor(rhoP * initialData.tsiCoe *
													dx / neff);
		nn = np;
		nf = RandomSampling(random).stochasticFloor(rhoF * initialData.tsiCoe *
													dx / neffF);

		cout << np << ' ' << nn << ' ' << nf << endl;

		// Particle1D3D * Sp = S_x->list('p');
		// Particle1D3D * Sn = S_x->list('n');
		// Particle1D3D * Sf = S_x->list('f');

		// update maxwellian part

		sX.rhoM = rhoF * initialData.tsiCoe;
		sX.u1M = 0.;
		sX.u2M = 0.;
		sX.u3M = 0.;
		sX.tprtM = tprt;

		// create F particles

		double vSq = 1.8;
		double maxF0 = exp(-vSq / 2) * (1 + 5 * vSq);
		std::vector<double> vp(3);
		double vmax = 6.;

		// int kf = 0;
		while (sX.size(ParticleKind::Full) < nf) {
			double v1 = (RandomSampling(random).uniform() - .5) * 2 * vmax;
			if (RandomSampling(random).uniform() <
				(exp(-v1 * v1 / 2) * (1 + 5 * v1 * v1) / maxF0)) {
				vp[0] = v1;
				for (int k = 1; k < 3; k++)
					vp[k] = RandomSampling(random).normal();

				Particle1D3D sOne(
					RandomSampling(random).uniform() * (x2 - x1) + x1, vp);
				sX.pushBack(sOne, ParticleKind::Full);
			}
		}

		// create P and N particles
		double coeM0 = rhoF / pow(sqrt(2. * pi * tprt), 3);
		double sqrtT = sqrt(tprt);
		while ((sX.size(ParticleKind::Positive) < np) ||
			   (sX.size(ParticleKind::Negative) < nn)) {
			// Velocity draws use the explicit RandomContext parameter.
			std::vector<double> vpSample{
				RandomSampling(random).normal() * sqrtT,
				RandomSampling(random).normal() * sqrtT,
				RandomSampling(random).normal() * sqrtT};
			double vsq = vpSample[0] * vpSample[0] + vpSample[1] * vpSample[1] +
						 vpSample[2] * vpSample[2];
			double f0 = exp(-vsq / 2) * (1 + 5 * vpSample[0] * vpSample[0]);
			double m0 = coeM0 * exp(-vsq / 2 / tprt);
			// cout << (f0/m0) /max_f_over_M << endl;
			if (RandomSampling(random).uniform() <
				abs((f0 - m0) / m0 / maxFOverM)) {
				if (f0 > m0) {
					if (sX.size(ParticleKind::Positive) < np) {
						Particle1D3D sOne(
							RandomSampling(random).uniform() * (x2 - x1) + x1,
							vpSample);
						sX.pushBack(sOne, ParticleKind::Positive);
					}
				} else {
					if (sX.size(ParticleKind::Negative) < nn) {
						Particle1D3D sOne(
							RandomSampling(random).uniform() * (x2 - x1) + x1,
							vpSample);
						sX.pushBack(sOne, ParticleKind::Negative);
					}
				}
			}
		}

		ParticleGroupOperations{}.assignPositions(sX, x1, x2, random);

		double m21 = initialData.tsiCoe * (initialData.tsiM21 - rhoF * tprt);
		double m22 = initialData.tsiCoe * (initialData.tsiM22 - rhoF * tprt);
		double m23 = initialData.tsiCoe * (initialData.tsiM23 - rhoF * tprt);
		// cout << "Now " <<  m21 << ' ' << m22 << ' ' << m23 << endl;
		ParticleConservation{}.enforce(0., 0., 0., 0., m21, m22, m23, sX, neff,
									   true, random);

	} else if (problemName == "BumpOnTail") {
		// decide the size
		int np =
			RandomSampling(random).stochasticFloor(initialData.rho * dx / neff);
		// int Nn = Np;
		int nf = RandomSampling(random).stochasticFloor(initialData.rho * dx /
														neffF);

		// Particle1D3D * Sp = S_x->list('p');
		// Particle1D3D * Sn = S_x->list('n');
		// Particle1D3D * Sf = S_x->list('f');
		// Particle1D3D * Sp = S_x->list('p');
		// Particle1D3D * Sn = S_x->list('n');

		// update maxwellian part

		sX.rhoM = initialData.rho;
		sX.u1M = initialData.velocity[0];
		sX.u2M = initialData.velocity[1];
		sX.u3M = initialData.velocity[2];
		sX.tprtM = initialData.tprt;

		double rho = initialData.rho;
		double rho1 = initialData.rho1;
		double rho2 = initialData.rho2;
		double u = initialData.velocity[0];
		double u1 = initialData.velocity1[0];
		double u2 = initialData.velocity2[0];
		double tprtBump = initialData.tprt;
		double tprt1 = initialData.tprt1;
		double tprt2 = initialData.tprt2;

		// create F particles

		std::vector<double> vf(3);
		double center[3] = {0, 0, 0}, sigma;
		for (int kf = 0; kf < nf; kf++) {
			if (RandomSampling(random).uniform() < rho1 / rho) {
				center[0] = u1;
				sigma = sqrt(tprt1);
			} else {
				center[0] = u2;
				sigma = sqrt(tprt2);
			}
			for (int k = 0; k < 3; k++)
				vf[k] = center[k] + sigma * RandomSampling(random).normal();
			Particle1D3D sOne(RandomSampling(random).uniform() * (x2 - x1) + x1,
							  vf);
			sX.pushBack(sOne, ParticleKind::Full);
		}

		// create P and N particles
		// int kp = 0, kn = 0;
		double coeT0T = sqrt(tprtBump / tprt1);
		coeT0T = coeT0T * coeT0T * coeT0T;
		for (int kf = 0; kf < np; kf++) {
			if (RandomSampling(random).uniform() < rho1 / rho) {
				center[0] = u;
				sigma = sqrt(tprtBump);
				for (int k = 0; k < 3; k++)
					vf[k] = center[k] + sigma * RandomSampling(random).normal();
				double rhoTemp =
					rho1 * coeT0T *
						exp(-(vf[0] * vf[0] + vf[1] * vf[1] + vf[2] * vf[2]) /
								2 / tprt1 +
							((vf[0] - u) * (vf[0] - u) + vf[1] * vf[1] +
							 vf[2] * vf[2]) /
								2 / tprtBump) -
					rho1 - rho2;
				if (RandomSampling(random).uniform() < (abs(rhoTemp) / rho1)) {
					if (rhoTemp > 0) {
						Particle1D3D sOne(
							RandomSampling(random).uniform() * (x2 - x1) + x1,
							vf);
						sX.pushBack(sOne, ParticleKind::Positive);
					} else {
						Particle1D3D sOne(
							RandomSampling(random).uniform() * (x2 - x1) + x1,
							vf);
						sX.pushBack(sOne, ParticleKind::Negative);
					}
				}

			} else {
				// generate one particle in the high bump
				center[0] = u2;
				sigma = sqrt(tprt2);
				for (int k = 0; k < 3; k++)
					vf[k] = center[k] + sigma * RandomSampling(random).normal();
				Particle1D3D sOne(
					RandomSampling(random).uniform() * (x2 - x1) + x1, vf);
				sX.pushBack(sOne, ParticleKind::Positive);
			}
		}

		// cout << rho1 << " " << rho2 << " " << u << " " << tprt1 << " " <<
		// tprt2
		// << " " << tprt << endl; cout << "[ "<< kp << ", " << kn << " ]" <<
		// endl;
	}

	// cout << "bounds " <<  S_x->getXMin() << ' ' << S_x->getXMax() << endl;

	sX.computeMoments();
}
} // namespace coulomb
