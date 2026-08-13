#include "ProjectionSampling.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "Grid.h"
#include "ParticleGroup.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"

namespace coulomb {
using std::abs;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::pow;
using std::sqrt;
using std::vector;

void sampleFromP3MCoeffVer3(NeParticleGroup& sX, double dt, double dx,
							double& a0, double& a11, double& a2, double& a21,
							double& a31) {
	double rhoM = sX.rhoM;
	double u1M = sX.u1M;
	double tprtM = sX.tprtM;

	// double dxRhoM = S_x . dxRhoM;
	double dxU1M = sX.dxU1M;
	double dxTprtM = sX.dxTprtM;

	const double dimen = 3.;

	// coefficients from -\Delta t (I-\Pi_M) (v\cdot\nabla_x  + E \cdot\nabla_v)
	// M

	double sqrtT = sqrt(tprtM);

	a0 = 0.;
	a11 = -dt * rhoM * (-(dimen + 2) / 2. / sqrtT * dxTprtM);
	a21 = -dt * rhoM * dxU1M;
	a2 = -dt * rhoM * (-dxU1M / dimen);
	a31 = -dt * rhoM * (dxTprtM / sqrtT / 2.);

	// coefficients from \Delta t \Pi_M (v\cdot\nabla_x  + E \cdot\nabla_v) (f_p
	// - f_n)
	sX.computeMoments();

	// inner product with 1
	double drho = sX.drhoG;
	double coe0 = drho;
	// inner product with (v_1 - u_1)
	double dm1 = sX.dm1G; // inner product with v_1
	double coe1 = dm1 - u1M * drho;
	// inner product with ( |v - u|^2/T - d )
	double dEnergy = sX.dEnergyG; // inner product with |v|^2
	double coe2 = 2. / tprtM * dEnergy - 2. * u1M / tprtM * dm1 +
				  (u1M * u1M / tprtM - dimen) * drho;

	// cout <<"a11: " << a11 << ' ' << coe_1 / sqrt(tprtM) << endl;
	// cout <<"a2: " << a2 << ' ' << coe_2 / 2. / dimen << endl;

	a0 += coe0 - .5 * coe2;
	a11 += coe1 / sqrt(tprtM);
	a2 += coe2 / 2. / dimen;

	// coefficients need to multiply the grid size dx

	a0 *= dx;
	a11 *= dx;
	a21 *= dx;
	a2 *= dx;
	a31 *= dx;

	// rhoM *= u1M; // non sense
}

//  Step 1, Determine the number of particles to be sampled

int sampleFromP3MGetsize(double a0, double a11, double a2, double a21,
						 double a31, double neff, RandomContext& random) {
	double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
					  (abs(a2) + abs(a21)) * 4 * exp(-1.) +
					  abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
	maxratio = maxratio * pow(sqrt(2), 3);
	return RandomSampling(random).stochasticFloor(maxratio / neff);
}

//  Step2, sample.

NeParticleGroup sampleFromP3MSample(double a0, double a11, double a2,
									double a21, double a31, int ntotal,
									RandomContext& random) {
	NeParticleGroup sNew;
	double maxratio = abs(a0) + abs(a11) * sqrt(2.) * exp(-0.5) +
					  (abs(a2) + abs(a21)) * 4 * exp(-1.) +
					  abs(a31) * (6 * sqrt(6.) + 4 * sqrt(2.)) * exp(-1.5);
	// double maxratio = abs(a0) + abs(b) * sqrt(2.)*exp(-0.5) +
	//                  abs(c) * 4*exp(-1.) + abs(d) * (6*sqrt(6.) +
	//                  4*sqrt(2.))*exp(-1.5);
	maxratio = maxratio * pow(sqrt(2), 3);

	std::array<double, 3> v{};
	double coeM0 = 1.0 / pow(sqrt(2. * pi), 3);
	double coeM1 = 1.0 / pow(sqrt(4. * pi), 3);

	Particle1D3D sOne;

	double sqrt2 = sqrt(2.0);

	for (int k = 0; k < ntotal; k++) {
		double vsq = 0.;
		for (int kv = 0; kv < 3; kv++) {
			v[kv] = sqrt2 * RandomSampling(random).normal();
			vsq += v[kv] * v[kv];
		}

		// double M0 = coe_M0 * (a + b * v[0] + c * vsq + d * v[0]*vsq) *
		// exp(-vsq/2.);
		double m0 = coeM0 *
					(a0 + a11 * v[0] + a2 * vsq + a21 * v[0] * v[0] +
					 a31 * v[0] * vsq) *
					exp(-vsq / 2.);
		double m1 = coeM1 * exp(-vsq / 4.);

		if (RandomSampling(random).uniform() < (abs(m0) / m1 / maxratio)) {
			if (m0 > 0) {
				sOne.setVelocity(v);
				sNew.pushBack(sOne, ParticleKind::Positive);
			} else {
				sOne.setVelocity(v);
				sNew.pushBack(sOne, ParticleKind::Negative);
			}
		}
	}

	return sNew;

	// cout << (kp+kn)/( (double) Ntotal) << endl;
}

//  Step3, enforce conservation
NeParticleGroup sampleFromP3MRescale(const NeParticleGroup& sNew, double u1,
									 double tprt) {
	NeParticleGroup sRescaled;
	const auto& sp = sNew.list(ParticleKind::Positive);
	const auto& sn = sNew.list(ParticleKind::Negative);

	std::vector<double> vRescale(3);

	std::vector<double> uCenter{u1, 0., 0.};
	double sqrtT = sqrt(tprt);

	for (const auto& sone : sp) {
		const auto& vNormalized = sone.velocity();
		for (int kv = 0; kv < 3; kv++)
			vRescale[kv] = uCenter[kv] + sqrtT * vNormalized[kv];
		sRescaled.pushBack(Particle1D3D(sone.position(), vRescale),
						   ParticleKind::Positive);
	}

	for (const auto& sone : sn) {
		const auto& vNormalized = sone.velocity();
		for (int kv = 0; kv < 3; kv++)
			vRescale[kv] = uCenter[kv] + sqrtT * vNormalized[kv];
		sRescaled.pushBack(Particle1D3D(sone.position(), vRescale),
						   ParticleKind::Negative);
	}

	return sRescaled;
}

/**
  The whole algorithm of sampling from the micro-macro projection of advection
  part Sample P and N particles from - (I-\Pi) T M + \Pi T g
*/

// in one grid
void ProjectionSampling::sampleHomogeneous(NeParticleGroup& sX,
										   const NumericGridClass& grid,
										   RandomContext& random) {
	double a0, a11, a2, a21, a31;
	int ntotal;

	sampleFromP3MCoeffVer3(sX, grid.dt, grid.dx, a0, a11, a2, a21, a31);
	// sample_from_P3M_coeff_nog(S_x, grid.dt, grid.neff, a0, a11, a2, a21,
	// a31);
	ntotal = sampleFromP3MGetsize(a0, a11, a2, a21, a31, grid.neff, random);

	if (sX.tprtM < 0) {
		cout << " (" << sX.rhoM << ' ' << sX.u1M << ' ' << sX.tprtM << ") ";
		cout << a0 << ' ' << a11 << ' ' << a2 << ' ' << a21 << ' ' << a31 << ' '
			 << ntotal << endl;
	}

	auto sXNew = sampleFromP3MSample(a0, a11, a2, a21, a31, ntotal, random);

	sXNew = sampleFromP3MRescale(sXNew, sX.u1M, sX.tprtM);

	ParticleGroupOperations{}.assignPositions(sXNew, sX.getXMin(), sX.getXMax(),
											  random);
	// cout << "( " << ptr_S_x_new.size('p') << ", " << ptr_S_x_new.size('n')
	// <<
	// ") ";

	ParticleGroupOperations{}.mergeSigned(sX, sXNew);

	if ((sX.size(ParticleKind::Positive) + sX.size(ParticleKind::Negative)) >
		200) {
		// ParticleConservation{}.enforceZero(S_x, grid.neff);
	}

	/*
	cout << "Np, Nn = " << S_x.size(ParticleKind::Positive) << ' '
		 << S_x.size(ParticleKind::Negative) << endl;
	cout << ", after cons 2d = " <<  S_x . positiveMoments.m0 - S_x .
	negativeMoments.m0 << endl; S_x . computeMoments(); cout << ", after cons
	2e = " <<  S_x . positiveMoments.m0 - S_x . negativeMoments.m0 << endl; if
	(abs(S_x . positiveMoments.m0 - S_x. negativeMoments.m0)>.5) { cout <<
	"conserv m0, out, " <<  S_x . positiveMoments.m0 << ' ' << S_x .
	negativeMoments.m0 << endl; exit(0);
	}
	*/
}

// over all grids
void ProjectionSampling::sample(std::vector<NeParticleGroup>& sX,
								const NumericGridClass& grid,
								RandomContext& random) {
	int nx = grid.nx;
	for (int kx = 0; kx < nx; kx++) {
		ProjectionSampling::sampleHomogeneous(sX[kx], grid, random);
	}
}

} // namespace coulomb
