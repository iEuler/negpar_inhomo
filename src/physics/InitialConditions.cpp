#include "InitialConditions.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "Grid.h"
#include "SimulationConfig.h"

namespace coulomb {
using std::abs;
using std::cos;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::pow;
using std::sin;
using std::sqrt;
using std::vector;

IniValClass InitialConditions::create(NumericGridClass& grid) {
	IniValClass initialData;

	initialData.problemName = "LandauDamping";
	// initialData.problemName = "TwoStreamInstab";
	// initialData.problemName = "test1";
	// initialData.problemName = "BumpOnTail";
	// initialData.problemName = "Analytic";
	// initialData.problemName = "Efficiency";

	if (initialData.problemName == "Efficiency")
		grid.lambdaPoisson = 0;
	if (initialData.problemName == "TwoStreamInstab")
		InitialConditions::configureTwoStream(initialData);

	cout << "Problem name = " << initialData.problemName << endl;
	initialData.totalMass = 0;
	return initialData;
}

void InitialConditions::configure(IniValClass& initialData,
								  const NumericGridClass& grid, int cell) {
	const double x = grid.x[cell];
	const double dx = grid.x[1] - grid.x[0];
	const double spatialVolume = grid.xmax - grid.xmin;
	const double spatialCenter = (grid.x[0] + grid.x[grid.nx - 1]) / 2.0;

	if ((initialData.problemName == "LandauDamping") ||
		(initialData.problemName == "Efficiency")) {
		initialData.rho = 1 + initialData.ldAlpha * sin(x);
		initialData.velocity[0] = 0.;
		initialData.velocity[1] = 0.;
		initialData.velocity[2] = 0.;
		initialData.tprt = 1.0;
		initialData.totalMass += initialData.rho * dx;
	} else if (initialData.problemName == "TwoStreamInstab") {
		constexpr double alpha = 0.5;
		constexpr double waveNumber = 0.5;
		initialData.tsiCoe =
			1. / 12 / pi * (1.0 + alpha * cos(waveNumber * x)) / 40.;
	} else if (initialData.problemName == "BumpOnTail") {
		const double rho0 =
			initialData.botRho0 / spatialVolume * initialData.botBeta;
		const double rho1 = initialData.botRho0 * (1 - initialData.botBeta) /
							spatialVolume *
							exp(-(x - spatialCenter) * (x - spatialCenter) /
								(2 * initialData.botTx));
		const double rhoTotal = rho0 + rho1;
		const double backgroundVelocity = 0.;
		const double bumpVelocity = initialData.botUb;
		const double momentum = rho1 * bumpVelocity;
		const double velocity = momentum / rhoTotal;
		const double backgroundTemperature = initialData.botTprt;
		const double bumpTemperature = initialData.botDTprt;
		const double energy =
			1.5 * rho0 * backgroundTemperature +
			rho1 * (.5 * bumpVelocity * bumpVelocity + 1.5 * bumpTemperature);
		const double temperature =
			(energy - .5 * rhoTotal * velocity * velocity) / 1.5 / rhoTotal;

		initialData.rho = rhoTotal;
		initialData.velocity[0] = velocity;
		initialData.velocity[1] = 0.;
		initialData.velocity[2] = 0.;
		initialData.tprt = temperature;

		initialData.rho1 = rho0;
		initialData.velocity1[0] = backgroundVelocity;
		initialData.velocity1[1] = 0.;
		initialData.velocity1[2] = 0.;
		initialData.tprt1 = backgroundTemperature;

		initialData.rho2 = rho1;
		initialData.velocity2[0] = bumpVelocity;
		initialData.velocity2[1] = 0.;
		initialData.velocity2[2] = 0.;
		initialData.tprt2 = bumpTemperature;
		initialData.totalMass += initialData.rho * dx;
	}
}

void InitialConditions::configureTwoStream(IniValClass& initialData) {
	double vmax = 6.;
	int nv = 200;
	double dv = 2.0 * vmax / nv;
	vector<double> v(nv);
	for (int kv = 0; kv < nv; kv++)
		v[kv] = (kv + 0.5) * dv - vmax;

	double v1, v2, v3, energyf = 0.;
	vector<double> m0(nv * nv * nv);
	vector<double> f0(nv * nv * nv);
	double rhoF = 0.;
	for (int kv1 = 0; kv1 < nv; kv1++) {
		v1 = v[kv1];
		for (int kv2 = 0; kv2 < nv; kv2++) {
			v2 = v[kv2];
			for (int kv3 = 0; kv3 < nv; kv3++) {
				v3 = v[kv3];
				int kk = kv3 + nv * (kv2 + nv * kv1);
				double vsq = v1 * v1 + v2 * v2 + v3 * v3;
				f0[kk] = exp(-vsq / 2.) * (1.0 + 5.0 * v1 * v1);
				rhoF += f0[kk];
				energyf += .5 * vsq * f0[kk];
			}
		}
	}
	rhoF *= dv * dv * dv;
	energyf *= dv * dv * dv;
	double tprt = energyf / (1.5 * rhoF);

	double rhoP = 0.;
	double maxFOverM = 0.;
	double m21 = 0., m22 = 0., m23 = 0.;
	for (int kv1 = 0; kv1 < nv; kv1++) {
		v1 = v[kv1];
		for (int kv2 = 0; kv2 < nv; kv2++) {
			v2 = v[kv2];
			for (int kv3 = 0; kv3 < nv; kv3++) {
				v3 = v[kv3];
				int kk = kv3 + nv * (kv2 + nv * kv1);
				double vsq = v1 * v1 + v2 * v2 + v3 * v3;
				m0[kk] =
					rhoF / pow(sqrt(2. * pi * tprt), 3) * exp(-vsq / 2. / tprt);
				maxFOverM = max(maxFOverM, abs(f0[kk] - m0[kk]) / m0[kk]);
				m21 += v1 * v1 * f0[kk];
				m22 += v2 * v2 * f0[kk];
				m23 += v3 * v3 * f0[kk];
				if (f0[kk] > m0[kk])
					rhoP += f0[kk] - m0[kk];
			}
		}
	}

	rhoP *= dv * dv * dv;
	m21 *= dv * dv * dv;
	m22 *= dv * dv * dv;
	m23 *= dv * dv * dv;

	initialData.tsiRhof = rhoF;
	initialData.tsiRhop = rhoP;
	initialData.tsiTprt = tprt;
	initialData.tsiMaxFOverM = maxFOverM;
	initialData.tsiM21 = m21;
	initialData.tsiM22 = m22;
	initialData.tsiM23 = m23;

	cout << "Initially " << m21 + m22 + m23 << " vs " << 3 * rhoF * tprt
		 << endl;
}

} // namespace coulomb
