#include "SimulationConfig.h"

namespace coulomb {

ParaClass::ParaClass() {
	method = SimulationMethod::HDP;
	FLAG_USE_OPENMP = true;
	dt = 0.01;
	coeff_binarycoll = 10.0;
	method_binarycoll = BinaryCollisionMethod::TA;
	resample_ratio = 1.2;
	Npickup_neg = 100;
	Nfreq = 30;
	useApproximation = true;
	collisionType = CollisionType::Coulomb;
	lambda_Poisson = 10.0;
	resample_spatial_ratio = 0.9;
	sync_time_interval = 0.5;
	resample_sync_ratio = 1.1;
}

IniValClass::IniValClass() {
	probname = "TwoPeaks";
	probname_ext = "BumpOnTail";
	totalmass = 0;
	rho1 = .9;
	rho2 = .1;
	velocity1[0] = 0.0;
	velocity1[1] = 0.0;
	velocity1[2] = 0.0;
	velocity2[0] = 5.0;
	velocity2[1] = 0.0;
	velocity2[2] = 0.0;
	Tprt1 = 1.0;
	Tprt2 = .01;
	LD_alpha = 0.4;
	BOT_beta = 0.9;
	BOT_rho0 = 1.0;
	BOT_Tprt = 1.0;
	BOT_dTprt = 0.01;
	BOT_ub = 5.0;
	BOT_Tx = 0.25;
}

} // namespace coulomb
