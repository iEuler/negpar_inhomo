#include "SimulationConfig.h"

namespace coulomb {

ParaClass::ParaClass() {
	method = SimulationMethod::HDP;
	flagUseOpenMp = true;
	dt = 0.01;
	coeffBinaryColl = 10.0;
	methodBinaryColl = BinaryCollisionMethod::TA;
	resampleRatio = 1.2;
	nPickupNeg = 100;
	nfreq = 30;
	useApproximation = true;
	collisionType = CollisionType::Coulomb;
	lambdaPoisson = 10.0;
	resampleSpatialRatio = 0.9;
	syncTimeInterval = 0.5;
	resampleSyncRatio = 1.1;
	hdpCouplingMode = HdpCouplingMode::Decoupled;
	effectiveWeightPolicy = EffectiveWeightPolicy::Fixed;
	collisionCoupling = CollisionCoupling::Standard;
	projectionMode = ProjectionMode::FullMicroMacro;
	deltaMMode = DeltaMMode::Enabled;
	weightedFourierResampling = false;
	partialResampling = false;
	conserveWeightedMoments = false;
	coulombBgkHybrid = false;
	partialResamplingCutoff = 3.0;
	signedWeightMin = 1e-7;
	signedWeightMax = 5e-3;
	fullWeightMin = 1e-7;
	fullWeightMax = 5e-3;
	cpuCostConstant = 0.205;
	cpuCostCollisionCoefficient = 3.277;
	bgkStrength = 0.0;
}

IniValClass::IniValClass() {
	problemName = "TwoPeaks";
	problemNameExt = "BumpOnTail";
	totalMass = 0;
	rho1 = .9;
	rho2 = .1;
	velocity1[0] = 0.0;
	velocity1[1] = 0.0;
	velocity1[2] = 0.0;
	velocity2[0] = 5.0;
	velocity2[1] = 0.0;
	velocity2[2] = 0.0;
	tprt1 = 1.0;
	tprt2 = .01;
	ldAlpha = 0.4;
	botBeta = 0.9;
	botRho0 = 1.0;
	botTprt = 1.0;
	botDTprt = 0.01;
	botUb = 5.0;
	botTx = 0.25;
}

} // namespace coulomb
