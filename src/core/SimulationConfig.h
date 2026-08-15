#pragma once

#include <string>

#include "SimulationTypes.h"

namespace coulomb {

class ParaClass {
  public:
	BinaryCollisionMethod methodBinaryColl;
	SimulationMethod method;
	CollisionType collisionType;
	double lambdaPoisson;
	double coeffBinaryColl;
	bool flagUseOpenMp;
	double dt;
	int nPickupNeg;
	double resampleRatio, resampleSpatialRatio, syncTimeInterval,
		resampleSyncRatio;
	int nfreq;
	bool useApproximation;
	HdpCouplingMode hdpCouplingMode;
	EffectiveWeightPolicy effectiveWeightPolicy;
	CollisionCoupling collisionCoupling;
	ProjectionMode projectionMode;
	DeltaMMode deltaMMode;
	bool weightedFourierResampling;
	bool partialResampling;
	bool conserveWeightedMoments;
	bool coulombBgkHybrid;
	double partialResamplingCutoff;
	double signedWeightMin, signedWeightMax;
	double fullWeightMin, fullWeightMax;
	double cpuCostConstant, cpuCostCollisionCoefficient;
	double bgkStrength;

	ParaClass();
};

class IniValClass {
  public:
	std::string problemName, problemNameExt;
	double totalMass{};
	double rho1{}, rho2{}, tprt1{}, tprt2{};
	double velocity1[3]{}, velocity2[3]{};
	double rho{}, tprt{}, dTprt{};
	double velocity[3]{};
	double tsiAlpha{}, tsiCoe{}, tsiRhof{}, tsiRhop{}, tsiTprt{},
		tsiMaxFOverM{}, tsiM21{}, tsiM22{}, tsiM23{};
	double a{}, b{}, c{};
	double ldAlpha{};
	double botBeta{}, botRho0{}, botTprt{}, botDTprt{}, botTx{}, botUb{};

	IniValClass();
};

} // namespace coulomb
