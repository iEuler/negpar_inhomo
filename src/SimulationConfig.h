#pragma once

#include <string>

#include "SimulationTypes.h"

namespace coulomb {

class ParaClass {
 public:
  BinaryCollisionMethod method_binarycoll;
  SimulationMethod method;
  CollisionType collisionType;
  double lambda_Poisson;
  double coeff_binarycoll;
  bool FLAG_USE_OPENMP;
  double dt;
  int Npickup_neg;
  double resample_ratio, resample_spatial_ratio, sync_time_interval,
      resample_sync_ratio;
  int Nfreq;
  bool useApproximation;

  ParaClass();
};

class IniValClass {
 public:
  std::string probname, probname_ext;
  double totalmass{};
  double rho1{}, rho2{}, Tprt1{}, Tprt2{};
  double velocity1[3]{}, velocity2[3]{};
  double rho{}, Tprt{}, dTprt{};
  double velocity[3]{};
  double TSI_alpha{}, TSI_coe{}, TSI_rhof{}, TSI_rhop{}, TSI_Tprt{},
      TSI_max_f_over_M{}, TSI_m21{}, TSI_m22{}, TSI_m23{};
  double a{}, b{}, c{};
  double LD_alpha{};
  double BOT_beta{}, BOT_rho0{}, BOT_Tprt{}, BOT_dTprt{}, BOT_Tx{}, BOT_ub{};

  IniValClass();
};

}  // namespace coulomb
