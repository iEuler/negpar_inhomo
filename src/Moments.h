#pragma once

#include <tuple>

#include "Classes.h"
#include "_global_variables.h"

namespace coulomb {

void uT2mE(const double& rho, double* u, const double& temperature,
           double* momentum, double& energy);
void mE2uT(const double& rho, double* momentum, const double& energy,
           double* u, double& temperature);

void uT2mE_x1v3(int grid_size, const std::vector<double>& rho,
                const std::vector<double>& u1,
                const std::vector<double>& temperature,
                std::vector<double>& m1, std::vector<double>& energy);
void mE2uT_x1v3(int grid_size, const std::vector<double>& rho,
                const std::vector<double>& m1,
                const std::vector<double>& energy,
                std::vector<double>& u1,
                std::vector<double>& temperature);

void particle2rhomE(const ParticleGroup& group, double effective_particles,
                    double& rho, double* momentum, double& energy);
void compute_rhouT(int grid_size, const std::vector<ParticleGroup>& groups,
                   double effective_particles, std::vector<double>& rho,
                   std::vector<double>& u1, std::vector<double>& u2,
                   std::vector<double>& u3,
                   std::vector<double>& temperature);

void update_rhouT(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& grid);
void update_rhouT_F(std::vector<NeParticleGroup>& groups,
                    const NumericGridClass& grid);
void update_dx_rhouT_M(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid);
void update_maxwellian(std::vector<NeParticleGroup>& groups,
                       const NumericGridClass& grid);
void compute_change_in_macro(std::vector<NeParticleGroup>& groups,
                             const NumericGridClass& grid,
                             SimulationState& state);
void compute_change_in_macro_onlykineitc(
    std::vector<NeParticleGroup>& groups, const NumericGridClass& grid);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
momentchange_g(NeParticleGroup* groups, const NumericGridClass& grid);
void momentchange_g_ver2(NeParticleGroup* groups, const NumericGridClass& grid,
                         std::vector<double>& drho,
                         std::vector<double>& dm1,
                         std::vector<double>& denergy);

}  // namespace coulomb
