#pragma once

#include <string>
#include <vector>

#include "Classes.h"
#include "_global_variables.h"

namespace coulomb {

void particleresample_inhomo(std::vector<NeParticleGroup>& groups,
                             NumericGridClass& grid, ParaClass& parameters,
                             SimulationState& state);
void sync_coarse(std::vector<NeParticleGroup>& groups, NumericGridClass& grid,
                 ParaClass& parameters, SimulationState& state);

void coulomb_collision_homo_PFNF(NeParticleGroup& groups,
                                  const ParaClass& parameters,
                                  RandomContext& random);
double evaluateM(const std::vector<double>& velocity,
                 const NeParticleGroup& groups);
double evaluateH(const std::vector<double>& velocity,
                 const std::vector<double>& source_velocity,
                 NeParticleGroup& groups, const ParaClass& parameters,
                 int mode = 0);
void finddeltambound(NeParticleGroup& groups, const ParaClass& parameters);
void finddeltambound_inhomo(std::vector<NeParticleGroup>& groups,
                            const NumericGridClass& grid,
                            const ParaClass& parameters);
int samplefromDeltamp_Npv(NeParticleGroup& groups, double effective_particles,
                          RandomContext& random);
void samplefromDeltam(NeParticleGroup& groups, NeParticleGroup& new_groups,
                      const ParaClass& parameters, double effective_particles,
                      RandomContext& random);
void merge_NeParticleGroup(NeParticleGroup& groups,
                           const NeParticleGroup& new_groups);
void mergeF_NeParticleGroup(NeParticleGroup& groups,
                            const NeParticleGroup& new_groups);
void mergeNeParticleGroup(NeParticleGroup& groups,
                          const NeParticleGroup& new_groups,
                          const std::string& particle_types);
void assign_positions(NeParticleGroup& groups, double xmin, double xmax,
                      RandomContext& random);

void NegPar_collision_homo(NeParticleGroup& groups,
                           const ParaClass& parameters, double effective_particles,
                           RandomContext& random);
void NegPar_collision(std::vector<NeParticleGroup>& groups,
                      const NumericGridClass& grid,
                      const ParaClass& parameters, RandomContext& random);
void NegPar_collision_openmp(std::vector<NeParticleGroup>& groups,
                             const NumericGridClass& grid,
                             const ParaClass& parameters, RandomContext& random);
void enforce_conservation(double m0, double m11, double m12, double m13,
                          double m21, double m22, double m23,
                          NeParticleGroup& groups, double effective_particles,
                          bool conserve_energy_vector, RandomContext& random);
void enforce_conservation_zero(NeParticleGroup& groups,
                               double effective_particles,
                               RandomContext& random);

void sample_from_P3M_coeff(NeParticleGroup& groups, double dt, double& a0,
                           double& a11, double& a2, double& a21, double& a31);
void sample_from_P3M_coeff_ver2(NeParticleGroup* groups, double dt,
                                double effective_particles, double& a0,
                                double& a11, double& a2, double& a21,
                                double& a31);
void sample_from_P3M_coeff_ver3(NeParticleGroup& groups, double dt,
                                double dx, double& a0, double& a11, double& a2,
                                double& a21, double& a31);
int sample_from_P3M_getsize(double a0, double a11, double a2, double a21,
                            double a31, double effective_particles,
                            RandomContext& random);
NeParticleGroup sample_from_P3M_sample(double a0, double a11, double a2,
                                       double a21, double a31, int total,
                                       RandomContext& random);
void sample_from_P3M_conserve(double a0, double a11, double a2, double a21,
                              double a31, NeParticleGroup& groups,
                              double effective_particles,
                              RandomContext& random);
NeParticleGroup sample_from_P3M_rescale(const NeParticleGroup& groups,
                                        double u1, double temperature);
void sample_from_MMprojection_homo(NeParticleGroup& groups,
                                   const NumericGridClass& grid,
                                   RandomContext& random);
void sample_from_MMprojection(std::vector<NeParticleGroup>& groups,
                              const NumericGridClass& grid,
                              RandomContext& random);

void update_macro(std::vector<NeParticleGroup>& groups,
                  const NumericGridClass& grid);
void NegPar_BGK_collision_homo(NeParticleGroup& groups, ParaClass& parameters,
                               RandomContext& random);
void NegPar_BGK_collision(std::vector<NeParticleGroup>& groups,
                          NumericGridClass& grid, ParaClass& parameters,
                          RandomContext& random);
void Negpar_inhomo_onestep(std::vector<NeParticleGroup>& groups,
                           NumericGridClass& grid, ParaClass& parameters,
                           SimulationState& state);
void Negpar_inhomo_onestep_ver2(std::vector<NeParticleGroup>& groups,
                                NumericGridClass& grid, ParaClass& parameters,
                                SimulationState& state);
void Negpar_inhomo_onestep_PIC(std::vector<NeParticleGroup>& groups,
                               NumericGridClass& grid, ParaClass& parameters,
                               SimulationState& state);
void Negpar_inhomo_onestep_stop(std::vector<NeParticleGroup>& groups,
                                NumericGridClass& grid, ParaClass& parameters,
                                int flag_stop, SimulationState& state);

}  // namespace coulomb
