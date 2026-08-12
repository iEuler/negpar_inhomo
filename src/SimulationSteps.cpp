#include "SimulationSteps.h"

#include <ctime>
#include <iostream>

#include "Advection.h"
#include "Collisions.h"
#include "Diagnostics.h"
#include "ElectricField.h"
#include "Grid.h"
#include "Moments.h"
#include "NegativeParticleCollisions.h"
#include "ParticleGroup.h"
#include "ParticleResampling.h"
#include "ProjectionSampling.h"
#include "SimulationConfig.h"
#include "SimulationState.h"

namespace coulomb {
using std::cout;
using std::endl;

/* ======================================================== *\
        Forward one step
\* ======================================================== */

/** One step
  Forward one step in time, with time splitting
*/

void SimulationSteps::advance_hdp(std::vector<NeParticleGroup> &S_x,
                                  NumericGridClass &grid, ParaClass &para,
                                  SimulationState &state) {
  // cout << "step start" << endl;

  // Step 1, collision.

  // Step 1.0 update all macro quantities
  MomentOperations::update_macro(S_x, grid);

  // Step 1.0 perform negative collisions

  state.t0Collision = clock();

  if (para.collisionType == CollisionType::Coulomb)
    // NegativeParticleCollisions::collide(S_x, grid, para);
    NegativeParticleCollisions::collide_parallel(S_x, grid, para, state.random);
  else if (para.collisionType == CollisionType::BGK)
    NegativeParticleCollisions::collide_bgk(S_x, grid, para, state.random);

  // cout << "step 1" << endl;

  state.t1Collision = clock();

  // step 2, advection

  state.t0Advection = state.t1Collision;

  // Step 2.0 update all macro quantities and electric field
  MomentOperations::update_macro(S_x, grid);
  ElectricFieldSolver::update(S_x, grid);

  for (int kx = 0; kx < grid.Nx; kx++)
    S_x[kx].copymoments();

  // cout << "step 2.0" << endl;

  // Switch 2.1 and 2.2

  // Step 2.1, compute moment change: S_x.drho, dm1, denergy
  MomentOperations::compute_macro_change(S_x, grid, state);
  // cout << "step 2.1" << endl;

  // Step 2.2, advect P N F particles.
  Advection::advance(S_x, grid, state);
  // cout << "step 2.2" << endl;

  // Step 2.3, Sample P and N particles from micro-macro projection
  ProjectionSampling::sample(S_x, grid, state.random);
  // cout << "step 2.3" << endl;

  // Step 2.4, update maxwellian part:S_x.rhoM, u1M, TprtM
  MomentOperations::update_maxwellian(S_x, grid);
  // cout << "step 2.4" << endl;

  // cout << "d(Np, Nn) = (" << Npcoll - Nplast << ", " << Nncoll - Nnlast
  //      << "), (" << Npadve - Npcoll << ", " << Nnadve - Nncoll << ")" <<
  //      endl;

  state.t1Advection = clock();

  // Step 3, resampling particles when needed
  // ParticleResampling::resample(S_x, grid, para, MLsol);

  state.t0Resampling = state.t1Advection;
  if (para.collisionType == CollisionType::Coulomb) {
    ParticleResampling::resample(S_x, grid, para, state);
  }
  state.t1Resampling = clock();

  ParticleResampling::synchronize_coarse(S_x, grid, para, state);

  // cout << "Np = " << Diagnostics::particle_count(S_x, grid.Nx, 'p')
  //      << "; Nn = " << Diagnostics::particle_count(S_x, grid.Nx, 'n')
  //      << "; Nf = " << Diagnostics::particle_count(S_x, grid.Nx, 'f') <<
  //      endl;
}

void SimulationSteps::advance_pic(std::vector<NeParticleGroup> &S_x,
                                  NumericGridClass &grid, ParaClass &para,
                                  SimulationState &state) {
  state.t0Collision = clock();

  for (int kx = 0; kx < grid.Nx; kx++) {
    auto &Sf = S_x[kx].list(ParticleKind::Full);
    CollisionOperator::collide_homogeneous(Sf, S_x[kx].size(ParticleKind::Full),
                                           para, state.random);
  }

  state.t1Collision = clock();

  ElectricFieldSolver::update_pic(S_x, grid);

  state.t0Advection = clock();
  Advection::advance(S_x, ParticleKind::Full, grid, state);
  state.t1Advection = clock();

  cout << "Np = "
       << Diagnostics::particle_count(S_x, grid.Nx, ParticleKind::Positive)
       << "; Nn = "
       << Diagnostics::particle_count(S_x, grid.Nx, ParticleKind::Negative)
       << "; Nf = "
       << Diagnostics::particle_count(S_x, grid.Nx, ParticleKind::Full) << endl;
}

} // namespace coulomb
