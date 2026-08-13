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

void SimulationSteps::advanceHdp(std::vector<NeParticleGroup>& sX) {
	auto& grid = gridRef;
	auto& para = parametersRef;
	auto& state = stateRef;
	// cout << "step start" << endl;

	// Step 1, collision.

	// Step 1.0 update all macro quantities
	MomentOperations{}.updateMacro(sX, grid);

	// Step 1.0 perform negative collisions

	state.t0Collision = clock();

	if (para.collisionType == CollisionType::Coulomb)
		// NegativeParticleCollisions{}.collide(S_x, grid, para);
		NegativeParticleCollisions(grid, para, state.random)
			.collideParallel(sX);
	else if (para.collisionType == CollisionType::BGK)
		NegativeParticleCollisions(grid, para, state.random).collideBgk(sX);

	// cout << "step 1" << endl;

	state.t1Collision = clock();

	// step 2, advection

	state.t0Advection = state.t1Collision;

	// Step 2.0 update all macro quantities and electric field
	MomentOperations{}.updateMacro(sX, grid);
	ElectricFieldSolver(grid).update(sX);

	for (int kx = 0; kx < grid.nx; kx++)
		sX[kx].copyMoments();

	// cout << "step 2.0" << endl;

	// Switch 2.1 and 2.2

	// Step 2.1, compute moment change: S_x.drho, dm1, dEnergy
	MomentOperations{}.computeMacroChange(sX, grid, state);
	// cout << "step 2.1" << endl;

	// Step 2.2, advect P N F particles.
	Advection(grid, state).advance(sX);
	// cout << "step 2.2" << endl;

	// Step 2.3, Sample P and N particles from micro-macro projection
	ProjectionSampling{}.sample(sX, grid, state.random);
	// cout << "step 2.3" << endl;

	// Step 2.4, update maxwellian part:S_x.rhoM, u1M, tprtM
	MomentOperations{}.updateMaxwellian(sX, grid);
	// cout << "step 2.4" << endl;

	// cout << "d(Np, Nn) = (" << Npcoll - Nplast << ", " << Nncoll - Nnlast
	//      << "), (" << Npadve - Npcoll << ", " << Nnadve - Nncoll << ")" <<
	//      endl;

	state.t1Advection = clock();

	// Step 3, resampling particles when needed
	// ParticleResampling{}.resample(S_x, grid, para, MLsol);

	state.t0Resampling = state.t1Advection;
	if (para.collisionType == CollisionType::Coulomb) {
		ParticleResampling(grid, para, state).resample(sX);
	}
	state.t1Resampling = clock();

	ParticleResampling(grid, para, state).synchronizeCoarse(sX);

	// cout << "Np = " << Diagnostics::particleCount(S_x, grid.nx, 'p')
	//      << "; Nn = " << Diagnostics::particleCount(S_x, grid.nx, 'n')
	//      << "; Nf = " << Diagnostics::particleCount(S_x, grid.nx, 'f') <<
	//      endl;
}

void SimulationSteps::advancePic(std::vector<NeParticleGroup>& sX) {
	auto& grid = gridRef;
	auto& para = parametersRef;
	auto& state = stateRef;
	state.t0Collision = clock();

	for (int kx = 0; kx < grid.nx; kx++) {
		auto& sf = sX[kx].list(ParticleKind::Full);
		CollisionOperator(para, state.random)
			.collideHomogeneous(sf, sX[kx].size(ParticleKind::Full));
	}

	state.t1Collision = clock();

	ElectricFieldSolver(grid).updatePic(sX);

	state.t0Advection = clock();
	Advection(grid, state).advance(sX, ParticleKind::Full);
	state.t1Advection = clock();

	cout << "Np = "
		 << Diagnostics(grid).particleCount(sX, grid.nx, ParticleKind::Positive)
		 << "; Nn = "
		 << Diagnostics(grid).particleCount(sX, grid.nx, ParticleKind::Negative)
		 << "; Nf = "
		 << Diagnostics(grid).particleCount(sX, grid.nx, ParticleKind::Full)
		 << endl;
}

} // namespace coulomb
