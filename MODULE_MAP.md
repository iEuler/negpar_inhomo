# Module Map

This is the current ownership map after the focused-module migration.

| Area | Current owner | Intended owner |
| --- | --- | --- |
| Simulation enums and mode names | `SimulationTypes.h/.cpp` | `simulation/types` |
| Physical and numerical parameters | `SimulationConfig.h/.cpp` | `simulation/configuration` |
| Numerical grid data | `Grid.h/.cpp` | `simulation/grid` |
| Particle value and kind | `Particle.h/.cpp` | `particles/particle` |
| Particle groups and moments | `ParticleGroup.h/.cpp` | `particles/group` |
| Typed particle-group merging and randomized placement | `ParticleGroupOperations.h/.cpp` | `particles/operations` |
| Shared 3-D tensor aliases | `TensorTypes.h` | `numerics/tensor_types` |
| FFT wrappers | `FFT.h/.cpp` | `numerics/fft` |
| 1-D finite-difference and kinetic Euler numerics | `Numerics.h/.cpp` | `numerics/utilities` |
| Macroscopic moments and Maxwellian updates | `Moments.h/.cpp` | `physics/moments` |
| Initial-condition definitions and Two-Stream preprocessing | `InitialConditions.h/.cpp` | `initialization/conditions` |
| Stochastic particle construction | `ParticleInitialization.h/.cpp` | `initialization/particles` |
| Grid initialization orchestration | `Initialization.h/.cpp` | `initialization/orchestration` |
| Advection | `Advection.h/.cpp` | `physics/advection` |
| Poisson solve and field updates | `ElectricField.h/.cpp` | `physics/electric_field` |
| Binary Coulomb collisions | `Collisions.h/.cpp` | `physics/collisions` |
| Negative-particle collision kernels and pipeline orchestration | `NegativeParticleCollisions.h/.cpp` | `physics/collisions` |
| Negative-particle collision-source sampling | `NegativeParticleSampling.h/.cpp` | `physics/collisions/sampling` |
| Signed particle conservation enforcement | `ParticleConservation.h/.cpp` | `physics/conservation` |
| HDP and PIC numerical-step sequencing | `SimulationSteps.h/.cpp` | `simulation/steps` |
| Maxwellian projection sampling | `ProjectionSampling.h/.cpp` | `physics/projection_sampling` |
| Signed Fourier resampling | `Resampler*.h/.cpp`, `ResamplingNumerics.*`, `ResamplingVelocity.*` | `resampling/` |
| Deterministic full-particle Fourier numerics | `FullParticleFourier.h/.cpp` | `resampling/full_particles/fourier` |
| Stochastic full-particle reconstruction | `FullParticleSampling.h/.cpp` | `resampling/full_particles/sampling` |
| Resampling policy and orchestration | `ParticleResampling.h/.cpp` | `resampling/policy` |
| Output path and numbering policy | `OutputPaths.h/.cpp` | `io/output/paths` |
| Macro, field, and moment output | `MacroOutput.h/.cpp` | `io/output/macro` |
| Distribution and particle output | `ParticleOutput.h/.cpp` | `io/output/particles` |
| Grid, parameter, and initial-condition metadata | `RunMetadataOutput.h/.cpp` | `io/output/metadata` |
| Particle-count diagnostics | `Diagnostics.h/.cpp` | `diagnostics` |
| Fixed-bin histogram construction | `Histogram.h/.cpp` | `diagnostics/histogram` |
| Legacy header compatibility | none retained in `src/` | callers must include the owning module header |
| Runtime options | `RunOptions.h/.cpp` | `simulation/configuration` |
| Mathematical constants | `Constants.h` | `numerics/constants` |
| Per-run RNG ownership and seeding | `RandomContext.h/.cpp` | `simulation/random_context` |
| Thread-local random draws and stochastic helpers | `RandomSampling.h/.cpp` | `simulation/random_sampling` |
| Mutable counters, flags, and timing | `SimulationState.h` | `simulation/state` |
| Macro-step and output orchestration | `Simulation.h/.cpp` | `simulation/` |

The extracted numerical modules and promoted `FourierResampler` are the active
numerical path. The former
`coulomb_*.h` compatibility umbrellas and forward-declaration header were
unused by production and test targets and have been removed. Public callers
should include the focused owning headers directly (`Initialization.h`,
`SimulationSteps.h`, `FullParticleSampling.h`, `MacroOutput.h`, and so on).
Those public headers forward-declare shared model types where possible and
include the focused owners when a complete value type is required. The former
`Classes.h/.cpp` monolith and umbrella are no longer present.

Focused modules expose their behavior through explicit class owners rather
than namespace-level functions. Examples include `Advection`,
`ElectricFieldSolver`, `CollisionOperator`, `MomentOperations`,
`ParticleResampling`, `MacroOutput`, and `RunMetadataOutput`. Stateless
operations are static members of their owning class; value objects such as
`RunOptions` and `RandomContext` own the behavior that acts on their state.
Translation-unit-only implementation helpers remain private in anonymous
namespaces.

Constants, random-context ownership, and mutable runtime state now have focused
headers; the historical `_global_variables.h` umbrella has been removed.
Mutable run state is owned by the caller as a `SimulationState` object.
Advection, collision-step, resampling, initialization, timing, and output
entry points receive that object explicitly. No process-wide aliases or
mutable state instances remain. Diagnostics are deliberately state-free pure
reductions because they do not consume or mutate run state.

The generic `utils.*` bucket has been retired. `RandomContext.*` owns seed and
generation state, `RandomSampling.*` owns thread-local engines and stochastic
draw helpers, and `Histogram.*` owns fixed-bin output aggregation. The former
vector-extrema helpers were private to collision-source bound construction and
were replaced there with standard-library algorithms. The uncalled
`momentchange_g_ver2` alternative was removed after a production-call audit.

Particle-group operations use the typed ParticleKind enum throughout active
code. The particle_kind_code and particle_kind_from_code helpers are kept only
for explicit compatibility and serialized character-code boundaries.

Grid simulation modes use SimulationMethod, and spatial boundary behavior uses
BoundaryCondition. Their character representations are produced only when
writing legacy-compatible parameter output.

Particle-cell moments are grouped in the initialized `Moments` value object.
`ParticleGroup` owns one record, while `NeParticleGroup` owns named current
positive, negative, and full records plus positive and negative snapshots from
before advection. All particle kinds share one accumulation implementation,
and aggregate macro-state refresh now belongs to `Moments.cpp`.

FFTW buffers and plans use the `FFTWBuffer` and `FFTWPlan` RAII handles,
including the legacy interpolation routines, so allocation or plan failures
cannot bypass cleanup.

Simulation orchestration now lives in `Simulation.cpp`. Its internal
`SimulationRunner` separates initialization, output scheduling, numerical-step
dispatch, checkpointing, and finalization, while `SimulationHistory` owns
diagnostic and timing records. The active numerical step variants are the HDP
and PIC paths; the unreferenced `ver2` and conditional-stop duplicates were
removed after a repository-wide usage check.

The historical
`inhomo_neg_coulomb.cpp` file is intentionally retained as an empty compatibility
translation unit for older project layouts. It is not part of any build target
and contains no implementation or global state.

Strict warning flags apply to all production code. The `sanitizer` preset runs
MSVC AddressSanitizer;
supported GCC/Clang builds additionally enable UndefinedBehaviorSanitizer.
`NEGPAR_ENABLE_THREAD_SANITIZER` provides a separate race-detection build for
supported GCC/Clang platforms.

OpenMP collision work partitions particle groups by spatial cell, so worker
threads do not mutate the same group. Each worker obtains a thread-local RNG
engine derived from the run-owned `RandomContext` seed and its OpenMP thread
identifier; reseeding is performed outside parallel regions. Regression
coverage checks fixed-schedule replay plus per-cell momentum and kinetic-energy
invariants with multiple threads. FFTW plans and buffers remain object-owned,
non-copyable, and are exercised through repeated forward/inverse lifecycles in
the sanitizer suite.

P/N-to-full Coulomb collisions and BGK relaxation are isolated in
`NegativeParticleCollisions.*`. Their fixed-seed characterization verifies
exact replay, particle-count behavior, finite rethermalized velocities, and
the established rule that the P/N-to-full path does not mutate full particles.
The module also owns the homogeneous and spatial collision pipelines; direct
pipeline characterization verifies fixed counts and positions, finite evolved
velocities, and exact fixed-seed replay across all particle kinds.

Collision-source Maxwellian evaluation, envelope construction, and Delta-M
sampling are isolated in `NegativeParticleSampling.*`. Direct characterization
checks the analytic Maxwellian value, finite positive sampling bounds,
deterministic virtual-candidate counts, exact fixed-seed replay, and finite
accepted particle velocities.

Signed mass, momentum, and energy enforcement is isolated in
`ParticleConservation.*`. Direct characterization perturbs a nondegenerate
particle group, verifies restoration of all seven requested signed moments,
and checks exact fixed-seed replay of the resulting particle velocities.

Cross-module particle-list merging and randomized spatial placement are owned
by `ParticleGroupOperations.*`. Their characterization verifies typed merge
selection and order, velocity preservation, position bounds, and exact
fixed-seed position replay.

The obsolete `coulomb2` implementation and its duplicate class/global files
are no longer present. A repository-wide dependency search found no remaining
`coulomb2` callers.

`resampling::FourierResampler` now owns signed P/N Fourier resampling in the
core application. Promotion followed exact same-seed particle equivalence with
the former legacy implementation plus end-to-end mass, energy, bounds, and
determinism characterization in exact and approximate modes. Shared frequency,
Taylor, and velocity-coordinate utilities live in `ResamplingNumerics.*` and
`ResamplingVelocity.*`. The duplicate legacy signed implementation and the
uninitialized `InhomoResampler` experiment were removed; full-particle
Maxwellian reconstruction is split between `FullParticleFourier.*` and
`FullParticleSampling.*`.

The full-particle public API was audited declaration by declaration before the
remaining implementation was split. `resample_F_from_MPN` is the only external
entry point and is called by `ParticleResampling.cpp`; its Fourier transform,
interpolation, Maxwellian, upper-bound, acceptance, and `filter_Fourier`
helpers are all live through that path. The separate `sampleF` /
`sampleF_inhomo` count-rescaling branch had no production or test caller and
was removed together with its grid, configuration, and diagnostics
dependencies.

Deterministic full-particle coefficient construction, inverse interpolation,
derivative grids, Maxwellian derivatives, filtering, and interpolation-cell
bounds now live in `FullParticleFourier.*`. Term-level transform and
Maxwellian helpers are private to that translation unit. The module contains no
random-number consumption or acceptance/rejection sampling; the remaining
stochastic reconstruction loop delegates all Fourier calculations through its
focused public API. The former `FullParticleResampling.*` name is retired;
`ParticleResampling.cpp` now delegates full-particle reconstruction to
`FullParticleSampling.*`, which owns the exact RNG-consuming acceptance loop.
Fixed-seed characterization checks exact output replay, particle-kind counts,
positions and velocities, finite restored coordinates, velocity bounds,
low-order mass and momentum, and preservation of the input particle lists.

Initialization is split into problem-specific macro definitions
(`InitialConditions.*`), RNG-driven particle construction
(`ParticleInitialization.*`), and grid orchestration (`Initialization.*`).
The selected Landau-damping defaults and all alternative problem formulas are
unchanged. The uncalled `initialize_distri_Negpar_test` branch and handwritten
local prototypes were removed. Focused tests characterize selected defaults,
exact fixed-seed Delta particles and moments, and Two-Stream preprocessing.

Output ownership is split across path policy, macro/field serialization,
particle/distribution serialization, and run metadata. The former
`Output.h/.cpp` umbrella is removed. Filenames, numbering, precision, and row
ordering remain characterized exactly, and `output_path` rejects absolute or
parent-traversal filenames so artifacts cannot escape the configured output
directory.

HDP and PIC step sequencing now belongs to `SimulationSteps.*`; the
misleading `NegativeParticle.*` module is removed. `Simulation.cpp` is its
only caller, and the collision, field, advection, projection, resampling, and
synchronization order is unchanged. The final dependency audit also removed
the remaining handwritten function prototype from a production `.cpp` file.
