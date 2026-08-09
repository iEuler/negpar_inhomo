# Refactoring Plan

This plan prioritizes preserving the numerical behavior of the simulation while
improving reproducibility, testability, structure, and portability.

## Refactor goals

The remaining work is tracked as sequential review checkpoints. Only one goal
is active at a time so structural changes stay attributable to a specific
validation result.

1. **Complete the module-ownership split.** Promote the characterized Fourier
   resampler, remove its superseded and uninitialized duplicates, finish the
   focused resampling modules, update both build systems, and pass Debug,
   Release, sanitizer, and short-reference validation.
2. **Separate configuration from runtime state.** Replace the historical
   `_global_variables.h` boundary with focused configuration, random-context,
   and simulation-state headers; reduce transitive includes and validate each
   migrated module.
3. **Harden resource and concurrency boundaries.** Finish warning cleanup,
   exercise FFT/resource ownership under sanitizers, and add supported
   multithreaded invariant or race-detection coverage without weakening the
   deterministic single-thread reference gate.
4. **Finish orchestration cleanup and handoff.** Keep the entry point thin,
   reduce the remaining broad `Classes` dependencies, confirm module ownership
   documentation matches the code, and run the complete validation protocol.
5. **Split the shared data-model monolith.** Extract simulation types and
   configuration, grid, particle, particle-group, and tensor ownership from
   `Classes.h/.cpp`; migrate all consumers and remove the obsolete umbrella.
6. **Extract negative-particle collision kernels.** Separate P/N-to-full
   Coulomb and BGK kernels from sampling, conservation, and time-step
   orchestration; characterize their fixed-seed behavior independently.
7. **Extract signed particle conservation.** Move signed mass, momentum, and
   energy enforcement behind a focused module boundary; migrate callers and
   build manifests, and directly characterize deterministic moment restoration.
8. **Extract negative-particle source sampling.** Isolate Maxwellian/source
   evaluation, sampling-bound construction, and Delta-M particle generation;
   characterize accepted samples and fixed-seed replay directly.
9. **Extract shared particle-group operations.** Move typed merge and random
   placement helpers behind a focused boundary, replace broad includes and
   handwritten declarations, and characterize data preservation and replay.
10. **Consolidate negative-particle collision orchestration.** Move the
    homogeneous, spatial, and OpenMP collision pipelines beside their kernels,
    leaving `SimulationSteps.*` focused on top-level time-step sequencing.

Goals 1-10 are complete. Each dependency pass remains separately validated
from the large resampling migration.

## Phase 1: Document and build the current program

- Document the current Windows and Linux build and run commands.
- Build the currently supported Debug and Release configurations.
- Record compiler and linker versions, dependencies, build flags, and warnings.
- Record the current default configuration, generated files, and output location.
- Record known limitations without changing behavior yet.

**Checkpoint:** the legacy program builds and its execution environment is
documented well enough to reproduce locally.

## Phase 2: Make reference runs reproducible

Add the minimum controls needed to create reliable characterization data:

- Allow the random seed to be specified and record it in the run metadata.
- Record RNG provenance in the run metadata, including the engine, seeding
  procedure, distribution algorithms or implementation details, standard
  library and version, and whether cross-platform numerical identity is
  expected.
- Allow the thread count to be specified.
- Validate the requested thread count; do not derive an invalid value from
  `omp_get_max_threads()` on small machines.
- Use one thread for reference runs.
- Investigate and remove the shared-RNG data race before treating OpenMP runs as
  reliable. Prefer explicitly owned or deterministically partitioned RNG state.
- Add a short-run configuration that exercises initialization, advection,
  collisions, resampling, diagnostics, and output where practical.
- Allow each run to use an isolated output directory.

Run the short reference configuration and save:

- The seed, thread count, configuration, compiler, and build type.
- Particle counts and conservation quantities at defined steps.
- Representative output files.
- NaN/Inf checks and generated-file inventory.
- Runtime as benchmark information, not as a correctness gate.

**Checkpoint:** the same single-thread configuration and seed reproduce the
reference results under the designated reference toolchain.

## Phase 3: Establish build and test infrastructure

- Add CMake support for Windows and Linux while retaining the Visual Studio
  project until CMake reaches parity.
- Define separate library, executable, and test targets so tests do not need to
  include the simulation entry point.
- Use Catch2 v3.x through CMake `FetchContent`, pinned to an exact tag or
  commit, with a documented offline/dependency-cache workflow.
- Add CMake presets for supported Debug and Release configurations.
- Register tests with CTest.
- Ensure test output is written to temporary or build-tree directories.
- Add CI only after the supported dependency installation process is clear.

**Checkpoint:** the application and an empty/smoke test suite build through
CMake on each supported platform.

## Phase 4: Add characterization and regression tests

Start with tests that capture current behavior before changing architecture.
Cover:

- Particle construction, insertion, deletion, and clearing.
- Particle and group moment calculations.
- FFT forward/inverse behavior, including the current normalization convention.
- Periodic and reflective boundary handling.
- Collision invariants intended by the algorithm.
- Resampling conservation properties.
- A short end-to-end reference run.

Classify each assertion as one of:

- Exact deterministic comparison.
- Absolute/relative numerical comparison.
- Invariant or bounded-drift check.
- Statistical/ensemble comparison for stochastic behavior.

For every numerical comparison, document the observable, absolute and relative
tolerances, reference scale, seed, thread count, and whether particle ordering
matters. Reject unexpected NaN and Inf values explicitly.

Do not silently turn a discovered legacy defect into the expected behavior.
Record it as a known defect, add a failing or quarantined reproducer where
appropriate, and fix it in a separate reviewed change.

**Checkpoint:** the reference run and focused characterization tests pass
without intentional numerical changes.

## Phase 5: Remove verified-dead duplicate code

The active implementation is `coulomb::` in the focused simulation, grid,
particle, and particle-group modules. The obsolete `coulomb2` implementation,
duplicate globals, and unused legacy header umbrellas have been removed after
dependency and compile checks. No production or test target depends on
`coulomb2`.

The only intentional compatibility boundary is the empty
`src/inhomo_neg_coulomb.cpp` translation unit, retained for older project
layouts; it is not compiled by the current CMake or Visual Studio targets.

**Checkpoint:** no production or test target depends on `coulomb2`, and behavior
matches the reference run.

## Phase 6: Inventory modules and existing refactors

Before creating a new directory hierarchy, map responsibilities and dependencies
across:

- Focused module headers (`Initialization.h`, `SimulationSteps.h`,
  `FullParticleSampling.h`, `ParticleResampling.h`, `MacroOutput.h`)
- `src/SimulationTypes.*`, `src/SimulationConfig.*`, `src/Grid.*`,
  `src/Particle.*`, `src/ParticleGroup.*`, and `src/TensorTypes.h`
- `src/FFT.h/.cpp`
- `src/Resampler.h/.cpp`
- `src/ResamplerHelper.h/.cpp`
- `src/ResamplingNumerics.h/.cpp`
- `src/ResamplingVelocity.h/.cpp`

Identify functions already superseded by newer `.cpp` components. Choose one
implementation for each responsibility and avoid introducing a third version.
Record the intended public interfaces and dependency direction.

**Checkpoint:** every major function has an identified owning module, and known
duplicate implementations have an explicit disposition.

## Phase 7: Introduce explicit run state while splitting headers

Separate constants, immutable configuration, mutable simulation state, timing,
and output policy. For example:

```cpp
struct SimulationState {
    int saveIndex{};
    int movedCount{};
    int resampleCount{};
    double syncTime{};
};
```

Keep mathematical constants such as `pi` out of mutable state. Give RNG state
clear ownership instead of leaving it as a process-wide global.

Move implementations from the large headers into cohesive `.cpp` modules while
passing the required state explicitly. Suggested boundaries, subject to the
Phase 6 inventory, are:

```text
src/
  simulation/
    simulation.h/.cpp
    configuration.h/.cpp
    state.h
  physics/
    advection.h/.cpp
    electric_field.h/.cpp
    collisions.h/.cpp
  particles/
    particle.h/.cpp
    particle_group.h/.cpp
  resampling/
    resampler.h/.cpp
    interpolation.h/.cpp
  io/
    initialization.h/.cpp
    output.h/.cpp
  numerics/
    fft.h/.cpp
    utilities.h/.cpp
```

Migrate one cohesive responsibility at a time and run the regression suite after
each migration. Avoid first moving code and then changing all interfaces in a
separate global-state phase.

**Checkpoint:** implementation-heavy legacy headers are reduced to declarations,
and mutable run state has explicit ownership.

## Phase 8: Improve type safety and data modeling

- Replace `'p'`, `'n'`, and `'f'` with `enum class ParticleKind`.
- Replace string modes such as `"HDP"` and `"PIC"` with enums.
- Validate grid and simulation parameters at construction time.
- Correct constructor patterns that create discarded temporaries instead of
  delegating construction.
- Group repeated moment fields into a `Moments` structure where doing so does
  not obscure the numerical formulas.
- Mark read-only methods and parameters `const`.
- Replace unchecked indexing only where profiling shows no unacceptable cost or
  where validation can occur at an API boundary.

Make these changes incrementally rather than as a repository-wide mechanical
rewrite.

## Phase 9: Resource, concurrency, and warning cleanup

- Immediately make owning FFT classes non-copyable; add move support only if it
  is needed and safe.
- Wrap FFTW allocations and plans in RAII types and handle plan-allocation
  failure.
- Enable OpenMP explicitly in supported builds after shared mutable state and RNG
  ownership are safe.
- Add deterministic single-thread tests and separate multithreaded invariant or
  stress tests.
- Raise compiler warning levels incrementally and fix newly enabled warnings.
- Enable AddressSanitizer and UndefinedBehaviorSanitizer where supported.
- Use platform-appropriate race-detection tooling where available.
- Remove bundled or stale dependency files only after the supported dependency
  strategy and license implications are documented.

## Phase 10: Refactor the simulation loop

The active entry point is now `src/main.cpp`; the historical
`src/inhomo_neg_coulomb.cpp` file is an empty, unbuilt compatibility unit.
`main.cpp` performs configuration handling and delegates execution to
`Simulation`:

```cpp
int main(int argc, char** argv) {
    const auto options = parse_run_options(argc, argv);
    SimulationState state;
    apply_run_options(options, state);
    return Simulation(options, state).run();
}
```

Move initialization, time stepping, diagnostics, output scheduling, and timing
into dedicated components with explicit dependencies. Preserve the established
ordering of numerical operations unless an intentional, separately validated
algorithm change is approved.

## Validation protocol

After every independently reviewable change:

1. Build the supported Debug and Release configurations.
2. Run unit, characterization, and applicable sanitizer tests.
3. Run the single-thread short reference simulation with its recorded seed.
4. Reject unexpected NaN or Inf values.
5. Compare particle counts, conservation quantities, and declared output
   artifacts using their documented comparison class and tolerances.
6. Confirm generated files stay within the selected output directory.
7. Record runtime separately as benchmark data; investigate large changes but do
   not make noisy wall-clock timing a correctness assertion.
8. Commit each coherent phase or migration separately.

Multithreaded and cross-platform runs are additional validation dimensions; they
must not replace the designated deterministic reference run. Bitwise equality
across compilers, standard libraries, FFTW builds, or thread schedules is not
assumed unless demonstrated and explicitly declared as an expected property of
the reference configuration.

## Review gates

- Changes to collision, advection, electric-field, FFT, resampling, or moment
  calculations require the relevant invariant results and reference-run
  comparison in the review description.
- Changes that intentionally alter numerical behavior require a stated rationale,
  updated tolerances or reference data, and focused review separate from
  mechanical restructuring.
- Dead-code deletion requires a recorded usage/dependency search and a successful
  full build and test run.
- Mechanical changes should remain separate from physics-sensitive changes so
  reviewers can attribute numerical differences.
