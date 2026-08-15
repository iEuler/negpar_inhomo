# Weighted HDP and Adaptive Resampling Plan

## Objective

Add the functionality found in `D:\Work\Research\NegPar_combine\code` that is
not present in the current repository, while preserving current numerical
behavior as the default until the new modes are characterized and validated.

The primary additions are:

- Variance-weighted coupling of coarse/full and deviational particle
  representations (weighted HDP, or WHDP).
- Frequency-dependent weighting during Fourier resampling.
- Partial core/tail resampling.
- Adaptive selection of signed and full effective particle weights.
- Optional collision-source, linearized-collision, projection, and hybrid BGK
  modes.

The comparison folder is a behavioral reference, not code to copy directly. It
contains invalid C++ in the unfinished Rosenbluth initialization, process-wide
mutable state, unchecked divisions, and numerical paths without focused
regression coverage.

## Guiding constraints

- Keep current decoupled HDP behavior as the default during development.
- Express modes with typed configuration rather than global booleans.
- Keep configuration, run state, numerical policy, and reconstruction in their
  existing focused modules.
- Make synchronization transactional: failure must not partially update cells
  or effective weights.
- Preserve deterministic replay for a fixed seed, thread count, build, and
  standard-library implementation.
- Reject non-finite weights, moments, Fourier coefficients, and particle
  coordinates explicitly.
- Separate intentional numerical changes from mechanical refactoring.

## Implementation protocol

Use this protocol when executing the plan:

1. Work on only one numbered phase at a time. Do not implement later phases
   merely because related code is nearby.
2. Begin each phase by identifying its owned files, current callers, baseline
   tests, and the smallest independently verifiable change.
3. Treat `D:\Work\Research\NegPar_combine\code` as read-only behavioral
   reference material. Port equations and intended behavior into current module
   boundaries; do not copy its global-state structure or broken code.
4. Preserve current defaults and fixed-seed behavior unless a phase explicitly
   authorizes and tests a numerical change.
5. Add or update tests in the same change as each behavior. Do not defer
   conservation, failure-path, or deterministic-replay coverage to a later
   cleanup phase.
6. Update CMake and Visual Studio project manifests together whenever source or
   test files are added, removed, or renamed.
7. Keep mechanical refactoring separate from physics-sensitive formula changes
   so any numerical difference has an attributable cause.
8. Run focused tests while iterating, followed by the complete applicable test
   suite, reference validation, and `git diff --check` before declaring a
   checkpoint complete.
9. At every checkpoint, report changed files, configuration/default changes,
   tests and validation results, remaining numerical uncertainties, and the
   next planned phase.
10. Stop and request direction when a legacy formula is ambiguous, required
    behavior conflicts with a current invariant, or proceeding would require a
    default change not authorized by this plan.

Do not port the unfinished Rosenbluth initialization or unrelated diagnostic
helpers. Begin implementation with Phase 1 and Phase 2; do not attempt all
phases in a single change.

## Phase 1: Specify and characterize the legacy algorithms

Document the equations and execution order for:

- Cell-level variance weights used to blend full and deviational estimates.
- Frequency-dependent weights used to blend Fourier coefficients.
- Combined density, momentum, energy, flux, and electric-field estimators.
- Maxwellian subtraction from the full-particle Fourier representation.
- Core/tail partitioning and tail-weight adjustment.
- Quadratic selection of new effective particle weights.
- The stability constraint applied after adaptive weight selection.
- Linearized Coulomb, Delta-M suppression, alternate projection, and
  Coulomb-plus-BGK paths.

For each algorithm, record its units, valid parameter ranges, degenerate-case
behavior, conservation invariants, and known defects or ambiguities in the
comparison implementation.

Create small deterministic characterization cases where the comparison code
can be isolated safely. Do not include the broken Rosenbluth branch in this
feature effort.

**Checkpoint:** every feature has a written formula, an owning current module,
defined degenerate-case behavior, and at least one planned test.

## Phase 2: Add typed configuration and metadata

Extend `SimulationTypes` with explicit modes:

```cpp
enum class HdpCouplingMode {
    Decoupled,
    VarianceWeighted,
};

enum class EffectiveWeightPolicy {
    Fixed,
    QuadraticAdaptive,
};

enum class CollisionCoupling {
    Standard,
    Linearized,
};

enum class ProjectionMode {
    FullMicroMacro,
    MaxwellianOnly,
};

enum class DeltaMMode {
    Enabled,
    Disabled,
};
```

Extend `ParaClass` with:

- HDP coupling mode.
- Effective-weight policy.
- Collision-coupling and projection modes.
- Delta-M mode.
- Partial-resampling enable flag and cutoff.
- Weighted-moment conservation flag.
- Minimum and maximum signed/full effective weights.
- CPU-cost model coefficients.
- Optional BGK strength following a Coulomb step.

Validate configuration at the boundary. Effective weights must be finite and
positive; cutoffs, cost coefficients, and BGK rates need documented ranges.

Record every setting in `RunMetadataOutput`. Add CLI options after internal
configuration and tests are stable.

**Checkpoint:** existing defaults produce unchanged fixed-seed results and
metadata unambiguously records every selected numerical mode.

### User-selectable configuration file

Add a checked-in `config/negpar.example.json` and support loading a run with:

```powershell
negpar_inhomo --config config/my-run.json
```

Use JSON for the public configuration format so the C++ executable and the
existing web UI can share the same representation. Use a maintained, pinned
JSON library rather than writing a partial JSON parser.

The example configuration should make every new feature visible and should
default to the current behavior:

```json
{
  "schema_version": 1,
  "features": {
    "weighted_hdp": false,
    "weighted_fourier_resampling": false,
    "partial_resampling": false,
    "adaptive_effective_weights": false,
    "linearized_coulomb": false,
    "delta_m_sampling": true,
    "projection_mode": "full_micro_macro",
    "coulomb_bgk_hybrid": false
  },
  "resampling": {
    "partial_cutoff_standard_deviations": 3.0,
    "conserve_weighted_moments": false,
    "signed_weight_min": 1e-7,
    "signed_weight_max": 5e-3,
    "full_weight_min": 1e-7,
    "full_weight_max": 5e-3,
    "cpu_cost_constant": 0.205,
    "cpu_cost_collision_coefficient": 3.277
  },
  "collisions": {
    "bgk_strength": 0.0
  }
}
```

The loader must:

- Require and validate `schema_version`.
- Reject unknown keys so misspelled feature names cannot be silently ignored.
- Apply compiled defaults first, configuration-file values second, and explicit
  CLI overrides last.
- Report validation errors with the JSON key path and invalid value.
- Validate feature dependencies. For example, partial or adaptive resampling
  must fail clearly if its required weighted reconstruction mode is disabled.
- Convert public strings and booleans into the typed internal configuration
  from this phase.
- Resolve relative paths consistently relative to the configuration file.
- Write the fully resolved effective configuration into the run output
  directory, including values inherited from defaults.

Add:

- `--validate-config <file>` to check a file without starting a simulation.
- `--print-effective-config` to show the fully resolved configuration.
- Loader tests for valid files, defaults, overrides, unknown keys, invalid
  types, unsupported schema versions, dependency failures, and non-finite or
  out-of-range values.
- UI support to load, edit, validate, save, and submit the same JSON structure.

**Checkpoint:** users can select any valid feature combination in one file,
validate it before a run, and reproduce the run from the effective
configuration saved with its outputs.

## Phase 3: Implement variance-weighted HDP coupling

Create a focused component, tentatively `WeightedHdpCoupling`, that owns
variance calculations and convex combinations. It should calculate:

- Deviational estimator variance.
- Full-particle estimator variance.
- A safe blend weight:

  ```text
  omega = variance_deviational /
          (variance_deviational + variance_full)
 ```

  In the reference blend, omega multiplies the full-particle estimator and
  1 - omega multiplies the signed/deviational estimator.

- Blended density, momentum, energy, and flux values.

Define degenerate cases explicitly:

- Both variances zero: use a documented deterministic fallback.
- Deviational variance zero: select the full-particle estimate.
- Full variance zero: select the signed/deviational estimate.
- Clamp small roundoff excursions to `[0, 1]`.
- Reject non-finite inputs or results.

Integrate the component into:

- `ElectricFieldSolver::update` for the WHDP electric field.
- `MomentOperations::computeMacroChange` for combined flux and source terms.
- Energy and moment diagnostics intended to describe the combined estimator.

Do not duplicate the weight formula inside individual physics modules.

**Acceptance criteria:**

- Decoupled mode remains numerically unchanged.
- Weighted results select the correct estimator at `omega = 0` and
  `omega = 1`.
- Intermediate results are convex combinations within floating-point
  tolerance.
- Empty cells and one-sided populations cannot produce NaN or Inf.

## Phase 4: Add weighted Fourier resampling

Extend the Fourier-resampling interface to consume positive, negative, and full
particles; Maxwellian moments; signed and full effective weights; spatial cell
width; and the partial-resampling cutoff.

Implement separately testable stages:

1. Normalize particle velocities to the Fourier domain.
2. Compute separate positive, negative, and full-particle coefficients.
3. Compute and subtract the Maxwellian Fourier contribution from the full
   representation.
4. Calculate a finite, frequency-dependent optimal weight.
5. Blend the full and signed coefficients.
6. Reconstruct signed or full particles at the requested target weight.
7. Restore velocity coordinates and assign valid spatial positions.
8. Enforce requested conservation constraints.
9. Return diagnostics, including the zero-frequency weight and sampling
   acceptance counts.

Use a result object rather than mutating the source group. Distinguish success,
a rejected reconstruction, and a sampling-attempt-budget failure.

**Acceptance criteria:**

- Exact fixed-seed replay under the reference configuration.
- Every frequency weight is finite and within `[0, 1]`.
- Reconstructed velocities and positions are finite and in bounds.
- Signed mass, momentum, and energy satisfy documented tolerances.
- Input particles remain unchanged after failure.
- Exact and approximate Fourier modes receive direct coverage.

## Phase 5: Add partial core/tail resampling

Create a typed particle-partition operation instead of embedding partition
logic directly in `ParticleResampling`.

- Define the core using local Maxwellian velocity, temperature, and the
  configured standard-deviation cutoff.
- Put particles outside the cutoff in typed positive, negative, and full tail
  groups.
- Resample only the core.
- Adjust tail multiplicity or weights consistently when an effective particle
  weight changes.
- Merge the untouched tail after successful reconstruction.
- Preserve particle kinds and deterministic ordering where required.

**Acceptance criteria:**

- Disabling partial resampling reproduces whole-cell resampling.
- Tail velocities are unchanged by core reconstruction.
- Tail positions remain inside the owning spatial cell.
- Conservation is checked before partitioning, after core reconstruction, and
  after the tail merge.
- Empty-core and empty-tail cases are supported.

## Phase 6: Implement adaptive effective-weight selection

Create an `EffectiveWeightSelector` responsible only for choosing target
signed and full weights.

Port the quadratic policy with:

- The three legacy operating points.
- Quadratic interpolation and tangent selection.
- Configured minimum and maximum clamps.
- The per-cell stability constraint.
- Explicit handling for zero signed particles, zero full particles, zero
  density, coincident interpolation points, and a degenerate quadratic.
- Rejection of non-finite or non-positive results.

Do not let this component resample particles or mutate the grid.

Update `ParticleResampling::synchronizeCoarse` to perform a transaction:

1. Collect per-cell estimator data required by the policy.
2. Calculate proposed signed and full effective weights.
3. Apply configured bounds and the stability constraint.
4. Reconstruct every cell into temporary particle groups.
5. Validate per-cell and global conservation, counts, and finite values.
6. Commit all reconstructed cells and both grid weights together.
7. Restore the pre-synchronization RNG state and leave simulation state
   unchanged if the operation fails, unless a documented retry policy is used.

**Acceptance criteria:**

- Fixed policy reproduces current synchronization behavior.
- Adaptive weights remain finite, positive, and within configured bounds.
- A failed cell cannot cause a partial grid update.
- Total full-particle mass is preserved when the full weight changes.
- Successful synchronization is deterministic for a fixed seed and thread
  count.
- Tests cover degenerate quadratic and zero-count cases.

## Phase 7: Add optional collision and projection modes

Implement these after weighted resampling and adaptive synchronization are
stable:

- Move the linearized positive/negative-to-Maxwellian collision alternative
  into `NegativeParticleCollisions`.
- Make Delta-M collision-source sampling configurable.
- Add the alternate projection-coefficient calculation to
  `ProjectionSampling`.
- Support optional BGK relaxation after a Coulomb step with an explicit,
  documented operator order.

Keep these choices independent of WHDP so each can be tested with decoupled and
variance-weighted coupling.

**Acceptance criteria:**

- Every mode has a direct unit or component test.
- Coulomb-only and BGK-only defaults remain unchanged.
- Hybrid mode runs Coulomb and BGK exactly once in the documented order.
- Disabling Delta-M changes only collision-source particle creation.
- Linearized and standard collision paths satisfy their intended conservation
  properties.

## Phase 8: Expose runtime controls and diagnostics

After the numerical APIs are stable, implement the Phase 2 JSON loader and add
CLI overrides for commonly changed modes and parameters. Update usage text,
validation errors, run metadata, and UI support. Keep the JSON file as the
complete interface; CLI flags are convenience overrides rather than a second,
independent configuration model.

Add diagnostic histories for:

- Signed and full effective particle weights.
- Minimum, maximum, and mean blend weight.
- Resampling attempts, accepted samples, and rejected reconstructions.
- Core and tail particle counts.
- Conservation error before and after synchronization.
- Synchronization time and particle-count change.

Keep all artifacts inside the configured output directory and use the existing
output-path validation.

**Checkpoint:** a run can be reproduced from its recorded options, and the
diagnostics explain every weight change or rejected synchronization.

## Phase 9: Verification matrix

### Unit tests

- Variance and blend-weight calculations.
- Weight boundaries and zero-variance behavior.
- Frequency-dependent coefficient blending.
- Quadratic interpolation and tangent selection.
- Effective-weight clamps and stability constraints.
- Core/tail partitioning and merging.
- Tail-weight adjustment.
- Configuration validation and metadata serialization.

### Component tests

- Weighted electric-field construction.
- Weighted macroscopic-change calculation.
- Weighted Fourier reconstruction in exact and approximate modes.
- Partial resampling with conservation enforcement.
- Adaptive synchronization success and rollback.
- Each collision and projection mode.

### End-to-end configurations

- Current decoupled HDP baseline with fixed weights.
- Variance-weighted HDP with fixed weights.
- Variance-weighted HDP with adaptive weights.
- Whole-cell and partial resampling variants.
- PIC regression.
- Coulomb-only, BGK-only, and Coulomb-plus-BGK runs.
- Single-thread deterministic reference runs.
- OpenMP fixed-schedule replay and multithreaded invariant tests.

For every numerical comparison, document whether it is exact, tolerance-based,
invariant-based, or statistical. Explicitly reject unexpected NaN and Inf.

Run Debug, Release, sanitizer, and applicable thread-sanitizer builds. Retain a
hard upper bound on rejection-sampling attempts and add a configurable guard
against unbounded particle growth.

## Phase 10: Numerical evaluation and rollout

Compare current and new algorithms using representative Landau-damping runs.
Record:

- Electric-field trajectory.
- Total and kinetic energy drift.
- Signed and full particle counts.
- Effective-weight history.
- Resampling frequency and acceptance rate.
- Conservation error at every synchronization.
- Runtime and peak memory as benchmark information.

Roll out in this order:

1. Merge configuration and inactive policy types with current defaults.
2. Merge weighted coupling behind an opt-in mode.
3. Merge weighted Fourier reconstruction behind an opt-in mode.
4. Merge partial resampling behind an opt-in flag.
5. Merge adaptive weight selection behind an opt-in policy.
6. Merge collision and projection alternatives individually.
7. Add CLI and UI exposure.
8. Review numerical evidence before considering any default change.

WHDP and adaptive effective weights remain opt-in until the full verification
matrix passes and their numerical behavior has been reviewed.

## Recommended implementation sequence

1. Behavioral specification and characterization fixtures.
2. Typed configuration and metadata.
3. Shared variance-weighted coupling component.
4. Electric-field and macroscopic-change integration.
5. Weighted Fourier coefficient construction and reconstruction.
6. Partial core/tail resampling.
7. Adaptive effective-weight selector.
8. Transactional synchronization.
9. Optional collision and projection modes.
10. Runtime controls, diagnostics, end-to-end evaluation, and rollout.

## Review gates

- Physics-sensitive changes require focused invariant results and a reference
  comparison in the review description.
- Mechanical ownership changes must be separate from numerical-formula changes.
- No adaptive synchronization change may merge without rollback tests.
- No stochastic path may merge without fixed-seed replay coverage and a finite
  sampling-attempt budget.
- No default may change until decoupled baseline tests and the complete
  weighted-mode verification matrix pass.
- Broken or unused comparison-folder code, including Rosenbluth initialization
  and diagnostic-only helpers, is outside scope unless separately specified and
  tested.
