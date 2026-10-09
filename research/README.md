# Research validation checkpoint

The original manuscript is preserved. Its revised counterpart is
`../../../paper/NegPar_mixing_revision.tex` (from this directory).
This checkpoint restores a runnable solver and creates reproducible evidence;
the paper still requires research validation before submission.

Run from the repository root after a Release build:

```powershell
python research/benchmark.py --executable build/release/Release/negpar_inhomo.exe --output research/runs/new_pilot --seeds 3 --steps 100 --cells 40
python research/report.py research/runs/new_pilot/summary.json
python research/refinement.py research/runs/new_pilot research/runs/new_refinement
python research/conservation_probe.py research/runs/new_pilot research/runs/new_refinement research/runs/new_conservation
```

The default pilot covers amplitudes 0.05 and 0.4, collision strengths 1 and
10, PIC, HDP, mixed and adaptive modes, two input weights and three seeds.
It measures error in the trajectory of integral E squared against independent
fine-PIC runs on the same grid and timestep. Those references have sampling
and discretization error. Equal input weights do not establish equal accuracy.
Wall times include initialization and output; run this sequentially on an
otherwise idle machine for publishable timing measurements.

Each experiment retains a frozen executable and DLLs, executable SHA-256,
source archive, inputs, metadata, console logs, numerical histories and
measurements. The output directory must be new. Generated runs are ignored by
Git but remain on disk. `pilot_v1` contains a failed diagnostic experiment;
`pilot_v2` predates the density normalization and diagnostic corrections.
Use `pilot_v3` for the current pilot.

Numerical corrections include separate positive/negative Fourier variances,
local collision-background density, degenerate pair scattering, Maxwellian
support during full reconstruction, realized background-count checks,
source/output weight separation, bounded inverse-weight allocation and
physical field-energy normalization. Histories include the final state and
explicit times, masses and weights.

Outstanding publication gates:

- Refine reference particle count, timestep, spatial grid and Fourier cutoff.
- Measure momentum and energy drift over longer simulations and synchronization.
- Exercise non-Maxwellian initial data and compare cost at overlapping errors.
- Validate the observable covariance and adaptive-weight bias assumptions.
- Independently validate the repaired two-term source sampler across broader
  parameter regimes; the exact split repair and its checks are described below.
- Certify source rejection envelopes and test adaptive failure rollback.
- Recover provenance for historical figures or replace them with new results.

The Release numerical test suite currently contains 82 cases. CTest runs the entire suite
as well as separately registered cases and seeded integration checks. Debug
and sanitizer validation must be recorded separately; Release success does not
imply those configurations were tested.

## Current measured outcome

All 108 `pilot_v3` runs completed with finite histories. PIC is faster in
these small tests; adaptive mixing does not consistently improve measured
field error. Mean maximum energy drift across pilot cases is approximately
0.0024%-0.038% for PIC and 0.54%-6.65% for the HDP variants. These results
require revising the historical uniform-efficiency claim.

The one-seed `refinement_v1` study includes exact baseline replay and separate
grid, timestep, Fourier and particle-count sensitivity tests. The default
adaptive 80-cell run failed at step 88 with invalid Maxwellian temperature.
`conservation_v1` enables signed-moment conservation: that case completes,
but maximum energy drift is still 2.55%. Neither the probe nor these
stochastic sensitivity tests establish convergence or justify promoting a
default. Numerical failures are retained alongside successful data.

Use `build/research-python/Scripts/python.exe` for plotting on this machine;
the global Matplotlib installation is incompatible with its NumPy version.
The benchmark and refinement runners themselves use only the standard library.

## Frozen-distribution mixing experiment

```powershell
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_frozen_mixing.py -v
build/research-python/Scripts/python.exe research/frozen_mixing.py --output research/runs/new_frozen
build/research-python/Scripts/python.exe research/frozen_mixing_plot.py research/runs/new_frozen/summary.json
```

This isolates scalar mixing for a positive Gaussian mixture with a known
Maxwellian baseline and signed Gaussian difference. It measures velocity,
second moment and a cosine mode over four mixture fractions. Exact means
and variances provide truth; 4,000 independent pilot replicas choose each
observable's weight and 20,000 fresh replicas evaluate it. Full and signed
sample streams are independent. The frozen baseline does not match the
mixture's moments; this is a controlled estimator test, not a macro-micro
closure test. NumPy and Matplotlib are installed in the isolated environment.

Two comparisons distinguish reuse of already-generated representations from
generating both components under an equal total Gaussian-draw budget. Raw
replicas, covariance, bias, MSE, measured runtime, exact predictions and the
executed source are retained. Mixing reduces variance for existing dual
representations, but does not beat the better single estimator at fixed draw
budget in this model. Pilot estimation costs are reported separately and also
charged in the from-scratch MSE-time efficiency metric. This experiment does
not validate the C++ solver's count-based weights or evolving dynamics.

## Homogeneous Coulomb mixing before resampling

Build the separate research target after configuring Release:

```powershell
cmake --build --preset release --target negpar_homogeneous --parallel 4
build/research-python/Scripts/python.exe research/homogeneous_experiment.py --executable build/release/Release/negpar_homogeneous.exe --output research/runs/new_homogeneous
build/research-python/Scripts/python.exe research/homogeneous_plot.py research/runs/new_homogeneous/summary.json
build/research-python/Scripts/python.exe research/homogeneous_proxy_analysis.py research/runs/new_homogeneous
build/research-python/Scripts/python.exe research/homogeneous_refinement.py research/runs/new_homogeneous research/runs/new_homogeneous_half_dt
```

`homogeneous_mixing.cpp` evolves a single homogeneous 3D velocity population
using the existing C++ Coulomb pair operator and signed source sampler. It
calls no spatial advection, projection, synchronization, adaptation or
resampling. The initial distribution is a mixture of the unit Maxwellian and
an anisotropic Gaussian with diagonal temperatures (2, 0.5, 0.5). Both have
the same exact invariant mass, mean and total energy, so the Maxwellian
baseline remains fixed. The full population remains constant; signed counts
can only grow. The harness fails if the collision-background constraint is
exhausted rather than invoking resampling.

The default experiment uses fractions 0.05 and 0.3, collision strength 5,
dt=0.01 and 20 steps. Each HDP trajectory starts with 2,048 full particles and
128 of each sign. Independent 256-replica pilot ensembles fit time- and
observable-specific weights; 1,024 fresh trajectories evaluate mixing.
Independent 64-replica PIC references use 32,768 particles, with a second
64-replica reference at half the timestep. Exact initial expectations and
zero-collision/same-seed replay controls verify initialization and the harness.

`homogeneous_v1` and `homogeneous_dt_half_v1` are superseded: the source
sampler read cached signed counts that the harness never populated, making
the rejection denominator zero. Their old drift and mixing measurements are
preserved with validity notices. The fix uses live particle-list sizes.
Corrected data are in `source_conservation_corrected_v1`. Mean mass and
total-second-moment drifts are consistent with zero at the tested precision
at both timesteps; individual trajectories still fluctuate. A held-out recheck
gives count-proxy final-MSE gains 1.090--1.265 (about 8%--21% error reduction),
with all six individual bootstrap intervals above one. Covariance-aware
mixing gains range 1.085--1.423. These are gains over the better component of
the same trajectories, conditional on the approximate PIC reference.

## Complete homogeneous cost at matched accuracy

```powershell
cmake --build --preset release --target negpar_homogeneous --parallel 4
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_matched_accuracy.py -v
build/research-python/Scripts/python.exe research/matched_accuracy.py --executable build/release/Release/negpar_homogeneous.exe --output research/runs/new_matched_accuracy
build/research-python/Scripts/python.exe research/matched_accuracy_plot.py research/runs/new_matched_accuracy/summary.json
```

The count sweep charges complete homogeneous initialization, bound setup,
signed-source sampling, signed/full collision evolution and observable
estimation/mixing. CSV output is excluded from compute timing; per-run wall
time including output is retained. The same eight diagnostic observables and
timing policy are used for PIC and HDP. No resampling or synchronization is
introduced. The reference uses known initial expectations as a control variate
to reduce its sampling uncertainty, with doubled-count and halved-dt checks.

The main metric is trajectory RMS over anisotropy, fourth moment and cosine
difference, normalized separately by their initial deviations from equilibrium.
Default fractions are 0.05, 0.01 and 0.3, with error thresholds 0.25 and 0.4.
Particle counts are searched separately for each method. Both ordinary
estimators and a temporal control variate based on exact initial moments are
tested; the stronger control variate is available to both PIC and HDP mixing.
Count-proxy mixing requires no fitted pilot weights.

Candidate selection uses the upper bootstrap accuracy bound. Selected counts
are rerun with fresh seeds, and unverified accuracy targets are reported
explicitly. Costs compare the cheapest tested configurations meeting a common
threshold, with achieved errors shown; they are not exact-equal-error or
globally optimal comparisons. Use `--resume` with identical settings after an
interruption only when existing run outputs completed successfully. Previous
experiments are preserved and new outputs must otherwise use a new directory.

The HDP trajectories in `research/runs/matched_accuracy_v1`, including its
confirmation, are superseded because of the source-cache bug. Pure-PIC
measurements are unaffected. Cross-regime count sweeps must be repeated;
their previous ratios are not current evidence. Selected epsilon=0.01 points
have been rechecked with a frozen corrected executable and fresh seeds:
ordinary PIC RMS 0.276 versus corrected HDP mixture RMS 0.283, compute times
2.444 versus 0.0515 seconds (about 47 times faster at closely matched errors).
With known initial moments available to both methods, PIC RMS is 0.172 at
0.575 seconds; corrected HDP RMS is 0.145 at 0.00997 seconds. This second
comparison achieves lower error for HDP rather than exactly equal error.
These are selected-point rechecks, not new count optimization. See
`source_conservation_corrected_v1/REPORT.md` for all uncertainties and counts.

## Signed-source conservation audit

`source_conservation.cpp` isolates source insertion and signed velocity
updates in paired fresh/stale-cache replicas. `source_conservation.py` freezes
the diagnostic executable and saves stage measurements and C++ kernel values.
`source_kernel_audit.py` cross-checks a vectorized kernel translation against
C++ and integrates around the source singularity. `source_conservation_corrected.py`
freezes the corrected homogeneous executable and reruns trajectories;
`source_conservation_report.py` summarizes them and marks old results.

The dominant previous drift came from the zero rejection denominator.
Kernel mass-integral defects are only around 1e-6 in the tested source cases.
A separate algebraic split discrepancy remains: the low branch uses
q*1(q<a) and the high branch max(q-a,0), whose sum omits a*1(q>=a).
Empirical envelopes, finite support and proposal-density normalization also
need an audit before a global consistency or long-time conservation claim.
No moment projection or resampling was introduced to conceal source errors.

## Exact signed-source split repair

The follow-up repair replaces the low/high legacy split with
`bounded = clamp(q, -a, a)` and `remainder = q - bounded`. The sum is exactly
q, including the omitted positive plateau and any negative spill beyond the
empirical cap. Reversing a background particle's sign flips the sampled
remainder sign as well. The bounded branch's proposal count and acceptance
envelope both use Np+Nn, giving a valid bound even when both populations
contribute with the same sign after subtraction. Both branches normalize by
the collision-background density rhoF; the previous remainder proposal count
used the signed reconstruction's density instead.

Acceptance ratios outside [0,1] now stop the run rather than silently clipping
an invalid envelope. Tiny rounding excesses within 1e-10 are allowed. The
bounded envelope is rigorous for the clamped target; the remainder envelope
is still empirical and finite radial support remains a truncation. This repair
does not certify every parameter regime or force sampled invariants exactly.

`source_split_experiment.py` archives 30,000 samples per source speed/sign,
checks mass and sign reversal, and runs fresh no-resampling trajectories at
epsilon 0.3 and 0.05 with dt 0.01/0.005, plus an epsilon 0.01 efficiency point.
The frozen protocol and measurements are under `research/runs/source_split_v1`.
Earlier cache-only corrected runs remain valid measurements of their frozen
version, rather than measurements of the new split.

The split-repair validation completed with bound checks enabled. Single-source
mean mass changes at speeds 1 and 2 are -0.000018 ± 0.000057 SE and
-0.000047 ± 0.000045 SE (30,000 replicas per sign). Sign reversal is exact
for sampled particle lists and mass, with summed velocity moments differing
only by roundoff. Some trajectory invariant changes are 2--2.5 SE from zero;
these finite ensembles do not establish exact conservation expectations.
The repaired epsilon=0.01 ordinary mixture gives RMS 0.279 versus PIC 0.276,
at 0.0700 versus 2.444 seconds, about 35 times lower measured compute cost.
The extra bounded-proposal work is included. The new report retains all drift
components, confidence intervals and frozen-source provenance.

## Repaired-source count sweep across perturbation sizes

`research/runs/matched_accuracy_split_v2` extends the repaired homogeneous
benchmark to epsilon 0.05 and 0.3. It tests 32/64/128 particles of each sign,
full-to-per-sign count ratios 8/16, and independent PIC count searches at a
common normalized trajectory RMS threshold of 0.4. Candidate ensembles use
96 replicas, selected allocations use 192 fresh replicas, and references use
128 replicas with doubled-count and halved-timestep checks. All runs completed
with rejection checks enabled; the frozen executable and source are retained.

The epsilon 0.05 ordinary HDP and CV PIC selections missed their fresh upper
accuracy bounds. These failures remain in the main report. The supplemental
`confirmation/REPORT.md` uses 384 new replicas per allocation. At epsilon 0.3,
PIC with 768 particles gives RMS 0.242 [0.229, 0.253], versus HDP with 512 full
and 32 of each sign at 0.244 [0.227, 0.259]. Mean complete compute times are
0.00754 and 0.01563 seconds: PIC is about 2.1 times faster at close errors.
At epsilon 0.05, conservative ordinary allocations give PIC RMS 0.319 at
0.06019 seconds and HDP RMS 0.264 at 0.08909 seconds. Both pass the threshold,
but their achieved errors differ, so this does not establish an equal-error
winner. Small timing differences need repeated measurements. No uniform HDP
efficiency claim is supported by the new sweep.

`split_sweep_analysis.py` runs these supplemental confirmations and quantifies
incremental mixing with paired trajectory/reference bootstrap intervals. In
the fresh ordinary confirmations, mixing reduces aggregate normalized MSE
relative to the better component of the same HDP populations by about 13%
at epsilon 0.05 (gain 1.149 [1.115, 1.184]) and 25% at epsilon 0.3
(gain 1.325 [1.203, 1.458]). These are individual intervals conditional on the
approximate reference, and compare already-generated components. They do not
measure the cost of a standalone signed-only solver. All five signed invariant
drift means and standard errors are retained in `confirmation/summary.json`.

Reproduction (the analysis reuses completed supplemental runs):

```powershell
build/research-python/Scripts/python.exe research/matched_accuracy.py --executable build/release/Release/negpar_homogeneous.exe --output research/runs/new_split_sweep --epsilons 0.05 0.3 --replicas 96 --large-replicas 64 --validation-replicas 192 --reference-replicas 128 --sign-counts 32 64 128 --full-ratios 8 16 --pic-counts 128 512 2048 8192 32768 --targets 0.4
build/research-python/Scripts/python.exe research/split_sweep_analysis.py research/runs/new_split_sweep
build/research-python/Scripts/python.exe research/matched_accuracy_plot.py research/runs/new_split_sweep/summary.json
```

## Logarithmic epsilon efficiency curve

`epsilon_efficiency.py` generates R(epsilon) = HDP-plus-mixing compute time
divided by ordinary PIC compute time, with epsilon sampled logarithmically
between 0.01 and 0.3 (seven points by default). The common normalized trajectory
RMS target is 0.4. The same homogeneous initial condition, dt=0.01, T=0.2 and
collision strength 5 are used at every epsilon. Known initial moments are
used only to improve the independent PIC reference, not the competitors.

Each method's count allocation is tuned using 48 pilot replicas; selection
requires its upper bootstrap error bound to meet 90% of the final target.
HDP searches signed counts 32/64/128/256 and full ratios 8/16. PIC counts are
proposed by inverse-square scaling and measured, with a larger fallback if
none passes. The selected allocation is evaluated with 192 fresh trajectories
in three sequential timing batches, alternating method order. Bootstrap ratio
intervals resample batches and trajectories. Only three timing batches are
available; the intervals do not cover every machine-load effect. Accuracy
intervals are conditional on the chosen allocations and approximate reference.

The two-panel `R_epsilon.png` and vector `R_epsilon.svg` show the ratio and
the achieved PIC/HDP errors, all against logarithmic epsilon. Unverified
targets are marked explicitly. This is a common-threshold comparison over a
finite count search, not exact-equal-error optimization. A proportional-to-epsilon
line is included only as a visual guide, not a fitted or proven scaling law.
Source/executable provenance, raw trajectories, reference checks and pilot
measurements are retained. No resampling or invariant correction is applied.

```powershell
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_epsilon_efficiency.py -v
build/research-python/Scripts/python.exe research/epsilon_efficiency.py --executable build/release/Release/negpar_homogeneous.exe --output research/runs/new_epsilon_efficiency
```

For interrupted completed subruns, repeat the identical command with `--resume`.
It rejects a changed protocol and preserves existing outputs. A failed partial
subrun is retained and must be investigated before resuming. Defaults can be
overridden with `--epsilon-min`, `--epsilon-max`, `--points`, `--target`,
`--pilot-replicas` and `--batch-replicas`.

The completed default experiment is `research/runs/epsilon_efficiency_v1`.
All seven fresh validations pass the common 0.4 threshold. In ascending epsilon
(0.01, 0.017627, 0.031072, 0.054772, 0.096549, 0.170190, 0.3), the measured
R values are 0.056, 0.110, 0.379, 0.593, 4.219, 2.809 and 3.693. This supports
the qualitative near-equilibrium advantage, with PIC favored by a few times
at larger perturbations in this count search. It does not establish a clean
linear or asymptotic law. The unequal achieved errors at some points and
coarse, statistically selected count allocations visibly affect the curve.
The final PNG was visually inspected; the SVG is available for reuse.

## Fixed-envelope resampling prototype

`FourierResamplerConfig::envelope` now supports an experimental
`ResamplingEnvelope::CertifiedQuadratic` mode. `LegacyAdaptive` remains the
default, including fixed-seed behavior. The new mode is currently available
through the C++ research API, not promoted to solver runtime configuration.
It requires quadratic reconstruction (`useApproximation=true`).

For normalized cell offsets |delta_i| <= h, the fixed envelope is

```text
B = |f| + h (|fx| + |fy| + |fz|)
        + h^2 [ (|fxx| + |fyy| + |fzz|)/2 + |fxy| + |fxz| + |fyz| ]
```

A small roundoff margin is added; nonfinite derivatives, overflow and envelope
violations are rejected. Each cell draws its proposal count once, with mean
B * cell_volume / output_weight, and accepts with probability |q|/B. Assigning
the sign of q therefore recovers the restricted piecewise quadratic source
in expectation. The existing spherical support restriction is retained.
This guarantee is not for the original empirical distribution or exact
Fourier interpolant, and does not enforce sampled invariants.

The legacy sampler uses neighboring grid values, multiplied by 1.5, as an
empirical envelope. When an interior sample exceeds it, the bound grows,
accepted particles are thinned, and the remaining proposal budget changes.
The new mode removes that adaptation and its need for a consistency argument.
Diagnostics record proposal attempts and envelope increases for both paths.

The frozen-population audit uses 512 particles per sign, epsilon=0.05,
output weight four times input weight, cutoffs 4/8/12, and 512 replicas per
mode/cutoff. Mode order alternates within each replica. No collisions or moment
projection are applied. `resampling_envelope_v1` preserves the initial prototype;
`resampling_envelope_v2` measures the final version after skipping an unused
legacy grid-envelope calculation. Scripts preserve executable/source provenance,
raw samples, all six moment changes and paired differences. This is an isolated
resampling study, not a repeated-resampling Landau damping validation.

```powershell
cmake --build --preset release --target negpar_tests negpar_resampling_probe --parallel 4
build/research-python/Scripts/python.exe research/resampling_experiment.py --executable build/release/Release/negpar_resampling_probe.exe --output research/runs/new_resampling
```

Validation: 82 Release numerical cases (7,264 assertions), all 54 Release CTest
checks including reference/synchronization, plus Debug and MSVC AddressSanitizer
probe runs at all three cutoffs (eight replicas per mode/cutoff, 48 finite rows
per configuration). Debug and sanitizer probes were configured with
`-DNEGPAR_BUILD_TESTS=OFF`; those are not full Debug/sanitizer unit suites.
Use `-DNEGPAR_BUILD_TESTS=ON` when configuring those suites later. Their raw
outputs and logs are under `research/runs/resampling_validation_v1`.

Final measured result (`resampling_envelope_v2`): 8,353 legacy envelope
increases across all replicas versus zero for the certified mode. Mean times
at cutoffs 4/8/12 are 7.282/16.975/46.310 ms for legacy and
7.030/17.570/45.183 ms for certified. These small, mixed differences do not
establish a uniform speed benefit. Anisotropy source-to-output RMSE is
0.011442/0.013775/0.014941 for legacy and 0.010398/0.013532/0.016383 for
certified, also mixed. Reliability of the envelope improved; physical
accuracy and conservation are not established by that improvement.

The next geometric audit was run as `resampling_geometry_v1`. Centered cells at
nodes 0, dx, ..., 2*pi-dx cover [-dx/2, 2*pi-dx/2] in each coordinate. The
legacy spherical mask therefore omits positive-side caps. Certified periodic
wrapping restores those caps; a 100,000-proposal geometric test matches the
analytic cap-volume difference and the wrapped sphere acceptance fraction.
Across 512 paired replicas at cutoffs 4/8/12, wrapping reduced the positive
v-squared bias by 0.00271, 0.000558, and 0.000543 (SE 0.000298, 0.000148,
0.000139). Mass moved closer to zero at cutoffs 4 and 8 but slightly farther
at 12. Anisotropy RMSE worsened at 4, was nearly unchanged at 8, and improved
slightly at 12; timings were mixed and wrapping was not uniformly faster.
This is evidence for a
support-geometry correction, not a general accuracy or efficiency guarantee.
Keep it experimental and opt-in. Tail/core preservation, moment projection,
and repeated-resampling dynamics remain untested. Results, raw rows, frozen
source and executable are under `research/runs/resampling_geometry_v1`.

## Core/tail resampling audit

`resampling_tail_probe.cpp` and `resampling_tail_experiment.py` compare full
certified wrapped reconstruction with core radii 2.5, 3, and 4 thermal units,
using the same frozen population and Fourier cutoffs 4/8/12 as the geometry
audit. Every candidate is applied through five consecutive calls to expose
accumulated reconstruction error; the production count-reduction gate is
recorded but not enforced. This is an isolated resampling experiment without
collisions, moment projection, position reassignment, or weighted Fourier
coupling. No production defaults are changed.

The unchanged-weight arm retains all tail particles exactly. A separate 4x
weight arm uniformly thins tails on the first call, matching the existing
equal-weight partial-resampling path, then retains tails exactly thereafter.
This distinction avoids representing unchanged tails with an incorrect new
weight. Every total count includes the tails. Timings exclude audit-only
moment evaluations and include partition, reconstruction, thinning, and merge.

```powershell
cmake --build --preset release --target negpar_resampling_tail_probe --parallel 4
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_resampling_tail_experiment.py
build/research-python/Scripts/python.exe research/resampling_tail_experiment.py --executable build/release/Release/negpar_resampling_tail_probe.exe --output research/runs/new_tail_audit --replicas 128 --rounds 5
```

The runner fails early if plotting dependencies are incompatible, requires a
fresh output directory, archives source/binary provenance and initial source
velocities, checks exact unchanged-weight tail moments, and pairs replicas by
ID. It records all 20 observables: mass, three momenta, total v-squared,
anisotropy, x-fourth and radial-fourth moments, and real/imaginary parts of six
fixed physical low Fourier modes. First/final-call RMSE ratios use 2,000 paired
bootstrap resamples. These are accuracy/count/cost tradeoffs; different
normalization boxes and retained counts prevent interpreting the comparisons
as matched-accuracy efficiency claims.

The completed `research/runs/resampling_tail_v1` study has 128 replicas and
15,360 finite resampling calls. At cutoff 8 with fourfold output weight, a
3-sigma core reduces one-call anisotropy/radial-fourth RMSE by 25%/48%, with
3% extra particles and similar runtime, but increases low Fourier-mode RMSE
by 6%. Repeated calls still accumulate distortion. See `ASSESSMENT.md`,
`REPORT.md`, `repeated.png`, and `tradeoffs.png` in that archive.

A deterministic source-support audit exposes an additional geometry problem:
normalization by separate axis extrema makes the spherical mask an ellipsoid
in physical velocities, clipping points already selected into the spherical
core. At a 3-sigma cutoff it excludes 10 positive and 21 negative initial
particles; their signed mass and energy closely predict the measured core
drift. The next correction should test fixed Maxwellian-centered core bounds
before adding moment projection or promoting this path. The archived
`audit_source_support.py` reproduces this diagnostic from initial velocities.

## Fixed physical core bounds

`FourierResamplerConfig::fixedVelocityBounds` is an opt-in six-value vector
`[xmin,xmax,ymin,ymax,zmin,zmax]`. Empty retains the existing extrema-based
normalization. Fixed bounds require certified periodic sampling, finite positive
spans, and all particles used by reconstruction inside the box. Bounds are
applied to the resampler's copy and used for both normalization and restoration.
For a physical spherical core, set each range to `u_i +/- cutoff*sqrt(T)`;
then its normalized spherical mask is exactly the partition's physical sphere.
No runtime JSON/CLI or production default was changed.

`resampling_domain_experiment.py` uses the extended tail probe's `aligned`
study option to run seven methods together: full wrapped control, three
extrema-based cores, and three aligned cores. It rotates timing order, pairs
replicas by ID, checks unchanged-weight tail retention, and retains all 20
observable statistics. Since changing bounds also changes physical grid
spacing, the experiment measures the whole domain correction, not masking in
isolation. All candidate transformations are applied; count-gate fractions
are diagnostics rather than production trajectory acceptance rates.

```powershell
cmake --build --preset release --target negpar_tests negpar_resampling_tail_probe --parallel 4
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_resampling*experiment.py
build/research-python/Scripts/python.exe research/resampling_domain_experiment.py --executable build/release/Release/negpar_resampling_tail_probe.exe --output research/runs/new_domain_audit --replicas 128 --rounds 5
```

For full numerical/CTest execution in a restricted sandbox, set both `TEMP`
and `TMP` to a writable workspace directory such as `build/test-temp` first.
The full Release numerical suite and all CTest checks passed with that setting.

The completed `research/runs/resampling_domain_v1` archive contains 26,880
calls (128 replicas, five calls, seven methods). All 15,360 extrema/full control
observations exactly reproduce prior non-timing outputs. At unchanged weight,
cutoff 8 and radius 3, aligned bounds reduce one-call anisotropy RMSE 39%,
radial-fourth RMSE 21%, and low-mode RMSE 15%, with 2% extra particles and
similar runtime. Mass/v-squared mean drifts become consistent with zero in
that one-call case. Repeated calls and coarser weights increase variance and
particle counts; the option does not ensure conservation or uniform efficiency.
See `ASSESSMENT.md` for the decision and limitations and `REPORT.md` for all
paired comparisons. Both `domain_comparison.png` and `rmse_ratios.png` were
visually inspected. Validation: 85 numerical cases / 7,728 assertions, 57/57
Release CTest checks without exclusions, and five Python analysis tests.


## Stratified envelope-proposal allocation

`FourierResamplerConfig::proposalAllocation` defaults to `IndependentRounding`.
The opt-in `Stratified` mode requires `CertifiedQuadratic`, so every cell's
envelope and proposal expectation are fixed before proposals in that cell.
`StratifiedProposalAllocator` places an independent uniform point in each
unit interval of cumulative expected proposals and assigns points to cells
in lexicographic grid order. Fully covered strata need no explicit draw;
the allocator retains the boundary stratum's draw across adjacent cells.
The count in each cell remains unbiased, and every cumulative prefix count
is floor/ceiling of its expectation. This is independent stratification,
not a single shared systematic offset. Uniform positions and rejection
sampling remain unchanged; weights and support are unchanged.

For cell expectation lambda_j = M_j * volume_j / outputWeight, unbiased
counts and independent within-cell rejection preserve the expected signed
quadratic reconstruction restricted to the existing sphere mask. This does
not establish unbiasedness relative to the original empirical source or
exact Fourier function, and does not guarantee lower variance for every
signed observable. Allocation correlations and rejection noise matter.

`resampling_stratified_experiment.py` invokes the tail probe's `stratified`
option: full independent control, three aligned independent cores, and three
aligned stratified cores. It compares all 20 observables, paired bootstrap
RMSE ratios, counts and their SD, proposal-count SD, cost, and count-gate
fractions over repeated forced calls. The default 128-replica/five-call run
checks all independent controls against `resampling_domain_v1`. Different
sizes need a matching preceding control archive via `--previous`.

```powershell
cmake --build --preset release --target negpar_tests negpar_resampling_tail_probe --parallel 4
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_resampling*experiment.py
build/research-python/Scripts/python.exe research/resampling_stratified_experiment.py --executable build/release/Release/negpar_resampling_tail_probe.exe --output research/runs/new_stratified_audit --replicas 128 --rounds 5
```

Production defaults and JSON/CLI are unchanged. These isolated forced calls
have no moment correction, collisions, or production rejection/rollback;
they are not matched-accuracy efficiency measurements.


Completed archive: `research/runs/resampling_stratified_v1`, 26,880 calls,
15,360 reproduced non-timing controls. At radius 3/cutoff 8/4x weight, one-call
mass, v-squared, anisotropy, radial-fourth and Fourier RMSE ratios are
0.827/0.771/0.790/0.758/0.937, each with a pointwise bootstrap interval below
one; particle counts are essentially equal. Unchanged-weight and repeated
results are mixed. Radius 3/cutoff 4/4x weight worsens one-call radial-fourth
RMSE, ratio 1.245 [1.055,1.477]. No universal improvement or solver-efficiency
claim follows. See ASSESSMENT.md for the decision, REPORT.md for the sweep,
and ALGORITHM.md for the expectation argument. Validation: 86 numerical
cases/8,034 assertions, 58 CTest checks plus two rebuilt-production reference
rechecks, seven Python analysis tests, verified hashes and inspected figures.


## Bounded signed-moment correction (research API)

`BoundedMomentCorrection::apply` is a separate opt-in correction for a signed
core, with explicit target mass, three momenta and three diagonal second
moments. It does not replace the legacy production conservation routine.
It copies the core and commits only on success; every rejection leaves the
input particles and metadata unchanged. Failure may advance the RNG through
attempted excess-sign deletion; this API does not promise RNG rollback.

Mass is computed from signed counts and must be representable at the output
weight within a roundoff-scale count tolerance. The correction deletes
uniformly selected particles of the excess sign, rejecting if insufficient
particles exist. It never creates particles or changes weights. Each velocity
step is the minimum-norm solution of linearized momentum/diagonal-second
constraints over currently free particles. Coordinates are centered and
scaled by the physical core radius. Support-blocking particles are frozen and
the step is recomputed; backtracking requires lower residual and spherical
support. This is a local iterative correction, not a globally optimal
velocity projection. Degenerate systems and blocked/failed convergence are
rejected safely. RMS movement is limited to 0.1 times the radius by default,
with 40 iterations and count-normalized residual tolerance 1e-10.

The `corrected` tail-probe study compares three aligned stratified cores with
and without correction plus a shared full independent control. The whole
input's low moments are the target. After optional tail coarsening, tail output
moments are subtracted to obtain the core target; these tail particles remain
fixed. Failed correction leaves the uncorrected reconstructed core candidate
in place. Every candidate/fallback is then applied to expose accumulated error;
this is deliberately not production rollback of the entire resampling call.
All unconditional error statistics include correction failures. Constrained
moment accuracy alone cannot establish distribution accuracy; radial-fourth
and Fourier errors are separate unconstrained measurements. Cross second
moments, fourth moments and nonzero Fourier modes are not constrained.

`resampling_correction_experiment.py` records status, attempted removals,
iterations and displacement, all seven per-call target/output moments, and
20 original-source observables. It independently checks every successful
correction's measured moments and compares uncorrected controls with the
preceding stratification archive. Alternate ensemble sizes require matching
control archives via `--previous`. Target extraction, correction and merge
are timed; audit-only moments are excluded. Production defaults are unchanged.

```powershell
cmake --build --preset release --target negpar_tests negpar_resampling_tail_probe negpar_inhomo --parallel 4
build/research-python/Scripts/python.exe -m unittest discover -s research -p test_resampling*experiment.py
build/research-python/Scripts/python.exe research/resampling_correction_experiment.py --executable build/release/Release/negpar_resampling_tail_probe.exe --output research/runs/new_correction_audit --replicas 128 --rounds 5
```


Completed archive: `research/runs/resampling_correction_v1`: 26,880 finite
calls and 15,360 exact non-timing control reproductions. At radius 3/cutoff 8
after five calls, unchanged-weight radial-fourth and Fourier RMSE ratios are
0.197 [0.168,0.231] and 0.633 [0.592,0.678], with 12% fewer particles. At
4x weight they are 0.284 [0.234,0.341] and 0.654 [0.611,0.699], with 19%
fewer particles. Isolated resampling time increases about 4%. Error statistics
include correction failures. Across the sweep, 9,096 of 11,520 attempts succeed;
2,422 hit displacement limits and two hit iteration limits. Maximum measured
successful physical moment error is 3.06e-10. Repeated count growth remains,
and one low-cutoff setting slightly worsens Fourier error. These are
conditional source-reconstruction results, not accepted solver trajectories
or uniform efficiency claims. See ASSESSMENT.md, REPORT.md, ALGORITHM.md,
repeated.png, rmse_ratios.png and validation.json. Validation: 88 numerical
cases/8,636 assertions, 60 CTest checks and ten resampling analysis tests.
