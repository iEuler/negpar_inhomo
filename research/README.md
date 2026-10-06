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
- Implement and independently validate the paper's proposed two-term source
  sampler. The executable currently uses the legacy source sampler.
- Certify source rejection envelopes and test adaptive failure rollback.
- Recover provenance for historical figures or replace them with new results.

The Release test suite currently contains 76 cases. CTest runs the entire suite
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
