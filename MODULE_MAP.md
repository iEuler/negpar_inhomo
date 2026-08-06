# Module Map

This is the current ownership map used before splitting the legacy headers.

| Area | Current owner | Intended owner |
| --- | --- | --- |
| Particle and grid data | `Classes.h/.cpp` | `particles/` |
| FFT wrappers | `FFT.h/.cpp` | `numerics/fft` |
| 1-D finite-difference and kinetic Euler numerics | `Numerics.h/.cpp` | `numerics/utilities` |
| Macroscopic moments and Maxwellian updates | `Moments.h/.cpp` | `physics/moments` |
| Grid initialization | `Initialization.h/.cpp` | `io/initialization` |
| Advection | `Advection.h/.cpp` | `physics/advection` |
| Poisson solve and field updates | `ElectricField.h/.cpp` | `physics/electric_field` |
| Binary Coulomb collisions | `Collisions.h/.cpp` | `physics/collisions` |
| Negative-particle collisions | `NegativeParticle.h/.cpp` | `physics/collisions` |
| Legacy Fourier resampling | `LegacyResampling.h/.cpp` | `resampling/legacy` initially |
| New resampling experiments | `Resampler*.h/.cpp`, `InhomoResampler.*` | `resampling/experimental` until dependencies are resolved |
| Generic macro, grid, field, distribution, particle, homogeneous, and parameter output writers | `Output.h/.cpp` | `io/output` |
| Particle-count diagnostics | `Diagnostics.h/.cpp` | `diagnostics` |
| Legacy header compatibility | none retained in `src/` | callers must include the owning module header |
| Runtime options | `RunOptions.h/.cpp` | `simulation/configuration` |
| Mutable counters, flags, and timing | `SimulationState` in `_global_variables.h` | `simulation/state` |
| Macro-step and output orchestration | `Simulation.h/.cpp` | `simulation/` |

The extracted `Initialization.cpp`, `NegativeParticle.cpp`, and
`LegacyResampling.cpp` files are the only active numerical path. The former
`coulomb_*.h` compatibility umbrellas and forward-declaration header were
unused by production and test targets and have been removed. Public callers
should include the focused owning headers directly (`Initialization.h`,
`NegativeParticle.h`, `LegacyResampling.h`, `Output.h`, and so on).

Mutable run state is now owned by the caller as a `SimulationState` object.
Advection, collision-step, resampling, initialization, timing, and output
entry points receive that object explicitly. No process-wide aliases or
mutable state instances remain. Diagnostics are deliberately state-free pure
reductions because they do not consume or mutate run state.

Simulation orchestration now lives in `Simulation.cpp`; the historical
`inhomo_neg_coulomb.cpp` file is intentionally retained as an empty compatibility
translation unit for older project layouts. It is not part of any build target
and contains no implementation or global state.

The obsolete `coulomb2` implementation and its duplicate class/global files
are no longer present. A repository-wide dependency search found no remaining
`coulomb2` callers.

The newer `Resampler*.cpp`, `InhomoResampler.cpp`, and `ResamplerHelper.cpp`
experiments remain outside the application target. Their public classes and
helpers now live in `coulomb::experimental`, while calls into the established
legacy implementation are explicitly qualified. CMake exposes them as the
opt-in `negpar_experimental_resampler` library and, when enabled, builds
`negpar_experimental_tests` to verify symbol isolation and explicit RNG
reproducibility. They remain opt-in until an end-to-end caller and full
conservation characterization are approved.
