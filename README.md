This project implements the negative particle method for collisional plasma described in the following work.

B. Yan, R. Caflisch, A Monte Carlo method with negative particles for Coulomb collisions, J. Comput. Phys., 298 (2015), pp. 711-740.

B. Yan, A hybrid method with deviational particles for spatial inhomogeneous plasma, J.Comput. Phys., 309 (2016), pp. 18-36.

## CMake build

The preferred cross-platform build uses CMake with the Visual Studio 2022
generator on Windows. It requires C++17,
OpenMP, and FFTW3. The repository contains Windows FFTW headers and libraries;
Linux builds should provide FFTW3 through the system package manager.

Configure and build a Debug tree:

```text
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

The Release equivalents use the `release` preset. Catch2 is pinned to v3.7.1
through CMake `FetchContent`. For offline builds, pre-populate CMake's
FetchContent source cache or configure with `-DNEGPAR_FETCH_TEST_DEPS=OFF` and
provide a locally installed Catch2 3.7.1 package. A local source checkout can
also be supplied directly with
`-DNEGPAR_CATCH2_SOURCE_DIR=<path-to-Catch2-v3.7.1>`; this avoids both Git and
HTTPS download requirements.

Run the sanitizer build with:

```text
cmake --preset sanitizer
cmake --build --preset sanitizer
ctest --preset sanitizer
```

This enables AddressSanitizer with MSVC and AddressSanitizer plus
UndefinedBehaviorSanitizer with supported GCC or Clang toolchains.

On a supported GCC or Clang platform, configure race detection separately:

```text
cmake -S . -B build/thread-sanitizer -DNEGPAR_ENABLE_THREAD_SANITIZER=ON
cmake --build build/thread-sanitizer
ctest --test-dir build/thread-sanitizer --output-on-failure
```

ThreadSanitizer cannot be combined with the address/undefined sanitizer option.
MSVC does not provide ThreadSanitizer, so this gate is intended for supported
Linux GCC/Clang environments.

Test output is written below the build tree. The regular Debug, Release, and
sanitizer configurations all include the Fourier resampler's deterministic,
conservation, bounds, and validation tests.

For a repeatable run, provide an explicit seed:

```text
build\\release\\Release\\negpar_inhomo.exe --seed 12345 --threads 1 --steps 100 --output-dir run_12345
```

The seed controls per-thread `std::mt19937` engines. Each engine is initialized
with `std::seed_seq{base_seed, OpenMP_thread_id}`. Distribution behavior still
depends on the selected standard-library implementation, so cross-platform
bitwise identity is not currently promised. Use `--threads 1` for the
designated deterministic reference run.
The output directory is created automatically. If omitted, the legacy `result`
directory is used.

Each run also writes `run_metadata.txt` with the effective seed, thread count,
RNG engine/distribution information, and reproducibility expectations.

### Feature configuration

The checked-in config/negpar.example.json exposes every opt-in numerical mode
while preserving the legacy defaults. Load it, or a copy of it, with:

~~~text
build\\release\\Release\\negpar_inhomo.exe --config config\\my-run.json
~~~

Configuration values are applied before explicit CLI overrides. Validate a
file without starting a simulation, or print the resolved configuration:

~~~text
build\\release\\Release\\negpar_inhomo.exe --validate-config config\\my-run.json
build\\release\\Release\\negpar_inhomo.exe --config config\\my-run.json --print-effective-config
~~~

Each simulation output contains effective_config.json alongside
run_metadata.txt, so a run can be reproduced from the recorded modes, weights,
collision/projection policy, and runtime values.

## Simulation Studio UI

The repository includes a dependency-free local web interface for configuring,
running, and visualizing the simulation. Build a Debug or Release preset first,
then start the Studio from the repository root with Python 3.10 or newer:

```text
python -m ui.server
```

The server binds to `127.0.0.1:8765` and opens the dashboard in the default
browser. Use `--no-open` when running headlessly, `--port <port>` to select a
different port, or `--executable <path>` to select a nonstandard simulation
binary. The Studio automatically discovers the regular CMake and standalone
Visual Studio output locations.

The UI maps directly to the supported CLI options and includes an editable
feature-configuration JSON panel. Load the checked-in example, validate it
through the simulation executable, save a copy under config/, or submit it
with a run; the submitted JSON is saved as ui_config.json in that run's
output. It streams console and step progress while the process runs, then
reads stable spatial snapshots and time series for density, velocity,
temperature, electric field, conservation, and particle population plots.
Output directories must stay inside the repository and must be empty before
launch so separate runs cannot silently overwrite each other.

Use the **Compare** mode in the command bar to select any two saved runs below
the repository's `result/` directory. The Studio synchronizes metric and
snapshot selection, renders both runs with identical plot scales, computes the
A-minus-B spatial delta, and summarizes relative L2 error, energy-drift
difference, CPU-runtime ratio, diagnostics, and run provenance. Use **Refresh**
after a new run finishes if the comparison selector is already open.

Run the UI backend and result-parser tests with:

```text
python -m unittest discover -s ui/tests -v
```

## Reference validation

The designated reference run is single-threaded and uses an explicit seed:

```text
build\\release\\Release\\negpar_inhomo.exe --steps 1 --seed 123 --threads 1 --output-dir validation_123
```

Validation should confirm that the process finishes successfully, writes
`run_metadata.txt`, and produces no unexpected `NaN` or `Inf` values in the
selected output directory. Keep the output directory isolated per run; remove
it after inspection when it is only a smoke test. Runtime is a benchmark
observation, not a correctness assertion.

Each CTest preset includes the one-step single-thread reference smoke test. For
a complete local check, run the Debug, Release, and sanitizer presets; use the
standalone command above only when inspecting the generated artifacts manually.

## Module ownership

Simulation orchestration lives in `Simulation.*`, while the ordered HDP and
PIC numerical steps live in `SimulationSteps.*`. Initialization is separated
into `InitialConditions.*`, `ParticleInitialization.*`, and the
`Initialization.*` grid orchestrator. Full-particle reconstruction is split
between deterministic `FullParticleFourier.*` numerics and stochastic
`FullParticleSampling.*`. Output callers include the focused
`MacroOutput.*`, `ParticleOutput.*`, `RunMetadataOutput.*`, or
`OutputPaths.*` owner directly; no compatibility umbrella is retained.

## Legacy Linux command

The original direct compilation command is retained for reference, but sources
are located under `src`, not the repository root:

```text
g++ src/*.cpp -Isrc -Isrc/fftw/include -o out -Lsrc/fftw/lib -lfftw3 -std=c++17 -fopenmp
./out
```
