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

## Legacy Linux command

The original direct compilation command is retained for reference, but sources
are located under `src`, not the repository root:

```text
g++ src/*.cpp -Isrc -Isrc/fftw/include -o out -Lsrc/fftw/lib -lfftw3 -std=c++17 -fopenmp
./out
```
