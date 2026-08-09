include(FetchContent)

set(CATCH_CONFIG_DISABLE_EXCEPTIONS OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "" FORCE)

set(NEGPAR_CATCH2_SOURCE_DIR "" CACHE PATH
  "Use a local Catch2 3.7.1 source tree instead of downloading it")

if(NEGPAR_CATCH2_SOURCE_DIR)
  if(NOT EXISTS "${NEGPAR_CATCH2_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "NEGPAR_CATCH2_SOURCE_DIR does not contain a Catch2 source tree: "
      "${NEGPAR_CATCH2_SOURCE_DIR}")
  endif()
  FetchContent_Declare(Catch2 SOURCE_DIR "${NEGPAR_CATCH2_SOURCE_DIR}")
else()
  FetchContent_Declare(
    Catch2
    # Use the release archive instead of Git so fresh Windows builds do not
    # depend on Git for Windows' POSIX signal-pipe support.
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.7.1.zip
  )
endif()

if(NEGPAR_FETCH_TEST_DEPS)
  FetchContent_MakeAvailable(Catch2)
else()
  find_package(Catch2 3.7.1 CONFIG REQUIRED)
endif()

add_executable(negpar_tests
  test_core.cpp
  test_full_particle_fourier.cpp
  test_full_particle_sampling.cpp
  test_resampler.cpp
  ../src/RunOptions.cpp
)
target_link_libraries(negpar_tests PRIVATE negpar_core Catch2::Catch2WithMain)
target_compile_features(negpar_tests PRIVATE cxx_std_17)
if(COMMAND negpar_copy_fftw_runtime)
  negpar_copy_fftw_runtime(negpar_tests)
endif()
if(NEGPAR_ENABLE_SANITIZERS AND MSVC)
  negpar_copy_sanitizer_runtime(negpar_tests)
endif()

# The repository carries the Windows FFTW import library but not its DLL.
# Delay-loading keeps non-FFT characterization tests runnable; FFT tests still
# require libfftw3-3.dll to be present at execution time.
if(MSVC)
  target_link_options(negpar_tests PRIVATE "/DELAYLOAD:libfftw3-3.dll")
  target_link_libraries(negpar_tests PRIVATE delayimp)
endif()

include(CTest)
include(Catch)
catch_discover_tests(negpar_tests)

add_test(
  NAME negpar_reference_smoke
  COMMAND ${CMAKE_COMMAND}
    -DREFERENCE_APP=$<TARGET_FILE:negpar_inhomo>
    -DREFERENCE_OUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/reference-smoke
    -P ${CMAKE_CURRENT_SOURCE_DIR}/validate_reference.cmake
)
set_tests_properties(negpar_reference_smoke PROPERTIES LABELS reference)
