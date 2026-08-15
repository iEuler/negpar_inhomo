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

add_executable(negpar_tests EXCLUDE_FROM_ALL
  test_core.cpp
  test_full_particle_fourier.cpp
  test_full_particle_sampling.cpp
  test_initialization.cpp
  test_resampler.cpp
  test_config.cpp
  test_weighted_hdp.cpp
  test_feature_modes.cpp
  ../src/app/ConfigFile.cpp
  ../src/app/RunOptions.cpp
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

set(PARSE_CATCH_TESTS_ADD_TARGET_IN_TEST_NAME OFF)
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(ParseAndAddCatchTests)
# These tests use literal TEST_CASE declarations, which the configure-time
# parser supports. Avoid build-time discovery because Visual Studio's external
# build adapter can enter a build/discovery feedback loop.
set(OptionalCatchTestLauncher "${CMAKE_COMMAND}" -E env --)
set(CMAKE_WARN_DEPRECATED OFF)
ParseAndAddCatchTests(negpar_tests)
unset(CMAKE_WARN_DEPRECATED)
unset(OptionalCatchTestLauncher)

add_test(
  NAME negpar.reference.smoke
  COMMAND ${CMAKE_COMMAND}
    -DREFERENCE_APP=$<TARGET_FILE:negpar_inhomo>
    -DREFERENCE_OUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/reference-smoke
    -P ${CMAKE_CURRENT_SOURCE_DIR}/validate_reference.cmake
)
set_tests_properties(negpar.reference.smoke PROPERTIES LABELS reference)
