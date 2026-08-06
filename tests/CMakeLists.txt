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
  ../src/RunOptions.cpp
)
target_link_libraries(negpar_tests PRIVATE negpar_core Catch2::Catch2WithMain)
target_compile_features(negpar_tests PRIVATE cxx_std_17)

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

if(TARGET negpar_experimental_resampler)
  add_executable(negpar_experimental_tests
    test_experimental.cpp
  )
  target_link_libraries(negpar_experimental_tests
    PRIVATE negpar_experimental_resampler Catch2::Catch2WithMain)
  target_compile_features(negpar_experimental_tests PRIVATE cxx_std_17)

  catch_discover_tests(negpar_experimental_tests)
endif()
