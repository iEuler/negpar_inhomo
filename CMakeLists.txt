cmake_minimum_required(VERSION 3.24)

project(negpar_inhomo VERSION 1.0 LANGUAGES CXX)

option(NEGPAR_BUILD_TESTS "Build regression tests" ON)
option(NEGPAR_FETCH_TEST_DEPS "Fetch pinned test dependencies when unavailable" ON)
option(NEGPAR_BUILD_EXPERIMENTAL_RESAMPLER
  "Compile the isolated experimental resampler library (not linked into the app)" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

find_path(FFTW3_INCLUDE_DIR fftw3.h
  PATHS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fftw/include"
)
find_library(FFTW3_LIBRARY NAMES fftw3-3 fftw3
  PATHS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/fftw/lib"
)

if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
  message(FATAL_ERROR "FFTW3 headers/library not found. Set FFTW3_INCLUDE_DIR and FFTW3_LIBRARY.")
endif()

add_library(negpar_core STATIC
  src/Classes.cpp
  src/Collisions.cpp
  src/Advection.cpp
  src/Diagnostics.cpp
  src/ElectricField.cpp
  src/FFT.cpp
  src/Initialization.cpp
  src/LegacyResampling.cpp
  src/Moments.cpp
  src/NegativeParticle.cpp
  src/Numerics.cpp
  src/Output.cpp
  src/utils.cpp
)

target_include_directories(negpar_core
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${FFTW3_INCLUDE_DIR}"
)
target_link_libraries(negpar_core PUBLIC OpenMP::OpenMP_CXX "${FFTW3_LIBRARY}")
target_compile_features(negpar_core PUBLIC cxx_std_17)

if(MSVC)
  target_compile_options(negpar_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(negpar_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(NEGPAR_BUILD_EXPERIMENTAL_RESAMPLER)
  add_library(negpar_experimental_resampler STATIC
    src/Resampler.cpp
    src/ResamplerHelper.cpp
    src/InhomoResampler.cpp
  )
  target_include_directories(negpar_experimental_resampler
    PUBLIC
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
      "${FFTW3_INCLUDE_DIR}"
  )
  target_link_libraries(negpar_experimental_resampler
    PRIVATE negpar_core OpenMP::OpenMP_CXX "${FFTW3_LIBRARY}")
  target_compile_features(negpar_experimental_resampler PUBLIC cxx_std_17)
  if(MSVC)
    target_compile_options(negpar_experimental_resampler PRIVATE /W4 /permissive-)
  else()
    target_compile_options(negpar_experimental_resampler PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

add_executable(negpar_inhomo
  src/main.cpp
  src/Simulation.cpp
  src/RunOptions.cpp
)
target_link_libraries(negpar_inhomo PRIVATE negpar_core)

if(WIN32 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/fftw/lib/libfftw3-3.dll")
  add_custom_command(TARGET negpar_inhomo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${CMAKE_CURRENT_SOURCE_DIR}/src/fftw/lib/libfftw3-3.dll"
      "$<TARGET_FILE_DIR:negpar_inhomo>/libfftw3-3.dll")
endif()

if(NEGPAR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
