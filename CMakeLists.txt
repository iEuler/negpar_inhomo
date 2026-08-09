cmake_minimum_required(VERSION 3.24)

project(negpar_inhomo VERSION 1.0 LANGUAGES CXX)

option(NEGPAR_BUILD_TESTS "Build regression tests" ON)
option(NEGPAR_FETCH_TEST_DEPS "Fetch pinned test dependencies when unavailable" ON)
option(NEGPAR_ENABLE_SANITIZERS
  "Enable AddressSanitizer and UndefinedBehaviorSanitizer where supported" OFF)
option(NEGPAR_ENABLE_THREAD_SANITIZER
  "Enable ThreadSanitizer with a supported GCC or Clang toolchain" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NEGPAR_ENABLE_SANITIZERS AND NEGPAR_ENABLE_THREAD_SANITIZER)
  message(FATAL_ERROR
    "Address/Undefined sanitizers cannot be combined with ThreadSanitizer")
endif()

if(NEGPAR_ENABLE_SANITIZERS)
  if(MSVC)
    if(MSVC_VERSION LESS 1929)
      message(FATAL_ERROR "AddressSanitizer requires Visual Studio 2019 16.9 or newer")
    endif()
    # MSVC's Debug runtime checks and incremental linking are incompatible
    # with AddressSanitizer.
    string(REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS_DEBUG
      "${CMAKE_CXX_FLAGS_DEBUG}")
    add_compile_options(/fsanitize=address /Zi)
    add_link_options(/INCREMENTAL:NO)
    get_filename_component(NEGPAR_COMPILER_RUNTIME_DIR
      "${CMAKE_CXX_COMPILER}" DIRECTORY)
    find_file(NEGPAR_ASAN_RUNTIME
      NAMES clang_rt.asan_dynamic-x86_64.dll
      HINTS "${NEGPAR_COMPILER_RUNTIME_DIR}"
      REQUIRED)

    function(negpar_copy_sanitizer_runtime target)
      add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${NEGPAR_ASAN_RUNTIME}"
          "$<TARGET_FILE_DIR:${target}>")
    endfunction()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
  else()
    message(FATAL_ERROR
      "NEGPAR_ENABLE_SANITIZERS is unsupported by ${CMAKE_CXX_COMPILER_ID}")
  endif()
endif()

if(NEGPAR_ENABLE_THREAD_SANITIZER)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
    add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
    add_link_options(-fsanitize=thread)
  else()
    message(FATAL_ERROR
      "ThreadSanitizer requires a supported GCC or Clang toolchain")
  endif()
endif()

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
  src/Grid.cpp
  src/Particle.cpp
  src/ParticleGroup.cpp
  src/ParticleGroupOperations.cpp
  src/SimulationConfig.cpp
  src/SimulationTypes.cpp
  src/Collisions.cpp
  src/Advection.cpp
  src/Diagnostics.cpp
  src/ElectricField.cpp
  src/FFT.cpp
  src/InitialConditions.cpp
  src/Initialization.cpp
  src/ParticleInitialization.cpp
  src/FullParticleFourier.cpp
  src/FullParticleSampling.cpp
  src/Moments.cpp
  src/SimulationSteps.cpp
  src/NegativeParticleCollisions.cpp
  src/NegativeParticleSampling.cpp
  src/ParticleConservation.cpp
  src/Numerics.cpp
  src/MacroOutput.cpp
  src/OutputPaths.cpp
  src/ParticleOutput.cpp
  src/RunMetadataOutput.cpp
  src/ParticleResampling.cpp
  src/ProjectionSampling.cpp
  src/ResamplingNumerics.cpp
  src/ResamplingVelocity.cpp
  src/Resampler.cpp
  src/ResamplerHelper.cpp
  src/Histogram.cpp
  src/RandomContext.cpp
  src/RandomSampling.cpp
)

target_include_directories(negpar_core
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${FFTW3_INCLUDE_DIR}"
)
target_link_libraries(negpar_core PUBLIC OpenMP::OpenMP_CXX "${FFTW3_LIBRARY}")
target_compile_features(negpar_core PUBLIC cxx_std_17)

function(negpar_enable_strict_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

negpar_enable_strict_warnings(negpar_core)

add_executable(negpar_inhomo
  src/main.cpp
  src/Simulation.cpp
  src/RunOptions.cpp
)
target_link_libraries(negpar_inhomo PRIVATE negpar_core)
negpar_enable_strict_warnings(negpar_inhomo)
if(NEGPAR_ENABLE_SANITIZERS AND MSVC)
  negpar_copy_sanitizer_runtime(negpar_inhomo)
endif()

if(WIN32 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/fftw/lib/libfftw3-3.dll")
  function(negpar_copy_fftw_runtime target)
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PROJECT_SOURCE_DIR}/src/fftw/lib/libfftw3-3.dll"
        "$<TARGET_FILE_DIR:${target}>/libfftw3-3.dll")
  endfunction()
  negpar_copy_fftw_runtime(negpar_inhomo)
endif()

if(NEGPAR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
