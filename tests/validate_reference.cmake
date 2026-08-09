if(NOT DEFINED REFERENCE_APP OR NOT DEFINED REFERENCE_OUTPUT_DIR)
  message(FATAL_ERROR "REFERENCE_APP and REFERENCE_OUTPUT_DIR are required")
endif()

file(REMOVE_RECURSE "${REFERENCE_OUTPUT_DIR}")
file(MAKE_DIRECTORY "${REFERENCE_OUTPUT_DIR}")

execute_process(
  COMMAND "${REFERENCE_APP}"
          --steps 1
          --seed 123
          --threads 1
          --output-dir "${REFERENCE_OUTPUT_DIR}"
  RESULT_VARIABLE run_result
  OUTPUT_VARIABLE run_stdout
  ERROR_VARIABLE run_stderr
)

if(NOT run_result EQUAL 0)
  message(FATAL_ERROR
    "Reference run failed with code ${run_result}\n"
    "stdout:\n${run_stdout}\n"
    "stderr:\n${run_stderr}")
endif()

set(metadata_file "${REFERENCE_OUTPUT_DIR}/run_metadata.txt")
if(NOT EXISTS "${metadata_file}")
  message(FATAL_ERROR "Reference run did not produce ${metadata_file}")
endif()

file(GLOB_RECURSE output_files LIST_DIRECTORIES false
     "${REFERENCE_OUTPUT_DIR}/*")
foreach(output_file IN LISTS output_files)
  file(READ "${output_file}" output_text)
  if(output_text MATCHES "(^|[^A-Za-z])(NaN|Inf|-Inf)([^A-Za-z]|$)")
    message(FATAL_ERROR "Non-finite value found in ${output_file}")
  endif()
endforeach()

message(STATUS "Reference smoke run passed: ${REFERENCE_OUTPUT_DIR}")
