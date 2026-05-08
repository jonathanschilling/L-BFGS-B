cmake_minimum_required(VERSION 3.10)

# set the project name
project(lbfgsb Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# gfortran warning hygiene:
#   -fimplicit-none  -- all variables must be declared. Catches missing
#                      `parameter (one=1.0d0)` style mistakes (one such bug
#                      shipped in src/bmv.f for ~5 years before being fixed).
#   -Wall -Wextra    -- standard warning sets; surface unused vars/params,
#                      uninitialised reads, real comparisons, etc.
#   -Wno-compare-reals -- L-BFGS-B legitimately uses == / /= on REAL values
#                      to detect exact algorithmic states (variable exactly
#                      at a bound, exact unit step, theta == 1). The compiler
#                      flag is a stylistic check; the comparisons are intentional.
#
# We deliberately do NOT enable -Wpedantic (which flags every `character*60`
# F77 declaration as obsolescent) or -Wimplicit-interface (the F77 source
# files call BLAS/LAPACK routines without explicit interface blocks, which
# is the standard pattern for legacy F77 numerical code).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
   add_compile_options(-fimplicit-none -Wall -Wextra -Wno-compare-reals)
endif()

# Optional coverage instrumentation. Enable with -DLBFGSB_COVERAGE=ON.
# Adds gcov/gcovr-compatible flags to liblbfgsb and to test executables.
# Used by the per-subroutine unit-test gate; see tests/unit/.
option(LBFGSB_COVERAGE "Build with code-coverage instrumentation (--coverage)" OFF)
if(LBFGSB_COVERAGE)
   if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      message(WARNING "LBFGSB_COVERAGE has only been validated with gfortran; gcov output for ${CMAKE_Fortran_COMPILER_ID} may be incomplete.")
   endif()
   add_compile_options(--coverage -O0 -g)
   add_link_options(--coverage)
endif()

# shared library for L-BFGS-B
set(liblbfgsb_src "")
add_subdirectory (src)
add_library(lbfgsb SHARED ${liblbfgsb_src})

if(BLAS_FOUND)
    target_link_libraries(lbfgsb PRIVATE ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

if(LAPACK_FOUND)
    target_link_libraries(lbfgsb PRIVATE ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

# JSON test instrumentation: linked into every driver so they honour the
# LBFGSB_JSON_OUTPUT environment variable. No-op at runtime if the env var
# is not set, so the example drivers behave as before. Built as a static
# library so the json_writer module is compiled once (avoids races in the
# generated .mod file when multiple targets build the same source).
add_library(json_writer STATIC tests/json_writer.f90)

# example executables
add_executable(driver1_f77 drivers/driver1.f)
add_executable(driver1_f90 drivers/driver1.f90)

add_executable(driver2_f77 drivers/driver2.f)
add_executable(driver2_f90 drivers/driver2.f90)

add_executable(driver3_f77 drivers/driver3.f)
add_executable(driver3_f90 drivers/driver3.f90)

target_link_libraries(driver1_f77 PRIVATE lbfgsb json_writer)
target_link_libraries(driver1_f90 PRIVATE lbfgsb json_writer)

target_link_libraries(driver2_f77 PRIVATE lbfgsb json_writer)
target_link_libraries(driver2_f90 PRIVATE lbfgsb json_writer)

target_link_libraries(driver3_f77 PRIVATE lbfgsb json_writer)
target_link_libraries(driver3_f90 PRIVATE lbfgsb json_writer)

# integration tests + per-subroutine unit tests + F77 conformance driver
enable_testing()
add_subdirectory(tests)
add_subdirectory(tests/unit)
add_subdirectory(tests/conformance)

# Coverage report target. Requires LBFGSB_COVERAGE=ON and gcovr installed.
# Usage: `cmake --build . --target coverage` after tests have been built.
if(LBFGSB_COVERAGE)
   find_program(GCOVR_EXECUTABLE gcovr)
   if(GCOVR_EXECUTABLE)
      add_custom_target(coverage
         COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage/html
         COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
         COMMAND ${GCOVR_EXECUTABLE}
            --root ${CMAKE_SOURCE_DIR}
            --filter ${CMAKE_SOURCE_DIR}/src/
            --html-details ${CMAKE_BINARY_DIR}/coverage/html/index.html
            --json ${CMAKE_BINARY_DIR}/coverage/coverage.json
            --print-summary
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
         COMMENT "Running tests and generating coverage report"
         VERBATIM
      )
   else()
      message(WARNING "gcovr not found in PATH; 'coverage' target unavailable. Install with: pip install --user gcovr")
   endif()
endif()
