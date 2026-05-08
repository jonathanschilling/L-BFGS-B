# Per-subroutine unit tests for L-BFGS-B src/.
#
# Each subroutine in src/ (except timer.f and prn*lb.f, which are out of
# scope) gets a hand-written test that exercises its branches. See the plan
# at docs (or in the planning artefact) for the per-subroutine contract.
#
# This file is the harness: it builds the assertion module and registers
# whatever test_*.f90 sources exist in this directory as ctest entries.
#
# Phase A delivers the harness + the coverage gate skeleton with no tests.
# Phases B–E add tests one tier at a time.

# Assertion helpers (compiled once to avoid .mod-file races under -j).
add_library(test_assert STATIC test_assert.f90)
set_target_properties(test_assert PROPERTIES Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Helper: register one test_<sub>.f90 file as a ctest entry.
function(lbfgsb_add_unit_test sub)
   set(target test_${sub})
   add_executable(${target} test_${sub}.f90)
   target_link_libraries(${target} PRIVATE lbfgsb test_assert)
   target_include_directories(${target} PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/modules)
   add_test(NAME ${target} COMMAND ${target})
endfunction()

# Per-subroutine tests, registered as ctest entries.
# Phase B (Tier 1+2): projgr, errclb, active, freev, hpsolb, cmprlb,
#                     bmv, formt, matupd, dcstep, formk
# Phase C (Tier 3):   dcsrch, lnsrlb
# Phase D (Tier 4):   cauchy, subsm
# Phase E (Tier 5):   setulb, mainlb
#
# Out of scope: timer, prn1lb, prn2lb, prn3lb (not core numerics).
lbfgsb_add_unit_test(projgr)
lbfgsb_add_unit_test(errclb)
lbfgsb_add_unit_test(active)
lbfgsb_add_unit_test(freev)
lbfgsb_add_unit_test(hpsolb)
lbfgsb_add_unit_test(cmprlb)
lbfgsb_add_unit_test(bmv)
lbfgsb_add_unit_test(formt)
lbfgsb_add_unit_test(matupd)
lbfgsb_add_unit_test(dcstep)
lbfgsb_add_unit_test(formk)
lbfgsb_add_unit_test(dcsrch)
lbfgsb_add_unit_test(lnsrlb)
lbfgsb_add_unit_test(cauchy)
lbfgsb_add_unit_test(subsm)
lbfgsb_add_unit_test(setulb)
lbfgsb_add_unit_test(mainlb)

# Branch coverage gate. Wired only when LBFGSB_COVERAGE=ON; otherwise there
# is no coverage data to check. In Phases A–E it runs in warn-only mode
# (always exits 0 with a printed report); Phase F flips it to fail-on-gap.
if(LBFGSB_COVERAGE)
   find_package(Python3 REQUIRED COMPONENTS Interpreter)
   add_test(
      NAME branch_coverage_gate
      COMMAND ${Python3_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/check_branch_coverage.py
         --root ${CMAKE_SOURCE_DIR}
         --build ${CMAKE_BINARY_DIR}
         --exclusions ${CMAKE_CURRENT_SOURCE_DIR}/COVERAGE_EXCLUSIONS.md
         --exclude-file timer.f
         --exclude-file prn1lb.f
         --exclude-file prn2lb.f
         --exclude-file prn3lb.f
         --warn-only
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   )
   # Run after all other tests so .gcda data exists.
   set_tests_properties(branch_coverage_gate PROPERTIES
      RUN_SERIAL TRUE)
endif()
