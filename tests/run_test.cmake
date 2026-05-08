# Per-test driver runner. Invoked by add_test() entries in tests/CMakeLists.txt.
#
# Required cache variables:
#   DRIVER          path to the driver executable to run
#   DRIVER_NAME     short name (driver1_f77, driver1_f90, ...) used by the checker
#   WORKDIR         per-test scratch directory (created here)
#   CHECK_SCRIPT    path to check_output.py
#   PYTHON          python interpreter

file(MAKE_DIRECTORY ${WORKDIR})
file(REMOVE
   ${WORKDIR}/output.json
   ${WORKDIR}/iterate.dat
)

# LBFGSB_TLIMIT and LBFGSB_NFG_LIMIT are set high so the algorithm reaches
# the gradient-based stopping criterion (driver2/driver3 otherwise abort
# at arbitrary eval-count or wallclock caps before convergence).
execute_process(
   COMMAND ${CMAKE_COMMAND} -E env
              LBFGSB_JSON_OUTPUT=${WORKDIR}/output.json
              LBFGSB_TLIMIT=86400
              LBFGSB_NFG_LIMIT=100000
              ${DRIVER}
   WORKING_DIRECTORY ${WORKDIR}
   OUTPUT_QUIET
   RESULT_VARIABLE drv_rc
)
if(NOT drv_rc EQUAL 0)
   message(FATAL_ERROR "Driver ${DRIVER} exited with code ${drv_rc}")
endif()

if(NOT EXISTS ${WORKDIR}/output.json)
   message(FATAL_ERROR "Driver did not produce ${WORKDIR}/output.json (env var not honored?)")
endif()

execute_process(
   COMMAND ${PYTHON} ${CHECK_SCRIPT} ${DRIVER_NAME} ${WORKDIR}/output.json
   RESULT_VARIABLE chk_rc
)
if(NOT chk_rc EQUAL 0)
   message(FATAL_ERROR "Output check failed for ${DRIVER_NAME}: ${WORKDIR}/output.json")
endif()
