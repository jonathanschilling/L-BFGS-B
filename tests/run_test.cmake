# Per-test driver runner. Invoked by add_test() entries in tests/CMakeLists.txt.
#
# Required cache variables:
#   DRIVER          path to the driver executable to run
#   EXPECTED        path to the committed reference JSON
#   WORKDIR         per-test scratch directory (created here)
#   COMPARE_SCRIPT  path to compare_json.py
#   PYTHON          python interpreter

file(MAKE_DIRECTORY ${WORKDIR})
file(REMOVE
   ${WORKDIR}/output.json
   ${WORKDIR}/iterate.dat
)

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
   COMMAND ${PYTHON} ${COMPARE_SCRIPT} ${EXPECTED} ${WORKDIR}/output.json
   RESULT_VARIABLE cmp_rc
)
if(NOT cmp_rc EQUAL 0)
   message(FATAL_ERROR
      "JSON mismatch.\n"
      "  expected: ${EXPECTED}\n"
      "  actual:   ${WORKDIR}/output.json\n"
      "If the algorithm change is intentional, run tests/regenerate_expected.sh and review the diff."
   )
endif()
