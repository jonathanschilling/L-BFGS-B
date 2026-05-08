# Integration tests: run each driver with JSON instrumentation and compare
# the output to a committed golden file using a relative-absolute tolerance.

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(LBFGSB_TEST_DRIVERS
   driver1_f77
   driver1_f90
   driver2_f77
   driver2_f90
   driver3_f77
   driver3_f90
)

foreach(drv IN LISTS LBFGSB_TEST_DRIVERS)
   add_test(
      NAME ${drv}
      COMMAND ${CMAKE_COMMAND}
         -DDRIVER=$<TARGET_FILE:${drv}>
         -DEXPECTED=${CMAKE_CURRENT_SOURCE_DIR}/expected/${drv}.json
         -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/${drv}
         -DCOMPARE_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/compare_json.py
         -DPYTHON=${Python3_EXECUTABLE}
         -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
   )
endforeach()
