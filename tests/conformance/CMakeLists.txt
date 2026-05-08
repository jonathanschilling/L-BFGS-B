# F77 conformance driver (Option 2 of the conformance scheme).
# Reads a structured text input file describing one JSON test vector,
# dispatches by subroutine name, calls the F77 routine, writes the
# actual outputs to a structured text output file. Used by
# docs/spec/runner/conformance.py with --engine fortran.

add_executable(conformance_driver
   conformance_io.f90
   conformance_driver.f90
)
target_link_libraries(conformance_driver PRIVATE lbfgsb)
