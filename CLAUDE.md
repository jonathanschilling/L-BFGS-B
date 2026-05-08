# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A repackaging of the original Fortran 77 **L-BFGS-B 3.0** code (Zhu, Byrd, Lu, Nocedal, Morales — limited-memory BFGS for bound-constrained optimization) obtained from Nocedal's Northwestern site. The maintainer (J. Schilling) made three substantive changes versus the upstream tarball:

1. Each subroutine is split into its own file under `src/` (upstream ships them concatenated).
2. Comments were converted to Doxygen markup (`c>` / `!>` prefixes).
3. The bundled BLAS/LINPACK was removed; the library now links against system BLAS/LAPACK. Specifically, `dtrsl` was replaced by LAPACK's `dtrsm` and `dpofa` by LAPACK's `dpotrf`. Keep this substitution in mind when comparing against any other L-BFGS-B fork.

There is a two-tier test suite under `tests/`:
- **Integration tests** (`tests/check_output.py`, 6 ctest entries): run each `drivers/driver{1,2,3}_f{77,90}` and verify final f / projg / x against bit-tight bounds derived from each driver's stopping criterion (see `tests/check_output.py:DRIVER_SPECS`). Drivers emit a JSON dump when `LBFGSB_JSON_OUTPUT` is set; the checker reads that.
- **Per-subroutine unit tests** (`tests/unit/test_<sub>.f90`, 17 ctest entries): hand-crafted test inputs for each in-scope subroutine (timer + prn{1,2,3}lb are out of scope per CLAUDE-tracked decision; everything else has its own `test_*.f90`). A coverage gate (`tests/unit/check_branch_coverage.py`) runs as the 24th ctest entry and fails if line coverage drops below 90% on the in-scope src/ files.

CI runs in `.github/workflows/test.yml` — two jobs, both required: `test` (Release build, integration tests) and `unit_tests` (coverage build, full ctest + gcovr report uploaded).

## Build

CMake-based, requires system BLAS and LAPACK:

```bash
mkdir build && cd build
cmake ..
make -j
ctest --output-on-failure        # run the suite
```

This produces `liblbfgsb.so` and six driver executables (`driver{1,2,3}_f{77,90}`) in `build/`. Each driver writes an `iterate.dat` file when run with `iprint > 0`, and a JSON summary when `LBFGSB_JSON_OUTPUT` is set; both are gitignored along with `x.*` outputs.

For coverage builds: `cmake -DLBFGSB_COVERAGE=ON ..; make; make coverage` (requires `gcovr`).

When adding a new source file: edit `src/CMakeLists.txt` to append it to `liblbfgsb_src` (the parent-scope variable pattern is intentional — `src/CMakeLists.txt` does not define a target itself).

## Code architecture

Public entry point is **`setulb`** (`src/setulb.f`). Users call `setulb` in a reverse-communication loop: the routine returns with `task` set to `FG` to request a function/gradient evaluation, `NEW_X` after a successful iteration, or a terminal status. The driver (the user's code) does the f/g evaluation and calls `setulb` again with the same workspace arrays.

`setulb` is a thin wrapper that partitions the user-provided `wa`/`iwa` workspaces into the named arrays the algorithm needs (`ws`, `wy`, `sy`, `ss`, `wt`, `wn`, `snd`, etc.) and forwards to **`mainlb`** (`src/mainlb.f`), which contains the actual algorithm state machine.

Below `mainlb`, the call graph is roughly:
- **`cauchy`** — generalized Cauchy point (identifies active bound set)
- **`subsm`** — subspace minimization on free variables
- **`lnsrlb`** → **`dcsrch`** → **`dcstep`** — More–Thuente line search
- **`formk`**, **`formt`**, **`bmv`**, **`matupd`** — compact L-BFGS matrix bookkeeping (the K, T = L L', and (B, theta) updates from Byrd/Nocedal/Schnabel 1994)
- **`active`**, **`freev`**, **`projgr`** — active-set / free-variable / projected-gradient utilities
- **`hpsolb`** — heap sort used by `cauchy` to order breakpoints
- **`prn1lb`/`prn2lb`/`prn3lb`** — diagnostic printing (controlled by `iprint`)
- **`errclb`** — input validation
- **`timer`** — wallclock helper

One subroutine per file, with the filename matching the subroutine name. Doxygen `@param` blocks at the top of each file are the canonical parameter documentation — keep them in sync if you change signatures.

## Documentation

Doxygen config is in `Doxyfile`. The GitHub Actions workflow `.github/workflows/doxygen.yml` runs doxygen on every push, builds the LaTeX manual, and on `master` deploys HTML + `L-BFGS-B.pdf` to the `gh-pages` branch. The other CI workflow is `.github/workflows/test.yml` (see Build section) which catches build/test regressions.

Reference PDFs (algorithm paper, ACM remark, original code listing) live in `docs/`.
