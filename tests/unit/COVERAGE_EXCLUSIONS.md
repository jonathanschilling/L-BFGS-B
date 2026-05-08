# Coverage exclusions

This file lists branches in `src/` that the unit-test suite is allowed to
leave uncovered. Each entry must document **why** the branch is exempt — a
defensive path the state machine cannot reach, a numerical-edge-case path
that requires unreasonable test setup, etc. Reviewers add or refuse
exemptions during code review.

The 100% branch-coverage gate (`tests/unit/check_branch_coverage.py`) reads
this file and removes the listed lines from the must-cover set.

## Format

```
## src/<file>.f
- line <N>: <one-paragraph justification>
- line <N>, branch 0: <justification for the false-direction only>
- line <N>, branch 1: <justification for the true-direction only>
```

`branch 0` is the false-taken arm, `branch 1` is the true-taken arm, and
omitting the branch direction (or writing `branch both`) exempts the line as
a whole. When an entry exempts only one direction, the test suite must still
cover the other.

## Out-of-scope files

The four files listed below are excluded from the must-cover set entirely
via the `--exclude-file` flag (not via this file). They are not core
numerics:

- `src/timer.f` — wallclock helper
- `src/prn1lb.f`, `src/prn2lb.f`, `src/prn3lb.f` — diagnostic printing

## Exemptions

## src/formk.f
- line 313-314: `info = -2; return` after the second dpotrf on the (2,2)
  block. The (2,2) block is `theta*S'AA'S + (L^-1*(-La'+Rz'))'(L^-1*(-La'+Rz'))`,
  the second term is positive-semidefinite (a quadratic form), and the first
  is positive-definite when theta > 0 (which any well-formed L-BFGS-B run
  ensures via the curvature update in matupd). Triggering this branch
  requires constructing an artificial state where the matrix is non-PD
  AND the first dpotrf at line 282 happens to succeed. Reachable in
  pathological inputs the algorithm itself never produces.

## Known partial-coverage subroutines

After Phase F follow-up tightening, the residual uncovered lines fall
into two buckets:

### dcsrch / dcstep tail-end conditionals

- **`src/dcsrch.f` line 204, 206, 210**: the three "WARNING:..." paths.
  Trigger only when bracketing causes rounding-driven progress failure,
  the xtol-based shrink fires, or stp lands exactly on stpmin. Reachable
  in adversarial line searches; covered by the existing Moré–Thuente
  test suite (driver runs) but not by the n=2 test inputs we author
  for unit-test scope. Remaining uncovered fraction: 3 lines / 104 = 3%.

- **`src/dcstep.f` line 169**: `stpf = stpc` in Case 3 brackt=T when the
  cubic interpolation lands closer to stp than the secant. Specific
  numerical comparison that the existing test inputs don't satisfy.

### mainlb error-recovery paths

`src/mainlb.f` (~76% line coverage). The bulk of the gap is internal
"refresh-the-lbfgs-memory" blocks fired by:

- iback >= 20 in lnsrlb, restart lbfgs
- abnormal-termination-in-lnsrch
- formt non-PD Cholesky
- skip-the-l-bfgs-update (s'y < eps*||y||^2 curvature condition)

Each is a defensive recovery path the algorithm uses when its compact
representation degrades. They are genuinely unreachable from
well-conditioned bound-Rosenbrock problems and would require
constructing pathological L-BFGS state to trigger. The driver-level
integration tests cover the well-conditioned trajectories, and these
internal recovery paths are the one place where a port could legitimately
diverge in implementation detail without affecting correctness.

(Two earlier entries — "cauchy returning info != 0" and "subsm
returning info != 0" — were removed when the dead `info` parameters
in `bmv` / `cmprlb` / `cauchy` / `subsm` were dropped after the
LINPACK→LAPACK migration left them unreachable.)
