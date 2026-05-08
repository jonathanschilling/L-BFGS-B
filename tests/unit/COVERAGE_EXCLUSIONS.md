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

## src/cmprlb.f
- line 61-63: `if (info .ne. 0)` after `call bmv(...)`. Unreachable in the
  current build because `bmv.f` was migrated from LINPACK's `dtrsl`
  (which set `info`) to LAPACK's `dtrsm` (which doesn't), and bmv now
  always sets `info = 0` explicitly before returning. Defensive code
  left from the migration; preserved as-is rather than removed because
  a future BLAS swap could resurrect the failure mode.

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

The following subroutines have residual uncovered lines that aren't
individually exempted but reflect specific algorithm-state combinations
not yet exercised. They are tracked here so the gate's line-coverage
threshold can stay snug rather than being relaxed to absorb them.

- **`src/cauchy.f`** (~91% line coverage): the iter==2 ibkmin-swap branch
  (line 335-336), some col>0 update sub-branches (428-433, 440-443),
  and `iprint >= 100` print blocks (350-352).
- **`src/dcsrch.f`** (~86% line coverage): the modified-function
  (psi-stage) branch (229-234, 239, 243-246) and the warning paths
  (204, 206, 210) need adversarial line-search inputs to fire.
- **`src/dcstep.f`** (~96% line coverage): three sub-branches (158, 169,
  176) of Case 3's r/stp safeguarding need targeted inputs where stp <
  stx with specific gamma signs.
- **`src/mainlb.f`** (~74% line coverage): inner-loop dispatch paths
  for `info != 0` returns from cauchy/subsm/formk (genuinely
  defensive given the rest of the algorithm), iprint=99/100 diagnostic
  blocks, and line-search re-entry replays on warning.
- **`src/subsm.f`** (~58% line coverage): the bound-clipping projection
  for nbd=1 and nbd=3 single-bound cases (247-249, 258-260) and the
  positive-directional-derivative backtracking branch (281-283 and
  the line-search clipping at 290-333). Worth closing in a dedicated
  follow-up since these are the largest remaining gap.
