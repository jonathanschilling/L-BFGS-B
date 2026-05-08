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

(none yet — added in Phase F as the gate transitions from warn-only to
fail-on-gap.)
