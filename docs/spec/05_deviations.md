# Deviations from the published algorithm

This document catalogues every place where the F77 implementation in
`src/` deviates from the published algorithms in `docs/algorithm.pdf`,
`docs/code.pdf`, and `docs/acm-remark.pdf`. Each deviation is
documented with the location, the published value, the implemented
value, and the rationale (or "unknown" if the source paper for the
specific constant was not located).

A port that wants to reproduce the F77 reference's numerical behavior
**must use the implemented values**, not the published ones; otherwise
the conformance suite (in `--strict` mode against the reference BLAS)
will fail.

## Confirmed deviations from upstream L-BFGS-B 3.0

This fork (jonathanschilling/L-BFGS-B) drops four `info` output
parameters from internal helper routines (`bmv`, `cmprlb`, `cauchy`,
`subsm`) that were left dead by the LINPACK->LAPACK migration. Upstream
3.0 still has them. The signatures of these four routines therefore
differ from any other L-BFGS-B fork that retains the LINPACK-era
parameters:

```
upstream:  call bmv(m, sy, wt, col, v, p, info)
this fork: call bmv(m, sy, wt, col, v, p)

upstream:  call cmprlb(..., cnstnd, info)
this fork: call cmprlb(..., cnstnd)

upstream:  call cauchy(..., sbgnrm, info, epsmch)
this fork: call cauchy(..., sbgnrm, epsmch)

upstream:  call subsm(..., wn, iprint, info)
this fork: call subsm(..., wn, iprint)
```

The `info` outputs were always `0` once `dtrsl` was replaced by
`dtrsm` (which cannot fail on a non-singular factor -- guaranteed by
`formt` / `formk` Cholesky checks). The dead-branch handling in
`mainlb` was removed at the same time.

`formt` and `formk` retain their `info` outputs because they perform
real Cholesky factorisations (`dpotrf`) that can legitimately fail.

## Confirmed deviations from the published algorithm

### 1. `ftol = 1.0e-3` (More-Thuente sufficient-decrease)

| Where | Value | Reference |
|-------|-------|-----------|
| Implemented | `ftol = 1.0d-3` | `src/lnsrlb.f:110` |
| Published (algorithm tech report) | `alpha = 1.0e-4` | More & Thuente 1994 sec.6 |

The algorithm.pdf tech report and the More-Thuente paper specify
`ftol = 1e-4` (the sufficient-decrease constant in the Wolfe
conditions). The L-BFGS-B implementation uses `1e-3`, a looser value.

This deviation is documented in the F77 source itself in the comment
block of `src/lnsrlb.f` (lines 99-105): *"the looser value 1.0d-3
matches the implementation that ships with Algorithm 778; neither the
ACM paper (docs/code.pdf) nor the 2011 remark documents the change
explicitly"*.

**Impact**: the Wolfe sufficient-decrease test
`f(stp) <= f(0) + ftol * stp * g(0)` is satisfied at slightly larger
function values with `ftol = 1e-3` than with `1e-4`. In practice the
line search accepts step lengths a bit earlier than the strict paper
algorithm would. The integration test bounds in `tests/check_output.py`
were tuned against the implemented value; changing `ftol` to `1e-4`
would require regenerating those bounds.

**For a port**: use `1e-3` for `--strict` conformance. A port that
prefers the published value can do so but cannot pass strict mode.

## Constants that match the papers (no deviation)

These were checked during the audit and do match the published
algorithm or its referenced sources.

### `gtol = 0.9` (curvature condition)

`src/lnsrlb.f:110`. Matches More & Thuente 1994 sec.6 and the
algorithm.pdf tech report (eq 2.6).

### `xtol = 0.1` (relative bracket-width tolerance)

`src/lnsrlb.f:110`. Matches More & Thuente 1994 sec.6.

### `p66 = 0.66` (bracket safeguarding factor)

`src/dcstep.f:81`, `src/dcsrch.f:111`. From More & Thuente 1994 sec.3,
the safeguarded-step heuristic. Not in algorithm.pdf (which doesn't
describe the line search in detail); present in code.pdf and the
MINPACK-2 reference. Treated as algorithmic, not arbitrary.

### `xtrapl = 1.1`, `xtrapu = 4.0` (dcsrch extrapolation factors)

`src/dcsrch.f:113`. From the MINPACK-2 line search (More-Thuente).
Inherited by L-BFGS-B; not in algorithm.pdf.

### `p5 = 0.5`, `two = 2.0`, `three = 3.0` (interpolation coefficients)

`src/dcstep.f`, `src/dcsrch.f`. Standard cubic / quadratic
interpolation; not algorithm-specific.

### `0.0d0`, `1.0d0`

Universal numeric literals; no deviation possible.

### `big = 1.0e10` (lnsrlb's no-bounds sentinel)

`src/lnsrlb.f:98`. Implementation detail, not algorithmic. Used
when `cnstnd = false` to set `stpmx` to a large finite value.

### `2.0` in `cauchy.f:420`

```fortran
f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
```

Coefficient in the Taylor expansion of `f1, f2` updates across a
breakpoint; matches the derivation in algorithm.pdf sec.4 and code.pdf
sec.3 (the `2` comes from the cross term in the quadratic).

## Numerical safeguards (not algorithmic deviations, but worth noting)

These are scaling clamps or division guards that don't appear in the
papers but are implementation necessities. None changes the algorithm
on well-conditioned inputs.

### `f2 = max(epsmch * f2_org, f2)` in `cauchy.f:423`

Floor for the segment-curvature variable to prevent
division-by-near-zero in `dtm = -f1 / f2`. Triggered only on
numerical noise; the floor is `epsmch ~= 2.22e-16` times the original
`f2`, so it is many orders of magnitude smaller than typical values
of `f2`.

### `dr <= epsmch * ddum` curvature-condition gate in `mainlb.f:566`

```fortran
if (dr .le. epsmch*ddum) then
   ! skip L-BFGS update; refresh
endif
```

Caller-side check (in `mainlb`) that gates calls to `matupd`. If the
new `(s, y)` pair has curvature `s^T y` smaller than `eps` times the
relevant scale, the pair is rejected to avoid corrupting the L-BFGS
cache. Standard refresh logic; the threshold `epsmch * ddum` keeps
the test scale-invariant.

## Constants from earlier suspicion that turned out NOT to exist

The original spec planning suggested suspect constants `1.0e-8` and
`1.0e-10` somewhere in `cauchy.f` / `subsm.f`. **The audit found no
such constants in the F77 source.** Cauchy's only safety floor uses
`epsmch * f2_org` (machine-eps-scaled, not a fixed value).

## Audit methodology

For each `src/*.f`:

1. `grep` for double-precision literals (`[0-9]\.[0-9]+(d|e)[+-]?[0-9]+`).
2. `grep` for `parameter (` blocks.
3. `grep` for `epsmch` usage to identify scaled tolerances.
4. Cross-reference each non-trivial constant against:
   - `algorithm.pdf` (Byrd/Lu/Nocedal/Zhu 1995) -- the original derivation.
   - `code.pdf` (Zhu/Byrd/Nocedal 1997) -- the algorithm-as-implemented paper.
   - `acm-remark.pdf` (Morales/Nocedal 2011) -- the `subsm` correction.
   - Cited references in the F77 source headers (e.g., More-Thuente
     1994 cited in `dcstep.f` and `dcsrch.f`).

If you spot a deviation not listed above, please open an issue or PR
against the spec pack -- the audit was performed once and may have
missed something subtle.

## Summary

**One** documented algorithm-relevant deviation (`ftol = 1e-3`).

All other numerical literals are either:
- algorithmic constants that match the cited reference (More-Thuente
  for the line-search ones; algorithm.pdf for the rest),
- universal numeric literals (`0`, `1`, `2`, `3`, `0.5`),
- implementation-detail sentinels (`big = 1e10`),
- or `epsmch`-scaled safety floors / gates.

The implementation is, by this measure, faithful to the published
algorithm modulo the documented `ftol` exception.
