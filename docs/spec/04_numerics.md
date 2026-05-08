# Numerical considerations

This document covers floating-point conventions, tolerances, and the
BLAS/LAPACK reference-build policy that ports must follow to achieve
bit-for-bit conformance.

The full magic-constants table (every numerical literal in the F77
source with its origin and whether ports may vary it) appears at the
end and is filled in during the deviations audit (see Phase D in the
plan).

## Precision

L-BFGS-B operates in **IEEE-754 double precision throughout**. All
real values are 64-bit binary doubles with normal/subnormal handling,
NaN propagation, and rounding mode = round-to-nearest-even (the IEEE
default).

Single-precision is **not supported**. The algorithm's
ill-conditioning behavior near convergence requires the headroom of
double precision; ports that wrap the algorithm in a single-precision
interface internally upcast to double for the algorithm itself.

Extended precision (`long double`, x87 80-bit) **must not be used**.
Some legacy compilers default to extended precision in registers; this
breaks bit-for-bit reproducibility against the F77 reference. Ports on
x86 must explicitly enforce `-msse2` (or equivalent SSE2 / FPU control
word setting) to keep all arithmetic in 64-bit IEEE-754.

## Machine epsilon

The algorithm uses machine epsilon (`eps`) in tolerance computations.

| Spec | F77 | Value |
|------|-----|-------|
| `eps` | `epsmch` | `~= 2.220446049250313e-16` (= `2^{-52}`) |

The F77 source recomputes `epsmch` at startup in `mainlb` via:

```
epsmch = 1
do while (1 + epsmch /= 1)
   epsmch = epsmch / 2
end do
epsmch = epsmch * 2
```

Ports may either:
- Use the same iterative computation (guaranteed to produce the same
  result on IEEE-754 hardware), or
- Use the language's built-in `DBL_EPSILON` / `numpy.finfo(float).eps`
  / equivalent, which on conforming IEEE-754 platforms equals
  `2.220446049250313e-16`.

Both choices produce identical tolerance arithmetic.

## Convergence tolerances

### Function-decrease tolerance: `factr`

The F77 user passes `factr` (default `1e7`) and the algorithm tests:

```
(f_old - f_new) / max(|f_old|, |f_new|, 1.0) <= factr * eps
```

Equivalent to: relative function decrease `<= factr * eps`. With
`factr = 1e7`, this is roughly `2.2e-9` -- about 9 digits of accuracy.
With `factr = 1e1`, it is roughly `2.2e-15` -- close to machine
precision.

Setting `factr = 0` disables this convergence test.

### Projected-gradient tolerance: `pgtol`

```
||g_proj||_Inf <= pgtol
```

where `g_proj` is the projected gradient defined in `subroutines/projgr.md`.
Default `pgtol = 1e-5`. Setting `pgtol = 0` disables this test.

### Both tests can fire

The algorithm tests both `factr` and `pgtol` after each iteration.
Whichever fires first determines the `info` code (see `02_api.md`).

## Line-search tolerances

The More-Thuente line search (`dcsrch`) uses three internal tolerances
plus step-length bounds.

| Symbol | F77 | Value | Meaning |
|--------|-----|-------|---------|
| `ftol` | `ftol` | `1.0e-3` | Sufficient-decrease (Armijo) coefficient. |
| `gtol` | `gtol` | `0.9` | Curvature condition coefficient. |
| `xtol` | `xtol` | `0.1` | Bracket-width relative tolerance. |
| `stpmin` | `stepmin` | `0.0` | Minimum allowed step. |
| `stpmax` | `stepmx` | computed | Max step (limited by box geometry). |

**`ftol = 1e-3` is a documented deviation from the published algorithm**,
which prescribes `alpha = 1e-4` in More-Thuente 1994. The L-BFGS-B
implementation uses `1e-3` and has been validated against the test
suite at this value. See `lnsrlb.md` for the comment in `src/lnsrlb.f`
explaining this. **Ports must use `ftol = 1e-3`** (not `1e-4`) to be
conformant. If we ever change this value, conformance test bounds will
need to be regenerated.

## Numerical safeguards

### Division guards

The algorithm performs several divisions where a denominator near zero
would be catastrophic. Each is guarded:

| Site | Guard | Behavior on guard fire |
|------|-------|------------------------|
| `cauchy`: breakpoint computation `(u - x) / d` | implicit (skip if `d = 0`) | Variable treated as having no breakpoint in this direction. |
| `cauchy`: `dt = -fp / fpp` | check `fpp > eps` (else use earlier breakpoint) | Algorithm exits the breakpoint loop. |
| `formt`: Cholesky pivot `T[i,i]^2 > 0` | LAPACK `dpotrf` returns `info > 0` | Caller (`mainlb`) signals refresh; L-BFGS history discarded. |
| `bmv`: `1 / sqrt(sy[i,i])` | precondition: `sy[i,i] > 0` enforced by `matupd` | If precondition violated, behavior undefined. |
| `matupd`: `s'y` curvature test | `s'y > eps * ||y||^2` (refresh threshold) | New pair rejected; `updatd = false`. |

### Overflow / underflow

The algorithm does not explicitly guard against IEEE-754 overflow to
infinity in arithmetic. Inputs that produce intermediate `Inf` or
`NaN` typically cause the line search to fail, returning
`INFO_ABNORMAL_LNSRLB`. The user should ensure `f_eval` and `g_eval`
return finite values for feasible inputs.

### Step-length floor

The line search has an implicit minimum step set by `stpmin = 0`.
After the bracket has shrunk so that `stx ~= sty`, `dcsrch` returns
with the best step found so far rather than continuing to shrink. See
`dcsrch.md`.

## BLAS / LAPACK policy

The algorithm calls these BLAS-level-1 and LAPACK routines:

| Call | Used in | Purpose |
|------|---------|---------|
| `dcopy` | many | Copy a vector. |
| `daxpy` | many | `y <- alpha*x + y`. |
| `dscal` | many | `y <- alpha*y`. |
| `ddot` | many | `x' y`. |
| `dnrm2` | line search, gradients | `||x||_2`. |
| `dgemv` | `formk`, `subsm`, `cauchy` | `y <- alpha*A*x + beta*y`. |
| `dtrsm` | `bmv`, `subsm` | Triangular solve `A x = b`. |
| `dpotrf` | `formt`, `formk` | Cholesky factorization. |

### Reference build for bit-for-bit conformance

Bit-for-bit reproducibility across BLAS implementations is **not
achievable**: OpenBLAS, MKL, Accelerate, and netlib-BLAS differ in
internal reduction order and SIMD behavior, and threaded versions add
nondeterminism through dynamic chunking.

For `--strict` conformance, ports and the F77 reference must link
against the **same reference BLAS/LAPACK build**:

- **Reference BLAS**: netlib BLAS, single-threaded, source build.
- **Reference LAPACK**: netlib LAPACK, single-threaded, source build.
- **Reference compiler**: `gfortran` (any version >= 9 in CI).
- **Reference flags**: `-O2 -fno-fast-math -fno-associative-math
  -frounding-math` to disable algebraic reordering and float-contract.

CI pins these versions and rebuilds them in the conformance job. The
exact version pin is in `.github/workflows/test.yml` (see Phase E of
the plan); the spec runner reads the pin from a metadata file in
`docs/spec/` so it stays in sync.

### Tolerance mode for non-reference BLAS

Ports that must use a non-reference BLAS (e.g., a
performance-optimized vendor implementation in production) can run the
conformance suite in `--tolerance` mode, which uses per-subroutine
numerical tolerances instead of bit-equality. See `07_conformance.md`
for the per-routine tolerances.

The `--tolerance` mode is for **validation**, not for "real"
conformance -- a port that only passes `--tolerance` may diverge from
the F77 reference under sufficiently long iteration sequences. For
production use of a non-reference BLAS, the port must accept this
risk; it is not an algorithmic bug, but a consequence of BLAS
nondeterminism.

## Order-of-operations constraints

For bit-for-bit conformance, **the order of additions and
subtractions matters** (floating-point addition is not associative).
Per-subroutine specs in `subroutines/` document the operation order
explicitly: which loops accumulate in which direction, which sums
group with which, etc.

In particular:
- **Reduction direction**: F77 BLAS-1 routines (`ddot`, `dnrm2`,
  `daxpy`) accumulate in increasing index order. Ports using their
  language's BLAS binding inherit this from the reference build.
- **Manual accumulations**: where the algorithm computes a sum without
  a BLAS call (e.g., the inner loops of `cauchy` and `subsm`), the
  spec specifies the loop order. Reversing is **not allowed** in
  `--strict` mode.
- **Fused multiply-add (FMA)**: not used in the F77 reference (the
  flags above disable it). Ports compiled with FMA enabled will
  produce different bit-level outputs even for `--strict` mode and are
  non-conformant.

## Magic-constants table

Populated by the Phase D deviations audit; see `05_deviations.md` for
the full audit notes.

### Algorithmic constants (require port to match)

| Constant | Site | Value | Origin | Tunable? |
|----------|------|-------|--------|----------|
| `ftol` | `lnsrlb.f`, `dcsrch.f` | `1.0e-3` | **Deviation from More-Thuente 1994 (which uses `1e-4`).** Code-only. See `05_deviations.md`. | No. Validated test bounds depend on this. |
| `gtol` | `lnsrlb.f`, `dcsrch.f` | `0.9` | More-Thuente 1994 sec.6 standard. | No. |
| `xtol` | `lnsrlb.f`, `dcsrch.f` | `0.1` | More-Thuente 1994 sec.6. | No. |
| `p66` | `dcstep.f`, `dcsrch.f` | `0.66` | More-Thuente 1994 sec.3 bracket safeguard. | No. |
| `xtrapl` | `dcsrch.f` | `1.1` | MINPACK-2 line search lower-extrapolation. | No. |
| `xtrapu` | `dcsrch.f` | `4.0` | MINPACK-2 line search upper-extrapolation. | No. |

### Universal numeric literals

| Constant | Sites | Value | Notes |
|----------|-------|-------|-------|
| `zero`, `one` | many | `0.0d0`, `1.0d0` | Universal. |
| `p5`, `two`, `three` | `dcstep.f`, `dcsrch.f` | `0.5`, `2.0`, `3.0` | Cubic / quadratic interpolation coefficients. |
| `2.0d0` | `cauchy.f:420` | `2.0` | Cross-term coefficient in the `f2` update across a breakpoint. From the Taylor expansion in `algorithm.pdf` sec.4. |

### Implementation-detail sentinels

| Constant | Site | Value | Notes |
|----------|------|-------|-------|
| `big` | `lnsrlb.f` | `1.0e10` | No-bounds-binding sentinel for `stpmx`. Not algorithmic. |

### `epsmch`-scaled safeguards

These don't appear as fixed numerical literals; they scale with the
machine precision of the running platform.

| Where | Form | Purpose |
|-------|------|---------|
| `cauchy.f:423` | `max(epsmch * f2_org, f2)` | Lower bound on segment curvature to prevent `0/0` in `dtm = -f1/f2`. |
| `mainlb.f:566` | `dr <= epsmch * ddum` | L-BFGS curvature-condition gate; rejects `(s, y)` pairs with too-small `s'y`. |
| `formt.f`, `formk.f` | LAPACK `dpotrf` failure | Cholesky pivot non-positivity -> caller refreshes L-BFGS history. |

### Constants from earlier suspicion that DO NOT exist

The original Phase D candidate list mentioned `1.0e-8` and `1.0e-10`
gates somewhere in `cauchy.f` / `subsm.f`. **The audit confirmed
neither exists in the F77 source.** Cauchy's only safety floor uses
`epsmch * f2_org` (machine-eps-scaled). See `05_deviations.md` for
the audit methodology.
