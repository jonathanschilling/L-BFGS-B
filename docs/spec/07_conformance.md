# Conformance criteria

A port is **conformant** if it passes `docs/spec/runner/conformance.py`
in **strict** mode against `docs/spec/data/`. Strict mode requires
exact IEEE-754 byte equality on every output, given the same input.

For ports that cannot use the reference BLAS / LAPACK (and therefore
cannot achieve bit-for-bit numerical agreement with the F77
reference; see `04_numerics.md`), the runner provides a **tolerance**
mode using per-subroutine numerical tolerances.

## Modes

### `--strict` (default)

```
python3 docs/spec/runner/conformance.py --strict
```

Each output value must match `np.float64(expected).tobytes() ==
np.float64(actual).tobytes()`. Integer and string values must match
exactly; boolean values must match exactly; arrays must have matching
length and pointwise bit-equal entries.

Strict mode is the gold-standard conformance test. It requires:

- The reference BLAS / LAPACK pin (see `04_numerics.md`).
- IEEE-754 round-to-nearest with no extended precision.
- Operation order matching the F77 source (documented per
  subroutine).

### `--tolerance`

```
python3 docs/spec/runner/conformance.py --tolerance
```

Output values are compared with absolute + relative tolerance per
subroutine (table below). Integer / string / boolean comparisons stay
exact.

This mode is for ports that:
- Use a non-reference BLAS for performance (MKL, OpenBLAS-threaded, ...).
- Use language-native linear-algebra primitives (Eigen, NumPy with
  default BLAS, LinearAlgebra.jl, ...).

A port that only passes `--tolerance` may diverge from the F77
reference under sufficiently long iteration sequences. This is not a
bug — it's the natural consequence of BLAS non-determinism. For most
applications this is acceptable; for bit-reproducible scientific
computing it is not.

## Per-subroutine tolerances

Used in `--tolerance` mode. Match constants in
`runner/conformance.py:TOLERANCES`.

| Subroutine | abs | rel | Justification |
|------------|-----|-----|---------------|
| `projgr` | 0 | 0 | Pure min/max/abs; integer arithmetic effectively. |
| `errclb` | 0 | 0 | Integer / string only. |
| `active` | 0 | 0 | Bound projection is exact. |
| `freev` | 0 | 0 | Integer / boolean. |
| `hpsolb` | 0 | 0 | Comparisons + assignments. |
| `bmv` | 1e-14 | 1e-14 | Triangular solves (small ULP from BLAS). |
| `cmprlb` | 1e-14 | 1e-14 | Calls `bmv`; inherits its tolerance. |
| `formt` | 1e-14 | 1e-14 | Cholesky factor (LAPACK). |
| `matupd` | 0 | 0 | Pure dot products at the BLAS level — F77 inner-product order is used. |
| `dcstep` | 1e-14 | 1e-14 | Cubic / secant interpolation. |
| `formk` | 1e-14 | 1e-14 | Two Choleskys + triangular-solve block. |
| `dcsrch` | 1e-14 | 1e-14 | Full line-search arithmetic. |
| `lnsrlb` | 1e-14 | 1e-14 | min/max + dot + projection. |
| `cauchy` | 1e-13 | 1e-13 | Accumulated `f1, f2` updates over breakpoints. |
| `subsm` | 1e-13 | 1e-13 | Subspace solve via `wn` factor. |
| `mainlb` | 1e-12 | 1e-12 | Full algorithm loop; tolerances accumulate. |

## How a port reports conformance

The port runs the conformance runner against its own implementation
(see `runner/README.md` for the I/O protocol). Output:

```
=== Conformance (strict) ===
  Pass: 101
  Fail: 0
  Skip: 0
```

For CI integration, the runner exits with status `0` on full pass,
`1` on any failure. A port that passes strict mode against the
reference BLAS pin earns the conformance badge; ports that pass
only tolerance mode are noted as such.

## Known caveats

- **`bmv_case_3` may fail strict mode** when the port uses a BLAS
  whose internal triangular-solve operation order differs from the
  reference. The expected value (`-1/7` etc.) is encoded as the
  canonical IEEE-754 representation of those rationals; some BLAS
  implementations produce a 1-ULP-off result on the input matrix
  for this specific case. Tolerance mode passes. This is documented
  as a known BLAS-reproducibility edge case rather than a defect in
  the spec or the F77 code.
- **Trajectory cases for `dcsrch` and `lnsrlb`** (full line searches)
  are not yet encoded as JSON. They require a multi-call replay
  protocol and are deferred to a future Phase C extension. Until
  then, line-search behavior is validated only at the single-call
  granularity (the 8 + 7 cases that exercise input validation and
  per-call branches).
- **`mainlb` integration cases** (full optimization runs) are
  similarly deferred. The Python reference impl is verified to
  converge on a quadratic and on bound-constrained Rosenbrock; full
  end-to-end JSON test vectors covering convergence trajectories
  are a future extension.
