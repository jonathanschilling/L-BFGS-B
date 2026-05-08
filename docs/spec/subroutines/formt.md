# formt

## Purpose

Form the upper triangle of the symmetric positive-definite matrix

```
T = theta * S'S + L * D^{-1} * L'
```

(in the notation of the compact L-BFGS representation; see
`01_glossary.md`), then Cholesky-factorize it in place to produce the
upper-triangular factor `T_chol` such that `T_chol^T * T_chol = T`.

The factor is stored in the upper triangle of `wt`, replacing the
unfactored matrix. It is consumed by `bmv` for the inverse-middle-
matrix solves.

If `T` turns out to be non-positive-definite (which can happen when
the L-BFGS history is degenerate), the Cholesky factorization fails;
`formt` reports `info = -3` and the caller (`mainlb`) refreshes the
L-BFGS history.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `m` | positive integer | Memory parameter (leading dimension of `wt`, `sy`, `ss`). |
| `col` | integer, `0 <= col <= m` | Active size. |
| `theta` | positive real | Hessian scaling. |
| `sy` | real matrix, `m * m` | Holds `D` on the diagonal and `L` (strict lower); see `01_glossary.md`. |
| `ss` | real matrix, `m * m` | Holds `S^TS` upper triangle (symmetric). |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `wt` | real matrix, `m * m` | Output: upper triangle holds the Cholesky factor `T_chol` of `theta*S^TS + L*D^{-1}*L^T`. The strict lower triangle is unspecified. |
| `info` | integer | `0` on success; `-3` if the Cholesky factorization fails (matrix not positive-definite). |

### Preconditions

- `sy[i, i] > 0` for `i = 1, ..., col` (positive curvature; ensured
  by `matupd`).
- `theta > 0`.
- `ss` upper triangle holds `s_i^T s_j` for `j >= i`.

### Postconditions

- On success (`info = 0`): `wt` upper triangle holds `T_chol` such
  that `T_chol^T * T_chol` reproduces the original `theta*S^TS + L*D^{-1}*L^T`.
- On failure (`info = -3`): `wt` contents are *implementation-defined*.
  Caller must not use them.

## Algorithm

### Phase 1: build upper triangle of `T = theta*S^TS + L*D^{-1}*L^T`

The diagonal of `T` is `T[i, i] = theta * s_i^T s_i + sum_{k<i} s_i^T y_k * (s_i^T y_k) / s_k^T y_k`.

The upper triangle for `j >= i` is `T[i, j] = theta * s_i^T s_j + sum_{k<i} s_i^T y_k * s_j^T y_k / s_k^T y_k`.

Special case: for `i = 1`, the inner sum is empty (no `k < 1`), so the
entire first row is just `theta * ss[1, :]`.

```
# First row: i = 1
for j = 1 to col:
    wt[1, j] = theta * ss[1, j]

# Remaining rows: i = 2..col, columns j = i..col
for i = 2 to col:
    for j = i to col:
        k1 = min(i, j) - 1                 # = i - 1 (since j >= i)
        sum = 0
        for k = 1 to k1:
            sum += sy[i, k] * sy[j, k] / sy[k, k]
        wt[i, j] = sum + theta * ss[i, j]
```

The `min(i, j) - 1` simplifies to `i - 1` for `j >= i`. The F77 source
keeps the `min` form for clarity.

### Phase 2: Cholesky factorization

Replace the upper-triangular content of `wt[1..col, 1..col]` with its
upper Cholesky factor:

```
call dpotrf('U', col, wt, m, info_lapack)
if info_lapack != 0:
    info = -3
else:
    info = 0
```

LAPACK's `dpotrf` with `'U'` consumes the upper triangle of the input
and writes the upper triangular Cholesky factor in place. On failure
(non-positive-definite or numerical breakdown), `dpotrf` returns
`info > 0`; `formt` overrides this with `info = -3`.

### Pseudocode

See above. Combined:

```
input:  m, col, theta, sy, ss
output: wt (upper triangle, [1..col, 1..col]), info

for j = 1 to col:
    wt[1, j] = theta * ss[1, j]

for i = 2 to col:
    for j = i to col:
        sum = 0
        for k = 1 to i - 1:
            sum += sy[i, k] * sy[j, k] / sy[k, k]
        wt[i, j] = sum + theta * ss[i, j]

cholesky_upper_in_place(wt[1..col, 1..col], info_lapack)
info = -3 if info_lapack != 0 else 0
```

### Magic constants

None.

### Numerical safeguards

- The Cholesky factorization is the only failure mode. Failure
  triggers the caller's L-BFGS refresh.
- `sy[k, k] > 0` is required (curvature condition). If violated,
  `bmv`/`formt` are undefined.

### Order-of-operations dependencies

- The first loop builds row 1 in ascending `j` order. Result is
  order-independent (each entry written once).
- The double loop fills `wt[i, j]` for `i >= 2`, `j >= i`. Outer in
  ascending `i`; inner in ascending `j`. The accumulation
  `sum += sy[i, k] * sy[j, k] / sy[k, k]` iterates `k = 1..i-1` in
  ascending order; `--strict` requires this.
- `dpotrf` uses the LAPACK reference build for `--strict`.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/formt_case_1.json` | `col = 0`: all loops empty, dpotrf with N=0 returns info=0 |
| 2 | `data/formt_case_2.json` | `col = 1`: only first-row loop fires |
| 3 | `data/formt_case_3.json` | `col = 2`: both loops; T is diagonal by construction |
| 4 | `data/formt_case_4.json` | `theta < 0`: dpotrf fails; `info = -3` |

## Reference implementation

`reference_impl/core/formt.py`.

## Cross-references

- **Paper**: Byrd/Nocedal/Schnabel 1994 sec.3 (compact representation),
  `code.pdf` sec.2.2.
- **Related subroutines**: called by `mainlb` once per L-BFGS update;
  output `wt` is consumed by `bmv`. Failure (`info = -3`) triggers
  history refresh.
- **F77 source**: `src/formt.f`.
- **Unit test**: `tests/unit/test_formt.f90` (4 case_* blocks).
