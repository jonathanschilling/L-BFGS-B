# bmv

## Purpose

Compute `p = M⁻¹ v`, where `M` is the `(2col × 2col)` middle matrix
of the compact L-BFGS representation:

```
M⁻¹ = [[-D,        L' ],
       [ L,    theta·S'S]]
```

The matrix is never materialized; the solve is done blockwise using
the diagonal `D = diag(s_i' y_i)`, the strict-lower L
(`L_{i,j} = s_i' y_j` for `i > j`), and the upper-triangular Cholesky
factor `T` (called `wt` in the F77 source) of
`theta · S'S + L · D⁻¹ · L'`.

This is the workhorse of the compact L-BFGS evaluation: every
Hessian-vector product (and inverse) for the L-BFGS approximation
goes through `bmv`.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `m` | positive integer | Memory parameter (used for the leading dimension of `sy` and `wt`). |
| `col` | integer, `0 ≤ col ≤ m` | Number of stored L-BFGS pairs. |
| `sy` | real matrix, `m × m`, lower-triangle-and-diagonal | `sy[i,i] = s_i' y_i = D[i]`; `sy[i,j] = s_i' y_j = L[i,j]` for `i > j`. Strict upper triangle unused. |
| `wt` | real matrix, `m × m`, upper-triangular | Cholesky factor `T` of `theta·S'S + L·D⁻¹·L'`, computed by `formt`. |
| `v` | real vector, length `2 col` | Right-hand side. Layout: `v[1..col]` is the "Y-block" component; `v[col+1..2col]` is the "S-block" component. |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `p` | real vector, length `2 col` | Result `p = M⁻¹ v`, same block layout as `v`. |

### Preconditions

- `sy[i,i] > 0` for `i = 1, ..., col` (positive curvature; ensured
  by `matupd`).
- `wt` is the upper Cholesky factor of a positive-definite matrix
  (ensured by `formt`).

### Postconditions

- `p` satisfies `M p = v` to within numerical precision.
- `v`, `sy`, `wt` are unchanged.

## Historical note

The F77 routine used to take an `info` output parameter for the
LINPACK `dtrsl` triangular-solve return status. Since LAPACK's
`dtrsm` cannot fail on a non-singular factor (and `formt` ensures
`wt` is non-singular), the parameter was always 0 and has been
removed.

## Algorithm

The factorization of `M⁻¹` used here:

```
M = [[-D,    L'],
     [ L,  thetaSS]]^(-1)

= [[D^(-1/2),         0],   [[-I,  D^(-1/2) L'   ],   [[D^(1/2),       0 ],
   [-L D^(-1/2),  J']]    *  [ 0,         I       ]] *  [-L D^(-1/2),   J ]]
                                                                              ^(-1)

(or the equivalent block-LDU; what matters is the
 sequence of triangular solves below.)
```

Equivalently, the F77 algorithm performs two block solves:

### Part I — solve the lower triangular block

```
[ D^{1/2}       0 ] [ p1 ]   [ v1 ]
[ -L D^{-1/2}   J ] [ p2 ] = [ v2 ]
```

In equations:
- `D^{1/2} p1 = v1` ⟹ `p1 = v1 / sqrt(D)`.
- `-L D^{-1/2} p1 + J p2 = v2` ⟹ `J p2 = v2 + L D^{-1/2} (v1 / sqrt(D)) = v2 + L D^{-1} v1`.

Pseudocode for Part I:

```
# Step 1: build (v2 + L D^{-1} v1) into p2 = p[col+1..2col].
p[col + 1] = v[col + 1]                               # nothing to add for i=1
for i = 2 to col:
    s = 0
    for k = 1 to i - 1:                               # accumulate L D^{-1} v1
        s = s + sy[i, k] * v[k] / sy[k, k]
    p[col + i] = v[col + i] + s

# Step 2: solve J p2 = (above) for p2 via triangular solve with J = wt.
solve J^T x = p[col+1..2col] for x                    # 'transpose' upper-triangular solve
p[col+1..2col] = x

# Step 3: solve D^{1/2} p1 = v1 for p1.
for i = 1 to col:
    p[i] = v[i] / sqrt(sy[i, i])
```

The "transpose upper" triangular solve is the F77 call
`dtrsm('l','u','t','n', col, 1, 1.0, wt, m, p(col+1), col)`.

### Part II — solve the upper triangular block

```
[ -D^{1/2}   D^{-1/2} L' ] [ p1 ]   [ p1 ]    (in-place; the operand
[ 0           J'         ] [ p2 ] = [ p2 ]     is the result of Part I)
```

In equations:
- `J' p2 = p2_pre` ⟹ solve via triangular solve with `J' = wt'`,
  no-transpose, upper.
- `-D^{1/2} p1 + D^{-1/2} L' p2 = p1_pre` ⟹
  `p1 = -D^{-1/2} (p1_pre - D^{-1/2} L' p2)`
       `= -D^{-1/2} p1_pre + D^{-1} L' p2`.

Pseudocode for Part II:

```
# Step 1: solve J^T p2 = p2 in place.
solve J x = p[col+1..2col] for x                     # 'no-transpose' upper-triangular solve
p[col+1..2col] = x

# Step 2: p1 = -D^{-1/2} p1_pre.
for i = 1 to col:
    p[i] = -p[i] / sqrt(sy[i, i])

# Step 3: add the D^{-1} L' p2 correction.
for i = 1 to col:
    s = 0
    for k = i + 1 to col:                            # L' = transpose of L; stored L_{k,i} for k > i
        s = s + sy[k, i] * p[col + k] / sy[i, i]
    p[i] = p[i] + s
```

The "no-transpose upper" triangular solve is the F77 call
`dtrsm('l','u','n','n', col, 1, 1.0, wt, m, p(col+1), col)`.

### Combined pseudocode

```
input:  m, col, sy, wt, v
output: p (length 2*col)

if col == 0:
    return    # p untouched

# Part I
p[col+1] = v[col+1]
for i = 2 to col:
    s = 0
    for k = 1 to i-1:
        s += sy[i, k] * v[k] / sy[k, k]
    p[col+i] = v[col+i] + s
solve_upper_triangular_T(wt[1..col, 1..col], p[col+1..2col])    # in-place
for i = 1 to col:
    p[i] = v[i] / sqrt(sy[i, i])

# Part II
solve_upper_triangular(wt[1..col, 1..col], p[col+1..2col])      # in-place
for i = 1 to col:
    p[i] = -p[i] / sqrt(sy[i, i])
for i = 1 to col:
    s = 0
    for k = i+1 to col:
        s += sy[k, i] * p[col+k] / sy[i, i]
    p[i] = p[i] + s

```

### Magic constants

None.

### Numerical safeguards

- `sqrt(sy[i, i])` requires `sy[i, i] > 0`. Caller (`matupd`) ensures
  this via the curvature condition `s' y > eps · ||y||²`. If a port
  receives `sy[i, i] ≤ 0` (corrupt state), `bmv` is undefined.
- The triangular solves require `wt[i, i] ≠ 0`. Caller (`formt`)
  ensures this; if the Cholesky factorization in `formt` fails, the
  caller refreshes the L-BFGS history rather than calling `bmv`.

### Order-of-operations dependencies

- The accumulation loops (`for k = 1 to i - 1` in Part I, `for k =
  i + 1 to col` in Part II) sum in ascending `k` order. Reversing
  changes the floating-point result. `--strict` conformance requires
  ascending `k`.
- The two triangular solves (`dtrsm` and equivalents) follow BLAS
  conventions; ports should use their language's BLAS binding for
  these (linked to the reference build for `--strict`).

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/bmv_case_1.json` | `col = 0` early return |
| 2 | `data/bmv_case_2.json` | `col = 1` diagonal-only path |
| 3 | `data/bmv_case_3.json` | `col = 2` with non-trivial L block |

## Reference implementation

`reference_impl/core/bmv.py`.

## Cross-references

- **Paper**: Byrd/Nocedal/Schnabel 1994 §3 (compact L-BFGS
  representation), `code.pdf` §2.2.
- **Related subroutines**: called by `cmprlb`, `cauchy`, `subsm`. Built
  on top of the Cholesky factor `T` produced by `formt`.
- **F77 source**: `src/bmv.f`.
- **Unit test**: `tests/unit/test_bmv.f90` (3 case_* blocks).
