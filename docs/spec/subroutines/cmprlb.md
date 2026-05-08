# cmprlb

## Purpose

Compute the reduced gradient `r` at the Cauchy point, restricted to
the currently-free variables. The reduced gradient is
`r = -B (xc - x) - g` (componentwise), evaluated only at the
indices in the free-variable index permutation.

`r` is the right-hand side for the subspace minimization problem in
`subsm`.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Number of variables. |
| `m` | positive integer | Memory parameter. |
| `col` | integer, `0 ≤ col ≤ m` | Number of stored L-BFGS pairs. |
| `head` | integer, `1 ≤ head ≤ m` | Index of oldest pair in the cyclic `S`/`Y` buffer. |
| `nfree` | integer, `0 ≤ nfree ≤ n` | Number of free variables. |
| `index` | integer vector, length `n` | Free/active permutation: `index[1..nfree]` are free indices. |
| `theta` | positive real | Hessian scaling. |
| `x`, `z` | real vectors, length `n` | Current iterate and Cauchy point. |
| `g` | real vector, length `n` | Gradient at `x`. |
| `ws`, `wy` | real matrices, `n × m` | L-BFGS history columns (`s`-vectors and `y`-vectors). |
| `sy`, `wt` | real matrices, `m × m` | Compact representation: `S'Y` packed and Cholesky factor `T`. |
| `wa` | real vector, length `4m` | **In/out**: on entry, `wa[2m+1..2m+2col]` holds `W'(z - x)` (filled by `cauchy`); on exit, `wa[1..2col]` holds `M⁻¹ W' (z - x)`. |
| `cnstnd` | boolean | True if any variable is bounded. |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `r` | real vector, length `n` | Filled in two regimes (see Algorithm): full `n`-vector when unconstrained with history, or `r[1..nfree]` for the free variables otherwise. |
| `wa` | (in-place) | `wa[1..2col]` filled with `M⁻¹ W'(z - x)`. |

### Preconditions

- `wa[2m+1..2m+2col]` contains `W' (z - x)` on entry. The first `col`
  entries are `Y' (z - x)`; the next `col` are `theta * S' (z - x)`.
  Caller (`cauchy`) sets these.
- `index[1..nfree]` are the indices of variables free at the Cauchy
  point.

### Postconditions

- `r` filled per the algorithm.
- `wa[1..2col]` filled.

## Historical note

The F77 routine used to take an `info` output parameter to forward
errors from the embedded `bmv` call. Since `bmv` cannot fail under
LAPACK `dtrsm`, the parameter was always 0 and has been removed.

## Algorithm

There are two paths.

### Path A: unconstrained with history (`!cnstnd` and `col > 0`)

```
for i = 1 to n:
    r[i] = -g[i]
return    # info unchanged
```

The reasoning: when there are no bounds, the Cauchy point equals the
current iterate (`z = x`), so `B (z - x) = 0` and `r = -g`. This
shortcut avoids the `bmv` call.

### Path B: constrained or `col = 0`

```
# 1. Build r on free variables: r[i] = -theta*(z[k] - x[k]) - g[k]
for i = 1 to nfree:
    k = index[i]
    r[i] = -theta * (z[k] - x[k]) - g[k]

# 2. Solve M^{-1} v for v = W'(z - x), already in wa[2m+1..2m+2col].
#    Result placed in wa[1..2col] by bmv.
call bmv(m, sy, wt, col, wa[2m+1..], wa[1..], info)
# 3. Apply the W*M^{-1}*W' correction column by column.
pointr = head
for j = 1 to col:
    a1 = wa[j]                         # j-th component of M^{-1} v on Y side
    a2 = theta * wa[col + j]           # theta-scaled j-th of M^{-1} v on S side
    for i = 1 to nfree:
        k = index[i]
        r[i] += wy[k, pointr] * a1 + ws[k, pointr] * a2
    pointr = (pointr mod m) + 1        # cyclic walk in chronological order
```

### Pseudocode (combined)

```
input:  n, m, col, head, nfree, index, theta, x, g, z, ws, wy, sy, wt,
        wa (in/out), cnstnd
output: r, wa[1..2col] mutated

if not cnstnd and col > 0:
    for i = 1 to n:
        r[i] = -g[i]
    return

# Unconstrained without history (col=0) takes Path B, since for col=0
# there is no Hessian information yet -- B = theta*I and the second
# loop just zeros out (since nothing is added). The result is still
# r[i] = -theta*(z[k]-x[k]) - g[k] for free k; with z = x in
# unconstrained problems, this reduces to -g[k].

for i = 1 to nfree:
    k = index[i]
    r[i] = -theta * (z[k] - x[k]) - g[k]

call bmv(m, sy, wt, col, wa[2m+1 ..], wa[1 ..])
pointr = head
for j = 1 to col:
    a1 = wa[j]
    a2 = theta * wa[col + j]
    for i = 1 to nfree:
        k = index[i]
        r[i] = r[i] + wy[k, pointr] * a1 + ws[k, pointr] * a2
    pointr = (pointr mod m) + 1
```

### Magic constants

None.

### Numerical safeguards

None. `bmv` cannot fail under the LAPACK `dtrsm` build, so the F77
routine no longer signals an error path.

### Order-of-operations dependencies

- The first loop (i = 1..nfree) is order-independent in its writes (each
  `r[i]` depends only on `i`).
- The accumulation loop `for j = 1..col` walks the columns in
  chronological order (`pointr = head, head+1, ..., wrapping`). The
  inner `for i = 1..nfree` accumulates into `r[i]` *across* `j`. The
  order of `j` matters for bit-for-bit summation; ports must walk
  `j = 1..col` and update `pointr` cyclically.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/cmprlb_case_1.json` | Unconstrained with history: shortcut `r = -g` |
| 2 | `data/cmprlb_case_2.json` | Constrained, `col = 0`: j-loop empty, only first loop runs |
| 3 | `data/cmprlb_case_3.json` | Constrained, `col = 1`, `z = x`: bmv sees zero input, j-loop adds zero |
| 4 | `data/cmprlb_case_4.json` | Constrained, non-trivial step: hand-derived `r = -B(z-x) - g` |

## Reference implementation

`reference_impl/core/cmprlb.py` (depends on `core.bmv`).

## Cross-references

- **Paper**: `code.pdf` §4 (subspace minimization setup). The reduced
  gradient is the residual for the subspace problem of
  Byrd/Lu/Nocedal/Zhu 1995 §5.
- **Related subroutines**: called by `mainlb` after `cauchy`. Calls
  `bmv` for the inverse-middle-matrix solve. Output `r` feeds
  `subsm`.
- **F77 source**: `src/cmprlb.f`.
- **Unit test**: `tests/unit/test_cmprlb.f90` (4 case_* blocks).
