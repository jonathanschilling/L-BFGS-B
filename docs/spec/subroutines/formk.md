# formk

## Purpose

Form the **`LEL'`-factorization** of the `2col * 2col` indefinite
matrix `K` used in subspace minimization (`subsm`):

```
K = [-D - Y'ZZ'Y/theta     L_a' - R_z'  ]
    [L_a - R_z         theta * S'AA'S   ]
```

where `Z` is the projection onto free variables, `A` onto active,
`L_a` is the strict-lower triangle of `S'AA'Y`, and `R_z` is the
upper triangle of `S'ZZ'Y`. (See `code.pdf` sec.4 and Byrd/Lu/Nocedal/Zhu
1995 sec.5.)

The factorization is incrementally maintained across iterations via
the helper matrix `wn1`, which holds (updated) inner-products that
feed the current `wn`. Two Cholesky factorizations occur (one per
diagonal block); failures are reported via `info = -1` (first block)
or `info = -2` (second block).

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Problem dimension. |
| `nsub` | integer >= 0 | Number of currently-free variables. |
| `ind` | integer vector, length `n` | `ind[1..nsub]` = currently-free indices, `ind[nsub+1..n]` = currently-active. |
| `nenter`, `ileave` | integers | From `freev`: change-detection counters. |
| `indx2` | integer vector, length `n` | From `freev`: entered (`[1..nenter]`) and leaving (`[ileave..n]`) variables. |
| `iupdat` | integer >= 0 | Total L-BFGS updates so far. |
| `updatd` | boolean | Whether L-BFGS history was just updated. |
| `m`, `col`, `head` | integers | Memory layout parameters. |
| `theta` | positive real | Hessian scaling. |
| `ws`, `wy` | real matrices, `n * m` | L-BFGS history columns. |
| `sy` | real matrix, `m * m` | Compact representation (diagonal `D`, lower `L`). |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `wn` | real matrix, `2m * 2m` | Upper triangle holds the `LEL'` factorization. |
| `wn1` (in/out) | real matrix, `2m * 2m` | Lower triangle of the auxiliary inner-product matrix (carried across iterations). |
| `info` | integer | `0` success; `-1` first Cholesky failed; `-2` second Cholesky failed. |

### Preconditions

- On the first iteration, `wn1` must be zeroed.
- `ind` and `indx2` are produced by `freev` for the current iteration.
- The first call must have `updatd = false` (no L-BFGS pair to add yet)
  unless `iupdat >= 1` is genuinely true.

### Postconditions

- On success: `wn[1..2col, 1..2col]` upper triangle is the `LEL'`
  factorization (Cholesky factor of (1,1), then triangular-solve
  block, then Cholesky of (2,2)).
- On failure: `wn` contents through the failed step are
  implementation-defined; `info` indicates which step failed.

## Algorithm

The algorithm has three phases.

### Phase 1: update `wn1` (the auxiliary inner-product cache)

#### 1a. If `updatd` and `iupdat > m`: shift `wn1`

The L-BFGS buffer wrapped, so the oldest pair is gone. Shift the
relevant blocks of `wn1` left and up by one to drop the oldest
column / row.

```
if updatd and iupdat > m:
    for jy = 1 to m - 1:
        js = m + jy
        # block (1,1): wn1[jy..m-1, jy] = wn1[jy+1..m, jy+1]
        for k = 0 to m - jy - 1:
            wn1[jy + k, jy] = wn1[jy + 1 + k, jy + 1]
        # block (2,2): wn1[js..2m-1, js] = wn1[js+1..2m, js+1]
        for k = 0 to m - jy - 1:
            wn1[js + k, js] = wn1[js + 1 + k, js + 1]
        # block (2,1) col `jy` from row m+2 down: shift
        for k = 0 to m - 2:
            wn1[m + 1 + k, jy] = wn1[m + 2 + k, jy + 1]
```

#### 1b. If `updatd`: fill new last row / column of `wn1`

The new L-BFGS pair contributes new inner-product values for the
last row/column of the three blocks.

```
if updatd:
    # New row `col` of (1,1) and (2,2) blocks via inner products on free / active sets.
    iy = col, is = m + col
    ipntr = head + col - 1   (cyclic, mod m)
    jpntr = head
    for jy = 1 to col:
        js = m + jy
        temp1 = sum over k in free indices of wy[k1, ipntr] * wy[k1, jpntr]
        temp2 = sum over k in active indices of ws[k1, ipntr] * ws[k1, jpntr]
        temp3 = sum over k in active indices of ws[k1, ipntr] * wy[k1, jpntr]
        wn1[iy, jy] = temp1
        wn1[is, js] = temp2
        wn1[is, jy] = temp3
        jpntr = next(jpntr)

    # New column `col` of (2,1) block, R_z part (S'ZZ'Y on free indices).
    jy = col, jpntr = head + col - 1 (cyclic)
    ipntr = head
    for i = 1 to col:
        is = m + i
        temp3 = sum over k in free indices of ws[k1, ipntr] * wy[k1, jpntr]
        wn1[is, jy] = temp3
        ipntr = next(ipntr)

    upcl = col - 1                 # don't re-modify the just-written column
else:
    upcl = col                     # all columns may have changed via the active-set drift
```

### Phase 2: modify `wn1` for active-set changes

Variables that *entered* the free set add to the (1,1) block (their
`Y'ZZ'Y` contribution wasn't in there before) and subtract from the
(2,2) block (their `S'AA'S` contribution was). Variables that *left*
have the opposite.

```
ipntr = head
for iy = 1 to upcl:
    is = m + iy
    jpntr = head
    for jy = 1 to iy:
        js = m + jy
        temp1, temp2, temp3, temp4 = 0
        for k = 1 to nenter:
            k1 = indx2[k]
            temp1 += wy[k1, ipntr] * wy[k1, jpntr]
            temp2 += ws[k1, ipntr] * ws[k1, jpntr]
        for k = ileave to n:
            k1 = indx2[k]
            temp3 += wy[k1, ipntr] * wy[k1, jpntr]
            temp4 += ws[k1, ipntr] * ws[k1, jpntr]
        wn1[iy, jy] += temp1 - temp3       # (1,1): entering add, leaving subtract
        wn1[is, js] += -temp2 + temp4      # (2,2): entering subtract, leaving add
        jpntr = next(jpntr)
    ipntr = next(ipntr)

# Modify (2,1) block similarly.
ipntr = head
for is = m+1 to m+upcl:
    jpntr = head
    for jy = 1 to upcl:
        temp1, temp3 = 0
        for k = 1 to nenter:
            temp1 += ws[k1, ipntr] * wy[k1, jpntr]
        for k = ileave to n:
            temp3 += ws[k1, ipntr] * wy[k1, jpntr]
        if is <= jy + m:
            wn1[is, jy] += temp1 - temp3
        else:
            wn1[is, jy] += -temp1 + temp3
        jpntr = next(jpntr)
    ipntr = next(ipntr)
```

### Phase 3: form `wn` and Cholesky-factorize

```
# Build wn upper triangle from wn1 / sy / theta.
for iy = 1 to col:
    is  = col + iy
    is1 = m + iy
    for jy = 1 to iy:
        js  = col + jy
        js1 = m + jy
        wn[jy, iy] = wn1[iy, jy] / theta              # (1,1) block scaled
        wn[js, is] = wn1[is1, js1] * theta            # (2,2) block scaled
    # (1,2) block: -L_a' for jy < iy, R_z for jy >= iy.
    for jy = 1 to iy - 1:
        wn[jy, is] = -wn1[is1, jy]
    for jy = iy to col:
        wn[jy, is] = wn1[is1, jy]
    wn[iy, iy] += sy[iy, iy]                          # add D diagonal

# First Cholesky on the (1,1) block.
dpotrf('U', col, wn, 2*m, info_lapack)
if info_lapack != 0: info = -1; return

# Triangular solve: form L^{-1}(-L_a' + R_z') in (1,2) block.
for js = col+1 to 2*col:
    dtrsm('l','u','t','n', col, 1, 1.0, wn, 2*m, wn[1, js], col)

# Form S'AA'S*theta + (...)' (...) in upper (2,2) block.
for is = col+1 to 2*col:
    for js = is to 2*col:
        wn[is, js] += dot_product(wn[1..col, is], wn[1..col, js])

# Second Cholesky on the (2,2) block.
dpotrf('U', col, wn[col+1, col+1], 2*m, info_lapack)
if info_lapack != 0: info = -2; return

info = 0
```

### Magic constants

None.

### Numerical safeguards

- Two Cholesky factorizations; either can fail (`info = -1` or `-2`).
  Caller refreshes the L-BFGS history on failure.
- The `theta`-scaled diagonal addition prevents the (1,1) block from
  collapsing in degenerate cases.

### Order-of-operations dependencies

- Many nested loops with reductions. **Each reduction loop iterates
  `k` in ascending order** (over free, active, or `1..nenter` /
  `ileave..n` indexed sets in the order specified). `--strict`
  conformance requires this.
- The Cholesky and triangular-solve calls follow BLAS reference
  conventions.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/formk_case_1.json` | `updatd = false`, no entering/leaving: phases 1a/1b/2 are no-ops |
| 2 | `data/formk_case_2.json` | `updatd = true`, `iupdat = 1`: phase 1b fills new row/col |
| 3 | `data/formk_case_3.json` | First Cholesky fails (non-PD `(1,1)` block) -> `info = -1` |

## Reference implementation

`reference_impl/core/formk.py`.

## Cross-references

- **Paper**: Byrd/Lu/Nocedal/Zhu 1995 sec.5 (subspace minimization),
  `code.pdf` sec.4.
- **Related subroutines**: called by `mainlb` per iteration when
  `wrk = true` (active set changed or L-BFGS updated). Output `wn`
  consumed by `subsm`.
- **F77 source**: `src/formk.f`.
- **Unit test**: `tests/unit/test_formk.f90` (3 case_* blocks).
