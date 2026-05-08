# subsm

## Purpose

Compute an approximate solution of the **subspace minimization
problem** at the generalized Cauchy point: minimize the quadratic
model `Q(x) = r&#39;(x - xc) + (1/2)(x - xc)&#39; B (x - xc)` over the free
variables, subject to bounds.

Uses the L-BFGS Hessian approximation `B = theta*I - W*M^{-1}*W&#39;`. The
Newton direction in the free subspace is computed via the **K**-matrix
factorization (`formk` output). Then the iterate is moved along this
direction, **projected onto the box**.

If projection causes the step to violate sufficient-descent (positive
directional derivative), `subsm` falls back to a **safeguarded
backtracking step** (the Morales/Nocedal 2011 fix to the original
algorithm).

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n`, `m` | positive integers | Dimensions. |
| `nsub` | integer >= 0 | Number of free variables. |
| `ind` | integer vector, length `nsub` | Coordinate indices of free variables. |
| `l`, `u`, `nbd` | (as in `02_api.md`) | Bounds. |
| `xx`, `gg` | real vectors, length `n` | **Current iterate** and gradient at it (NOT `xc`). Used by the safeguarding check. |
| `ws`, `wy` | real matrices, `n * m` | L-BFGS history. |
| `wn` | real matrix, `2m * 2m` | LEL' factorization from `formk`. Upper triangle. |
| `theta` | positive real | Hessian scaling. |
| `col`, `head` | integers | L-BFGS layout. |
| `iprint` | integer | Diagnostic verbosity. |

### Logical outputs (in/out)

| Name | Type | Description |
|------|------|-------------|
| `x` (in/out) | real vector, length `n` | On entry: Cauchy point `xc`. On exit: subspace minimizer. |
| `d` (in/out) | real vector, length `nsub` | On entry: reduced gradient. On exit: Newton direction (then scaled). |
| `xp` (out) | real vector, length `n` | Saved copy of `xc` for backtracking restore. |
| `wv` (out) | real vector, length `2m` | Workspace: `K^{-1}*W&#39;*Z*d`. |
| `iword` | integer | `0` if minimizer is in the box; `1` if some bound was encountered. |

### Preconditions

- `nsub >= 0`. If `0`, immediate return (no free variables means
  nothing to minimize).
- `wn[:2col, :2col]` upper triangle holds the LEL' factorization
  from `formk`.
- `x` on entry is the GCP from `cauchy`.

### Postconditions

- `x` is the (possibly safeguarded) subspace minimizer, feasible.
- `iword = 1` if any bound was hit during the projection.

## Historical note

The F77 routine used to take an `info` output parameter for an
"ill-conditioned K" status. Since K's conditioning is checked in
`formk` / `formt` and the two `dtrsm` calls cannot fail on a
non-singular factor, the parameter was always 0 and has been removed.

## Algorithm

### Phase 1: Compute the Newton direction

```
# Step 1a: wv = W' Z d  (Z is selection of free variables in ind[])
pointr = head
for i = 1 to col:
    temp1 = sum over j=1..nsub: wy[ind[j], pointr] * d[j]
    temp2 = sum over j=1..nsub: ws[ind[j], pointr] * d[j]
    wv[i] = temp1
    wv[col + i] = theta * temp2
    pointr = next(pointr)

# Step 1b: wv := K^{-1} wv via the LEL' factorization stored in wn.
solve_upper_triangular_T(wn[:2col, :2col], wv[:2col])
for i = 1 to col:
    wv[i] = -wv[i]                        # the E-block sign flip
solve_upper_triangular_N(wn[:2col, :2col], wv[:2col])

# Step 1c: d = (1/theta) d + (1/theta^2) Z' W wv
pointr = head
for jy = 1 to col:
    js = col + jy
    for i = 1 to nsub:
        k = ind[i]
        d[i] = d[i] + wy[k, pointr] * wv[jy] / theta + ws[k, pointr] * wv[js]
    pointr = next(pointr)
d = d / theta                             # final division
```

### Phase 2: Project x = x + d onto bounds

```
xp = x       # save for potential restore
iword = 0
for i = 1 to nsub:
    k = ind[i]
    dk = d[i]; xk = x[k]
    if nbd[k] == 0:
        x[k] = xk + dk                    # free, take full step
    elif nbd[k] == 1:                     # lower only
        x[k] = max(l[k], xk + dk)
        if x[k] == l[k]: iword = 1
    elif nbd[k] == 2:                     # both bounds
        x[k] = min(u[k], max(l[k], xk + dk))
        if x[k] == l[k] or x[k] == u[k]: iword = 1
    elif nbd[k] == 3:                     # upper only
        x[k] = min(u[k], xk + dk)
        if x[k] == u[k]: iword = 1
```

### Phase 3: If `iword = 0`, done

```
if iword == 0: return
```

### Phase 4: Safeguarding (Morales/Nocedal 2011)

If a bound was hit, check the directional derivative against the
*current iterate's* gradient `gg`:

```
dd_p = sum over i=1..n: (x[i] - xx[i]) * gg[i]
if dd_p > 0:
    # Positive directional derivative -> bad step; restore and backtrack.
    x = xp
    # (write diagnostic warning to log)
    # fall through to alpha-scaling
else:
    return                                # accept the projected step
```

### Phase 5: Backtracking step (only if `dd_p > 0`)

Compute the largest `alpha in (0, 1]` such that `x_proj + alpha * d`
stays in the box:

```
alpha = 1.0
ibd = 0                               # index of binding constraint
for i = 1 to nsub:
    k = ind[i]; dk = d[i]
    if nbd[k] != 0:
        if dk < 0 and nbd[k] <= 2:     # toward lower
            temp2 = l[k] - x[k]
            if temp2 >= 0: temp1 = 0
            elif dk * alpha < temp2: temp1 = temp2 / dk
        elif dk > 0 and nbd[k] >= 2:   # toward upper
            temp2 = u[k] - x[k]
            if temp2 <= 0: temp1 = 0
            elif dk * alpha > temp2: temp1 = temp2 / dk
        if temp1 < alpha:
            alpha = temp1
            ibd = i

# Pin the binding variable at its bound; zero d[ibd] so it doesn't move further.
if alpha < 1.0:
    dk = d[ibd]; k = ind[ibd]
    if dk > 0:
        x[k] = u[k]; d[ibd] = 0
    elif dk < 0:
        x[k] = l[k]; d[ibd] = 0

# Take the safeguarded step.
for i = 1 to nsub:
    k = ind[i]
    x[k] = x[k] + alpha * d[i]
```

### Magic constants

None. The Morales/Nocedal fix introduces no thresholds; it is purely
algorithmic.

### Numerical safeguards

- `nsub == 0` early return.
- The Morales/Nocedal 2011 fix prevents the algorithm from accepting
  a projected step that doesn't actually decrease `f` (a defect in
  the original 1997 implementation).
- The exact-equality test `x[k] == l[k]` (line 245 in `subsm.f`) is
  intentional: it detects when the projection has *exactly* clamped
  to a bound, and only then sets `iword = 1`. Using a tolerance
  here would slightly change which iterates trigger the backtrack.

### Order-of-operations dependencies

- The reductions in Phase 1 (computing `wv` and the `d` update) sum
  in ascending index order. `--strict` requires this.
- Phase 2 iterates `i = 1..nsub` in ascending order; per-variable
  updates are independent so order doesn't affect numerics.
- Phase 5's `alpha`-shrinking iterates ascending; binds the *first*
  binding variable when ties (rare but possible).

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/subsm_case_1.json` | `nsub = 0`: immediate return |
| 2 | `data/subsm_case_2.json` | All free, no clip: `iword = 0` early-exit |
| 3 | `data/subsm_case_3.json` | `nbd = 2` clip to upper: `iword = 1`, `dd_p < 0`, accept |
| 4 | `data/subsm_case_4.json` | `nbd = 1` clip to lower |
| 5 | `data/subsm_case_5.json` | `nbd = 3` clip to upper only |
| 6 | `data/subsm_case_6.json` | `dd_p > 0`: triggers backtracking (Morales/Nocedal 2011 fix) |
| 7 | `data/subsm_case_7.json` | Backtrack with lower-bound clip |
| 8 | `data/subsm_case_8.json` | Backtrack: variable already at lower bound (`a2 >= 0`) |
| 9 | `data/subsm_case_9.json` | Backtrack: variable already at upper bound (`a2 <= 0`) |

## Reference implementation

`reference_impl/core/subsm.py`.

## Cross-references

- **Paper**: `acm-remark.pdf` (Morales/Nocedal 2011) for the
  safeguarding fix; original derivation in
  Byrd/Lu/Nocedal/Zhu 1995 sec.5.
- **Related subroutines**: called by `mainlb` after `cauchy`. Uses
  `formk`'s output `wn`. Does *not* call `bmv` directly (the K-matrix
  factorization handles the inverse).
- **F77 source**: `src/subsm.f`.
- **Unit test**: `tests/unit/test_subsm.f90` (9 case_* blocks).
