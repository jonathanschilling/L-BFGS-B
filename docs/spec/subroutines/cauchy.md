# cauchy

## Purpose

Compute the **Generalized Cauchy Point** (GCP): the first local
minimizer of the quadratic model `Q(x + s) = g^T*s + (1/2)*s^T*B*s`
along the projected steepest-descent path `P(x - t*g, l, u)` for
`t >= 0`.

The GCP determines which bounds become active in this iteration; the
free variables at the GCP form the subspace for `subsm`.

The algorithm walks **breakpoints** along the projected gradient
direction in order of increasing `t`. At each breakpoint, one
variable hits a bound and is "fixed"; between breakpoints, the
quadratic is minimized analytically. The minimizer in the current
piecewise-linear segment is the GCP. (See Byrd/Lu/Nocedal/Zhu 1995
sec.4 for the derivation.)

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Dimension. |
| `x`, `l`, `u`, `nbd` | (as in `02_api.md`) | Iterate, bounds, types. |
| `g` | real vector, length `n` | Gradient (must be nonzero; checked via `sbgnrm > 0`). |
| `m`, `col`, `head` | integers | L-BFGS memory layout. |
| `theta` | positive real | Hessian scaling. |
| `ws`, `wy` | real matrices, `n * col` | L-BFGS history. |
| `sy`, `wt` | real matrices, `m * m` | Compact representation. |
| `iwhere` (in/out) | integer vector, length `n` | Coarse bound status from `active`. Refined by `cauchy`. |
| `sbgnrm` | nonneg real | Projected-gradient infinity norm (precomputed by `projgr`). |
| `epsmch` | real | Machine epsilon. |
| `iprint` | integer | Diagnostic verbosity. |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `xcp` | real vector, length `n` | Generalized Cauchy point. |
| `d` | real vector, length `n` | Cauchy direction `P(x - t*g, l, u) - x` at the GCP. |
| `iwhere` (in/out) | integer vector | Refined: see `01_glossary.md`. |
| `iorder` | integer vector, length `n` | Permutation of variables: free, bound-encountered, bound-not-encountered. |
| `t` | real vector, length `n` | Workspace: breakpoints (used for `hpsolb`). |
| `p` | real vector, length `2col` | `W^T d`. |
| `c` | real vector, length `2col` | `W^T (xcp - x)`. |
| `wbp`, `v` | real vectors, length `2m` | Workspace. |
| `nseg` | integer | Number of piecewise-quadratic segments traversed. |

### Preconditions

- `sbgnrm >= 0`. If `0`, the routine returns immediately with `xcp = x`.
- `iwhere` set by `active` initially (only `-1`, `0`, `3`).
- L-BFGS state (`ws`, `wy`, `sy`, `wt`, `theta`, `col`, `head`)
  consistent (verified by `formt`/`bmv`).

### Postconditions

- `xcp` is the generalized Cauchy point (within numerical precision).
- `iwhere` updated: variables fixed at bounds during the walk get
  `iwhere = 1` (lower) or `iwhere = 2` (upper); free movers get
  `iwhere = 0`; free non-movers (zero gradient component) get
  `iwhere = -3`.

## Historical note

The F77 routine used to take an `info` output parameter to forward
errors from the embedded `bmv` calls. Since `bmv` cannot fail under
LAPACK `dtrsm`, the parameter was always 0 and has been removed.

## Algorithm

### Phase 1: setup -- bound status and breakpoints

For each variable `i`:

```
neggi = -g[i]
if iwhere[i] not in {3, -1}:                        # has bounds and not fixed
    tl = x[i] - l[i] (if has lower)
    tu = u[i] - x[i] (if has upper)
    xlower = (has lower) and tl <= 0                # at lower bound
    xupper = (has upper) and tu <= 0                # at upper bound

    iwhere[i] = 0
    if xlower: if neggi <= 0: iwhere[i] = 1         # at-lower, can't move down
    elif xupper: if neggi >= 0: iwhere[i] = 2       # at-upper, can't move up
    else: if |neggi| <= 0: iwhere[i] = -3           # zero gradient, can't move

if iwhere[i] in {0, -1}:                            # variable will move
    d[i] = neggi
    f1 -= neggi * neggi
    p += -W[i] * g[i]                                # accumulate W' d
    if has bounded breakpoint along d:
        compute t[i] = (l[i] - x[i]) / d[i] or (u[i] - x[i]) / d[i]
        track smallest breakpoint
    else:
        nfree -= 1
        iorder[nfree] = i
        bnded = false (unbounded direction available)
else:
    d[i] = 0
```

After the loop, `iorder[1..nbreak]` lists variables with breakpoints
in *unsorted* order; `iorder[nfree..n]` lists variables with no
breakpoint along `d`.

### Phase 2: initial GCP candidate

Initialize:

```
xcp = x
c = 0
f2 = -theta * f1                                    # before W'd correction
if col > 0:
    v = M^{-1} p  (via bmv)
    f2 -= p . v                                     # final initial curvature
dtm = -f1 / f2                                      # candidate step

if nbreak == 0:                                     # no breakpoints
    if nfree == n+1:                                # nothing moves
        return with xcp = x
    else:
        skip to phase 4 (locate GCP in this segment)
```

### Phase 3: piecewise-quadratic walk

Loop over breakpoints, smallest first:

```
nleft = nbreak
tj = 0
iter = 1
loop:
    tj0 = tj
    if iter == 1:                                   # use the precomputed smallest
        tj = bkmin
        ibp = iorder[ibkmin]
    else:
        if iter == 2: replace iorder[ibkmin] with iorder[nbreak] (heap setup)
        call hpsolb(nleft, t, iorder, iter - 2)     # extract next-smallest
        tj = t[nleft]
        ibp = iorder[nleft]

    dt = tj - tj0
    if dtm < dt:                                    # GCP lies in current segment
        break

    # Otherwise fix variable ibp at its bound.
    tsum += dt
    nleft -= 1
    iter += 1
    dibp = d[ibp]
    d[ibp] = 0
    if dibp > 0:
        zibp = u[ibp] - x[ibp]; xcp[ibp] = u[ibp]; iwhere[ibp] = 2
    else:
        zibp = l[ibp] - x[ibp]; xcp[ibp] = l[ibp]; iwhere[ibp] = 1

    if nleft == 0 and nbreak == n:                  # all fixed; stop
        dtm = dt; goto exit

    # Update derivatives f1, f2 across the segment.
    f1 += dt * f2 + dibp**2 - theta * dibp * zibp
    f2 -= theta * dibp**2
    if col > 0:
        c += dt * p
        wbp = W[ibp] (with theta scaling on the S half)
        v = M^{-1} wbp  (via bmv)
        wmc = c . v
        wmp = p . v
        wmw = wbp . v
        p -= dibp * wbp
        f1 += dibp * wmc
        f2 += 2 * dibp * wmp - dibp**2 * wmw

    f2 = max(epsmch * f2_org, f2)                   # numerical guard
    if nleft > 0:
        dtm = -f1 / f2
        continue loop
    elif bnded:                                     # all bounded, no unbounded direction
        f1 = 0; f2 = 0; dtm = 0
    else:
        dtm = -f1 / f2

# Phase 4: locate GCP in the final (or only) segment.
if dtm <= 0: dtm = 0
tsum += dtm
xcp += tsum * d                                     # apply remaining motion to free vars
if col > 0: c += dtm * p                            # update W^T (xcp - x)
return
```

### Magic constants

| Constant | F77 | Value | Meaning |
|----------|-----|-------|---------|
| `epsmch * f2_org` floor | (loop) | `eps * f2_org` | Lower bound on `f2` to avoid division by tiny positive number when curvature collapses. |

### Numerical safeguards

- `sbgnrm = 0` early return: the projected gradient is zero, so `x` is
  already a Cauchy point.
- The `f2 = max(epsmch * f2_org, f2)` guard prevents `dtm` from
  blowing up when the quadratic curvature `f2` collapses near zero
  due to numerical noise.
- The `bnded` flag detects an unbounded descent direction (variable
  with no bound along `d` and nonzero `d`). When `bnded = false` and
  all breakpoints are exhausted, the GCP lies along the unbounded
  ray; `dtm = -f1/f2` proceeds normally.

### Order-of-operations dependencies

- The setup loop (`for i = 1..n`) iterates in ascending order. `f1`
  accumulates and `p` accumulates (both reductions); `--strict`
  conformance requires this order.
- The main loop processes breakpoints in increasing `t` order via
  `hpsolb`. The first breakpoint is the precomputed `bkmin`; subsequent
  ones come from the heap. **Note**: `hpsolb` ties are broken by
  insertion order (see `hpsolb.md`); for `--strict` ports the heap
  must match.
- Inside the breakpoint update, `f1` and `f2` are updated by
  multi-term expressions; the F77 source orders the additions
  specifically. Ports must follow.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/cauchy_case_1.json` | `sbgnrm = 0`: early return, `xcp = x` |
| 2 | `data/cauchy_case_2.json` | All variables fixed (`iwhere = 3`): `nbreak = 0`, second early return |
| 3 | `data/cauchy_case_3.json` | All unbounded: `nbreak = 0`, `nfree < n+1`; `xcp = x - g/theta` |
| 4 | `data/cauchy_case_4.json` | One breakpoint, `col = 0`: simple bounded case |
| 5 | `data/cauchy_case_5.json` | With L-BFGS history (`col = 1`) |
| 6 | `data/cauchy_case_6.json` | Variable at upper bound with positive `neggi`: stays at upper |
| 7 | `data/cauchy_case_7.json` | Mixed bounded/unbounded variables |
| 8 | `data/cauchy_case_8.json` | All bounded with one fixed variable |
| 9 | `data/cauchy_case_9.json` | `iprint = 100`: diagnostic-print path |

## Reference implementation

`reference_impl/core/cauchy.py` (depends on `core.bmv`, `core.hpsolb`).

## Cross-references

- **Paper**: Byrd/Lu/Nocedal/Zhu 1995 sec.4 (the original GCP derivation),
  `code.pdf` sec.3.
- **Related subroutines**: called by `mainlb` per iteration. Calls
  `bmv` (twice per breakpoint when `col > 0`) and `hpsolb` (heap-sort
  the breakpoints).
- **F77 source**: `src/cauchy.f`.
- **Unit test**: `tests/unit/test_cauchy.f90` (9 case_* blocks).
