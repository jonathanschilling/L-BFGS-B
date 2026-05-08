# freev

## Purpose

Identify the free / active variable partition at the generalized
Cauchy point (GCP) and detect which variables changed status since
the previous iteration. The result feeds the subspace minimization
(`subsm`) and signals whether the limited-memory bookkeeping
(`wn` matrix) needs to be rebuilt.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Number of variables. |
| `nfree` (in) | integer, `0 <= nfree <= n` | Previous-iteration count of free variables. Only consulted when `iter > 0` and `cnstnd`; otherwise unused. |
| `index` (in) | integer vector, length `n` | Previous-iteration partition: `index[1..nfree]` are previously-free indices, `index[nfree+1..n]` are previously-active. Only consulted when `iter > 0` and `cnstnd`. |
| `iwhere` | integer vector, length `n` | Current bound status per variable (set by `cauchy`): `<= 0` = free at GCP, `> 0` = at a bound. |
| `updatd` | boolean | `true` if the L-BFGS matrix was updated in the previous iteration. |
| `cnstnd` | boolean | `true` if at least one variable has a bound. |
| `iter` | integer >= 0 | Outer iteration counter. |
| `iprint` | integer | Diagnostic verbosity. |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `nfree` (out) | integer | New count of free variables (number of `i` with `iwhere[i] <= 0`). |
| `index` (out) | integer vector, length `n` | New partition: `index[1..nfree]` lists currently-free indices in ascending order, `index[nfree+1..n]` lists currently-active indices in descending order (see Algorithm). |
| `nenter` | integer >= 0 | Number of variables that *entered* the free set this iteration (were active, now free). |
| `ileave` | integer, `1 <= ileave <= n+1` | Sentinel: `indx2[ileave..n]` lists the variables that *left* the free set this iteration. `ileave = n+1` means none left. |
| `indx2` | integer vector, length `n` | Change record: `indx2[1..nenter]` are entering variables (in ascending order of their original index `i`); `indx2[ileave..n]` are leaving variables (in descending order). |
| `wrk` | boolean | `true` if the active set changed *or* `updatd` was `true` -- signals that the `wn` matrix must be rebuilt before subspace minimization. |

### Preconditions

- `iwhere` is consistent with `cauchy`'s output: each `iwhere[i] in {-1, 0, 1, 2, 3}`.
  Values `<= 0` (`-1`, `0`) indicate "free"; values `> 0` (`1`, `2`, `3`) indicate "at a bound".
- When `iter > 0` and `cnstnd`, the input `(nfree, index)` reflects the
  previous iteration's partition.
- `iter == 0` ==> no previous partition; the function uses only `iwhere`
  and skips change detection.

### Postconditions

- `nfree = |{i : iwhere[i] <= 0}|`.
- `index` partitions `{1, ..., n}` into free and active.
- `nenter + (n + 1 - ileave) <= n` (no variable both entered and left).
- `wrk = (ileave < n+1) or (nenter > 0) or updatd`.

## Algorithm

### Phase 1: detect entering / leaving variables (only if `iter > 0` and `cnstnd`)

```
nenter = 0
ileave = n + 1

if iter > 0 and cnstnd:
    # Variables previously free, now at a bound -> "leave"
    for i = 1 to nfree (input):
        k = index[i]                                   # previously-free index
        if iwhere[k] > 0:
            ileave = ileave - 1
            indx2[ileave] = k

    # Variables previously active, now free -> "enter"
    for i = (nfree + 1) to n (input):
        k = index[i]                                   # previously-active index
        if iwhere[k] <= 0:
            nenter = nenter + 1
            indx2[nenter] = k
```

When `iter == 0` or `not cnstnd`, the change detection is skipped; the
sentinels stay at their initial values (`nenter = 0`, `ileave = n+1`).

### Phase 2: set the work flag

```
wrk = (ileave < n+1) or (nenter > 0) or updatd
```

### Phase 3: rebuild `index` partition from the current `iwhere`

This phase always runs, regardless of `iter` and `cnstnd`. **It
overwrites the input `index` and `nfree`.**

```
nfree = 0
iact  = n + 1
for i = 1 to n:                                        # ascending order
    if iwhere[i] <= 0:
        nfree = nfree + 1
        index[nfree] = i
    else:
        iact = iact - 1
        index[iact] = i
```

After this loop:
- `index[1..nfree]` = currently-free indices, ordered ascending by `i`.
- `index[nfree+1..n]` = currently-active indices, ordered descending by `i`
  (because the active part is built from the back).

### Pseudocode (combined)

```
input:  n, nfree (in), index (in/out), iwhere, updatd, cnstnd, iter, iprint
output: nfree (out), index (out), nenter, ileave, indx2, wrk

# Phase 1: change detection
nenter = 0
ileave = n + 1
if iter > 0 and cnstnd:
    for i = 1 to nfree:                        # over previously-free
        k = index[i]
        if iwhere[k] > 0:
            ileave -= 1
            indx2[ileave] = k
    for i = (nfree + 1) to n:                  # over previously-active
        k = index[i]
        if iwhere[k] <= 0:
            nenter += 1
            indx2[nenter] = k

# Phase 2: rebuild flag
wrk = (ileave < n+1) or (nenter > 0) or updatd

# Phase 3: rebuild partition from iwhere
nfree = 0
iact  = n + 1
for i = 1 to n:
    if iwhere[i] <= 0:
        nfree += 1
        index[nfree] = i
    else:
        iact -= 1
        index[iact] = i

# Optional diagnostic prints (port-defined)
```

### Magic constants

None.

### Numerical safeguards

None (entirely integer arithmetic and array indexing).

### Order-of-operations dependencies

- The change-detection loops iterate `i = 1..nfree` (then
  `nfree+1..n`) in ascending order. Order matters because `indx2`
  records *positional* slots: with `iprint >= 100`, ports that match
  the F77 diagnostic output must visit the variables in the same
  order.
- The partition loop iterates `i = 1..n` in ascending order;
  changing the order would change which active variables land in
  which slot of `index[nfree+1..n]`. `--strict` conformance requires
  ascending.

## Special cases

- **`iter == 0` (first call)**: no previous partition exists. Phase 1
  is skipped; sentinels remain `nenter = 0`, `ileave = n+1`.
- **Unconstrained problem (`cnstnd = false`)**: every `iwhere[i] = -1`,
  so phase 3 produces `nfree = n` and `index = [1, 2, ..., n]`. Phase 1
  is also skipped.
- **Fixed variables (`iwhere[i] = 3`)**: treated as "at bound" (`> 0`),
  so they fall in the active partition. They never re-enter the free
  set across iterations.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/freev_case_1.json` | `iter = 0`, mixed `iwhere`: phase 1 skipped, partition built |
| 2 | `data/freev_case_2.json` | `iter > 0` but `cnstnd = false`: phase 1 skipped |
| 3 | `data/freev_case_3.json` | `iter > 0` cnstnd, no changes, `updatd = true`: `wrk` via OR clause |
| 4 | `data/freev_case_4.json` | `iter > 0`: one variable leaves the free set |
| 5 | `data/freev_case_5.json` | `iter > 0`: one variable enters the free set |
| 6 | `data/freev_case_6.json` | `iter > 0`: both events; `iprint = 100` exercises diagnostic prints |

## Reference implementation

`reference_impl/core/freev.py`.

## Cross-references

- **Paper**: `code.pdf` sec.3 (active-set bookkeeping). The
  enter/leave detection is needed to know when the reduced Hessian
  must be refactored.
- **Related subroutines**: called by `mainlb` after `cauchy`. The
  outputs feed `formk` (uses `wrk`) and `subsm` (uses `index`,
  `nfree`).
- **F77 source**: `src/freev.f`.
- **Unit test**: `tests/unit/test_freev.f90` (6 case_* blocks).
