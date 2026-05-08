# lnsrlb

## Purpose

Drive the More-Thuente line search (`dcsrch`) along a search direction
`d`, with safeguarding so that all trial points stay inside the
feasible box. Computes the actual `x = t + stp * d` at each trial
(or `x = z` exactly when `stp = 1`), and clamps the maximum step
`stpmx` to keep `x` feasible.

This is the bridge between `mainlb`'s outer iteration and `dcsrch`'s
1D search. In F77 it is reverse-comm; the wrapper accepts
`task = 'FG_LN'` (continuation) vs `'START'` (new search).

In a callback-based port, `lnsrlb` is just an internal function in
the `mainlb` loop that calls the user's `f_eval` / `g_eval` between
its own logic and `dcsrch`.

## Mathematical contract

### F77 reverse-comm signature (in/out arguments are mostly all stateful)

| Name | Kind | Description |
|------|------|-------------|
| `n` | in | Dimension. |
| `l, u, nbd` | in | Bounds. |
| `x, f, g` | in/out | Iterate, value, gradient. |
| `fold` | in/out | Saved `f` at start of search (for restore). |
| `gd, gdold` | in/out | `g&#39;d` at current trial / at `stp=0`. |
| `d` | in | Search direction. |
| `r, t` | in/out | Saved `g` and `x` at start of search. |
| `z` | in | The Cauchy/subsm output point (used when `stp = 1`). |
| `stp, dnorm, dtd, xstep` | in/out | Step length, `||d||`, `d&#39;d`, `stp * ||d||`. |
| `stpmx` | in/out | Maximum allowed step (computed from box geometry). |
| `iter` | in | Outer iteration counter. |
| `ifun, iback, nfgv` | in/out | f/g eval counters. |
| `info` | out | `0` on success; `-4` if not a descent direction. |
| `task` | in/out | `&#39;START&#39;` (new search), `&#39;FG_LN&#39;` (continuation), `&#39;FG_LNSRCH&#39;` (out: more f/g needed), `&#39;NEW_X&#39;` (out: search done). |
| `boxed, cnstnd` | in | Box-constraint flags. |
| `csave, isave, dsave` | in/out | Inner state for `dcsrch`. |

### Algorithm

```
if task starts with 'FG_LN':
    skip setup (continuation: jump to step 4)

# Step 1: setup (first call only)
dtd = d'd; dnorm = sqrt(dtd)

# Step 2: compute stpmx from box geometry (only if cnstnd)
stpmx = big = 1e10
if cnstnd:
    if iter == 0:
        stpmx = 1
    else:
        for i = 1 to n:
            a1 = d[i]
            if nbd[i] != 0:
                if a1 < 0 and nbd[i] <= 2:        # would push x toward lower bound
                    a2 = l[i] - x[i]
                    if a2 >= 0: stpmx = 0          # already at or past bound
                    elif a1 * stpmx < a2: stpmx = a2 / a1
                elif a1 > 0 and nbd[i] >= 2:      # would push toward upper
                    a2 = u[i] - x[i]
                    if a2 <= 0: stpmx = 0
                    elif a1 * stpmx > a2: stpmx = a2 / a1

# Step 3: pick initial stp
if iter == 0 and not boxed:
    stp = min(1 / dnorm, stpmx)
else:
    stp = 1                                    # default unit step

# Step 3a: save state for restore on failure
t = x; r = g; fold = f; ifun = 0; iback = 0; csave = 'START'

# Step 4: enter dcsrch
gd = g'd
if ifun == 0:
    gdold = gd
    if gd >= 0:                                # not a descent direction
        info = -4
        return

call dcsrch(f, gd, stp, ftol=1e-3, gtol=0.9, xtol=0.1, 0, stpmx, csave, isave, dsave)

xstep = stp * dnorm

# Step 5: dispatch on dcsrch's reply
if csave is not 'CONV' or 'WARN':
    task = 'FG_LNSRCH'                         # need another f/g
    ifun += 1; nfgv += 1; iback = ifun - 1
    if stp == 1:
        x = z                                  # exact unit step: use precomputed Cauchy/subsm point
    else:
        for i = 1 to n: x[i] = stp * d[i] + t[i]
    return
else:
    task = 'NEW_X'                             # line search done
    return
```

### Magic constants

| Constant | F77 | Value | Meaning |
|----------|-----|-------|---------|
| `big` | `big` | `1.0d10` | Huge fallback for `stpmx` when no bound binds. |
| `ftol` | `ftol` | `1.0d-3` | Sufficient-decrease tolerance for `dcsrch`. **DEVIATES from More-Thuente's `1e-4`**; see `04_numerics.md`. |
| `gtol` | `gtol` | `0.9d0` | Curvature tolerance. |
| `xtol` | `xtol` | `0.1d0` | Bracket-width tolerance. |

### Numerical safeguards

- `stpmx = 0` when a variable is already at the binding bound (`a2 = 0`):
  no step in that direction is allowed. Caller (`mainlb`) detects this
  and either refreshes the L-BFGS history or terminates.
- The exact `stp == 1` test (line 180) is intentional: when the
  subspace minimizer happens to land at the unit step, `lnsrlb` uses
  the precomputed `z` directly (avoids redundant arithmetic). The
  numerical result is the same as computing `t + stp*d`; the branch is
  for efficiency.
- The descent-direction check (`gd >= 0` ==> `info = -4`) catches cases
  where the subspace minimizer produces a non-descent `d` (rare but
  possible with degenerate Hessian approximations).

### Order-of-operations dependencies

- The `stpmx`-computation loop iterates `i = 1..n` in ascending order.
  `stpmx` is updated multiplicatively and via min/max; the result is
  order-independent in the limit, but **bit-exact reproducibility
  requires ascending order** (the `stpmx` variable's intermediate
  values matter for the `a1 * stpmx < a2` comparison at each step).
- The `x[i] = stp * d[i] + t[i]` loop is a `daxpy`-style operation,
  order-independent.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/lnsrlb_case_1.json` | `cnstnd = false`, `iter = 0`, not boxed: `stp = min(1/dnorm, big)` |
| 2 | `data/lnsrlb_case_2.json` | `cnstnd = true`, `iter = 0`: `stpmx = 1` |
| 3 | `data/lnsrlb_case_3.json` | `cnstnd = true`, `iter > 0`, `d < 0` lower-bound: `stpmx = (l-x)/d` |
| 4 | `data/lnsrlb_case_4.json` | `cnstnd = true`, `iter > 0`, `d > 0` upper-bound: `stpmx = (u-x)/d` |
| 5 | `data/lnsrlb_case_5.json` | `gd >= 0` ascent direction: `info = -4` |
| 6 | `data/lnsrlb_case_6.json` | Continuation path: `task = &#39;FG_LN&#39;` skips setup |
| 7 | `data/lnsrlb_case_7.json` | Variable already at lower bound: `a2 = 0 -> stpmx = 0` |
| 8 | `data/lnsrlb_case_8.json` | Variable already at upper bound: `a2 = 0 -> stpmx = 0` |

## Reference implementation

`reference_impl/core/lnsrlb.py` (depends on `core.dcsrch`).

## Cross-references

- **Paper**: Wraps More-Thuente 1994; the bound-projection logic is
  in `code.pdf` sec.3.
- **Related subroutines**: called by `mainlb` per iteration. Calls
  `dcsrch` for the 1D search.
- **F77 source**: `src/lnsrlb.f`.
- **Unit test**: `tests/unit/test_lnsrlb.f90` (8 case_* blocks).
