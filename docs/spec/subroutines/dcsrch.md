# dcsrch

## Purpose

Drive the More-Thuente line search to find a step `stp` along a
descent direction that satisfies both the **sufficient decrease**
condition

```
phi(stp) <= phi(0) + ftol * stp * (dphi/dstp)(0)
```

and the **curvature** condition

```
|dphi/dstp| <= gtol * |(dphi/dstp)(0)|
```

(strong Wolfe conditions). Internally, the algorithm maintains an
interval `[stx, sty]` that contains a minimizer, and uses `dcstep`
to compute each next trial step.

The F77 implementation uses **reverse-communication**: the caller
supplies `phi(stp_trial)` and `(dphi/dstp)(stp_trial)` between calls. Ports
targeting modern languages can wrap this in a callback style:

```
def dcsrch(phi, ftol, gtol, xtol, stpmin, stpmax, stp_initial) -> (stp_final, status, n_evals):
    state = init(...)
    while True:
        stp = compute_next_step(state)
        if state.task != 'FG':
            return stp, state.task, state.n_evals
        f, g = phi(stp)
        update(state, stp, f, g)
```

The two styles are equivalent; this spec describes both.

## Mathematical contract

### Reverse-communication interface (F77 style)

| Name | Type | Description |
|------|------|-------------|
| `f` (in/out) | real | On `START`: `phi(0)`. On subsequent: `phi(stp_trial)`. On exit: `phi(stp_final)`. |
| `g` (in/out) | real | On `START`: `(dphi/dstp)(0)` (must be `< 0`). On subsequent: `(dphi/dstp)(stp_trial)`. On exit: `(dphi/dstp)(stp_final)`. |
| `stp` (in/out) | real | On `START`: initial guess. On subsequent: trial step from previous call. On exit: next trial / final step. |
| `ftol`, `gtol`, `xtol` | real, >= 0 | Tolerances. |
| `stpmin`, `stpmax` | real, `0 <= stpmin <= stpmax` | Step bounds. |
| `task` (in/out) | string | Reverse-comm state code. |
| `isave[2]`, `dsave[13]` | (state arrays) | Carry inner state across calls. |

### Callback-style interface (recommended for ports)

```
result dcsrch(
    phi: function(stp) -> (f, g),       # phi(stp), dphi/dstp
    f0: real,                            # phi(0) = caller-supplied
    g0: real,                            # (dphi/dstp)(0), must be < 0
    stp_initial: real,                   # first trial step
    ftol, gtol, xtol: real >= 0,
    stpmin, stpmax: real,
)
```

returning:

```
result {
    stp_final: real,
    f_final: real,
    g_final: real,
    status: enum('CONV', 'WARN_*', 'ERROR_*'),
    n_evals: int,
}
```

### Task / status codes

| F77 `task` prefix | callback `status` | Meaning |
|-------------------|-------------------|---------|
| `START` | (input only) | First call. |
| `FG` | (internal) | Reverse-comm: caller must compute `f` and `g`. |
| `CONV` | `CONV` | Both Wolfe conditions met. |
| `WARN: ROUNDING ERRORS PREVENT PROGRESS` | `WARN_ROUNDING` | Bracket too tight to advance. |
| `WARN: XTOL TEST SATISFIED` | `WARN_XTOL` | `(stmax - stmin) <= xtol * stmax`. |
| `WARN: STP = STPMAX` | `WARN_STPMAX` | Step at upper bound, function still decreasing. |
| `WARN: STP = STPMIN` | `WARN_STPMIN` | Step at lower bound, no improvement. |
| `ERROR: STP .LT. STPMIN` | `ERROR_STP_LT_STPMIN` | Initial `stp` below `stpmin`. |
| `ERROR: STP .GT. STPMAX` | `ERROR_STP_GT_STPMAX` | Initial `stp` above `stpmax`. |
| `ERROR: INITIAL G .GE. ZERO` | `ERROR_G_NONNEG` | Initial gradient not in descent direction. |
| `ERROR: FTOL .LT. ZERO` | `ERROR_FTOL_NEG` | Negative `ftol`. |
| `ERROR: GTOL .LT. ZERO` | `ERROR_GTOL_NEG` | Negative `gtol`. |
| `ERROR: XTOL .LT. ZERO` | `ERROR_XTOL_NEG` | Negative `xtol`. |
| `ERROR: STPMIN .LT. ZERO` | `ERROR_STPMIN_NEG` | Negative `stpmin`. |
| `ERROR: STPMAX .LT. STPMIN` | `ERROR_STPMAX_LT_STPMIN` | Inverted step bounds. |

### Preconditions

- On `'START'`: `f = phi(0)`, `g = (dphi/dstp)(0) < 0`, `stp` in
  `[stpmin, stpmax]`.
- The phi function (or the caller-supplied trial f, g) must be the
  derivative-and-value of a 1D function along the search direction.

### Postconditions

- `stp` is the trial step (`task = 'FG'`) or final step (`task = 'CONV'/'WARN'`).
- `f` and `g` reflect `phi(stp)` and `dphi/dstp`.

## Algorithm

The algorithm has two stages:

- **Stage 1**: search uses a "modified function" `psi(stp) = phi(stp) - phi(0) - ftol * stp * (dphi/dstp)(0)`
  -- a curve translated so the sufficient-decrease line is the
  zero-axis. Once `psi(stp) <= 0` and `dphi/dstp >= 0`, switch to
  stage 2.
- **Stage 2**: search uses the original `phi` directly.

Per call, `dcsrch`:

1. Restore local state from `isave`/`dsave` (or initialize on
   `'START'`).
2. Update `stage` if the stage-2 entry condition is met.
3. Test for warnings (bracket too small, `stp` at boundary).
4. Test for convergence (Wolfe conditions).
5. If terminated, save state and return.
6. Otherwise, call `dcstep` (with `phi` or `psi` depending on stage)
   to compute the next trial step.
7. Apply bisection if the bracket isn't shrinking fast enough.
8. Update step bounds `[stmin, stmax]`.
9. Clamp `stp` to `[stpmin, stpmax]`.
10. Save state, return with `task = 'FG'`.

### Magic constants

| Constant | F77 | Value | Meaning |
|----------|-----|-------|---------|
| `p5` | `p5` | `0.5d0` | Bisection step factor. |
| `p66` | `p66` | `0.66d0` | Width-shrinkage threshold. |
| `xtrapl` | `xtrapl` | `1.1d0` | Lower-extrapolation factor when no bracket. |
| `xtrapu` | `xtrapu` | `4.0d0` | Upper-extrapolation factor when no bracket. |

### Initial-stage tolerances

`ftol = 1.0e-3` (NOT More-Thuente's `1e-4`; see `04_numerics.md`).
`gtol = 0.9` (curvature). `xtol = 0.1` (bracket width).

### Numerical safeguards

- The "modified function" `psi` avoids numerical issues when `phi`
  decreases very slowly near the minimum.
- Bisection step (when `|sty - stx| >= 0.66 * width1`) ensures the
  bracket halves regularly.
- Forced `stp in [stpmin, stpmax]` after each computation.
- If the bracket interval gets stuck (not advancing past `stmin` or
  `stmax`), `dcsrch` falls back to `stx` (the best step seen).

## Test vectors

Test vectors split into **single-call** (input validation) and
**multi-call trajectory** cases. The latter exercise the full
state-machine; ports validate by replaying the trajectory.

| Case | File | Type | Branch exercised |
|------|------|------|------------------|
| 1 | `data/dcsrch_case_1.json` | single | `stp < stpmin` |
| 2 | `data/dcsrch_case_2.json` | single | `stp > stpmax` |
| 3 | `data/dcsrch_case_3.json` | single | `g >= 0` |
| 4 | `data/dcsrch_case_4.json` | single | `ftol < 0` |
| 5 | `data/dcsrch_case_5.json` | single | `gtol < 0` |
| 6 | `data/dcsrch_case_6.json` | single | `xtol < 0` |
| 7 | `data/dcsrch_case_7.json` | single | `stpmin < 0` |
| 8 | `data/dcsrch_case_8.json` | single | `stpmax < stpmin` |

Multi-call trajectory cases (`*_traj_*.json`) are deferred to the
conformance runner (Phase C) -- they require a phi-function dispatch
table plus a step-by-step replay protocol. The four candidate
trajectories from the F77 unit test:

- `phi(t) = (t-1)^2 - 1` (quadratic, converges to `stp = 1`).
- `phi(t) = -t` (monotone descent -> `WARN_STPMAX`).
- `phi(t) = -t / (t^2 + 2)` (More-Thuente phi1, brackets across overshoot).
- `phi(t) = -t + 0.5*t^2` with `stp_initial = 2` (modified-function path).

## Reference implementation

`reference_impl/core/dcsrch.py` (mirrors the F77 reverse-comm
protocol for ease of validation; ports may rewrap as callbacks).

## Cross-references

- **Paper**: More & Thuente 1994 (the algorithm).
- **Related subroutines**: called by `lnsrlb` (which provides
  the reverse-comm bridge to `mainlb`'s outer `f`/`g` evaluation).
  Calls `dcstep` to compute each trial step.
- **F77 source**: `src/dcsrch.f`.
- **Unit test**: `tests/unit/test_dcsrch.f90` (12 case_* blocks; 8
  validation, 4 trajectory).
