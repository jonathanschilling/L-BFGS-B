# dcsrch

## Purpose

Drive the More-Thuente line search to find a step `stp` along a
descent direction that satisfies both the **sufficient decrease**
condition

```
phi(stp) <= phi(0) + ftol * stp * phi'(0)
```

and the **curvature** condition

```
|phi'(stp)| <= gtol * |phi'(0)|
```

(strong Wolfe conditions). Internally, the algorithm maintains an
interval `[stx, sty]` that contains a minimizer, and uses `dcstep`
to compute each next trial step.

The F77 implementation uses **reverse-communication**: the caller
supplies `phi(stp_trial)` and `phi&#39;(stp_trial)` between calls. Ports
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
| `f` (in/out) | real | On `&#39;START&#39;`: `phi(0)`. On subsequent: `phi(stp_trial)`. On exit: `phi(stp_final)`. |
| `g` (in/out) | real | On `&#39;START&#39;`: `phi&#39;(0)` (must be `< 0`). On subsequent: `phi&#39;(stp_trial)`. On exit: `phi&#39;(stp_final)`. |
| `stp` (in/out) | real | On `&#39;START&#39;`: initial guess. On subsequent: trial step from previous call. On exit: next trial / final step. |
| `ftol`, `gtol`, `xtol` | real, >= 0 | Tolerances. |
| `stpmin`, `stpmax` | real, `0 <= stpmin <= stpmax` | Step bounds. |
| `task` (in/out) | string | Reverse-comm state code. |
| `isave[2]`, `dsave[13]` | (state arrays) | Carry inner state across calls. |

### Callback-style interface (recommended for ports)

```
result dcsrch(
    phi: function(stp) -> (f, g),       # phi(stp), phi'(stp)
    f0: real,                            # phi(0) = caller-supplied
    g0: real,                            # phi'(0), must be < 0
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
| `&#39;START&#39;` | (input only) | First call. |
| `&#39;FG&#39;` | (internal) | Reverse-comm: caller must compute `f` and `g`. |
| `&#39;CONV&#39;` | `CONV` | Both Wolfe conditions met. |
| `&#39;WARN: ROUNDING ERRORS PREVENT PROGRESS&#39;` | `WARN_ROUNDING` | Bracket too tight to advance. |
| `&#39;WARN: XTOL TEST SATISFIED&#39;` | `WARN_XTOL` | `(stmax - stmin) <= xtol * stmax`. |
| `&#39;WARN: STP = STPMAX&#39;` | `WARN_STPMAX` | Step at upper bound, function still decreasing. |
| `&#39;WARN: STP = STPMIN&#39;` | `WARN_STPMIN` | Step at lower bound, no improvement. |
| `&#39;ERROR: STP .LT. STPMIN&#39;` | `ERROR_STP_LT_STPMIN` | Initial `stp` below `stpmin`. |
| `&#39;ERROR: STP .GT. STPMAX&#39;` | `ERROR_STP_GT_STPMAX` | Initial `stp` above `stpmax`. |
| `&#39;ERROR: INITIAL G .GE. ZERO&#39;` | `ERROR_G_NONNEG` | Initial gradient not in descent direction. |
| `&#39;ERROR: FTOL .LT. ZERO&#39;` | `ERROR_FTOL_NEG` | Negative `ftol`. |
| `&#39;ERROR: GTOL .LT. ZERO&#39;` | `ERROR_GTOL_NEG` | Negative `gtol`. |
| `&#39;ERROR: XTOL .LT. ZERO&#39;` | `ERROR_XTOL_NEG` | Negative `xtol`. |
| `&#39;ERROR: STPMIN .LT. ZERO&#39;` | `ERROR_STPMIN_NEG` | Negative `stpmin`. |
| `&#39;ERROR: STPMAX .LT. STPMIN&#39;` | `ERROR_STPMAX_LT_STPMIN` | Inverted step bounds. |

### Preconditions

- On `&#39;START&#39;`: `f = phi(0)`, `g = phi&#39;(0) < 0`, `stp` in
  `[stpmin, stpmax]`.
- The phi function (or the caller-supplied trial f, g) must be the
  derivative-and-value of a 1D function along the search direction.

### Postconditions

- `stp` is the trial step (`task = &#39;FG&#39;`) or final step (`task = &#39;CONV&#39;/&#39;WARN&#39;`).
- `f` and `g` reflect `phi(stp)` and `phi&#39;(stp)`.

## Algorithm

The algorithm has two stages:

- **Stage 1**: search uses a "modified function" `psi(stp) = phi(stp) - phi(0) - ftol * stp * phi&#39;(0)`
  -- a curve translated so the sufficient-decrease line is the
  zero-axis. Once `psi(stp) <= 0` and `phi&#39;(stp) >= 0`, switch to
  stage 2.
- **Stage 2**: search uses the original `phi` directly.

Per call, `dcsrch`:

1. Restore local state from `isave`/`dsave` (or initialize on
   `&#39;START&#39;`).
2. Update `stage` if the stage-2 entry condition is met.
3. Test for warnings (bracket too small, `stp` at boundary).
4. Test for convergence (Wolfe conditions).
5. If terminated, save state and return.
6. Otherwise, call `dcstep` (with `phi` or `psi` depending on stage)
   to compute the next trial step.
7. Apply bisection if the bracket isn't shrinking fast enough.
8. Update step bounds `[stmin, stmax]`.
9. Clamp `stp` to `[stpmin, stpmax]`.
10. Save state, return with `task = &#39;FG&#39;`.

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
