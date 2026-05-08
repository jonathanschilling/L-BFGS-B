# Algorithm interface

The canonical L-BFGS-B interface is a callback-based function. The
caller passes the initial point, the bounds, the algorithm parameters,
and two closures that evaluate the objective function and its gradient
on demand. The optimizer drives the callbacks internally and returns
the final iterate.

This is the natural interface for languages with first-class function
values. Ports are free to expose their own user-facing API on top of
this -- generators, async iterators, builder patterns, observer
callbacks, or even the F77 reverse-comm protocol described in
`08_legacy_reverse_comm.md` -- but the spec describes the algorithm in
callback terms.

## Signature

In language-neutral pseudocode:

```
result minimize(
    n: positive integer,                    # dimension
    x0: real vector of length n,            # initial point
    l: real vector of length n,             # lower bounds
    u: real vector of length n,             # upper bounds
    nbd: integer vector of length n,        # bound type per variable
    m: positive integer,                    # max L-BFGS pairs (typical 5-17)
    f_eval: function(x: real vector) -> real,
    g_eval: function(x: real vector) -> real vector of length n,
    factr: nonneg real = 1e7,               # function-decrease tolerance
    pgtol: nonneg real = 1e-5,              # projected-gradient tolerance
    max_iter: positive integer = none,      # iteration cap (none = unlimited)
    max_fg: positive integer = none,        # f/g evaluation cap
    iprint: integer = -1,                   # diagnostic verbosity
)
```

Returns:

```
result {
    x: real vector of length n,        # final iterate
    f: real,                           # f(x) at final iterate
    g: real vector of length n,        # gradient at final iterate
    info: integer,                     # termination code (see below)
    n_iter: integer,                   # iterations performed
    n_fg: integer,                     # f/g callback invocations
    final_projg: real,                 # final ||g_proj||_inf
    message: string,                   # human-readable termination reason
}
```

## Inputs

### `n` -- dimension

Positive integer. Must be `>= 1`.

### `x0` -- initial point

Real vector of length `n`. May be infeasible: the optimizer will
project onto the box `l <= x <= u` before iterating; the caller is
notified via `info` if projection occurred.

### `l`, `u` -- bounds

Real vectors of length `n`. Entries are **only consulted** where
`nbd[i]` indicates the corresponding bound is active (see below).
Entries with `nbd[i] = 0` may be any value (commonly `0`).

### `nbd` -- bound types

Integer vector of length `n`. Per variable:

| `nbd[i]` | Meaning |
|----------|---------|
| `0` | Variable is unbounded. `l[i]`, `u[i]` ignored. |
| `1` | Lower bound active: `x[i] >= l[i]`. |
| `2` | Both bounds active: `l[i] <= x[i] <= u[i]`. Allowed `l[i] = u[i]` (variable fixed). |
| `3` | Upper bound active: `x[i] <= u[i]`. |

Validation: `nbd[i] in {0,1,2,3}` is required; if `nbd[i] = 2` then
`l[i] <= u[i]` is required. Violations cause an immediate return with
`info` set to an error code (see `errclb.md`).

### `m` -- memory parameter

Positive integer. Number of L-BFGS pairs `(s_i, y_i)` retained.
Typical values: `5` (default), up to `17`. Larger `m` improves the
Hessian approximation at the cost of `O(mn)` extra storage and
`O(m^2n)` per-iteration work.

### `f_eval`, `g_eval` -- callbacks

Pure functions called by the optimizer. Contract:

- **Pure**: a port may rely on `f_eval(x) = f_eval(x)` for the same `x`
  (the optimizer may call them with the same `x` more than once,
  although it tries to avoid this).
- **Inputs**: `f_eval` and `g_eval` receive the current trial iterate
  `x`. The optimizer guarantees `l <= x <= u` (the trial point is always
  feasible after the initial projection).
- **Outputs**: `f_eval` returns a finite real (`NaN` / `Inf` cause
  abnormal termination); `g_eval` returns a real vector of length `n`.
- **Cost**: the optimizer minimizes calls but cannot avoid them -- for
  expensive objectives, the user controls cost via `max_fg` and
  `factr`/`pgtol` tolerances.

A port may merge `f_eval` and `g_eval` into a single callback returning
both `(f, g)`; this is more efficient when the gradient is computed as
a byproduct of the function evaluation. Both styles are valid; this
spec describes them as separate callbacks for clarity.

### `factr` -- relative function-decrease tolerance

Convergence triggers when

```
(f_old - f_new) / max(|f_old|, |f_new|, 1) <= factr * eps
```

where `eps` is the machine epsilon (`~= 2.22e-16`). Recommended values:

- `factr = 1e12`: low accuracy (~6 digits).
- `factr = 1e7`: medium accuracy (~10 digits, default).
- `factr = 1e1`: machine accuracy (~15 digits).

Setting `factr = 0` disables the test (only `pgtol` and limits remain).

### `pgtol` -- projected-gradient tolerance

Convergence triggers when `||g_proj||_Inf <= pgtol`, where `g_proj` is the
projected gradient (see `subroutines/projgr.md` for the exact
definition).

Setting `pgtol = 0` disables the test.

### `max_iter`, `max_fg` -- limits

Optional caps. If reached, the optimizer terminates with
`info = INFO_LIMIT` (see below). The F77 default in the drivers is
`max_iter = 99` and unlimited `max_fg`; ports may choose their own
defaults.

### `iprint` -- diagnostic verbosity

| Value | Behavior |
|-------|----------|
| `< 0` | Silent (default for libraries). |
| `0` | Print only at start and end. |
| `1` to `98` | Print every `iprint`-th iteration; final summary. |
| `99` | Print details at every iteration. |
| `100` | Plus print Cauchy point variables and active-set changes. |
| `> 100` | Plus print working vectors at each iteration. |

Diagnostic format is **port-defined**. The F77 implementation writes
to unit 6 and an `iterate.dat` file; ports may use any logging
mechanism. Conformance tests do not check diagnostic output.

## Termination codes (`info`)

| Code | Symbolic name | Meaning |
|------|---------------|---------|
| `0` | `INFO_CONVERGED_PGTOL` | `||g_proj||_Inf <= pgtol`. |
| `1` | `INFO_CONVERGED_FACTR` | Relative function decrease below `factr * eps`. |
| `2` | `INFO_LIMIT_ITER` | Reached `max_iter` without converging. |
| `3` | `INFO_LIMIT_FG` | Reached `max_fg` without converging. |
| `4` | `INFO_ABNORMAL_LNSRLB` | Line search failed (function may be unbounded below, or gradient is inconsistent with `f`). |
| `5` | `INFO_USER_STOP` | User requested early termination. |
| `<= -1` | `INFO_INPUT_ERROR_*` | Input validation failed. See `errclb.md` for the specific code per validation. |

Ports may use enums, exceptions, or sentinel values to expose this
information; the codes above are the spec convention.

The `message` string in the result is a human-readable description
suitable for logging. The F77 implementation uses the `task` string;
ports may produce any text but should match the meaning of the
corresponding `info` code.

## Behavioral guarantees

- **Initial projection**: if any `x0[i]` is outside its bounds, the
  optimizer projects to the nearest feasible point before any callback
  is invoked. The callbacks always receive feasible points.
- **Monotonic decrease (within tolerance)**: each iteration either
  reduces `f` or terminates. The line search (More-Thuente Wolfe)
  enforces sufficient decrease and curvature conditions; rejected
  steps trigger refresh of the L-BFGS history rather than accepting an
  ascent step.
- **Bounded callback count**: the algorithm performs at most
  `max_iter` outer iterations and at most `max_fg` total callback
  evaluations.
- **Determinism**: given identical inputs, identical callbacks, and a
  reference BLAS/LAPACK pin (see `04_numerics.md`), two runs produce
  bit-identical outputs. The algorithm has no randomness.

## Behavior on invalid input

The optimizer validates inputs before the first callback. On
validation failure it returns immediately with `info <= -1` and
**no callbacks have been invoked**. See `errclb.md` for the full list
of error cases and codes.

NaN / Inf in the user-supplied `x0`, `l`, or `u` is allowed only where
the bound is inactive (e.g., `u[i]` may be `+Inf` if `nbd[i] in {0, 1}`).
NaN / Inf returned from `f_eval` or `g_eval` causes abnormal
termination with `info = INFO_ABNORMAL_LNSRLB`.

## Optional: stopping callback

Some ports may want to allow the user to interrupt the optimization
between iterations (e.g., timeout, external signal). The F77 `setulb`
exposes this via `task = &#39;STOP&#39;` mid-iteration.

For a callback-based port, the recommended idiom is an optional
`should_stop: function(state) -> boolean` callback invoked at each
iteration boundary; if it returns `true`, the optimizer terminates
with `info = INFO_USER_STOP`. The state passed to `should_stop` may
include `iter`, `n_fg`, current `f`, current `||g_proj||_Inf`, and elapsed
wall time -- port's choice.

This is **not part of the algorithm spec** (the algorithm itself does
nothing different); it is an interface convenience that ports may
or may not provide.

## Mapping to F77 `setulb`

For readers comparing to the F77 source:

| Spec parameter | F77 equivalent |
|----------------|----------------|
| `n`, `m`, `x0`, `l`, `u`, `nbd` | Same names. |
| `factr`, `pgtol`, `iprint` | Same names. |
| `f_eval`, `g_eval` | Caller computes `f`, `g` between `setulb` calls when `task = &#39;FG&#39;`. |
| `max_iter`, `max_fg` | F77 has no caps in `setulb` itself; the drivers loop and break on `task != &#39;FG&#39; && task != &#39;NEW_X&#39;`. |
| Result `x`, `f`, `g`, `info` | F77 leaves these in-place in the user's arrays. |
| Result `n_iter`, `n_fg` | F77 stores in `isave(30)`, `isave(34)`. |
| Result `message` | F77 returns in `task`. |

The full F77 protocol is in `08_legacy_reverse_comm.md`.
