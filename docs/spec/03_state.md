# Logical state model

This document describes the algorithm's state in terms of logical
objects: what L-BFGS-B holds in memory, what changes per iteration,
and how state is initialized and refreshed.

It does **not** describe the F77 `wa` / `iwa` workspace packing -- that
is an implementation choice of the F77 source. Ports may use any
storage layout.

## Overview

The state divides into three groups:

1. **Persistent state** -- survives across iterations. Updated each
   successful iteration.
2. **Transient state** -- built and discarded within a single iteration.
3. **Reverse-comm state** *(legacy only)* -- only meaningful for ports
   that mimic the F77 `setulb` reverse-comm interface. Encodes how to
   resume mid-iteration. See `08_legacy_reverse_comm.md`.

A callback-based port (the canonical interface; see `02_api.md`) needs
groups 1 and 2 only. The "state" between callbacks is the natural
local-variable state of the running optimization function -- no
packing required.

## Persistent state

Held across iteration boundaries. After an iteration completes, this
is what the optimizer needs to start the next one.

### Iterate

| Object | Type | Description |
|--------|------|-------------|
| `x` | real vector, length `n` | Current iterate. Always feasible (`l <= x <= u`). |
| `f` | real | `f(x)` at current iterate. |
| `g` | real vector, length `n` | `gradf(x)` at current iterate. |
| `xp` | real vector, length `n` | Previous iterate (used for refresh decisions). |

### Counters

| Object | Type | Description |
|--------|------|-------------|
| `iter` | integer >= 0 | Iteration counter; `iter = 0` is the initial point. |
| `n_fg` | integer >= 0 | Total callback invocations so far. |

### L-BFGS history

The history is a **logical pair list**, ordered from oldest to newest:
`[(s_1, y_1), (s_2, y_2), ..., (s_col, y_col)]`. After each successful
iteration with a sufficiently positive curvature, a new pair is
appended; if `col = m` is already at the limit, the oldest pair is
dropped first.

| Object | Type | Description |
|--------|------|-------------|
| `col` | integer, `0 <= col <= m` | Current pair count. |
| `S` | real matrix, `n * col` | Columns are `s_i = x_{k_i+1} - x_{k_i}` in insertion order. |
| `Y` | real matrix, `n * col` | Columns are `y_i = g_{k_i+1} - g_{k_i}` in insertion order. |
| `theta` | real, `> 0` | Hessian scaling (`y_k&#39; y_k / s_k&#39; y_k` after each update; initialized to `1`). |

Ports may store `S` and `Y` as column-major matrices, row-major
matrices, lists of vectors, or any other structure. The F77 source uses
column-major `n * m` arrays with `head`/`itail` indices for cyclic
storage when `col = m`; this is one valid choice but not required.

### Cached quantities for the compact representation

Maintained alongside `S`, `Y` to avoid recomputing `O(col^2 * n)` work
per iteration. **All three are derived deterministically from `S`,
`Y`, and `theta`.** A port may recompute on demand instead of caching;
the only consequence is more `O(col^2 n)` work per iteration. If
cached, they must be updated when `(S, Y, theta)` changes.

| Object | Type | Defined as |
|--------|------|------------|
| `sy` | real matrix, `col * col` | `sy[i,j] = s_i&#39; * y_j`. The diagonal is `D = diag(s_i&#39; y_i)`. The strict lower triangle is `L_{i,j} = s_i&#39; y_j` for `i > j`. The strict upper triangle is unused (in the F77 packing). |
| `ss` | real matrix, `col * col` | `ss[i,j] = s_i&#39; * s_j`. Symmetric; the F77 source stores only the upper triangle. |
| `T` | real matrix, `col * col`, upper-triangular | Cholesky factor: `T&#39; T = theta * S&#39;S + L * D^{-1} * L&#39;` (the Cholesky factor of the middle matrix `M_inv`'s lower-right block). Computed by `formt` from `sy`, `ss`, `theta`. |

### Bound state

| Object | Type | Description |
|--------|------|-------------|
| `iwhere` | integer vector, length `n` | Per-variable bound status: `-1` (unbounded), `0` (free, has bounds), `1` (at lower), `2` (at upper), `3` (fixed). See `01_glossary.md`. Initialized by `active`; refined by `cauchy`. |
| `prjctd` | boolean | `true` if the initial `x0` was infeasible and projected. Set once by `active` at startup. |
| `cnstnd` | boolean | `true` if any variable has at least one bound. Set once at startup. |
| `boxed` | boolean | `true` if every variable has both lower and upper bounds. Set once at startup. |

### Refresh tracking

| Object | Type | Description |
|--------|------|-------------|
| `updatd` | boolean | `true` if the L-BFGS history was updated in the *previous* iteration. Used by `cauchy` and `subsm` to decide whether `T` and `K` factorizations need recomputation. |
| `head` | integer | F77 implementation detail: index of oldest column in the cyclic `S`/`Y` buffer. Ports using non-cyclic storage do not need this. |

## Transient state (per iteration)

Built up during one iteration and discarded at the iteration boundary.
A callback-based port can hold these as local variables in the iteration
loop.

### Cauchy step

| Object | Type | Description |
|--------|------|-------------|
| `xc` | real vector, length `n` | Generalized Cauchy point. Computed by `cauchy`. |
| `c` | real vector, length `2 col` | `c = W&#39; (xc - x)`; built by `cauchy`, used by `subsm`. |

### Active set after Cauchy

| Object | Type | Description |
|--------|------|-------------|
| `nfree` | integer | Number of currently free variables. |
| `index` | integer vector, length `n` | Permutation: `index[1..nfree]` are the free-variable indices; `index[nfree+1..n]` are the active (bound) ones. |
| `indx2` | integer vector, length `n` | Indices of variables that *changed* status this iteration: `indx2[1..nenter]` entered the free set; `indx2[ileave..n]` leaved. |
| `nenter` | integer | Number of variables that became free this iteration. |
| `ileave` | integer <= `n+1` | Index in `indx2` where leaving variables start (`n+1` means none left). |

### Subspace step

| Object | Type | Description |
|--------|------|-------------|
| `r` | real vector, length `n` | Reduced gradient at the Cauchy point on the free variables. |
| `d` | real vector, length `n` | Search direction for the line search. |

### Line search

| Object | Type | Description |
|--------|------|-------------|
| `stp` | real, `> 0` | Step length along `d`. |
| `xnew` | real vector, length `n` | Trial point `x + stp * d`. |
| `f_new`, `g_new` | real, real vector | Trial function value and gradient. |
| (line-search internal state) | -- | Bracket interval `[stx, sty]`, plus best/end-point function and slope values. See `dcsrch.md` for the full set. |

### Diagnostics and counters (per iteration)

| Object | Type | Description |
|--------|------|-------------|
| `g_proj_inf` | real | `||g_proj||_Inf` at the start of the iteration. Used for the convergence test. |
| `f_old` | real | `f` at the start of the iteration; used by the relative-decrease convergence test. |
| `nfgv` | integer | Callbacks consumed within this iteration's line search. |
| `nint` | integer | Internal step-count counter (cumulative across iterations). |

## Initialization

Before the first iteration:

1. **Validate inputs** (`errclb`). Failure -> return with `info < 0`.
2. **Project `x0`** to the feasible set if needed (`active`); set
   `iwhere`, `prjctd`, `cnstnd`, `boxed`.
3. **Allocate** `S`, `Y` empty (`col = 0`); set `theta = 1`.
4. **Set counters**: `iter = 0`, `n_fg = 0`, `updatd = false`.
5. **First evaluation**: invoke `f_eval(x0)` and `g_eval(x0)` to get
   `f` and `g`. Increment `n_fg`.
6. **Initial convergence check**: compute `g_proj_inf`; if it is
   `<= pgtol`, return with `info = INFO_CONVERGED_PGTOL` immediately.

The legacy reverse-comm version splits steps 1-4 into the `START` task
and yields after step 5 with `task = &#39;FG&#39;`; ports using callbacks
collapse this into straight code.

## Per-iteration update

After computing the next iterate `x_new` with function `f_new` and
gradient `g_new`:

1. **Check curvature**: compute `s = x_new - x` and `y = g_new - g`. If
   `s&#39; y / |s&#39; s| > eps * ||y||^2` (a refresh threshold; see `mainlb.md`
   for the exact formula), accept the pair.
2. **If accepted**: append `(s, y)` to `S`, `Y`. If `col = m` already,
   drop the oldest column first. Update `theta = (y&#39; y) / (s&#39; y)`.
   Update cached `sy`, `ss`. Recompute `T = formt(sy, ss, theta)`. Set
   `updatd = true`.
3. **If rejected** (curvature too small or `formt` Cholesky failed):
   discard the proposed pair; set `updatd = false`. The `cauchy` and
   `subsm` calls in the next iteration will reuse the existing
   factorizations.
4. **Set** `xp = x`, `x = x_new`, `f = f_new`, `g = g_new`. Increment
   `iter`.

## Refresh / restart

In adverse cases (e.g., the line search fails, or repeatedly rejects
steps), the algorithm "refreshes" by discarding all L-BFGS history and
starting from the current point with `col = 0`, `theta = 1`. This
effectively resets the Hessian approximation to the identity. See
`mainlb.md` and `lnsrlb.md` for the exact triggers.

A refresh is a complete reset of the L-BFGS history; persistent state
*outside* the history (`x`, `f`, `g`, `iwhere`, counters) survives.

## State at termination

When the optimizer returns, the persistent state contains:
- The final iterate `x` and its `f`, `g`.
- Counters `iter`, `n_fg`.
- Termination reason in `info` and `message` (see `02_api.md`).

Transient state and L-BFGS history are no longer needed and may be
discarded. Ports that wrap the optimizer in an object are free to
retain history if useful (e.g., for warm-starting a follow-up
optimization), but the spec does not require this.

## Cardinality summary

For dimension `n` and memory `m`:

| Storage class | Approximate size (doubles) |
|---------------|----------------------------|
| Persistent -- iterate (`x`, `g`, `xp`) | `3n` |
| Persistent -- L-BFGS history (`S`, `Y`) | `2 * m * n` |
| Persistent -- cached (`sy`, `ss`, `T`) | `3 * m^2` |
| Transient -- Cauchy/subsm (`xc`, `c`, `r`, `d`, `K`, `K_chol`) | `~5n + 8m^2` |

Total typical: `5n + 2mn + 11m^2` doubles, plus `O(n)` integer state.
This matches the F77 `wa` (`2mn + 5n + 11m^2 + 8m`) and `iwa` (`3n`)
workspace sizes; the extra `8m` is small scratch in `cauchy`/`bmv`
that ports can also hold as locals.
