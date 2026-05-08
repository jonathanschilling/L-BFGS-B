# mainlb

## Purpose

The L-BFGS-B algorithm loop. Drives the iteration: per outer step,
identifies the active set, computes the generalized Cauchy point,
performs subspace minimization, runs a projected line search, updates
the L-BFGS history, and tests for convergence.

In the F77 source `mainlb` is reverse-comm: it returns to `setulb`
(and thence to the user) when `f`/`g` are needed, with state preserved
in `isave`/`dsave`. **In this spec, `mainlb` is described in
callback form**: it takes `f_eval` and `g_eval` callbacks and runs the
algorithm as a straight loop. Ports may implement either form; the
state-machine F77 form is captured separately in
`08_legacy_reverse_comm.md`.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n`, `m` | positive integers | Dimension and memory parameter. |
| `x0` | real vector, length `n` | Initial point (will be projected feasible). |
| `l`, `u`, `nbd` | (as in `02_api.md`) | Bounds. |
| `f_eval(x) -> f` | callback | Function evaluator. |
| `g_eval(x) -> g` | callback | Gradient evaluator. |
| `factr`, `pgtol` | nonneg reals | Convergence tolerances. |
| `max_iter`, `max_fg` | optional integer | Iteration / evaluation caps. |
| `iprint` | integer | Diagnostic verbosity. |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `x` | real vector, length `n` | Final iterate. |
| `f` | real | `f(x)` at final iterate. |
| `g` | real vector, length `n` | `gradf(x)` at final iterate. |
| `info` | integer | Termination code (see `02_api.md`). |
| `n_iter`, `n_fg` | integers | Counters. |
| `message` | string | Human-readable termination reason. |

## Algorithm (callback form)

```
def mainlb(n, m, x0, l, u, nbd, f_eval, g_eval, factr, pgtol, max_iter, max_fg, iprint):

    # ============================================================
    # Phase 0: validation and initialization
    # ============================================================
    info, message = errclb(n, m, factr, l, u, nbd)
    if info != 0:
        return result(info=info, message=message)

    x = copy(x0)
    iwhere = zeros(n, int)
    prjctd, cnstnd, boxed = active(n, l, u, nbd, x, iwhere, iprint)

    eps = compute_epsmch()                    # IEEE: 2.22e-16

    # L-BFGS state (initially empty).
    S = empty(n, m); Y = empty(n, m); col = 0; head = 1; itail = 0; theta = 1.0
    sy = zeros(m, m); ss = zeros(m, m); T = zeros(m, m)
    wn = zeros(2m, 2m); wn1 = zeros(2m, 2m)

    # Counters and flags.
    n_iter = 0; n_fg = 0; iupdat = 0
    updatd = false                            # L-BFGS updated last iter?

    # First f, g.
    f = f_eval(x); g = g_eval(x); n_fg += 1

    # Initial convergence check: |g_proj|_inf <= pgtol.
    sbgnrm = projgr(n, l, u, nbd, x, g)
    if sbgnrm <= pgtol:
        return result(x, f, g, info=INFO_CONVERGED_PGTOL, n_iter=0, n_fg=1)

    # ============================================================
    # Phase 1: outer iteration loop
    # ============================================================
    while true:
        n_iter += 1
        if max_iter and n_iter > max_iter:
            return result(info=INFO_LIMIT_ITER)

        f_old = f
        x_prev = copy(x)

        # 1a. Generalized Cauchy point.
        xc, d_cauchy, c, info = cauchy(n, x, l, u, nbd, g, iwhere, ws=S, wy=Y,
                                        sy=sy, wt=T, theta=theta, col=col, head=head,
                                        sbgnrm=sbgnrm, eps=eps)
        if info != 0:                          # bmv failure inside cauchy
            handle_bmv_failure(); continue     # refresh L-BFGS history

        # 1b. Identify free / active partition; detect entering/leaving.
        nfree, index, nenter, ileave, indx2, wrk = freev(
            n, prev_nfree=col, prev_index=index, iwhere=iwhere,
            updatd=updatd, cnstnd=cnstnd, iter=n_iter, iprint=iprint,
        )

        # 1c. Compute reduced gradient r at xc on free variables.
        r = empty(n)
        info = cmprlb(n, m, x, g, S, Y, sy, T, xc, r, wa, index,
                      theta, col, head, nfree, cnstnd)
        if info != 0:
            handle_bmv_failure(); continue

        # 1d. Subspace minimization (only if nfree > 0 and col > 0).
        if nfree > 0:
            if wrk:
                info = formk(n, nfree, index, nenter, ileave, indx2,
                             iupdat, updatd, wn, wn1, m, S, Y, sy,
                             theta, col, head)
                if info != 0:
                    handle_formk_failure(); continue

            x_subsm = copy(x)
            iword, info = subsm(n, m, nfree, index, l, u, nbd, x_subsm,
                                d_subsm=r[:nfree], xp=empty(n),
                                ws=S, wy=Y, theta=theta, xx=x_prev, gg=g,
                                col=col, head=head, wv=zeros(2m), wn=wn)
        else:
            x_subsm = xc                       # no free vars: stay at GCP

        # 1e. Projected line search.
        d = x_subsm - x_prev                   # search direction
        f, x_new, g_new, lnsrlb_status = lnsrlb(
            n, l, u, nbd, x_prev, f, g, d, z=x_subsm,
            iter=n_iter, boxed=boxed, cnstnd=cnstnd,
            f_eval=f_eval, g_eval=g_eval,
        )
        if lnsrlb_status == ABNORMAL:
            return result(info=INFO_ABNORMAL_LNSRLB)
        n_fg += lnsrlb.n_evals

        x = x_new; g = g_new
        sbgnrm = projgr(n, l, u, nbd, x, g)

        # 1f. Convergence tests.
        if sbgnrm <= pgtol:
            return result(x, f, g, info=INFO_CONVERGED_PGTOL, n_iter=n_iter, n_fg=n_fg)
        rel_reduction = (f_old - f) / max(abs(f_old), abs(f), 1.0)
        if rel_reduction <= factr * eps:
            return result(x, f, g, info=INFO_CONVERGED_FACTR, n_iter=n_iter, n_fg=n_fg)

        if max_fg and n_fg >= max_fg:
            return result(info=INFO_LIMIT_FG)

        # 1g. Update L-BFGS history.
        s = x - x_prev
        y = g - g_old
        ddum = max(abs(f_old), abs(f), 1.0)    # for relative test
        sy_new = dot(s, y)
        if sy_new > eps * ||y||**2:              # curvature condition
            iupdat += 1
            matupd(n, m, S, Y, sy, ss, d=s, r=y, iupdat=iupdat,
                   itail=itail, col=col, head=head, theta=theta,
                   rr=||y||**2, dr=sy_new, stp=1.0, dtd=||s||**2)
            updatd = true
            info = formt(m, T, sy, ss, col, theta)
            if info != 0:
                handle_formt_failure(); continue
        else:
            updatd = false                     # skip update this iter

        # Loop back.
```

### Magic constants

| Constant | F77 | Value | Meaning |
|----------|-----|-------|---------|
| `eps * ||y||^2` curvature gate | (matupd caller) | machine eps | Refresh threshold. |
| All others | (in callees) | -- | See `04_numerics.md`. |

### Numerical safeguards

- `errclb` validates inputs; failure -> immediate return.
- Sub-calls (`cauchy`, `cmprlb`, `formk`, `formt`) report Cholesky /
  matrix-solve failures via nonzero `info`. The handler refreshes
  the L-BFGS history (`col = 0`, `theta = 1`) and continues.
- The curvature gate prevents bad pairs from entering the L-BFGS
  cache; affected iterations have `updatd = false`.
- `max_iter` / `max_fg` caps prevent infinite loops.

### Order-of-operations dependencies

Within `mainlb` itself, the per-iteration sub-routine call sequence
is fixed:
**`projgr -> cauchy -> freev -> cmprlb -> formk(if wrk) -> subsm -> lnsrlb ->
projgr (again) -> matupd / formt`**.

Reordering would change the algorithm. The order matches the F77
source's state-machine flow.

The reductions internal to each sub-routine are documented in their
respective spec files.

## Reverse-communication form

The F77 implementation interleaves the loop with the user-level
`f_eval`/`g_eval` calls via reverse-comm. Each spot where the
callback would be invoked above corresponds to a `task = 'FG'` return
in F77. The user supplies `f`/`g`, then re-enters; `mainlb` resumes
from the saved program-counter (encoded in `isave`).

For ports that need this protocol, see `08_legacy_reverse_comm.md`.
The Python reference impl uses callback form.

## Test vectors

The unit test `test_mainlb.f90` exercises **integration-style** flows
(driving full optimizations on small quadratics, with stopping
criteria triggered on `pgtol`, `factr`, user `'STOP'`, etc.). These
do not lend themselves to single-call JSON test vectors; they are
covered by the **conformance runner** (Phase C) in
`docs/spec/runner/conformance.py`, which drives the Python ref impl
through full optimizations and compares to the F77 reference.

The 7 F77 test cases:

| F77 case | Conformance scenario |
|---|---|
| `case_invalid_input_returns_error` | `n = 0` -> `INFO_INPUT_ERROR_*`. |
| `case_pgtol_convergence` | Quadratic, `pgtol = 1e-3` -> converges via projected gradient. |
| `case_user_signals_stop` | User sets `task = STOP` mid-iteration. |
| `case_factr_convergence` | Slow `pgtol`; `factr` fires first. |
| `case_immediate_pgtol_convergence` | `x0` already at minimum. |
| `case_user_stop_cpu_restores_iterate` | User stop after a partial iteration restores last iterate. |
| `case_iprint_99_diagnostics` | High-verbosity printing path. |

## Reference implementation

`reference_impl/core/mainlb.py` (callback form; depends on every
other `core` module).

## Cross-references

- **Paper**: `code.pdf` sec.3 (the algorithm flow). `algorithm.pdf`
  sec.1-sec.5 for derivation. `acm-remark.pdf` for the `subsm`
  safeguarding.
- **Related subroutines**: every other in-scope routine.
- **F77 source**: `src/mainlb.f`.
- **Unit test**: `tests/unit/test_mainlb.f90` (7 cases, integration-level).
