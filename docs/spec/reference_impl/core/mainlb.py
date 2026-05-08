"""mainlb -- L-BFGS-B algorithm loop in callback form.

Reference implementation following ``docs/spec/subroutines/mainlb.md``.
This is a callback-based reimplementation, NOT a literal port of
``src/mainlb.f``'s reverse-comm state machine. It calls the other
core/* modules (projgr, errclb, active, cauchy, freev, cmprlb, formk,
subsm, lnsrlb, matupd, formt) directly.

Functional correctness on small problems is validated by the
conformance runner (see docs/spec/runner/, Phase C). Bit-for-bit
agreement with the F77 reference under the same BLAS pin is the
target for ``--strict`` conformance mode; the conformance runner
exercises that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from .active import active
from .bmv import bmv
from .cauchy import cauchy
from .cmprlb import cmprlb
from .dcsrch import DcsrchState
from .errclb import ErrclbState, errclb
from .formk import formk
from .formt import formt
from .freev import freev
from .lnsrlb import LnsrlbState, lnsrlb
from .matupd import MatupdState, matupd
from .projgr import projgr


INFO_CONVERGED_PGTOL = 0
INFO_CONVERGED_FACTR = 1
INFO_LIMIT_ITER = 2
INFO_LIMIT_FG = 3
INFO_ABNORMAL_LNSRLB = 4
INFO_USER_STOP = 5


@dataclass
class MinimizeResult:
    x: np.ndarray
    f: float
    g: np.ndarray
    info: int
    n_iter: int
    n_fg: int
    final_projg: float
    message: str


def minimize(
    n: int,
    x0: np.ndarray,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    m: int,
    f_eval: Callable[[np.ndarray], float],
    g_eval: Callable[[np.ndarray], np.ndarray],
    factr: float = 1.0e7,
    pgtol: float = 1.0e-5,
    max_iter: Optional[int] = None,
    max_fg: Optional[int] = None,
    iprint: int = -1,
) -> MinimizeResult:
    """Minimize ``f`` over the box ``l <= x <= u`` via L-BFGS-B.

    See ``docs/spec/02_api.md`` for the contract.
    """
    # ---- Phase 0: validation, projection, initial state.
    err = ErrclbState(task="START", info=0, k=0)
    errclb(n, m, factr, l, u, nbd, err)
    if err.task[:5] == "ERROR":
        return MinimizeResult(
            x=np.array(x0, copy=True), f=0.0, g=np.zeros(n),
            info=-abs(err.info) if err.info < 0 else -1, n_iter=0, n_fg=0,
            final_projg=0.0, message=err.task,
        )

    x = np.array(x0, dtype=np.float64).copy()
    iwhere = np.zeros(n, dtype=np.int32)
    a = active(n, l, u, nbd, x, iwhere, iprint)
    cnstnd, boxed = a.cnstnd, a.boxed

    eps = np.finfo(np.float64).eps

    # L-BFGS history.
    ws = np.zeros((n, m))
    wy = np.zeros((n, m))
    sy = np.zeros((m, m))
    ss = np.zeros((m, m))
    wt = np.zeros((m, m))
    wn = np.zeros((2 * m, 2 * m))
    wn1 = np.zeros((2 * m, 2 * m))

    head = 1
    itail = 0
    col = 0
    theta = 1.0
    iupdat = 0
    updatd = False
    nfree_prev = 0
    index = np.arange(1, n + 1, dtype=np.int32)
    indx2 = np.zeros(n, dtype=np.int32)

    # Initial f, g.
    f = float(f_eval(x))
    g = np.array(g_eval(x), dtype=np.float64).copy()
    n_fg = 1

    sbgnrm = projgr(n, x, l, u, nbd, g)
    if sbgnrm <= pgtol:
        return MinimizeResult(
            x=x, f=f, g=g, info=INFO_CONVERGED_PGTOL,
            n_iter=0, n_fg=n_fg, final_projg=sbgnrm,
            message="CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",
        )

    # Workspaces sized for the largest col (= m).
    p = np.zeros(2 * m)
    c = np.zeros(2 * m)
    wbp = np.zeros(2 * m)
    v = np.zeros(2 * m)
    iorder = np.zeros(n, dtype=np.int32)
    t = np.zeros(n)
    d_dir = np.zeros(n)
    xc = np.zeros(n)
    xp = np.zeros(n)
    r = np.zeros(n)
    wa = np.zeros(4 * m)

    n_iter = 0
    while True:
        n_iter += 1
        if max_iter is not None and n_iter > max_iter:
            return MinimizeResult(
                x=x, f=f, g=g, info=INFO_LIMIT_ITER,
                n_iter=n_iter - 1, n_fg=n_fg, final_projg=sbgnrm,
                message=f"ITERATION LIMIT ({max_iter}) REACHED",
            )

        f_old = f
        x_prev = x.copy()
        g_old = g.copy()

        # Cauchy point. (No info return: bmv cannot fail under LAPACK.)
        nseg = cauchy(
            n, x, l, u, nbd, g, iorder, iwhere, t, d_dir, xc,
            m, wy, ws, sy, wt, theta, col, head,
            p, c, wbp, v, sbgnrm, eps, iprint,
        )

        # Free / active partition + change detection.
        fr = freev(n, nfree_prev, index, indx2, iwhere,
                   updatd, cnstnd, n_iter - 1, iprint)
        nfree = fr.nfree
        nenter = fr.nenter
        ileave = fr.ileave
        wrk = fr.wrk

        # Reduced gradient at Cauchy point on free variables.
        # cmprlb expects wa[2m+1..2m+2col] to hold W'(z-x); cauchy filled c with that.
        wa_for_cmprlb = wa.copy()
        wa_for_cmprlb[2 * m : 2 * m + 2 * col] = c[: 2 * col]
        cmprlb(n, m, x_prev, g, ws, wy, sy, wt, xc, r, wa_for_cmprlb,
               index, theta, col, head, nfree, cnstnd)

        # Subspace minimization (only if there are free variables).
        x_subsm = xc.copy()
        if nfree > 0 and col > 0:
            if wrk:
                info_fk = formk(n, nfree, index, nenter, ileave, indx2,
                                iupdat, updatd, wn, wn1, m, ws, wy, sy,
                                theta, col, head)
                if info_fk != 0:
                    col = 0; theta = 1.0; iupdat = 0; updatd = False
                    continue
            d_sub = r[:nfree].copy()
            iword = subsm_call(n, m, nfree, index, l, u, nbd,
                               x_subsm, d_sub, xp, ws, wy,
                               theta, x_prev, g, col, head, wn, iprint)

        # Projected line search.
        d = x_subsm - x_prev
        ls = LnsrlbState(
            fold=0.0, gd=0.0, gdold=0.0, stp=0.0,
            dnorm=0.0, dtd=0.0, xstep=0.0, stpmx=0.0,
            ifun=0, iback=0, nfgv=n_fg, info=0,
            task="START", csave="",
            isave=np.zeros(2, dtype=np.int32),
            dsave=np.zeros(13, dtype=np.float64),
        )
        x_search = x_prev.copy()
        f_cur = f_old
        g_cur = g_old.copy()
        r_buf = np.zeros(n)
        t_buf = np.zeros(n)

        # Drive the line-search reverse-comm loop.
        while True:
            f_cur = lnsrlb(
                n, l, u, nbd, x_search, f_cur, g_cur, d, r_buf, t_buf, x_subsm,
                n_iter - 1, boxed, cnstnd, ls,
            )
            if ls.info == -4:
                # Ascent direction: refresh and retry.
                col = 0; theta = 1.0; iupdat = 0; updatd = False
                break
            if ls.task[:9] == "FG_LNSRCH":
                f_cur = float(f_eval(x_search))
                g_cur = np.array(g_eval(x_search), dtype=np.float64).copy()
                n_fg += 1
                if max_fg is not None and n_fg >= max_fg:
                    return MinimizeResult(
                        x=x_search, f=f_cur, g=g_cur, info=INFO_LIMIT_FG,
                        n_iter=n_iter, n_fg=n_fg, final_projg=sbgnrm,
                        message=f"F/G EVAL LIMIT ({max_fg}) REACHED",
                    )
            elif ls.task[:5] == "NEW_X":
                f = f_cur
                x = x_search.copy()
                g = g_cur.copy()
                break
            else:
                # CONV / WARN / ABNORMAL from dcsrch.
                if ls.csave[:4] == "WARN" or ls.csave[:5] == "ERROR":
                    return MinimizeResult(
                        x=x_prev, f=f_old, g=g_old, info=INFO_ABNORMAL_LNSRLB,
                        n_iter=n_iter, n_fg=n_fg, final_projg=sbgnrm,
                        message="ABNORMAL TERMINATION IN LNSRLB",
                    )
                f = f_cur; x = x_search.copy(); g = g_cur.copy()
                break

        if ls.info == -4:
            continue                            # refreshed; redo iteration

        # Convergence checks.
        sbgnrm = projgr(n, x, l, u, nbd, g)
        if sbgnrm <= pgtol:
            return MinimizeResult(
                x=x, f=f, g=g, info=INFO_CONVERGED_PGTOL,
                n_iter=n_iter, n_fg=n_fg, final_projg=sbgnrm,
                message="CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",
            )
        ddum = max(abs(f_old), abs(f), 1.0)
        if (f_old - f) / ddum <= factr * eps:
            return MinimizeResult(
                x=x, f=f, g=g, info=INFO_CONVERGED_FACTR,
                n_iter=n_iter, n_fg=n_fg, final_projg=sbgnrm,
                message="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",
            )

        # L-BFGS update.
        s = x - x_prev
        y = g - g_old
        sy_new = float(np.dot(s, y))
        yy = float(np.dot(y, y))
        if sy_new > eps * yy:
            iupdat += 1
            mu = MatupdState(itail=itail, col=col, head=head, theta=theta)
            matupd(n, m, ws, wy, sy, ss, s, y, iupdat, mu, yy, sy_new, 1.0,
                   float(np.dot(s, s)))
            itail = mu.itail; col = mu.col; head = mu.head; theta = mu.theta
            info_t = formt(m, wt, sy, ss, col, theta)
            if info_t != 0:
                # Cholesky failure -> refresh.
                col = 0; theta = 1.0; iupdat = 0; updatd = False
                continue
            updatd = True
        else:
            updatd = False

        nfree_prev = nfree


def subsm_call(n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, theta, xx, gg,
               col, head, wn, iprint):
    """Wrapper: subsm needs full-length d but uses only the first nsub entries."""
    from .subsm import subsm
    wv = np.zeros(2 * m)
    return subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, theta, xx, gg,
                 col, head, wv, wn, iprint)
