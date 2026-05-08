"""lnsrlb -- line-search wrapper around dcsrch with bound projection.

Reference implementation following ``docs/spec/subroutines/lnsrlb.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .dcsrch import dcsrch, DcsrchState


FTOL = 1.0e-3
GTOL = 0.9
XTOL = 0.1
BIG = 1.0e10


@dataclass
class LnsrlbState:
    """State carried across F77-style reverse-comm calls."""

    fold: float
    gd: float
    gdold: float
    stp: float
    dnorm: float
    dtd: float
    xstep: float
    stpmx: float
    ifun: int
    iback: int
    nfgv: int
    info: int
    task: str
    csave: str
    isave: np.ndarray              # length 2
    dsave: np.ndarray              # length 13


def lnsrlb(
    n: int,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    x: np.ndarray,
    f: float,
    g: np.ndarray,
    d: np.ndarray,
    r: np.ndarray,
    t: np.ndarray,
    z: np.ndarray,
    iter: int,
    boxed: bool,
    cnstnd: bool,
    state: LnsrlbState,
) -> float:
    """Single reverse-comm call. Mutates state and x; returns updated f."""
    if state.task[:5] != "FG_LN":
        # Setup path.
        state.dtd = float(np.dot(d, d))
        state.dnorm = float(np.sqrt(state.dtd))

        # Determine stpmx.
        state.stpmx = BIG
        if cnstnd:
            if iter == 0:
                state.stpmx = 1.0
            else:
                for i in range(n):
                    a1 = d[i]
                    if nbd[i] != 0:
                        if a1 < 0.0 and nbd[i] <= 2:
                            a2 = l[i] - x[i]
                            if a2 >= 0.0:
                                state.stpmx = 0.0
                            elif a1 * state.stpmx < a2:
                                state.stpmx = a2 / a1
                        elif a1 > 0.0 and nbd[i] >= 2:
                            a2 = u[i] - x[i]
                            if a2 <= 0.0:
                                state.stpmx = 0.0
                            elif a1 * state.stpmx > a2:
                                state.stpmx = a2 / a1

        if iter == 0 and not boxed:
            state.stp = min(1.0 / state.dnorm, state.stpmx)
        else:
            state.stp = 1.0

        t[:] = x
        r[:] = g
        state.fold = f
        state.ifun = 0
        state.iback = 0
        state.csave = "START"

    # Continuation path (and end of setup path).
    state.gd = float(np.dot(g, d))
    if state.ifun == 0:
        state.gdold = state.gd
        if state.gd >= 0.0:
            # Ascent direction.
            state.info = -4
            return f

    # Call dcsrch.
    ds = DcsrchState(task=state.csave, isave=state.isave, dsave=state.dsave)
    f_new, gd_new, stp_new = dcsrch(
        f, state.gd, state.stp, FTOL, GTOL, XTOL, 0.0, state.stpmx, ds,
    )
    state.csave = ds.task
    state.isave[:] = ds.isave
    state.dsave[:] = ds.dsave
    state.gd = gd_new
    state.stp = stp_new

    state.xstep = state.stp * state.dnorm
    if state.csave[:4] not in ("CONV", "WARN"):
        state.task = "FG_LNSRCH"
        state.ifun += 1
        state.nfgv += 1
        state.iback = state.ifun - 1
        if state.stp == 1.0:
            x[:] = z
        else:
            for i in range(n):
                x[i] = state.stp * d[i] + t[i]
    else:
        state.task = "NEW_X"

    return f
