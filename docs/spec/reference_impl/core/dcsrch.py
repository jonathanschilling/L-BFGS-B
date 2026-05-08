"""dcsrch -- More-Thuente line search reverse-comm driver.

Reference implementation following ``docs/spec/subroutines/dcsrch.md``.
Mirrors F77's reverse-comm style for cross-validation; modern ports
should wrap this in a callback API.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .dcstep import dcstep, DcstepState


P5 = 0.5
P66 = 0.66
XTRAPL = 1.1
XTRAPU = 4.0


@dataclass
class DcsrchState:
    """Carries dcsrch's working state across reverse-comm calls."""

    task: str
    isave: np.ndarray              # length 2
    dsave: np.ndarray              # length 13


def _save(brackt: bool, stage: int, ginit: float, gtest: float,
          gx: float, gy: float, finit: float, fx: float, fy: float,
          stx: float, sty: float, stmin: float, stmax: float,
          width: float, width1: float, isave: np.ndarray, dsave: np.ndarray) -> None:
    isave[0] = 1 if brackt else 0
    isave[1] = stage
    dsave[0] = ginit
    dsave[1] = gtest
    dsave[2] = gx
    dsave[3] = gy
    dsave[4] = finit
    dsave[5] = fx
    dsave[6] = fy
    dsave[7] = stx
    dsave[8] = sty
    dsave[9] = stmin
    dsave[10] = stmax
    dsave[11] = width
    dsave[12] = width1


def _restore(isave: np.ndarray, dsave: np.ndarray):
    brackt = (isave[0] == 1)
    stage = int(isave[1])
    return (brackt, stage,
            dsave[0], dsave[1], dsave[2], dsave[3],
            dsave[4], dsave[5], dsave[6],
            dsave[7], dsave[8], dsave[9], dsave[10],
            dsave[11], dsave[12])


def dcsrch(
    f: float,
    g: float,
    stp: float,
    ftol: float,
    gtol: float,
    xtol: float,
    stpmin: float,
    stpmax: float,
    state: DcsrchState,
) -> tuple[float, float, float]:
    """Single reverse-comm call. Returns updated (f, g, stp) and updates state."""
    if state.task[:5] == "START":
        # Input validation: F77 sequential overwrite.
        if stp < stpmin:
            state.task = "ERROR: STP .LT. STPMIN"
        if stp > stpmax:
            state.task = "ERROR: STP .GT. STPMAX"
        if g >= 0.0:
            state.task = "ERROR: INITIAL G .GE. ZERO"
        if ftol < 0.0:
            state.task = "ERROR: FTOL .LT. ZERO"
        if gtol < 0.0:
            state.task = "ERROR: GTOL .LT. ZERO"
        if xtol < 0.0:
            state.task = "ERROR: XTOL .LT. ZERO"
        if stpmin < 0.0:
            state.task = "ERROR: STPMIN .LT. ZERO"
        if stpmax < stpmin:
            state.task = "ERROR: STPMAX .LT. STPMIN"

        if state.task[:5] == "ERROR":
            return f, g, stp

        # Initialize.
        brackt = False
        stage = 1
        finit = f
        ginit = g
        gtest = ftol * ginit
        width = stpmax - stpmin
        width1 = width / P5
        stx = 0.0
        fx = finit
        gx = ginit
        sty = 0.0
        fy = finit
        gy = ginit
        stmin = 0.0
        stmax = stp + XTRAPU * stp
        state.task = "FG"
        _save(brackt, stage, ginit, gtest, gx, gy, finit, fx, fy,
              stx, sty, stmin, stmax, width, width1,
              state.isave, state.dsave)
        return f, g, stp

    # Resume.
    (brackt, stage, ginit, gtest, gx, gy, finit, fx, fy,
     stx, sty, stmin, stmax, width, width1) = _restore(state.isave, state.dsave)

    ftest = finit + stp * gtest
    if stage == 1 and f <= ftest and g >= 0.0:
        stage = 2

    # Warnings.
    if brackt and (stp <= stmin or stp >= stmax):
        state.task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS"
    if brackt and (stmax - stmin) <= xtol * stmax:
        state.task = "WARNING: XTOL TEST SATISFIED"
    if stp == stpmax and f <= ftest and g <= gtest:
        state.task = "WARNING: STP = STPMAX"
    if stp == stpmin and (f > ftest or g >= gtest):
        state.task = "WARNING: STP = STPMIN"

    # Convergence.
    if f <= ftest and abs(g) <= gtol * (-ginit):
        state.task = "CONVERGENCE"

    # Termination check.
    if state.task[:4] in ("WARN", "CONV"):
        _save(brackt, stage, ginit, gtest, gx, gy, finit, fx, fy,
              stx, sty, stmin, stmax, width, width1,
              state.isave, state.dsave)
        return f, g, stp

    # dcstep call (modified function in stage 1 if appropriate).
    if stage == 1 and f <= fx and f > ftest:
        fm = f - stp * gtest
        fxm = fx - stx * gtest
        fym = fy - sty * gtest
        gm = g - gtest
        gxm = gx - gtest
        gym = gy - gtest

        ds = DcstepState(stx=stx, fx=fxm, dx=gxm, sty=sty, fy=fym, dy=gym,
                         stp=stp, brackt=brackt)
        dcstep(ds, fm, gm, stmin, stmax)
        stx, fxm, gxm = ds.stx, ds.fx, ds.dx
        sty, fym, gym = ds.sty, ds.fy, ds.dy
        stp, brackt = ds.stp, ds.brackt

        fx = fxm + stx * gtest
        fy = fym + sty * gtest
        gx = gxm + gtest
        gy = gym + gtest
    else:
        ds = DcstepState(stx=stx, fx=fx, dx=gx, sty=sty, fy=fy, dy=gy,
                         stp=stp, brackt=brackt)
        dcstep(ds, f, g, stmin, stmax)
        stx, fx, gx = ds.stx, ds.fx, ds.dx
        sty, fy, gy = ds.sty, ds.fy, ds.dy
        stp, brackt = ds.stp, ds.brackt

    # Bisection if bracket isn't shrinking fast enough.
    if brackt:
        if abs(sty - stx) >= P66 * width1:
            stp = stx + P5 * (sty - stx)
        width1 = width
        width = abs(sty - stx)

    # Step bounds.
    if brackt:
        stmin = min(stx, sty)
        stmax = max(stx, sty)
    else:
        stmin = stp + XTRAPL * (stp - stx)
        stmax = stp + XTRAPU * (stp - stx)

    # Force stp in [stpmin, stpmax].
    stp = max(stp, stpmin)
    stp = min(stp, stpmax)

    # If progress impossible, fall back to best step.
    if (brackt and (stp <= stmin or stp >= stmax)) or \
       (brackt and (stmax - stmin) <= xtol * stmax):
        stp = stx

    state.task = "FG"
    _save(brackt, stage, ginit, gtest, gx, gy, finit, fx, fy,
          stx, sty, stmin, stmax, width, width1,
          state.isave, state.dsave)
    return f, g, stp
