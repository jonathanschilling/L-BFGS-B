"""dcstep -- safeguarded More-Thuente cubic-step generator.

Reference implementation following ``docs/spec/subroutines/dcstep.md``.
Literal port of ``src/dcstep.f`` to preserve bit-for-bit behavior on
the four-case dispatch and the cubic / quadratic / secant arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt


@dataclass
class DcstepState:
    """Bracket-and-trial state passed in/out of dcstep."""

    stx: float
    fx: float
    dx: float
    sty: float
    fy: float
    dy: float
    stp: float
    brackt: bool


def dcstep(
    state: DcstepState,
    fp: float,
    dp: float,
    stpmin: float,
    stpmax: float,
) -> None:
    """Compute next safeguarded step; mutate ``state`` in place.

    Parameters
    ----------
    state : DcstepState
        Carries (stx, fx, dx, sty, fy, dy, stp, brackt).
    fp, dp : float
        Function value and derivative at the current trial ``state.stp``.
    stpmin, stpmax : float
        Step-length bounds.
    """
    stx, fx, dx = state.stx, state.fx, state.dx
    sty, fy, dy = state.sty, state.fy, state.dy
    stp, brackt = state.stp, state.brackt

    p66 = 0.66
    two = 2.0
    three = 3.0

    sgnd = dp * (dx / abs(dx))

    if fp > fx:
        # Case 1: higher function value.
        theta = three * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), abs(dx), abs(dp))
        gamma = s * sqrt((theta / s) ** 2 - (dx / s) * (dp / s))
        if stp < stx:
            gamma = -gamma
        p = (gamma - dx) + theta
        q = ((gamma - dx) + gamma) + dp
        r = p / q
        stpc = stx + r * (stp - stx)
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx)
        if abs(stpc - stx) < abs(stpq - stx):
            stpf = stpc
        else:
            stpf = stpc + (stpq - stpc) / two
        brackt = True

    elif sgnd < 0.0:
        # Case 2: opposite-sign derivatives.
        theta = three * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), abs(dx), abs(dp))
        gamma = s * sqrt((theta / s) ** 2 - (dx / s) * (dp / s))
        if stp > stx:
            gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dx
        r = p / q
        stpc = stp + r * (stx - stp)
        stpq = stp + (dp / (dp - dx)) * (stx - stp)
        if abs(stpc - stp) > abs(stpq - stp):
            stpf = stpc
        else:
            stpf = stpq
        brackt = True

    elif abs(dp) < abs(dx):
        # Case 3: decreasing-magnitude same-sign derivative.
        theta = three * (fx - fp) / (stp - stx) + dx + dp
        s = max(abs(theta), abs(dx), abs(dp))
        gamma = s * sqrt(max(0.0, (theta / s) ** 2 - (dx / s) * (dp / s)))
        if stp > stx:
            gamma = -gamma
        p = (gamma - dp) + theta
        q = (gamma + (dx - dp)) + gamma
        r = p / q
        if r < 0.0 and gamma != 0.0:
            stpc = stp + r * (stx - stp)
        elif stp > stx:
            stpc = stpmax
        else:
            stpc = stpmin
        stpq = stp + (dp / (dp - dx)) * (stx - stp)

        if brackt:
            if abs(stpc - stp) < abs(stpq - stp):
                stpf = stpc
            else:
                stpf = stpq
            if stp > stx:
                stpf = min(stp + p66 * (sty - stp), stpf)
            else:
                stpf = max(stp + p66 * (sty - stp), stpf)
        else:
            if abs(stpc - stp) > abs(stpq - stp):
                stpf = stpc
            else:
                stpf = stpq
            stpf = min(stpmax, stpf)
            stpf = max(stpmin, stpf)

    else:
        # Case 4: same-sign, magnitude not decreasing.
        if brackt:
            theta = three * (fp - fy) / (sty - stp) + dy + dp
            s = max(abs(theta), abs(dy), abs(dp))
            gamma = s * sqrt((theta / s) ** 2 - (dy / s) * (dp / s))
            if stp > sty:
                gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p / q
            stpc = stp + r * (sty - stp)
            stpf = stpc
        elif stp > stx:
            stpf = stpmax
        else:
            stpf = stpmin

    # Update bracket data.
    if fp > fx:
        sty = stp
        fy = fp
        dy = dp
    else:
        if sgnd < 0.0:
            sty = stx
            fy = fx
            dy = dx
        stx = stp
        fx = fp
        dx = dp

    # New step.
    stp = stpf

    # Write back.
    state.stx, state.fx, state.dx = stx, fx, dx
    state.sty, state.fy, state.dy = sty, fy, dy
    state.stp, state.brackt = stp, brackt
