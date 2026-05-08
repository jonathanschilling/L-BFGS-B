"""cauchy -- Generalized Cauchy Point computation.

Reference implementation following ``docs/spec/subroutines/cauchy.md``.
Literal port of ``src/cauchy.f`` to preserve the breakpoint-walk order
and the exact F77 arithmetic for f1, f2, p, c updates.

Depends on ``core.bmv`` and ``core.hpsolb``.
"""

from __future__ import annotations

import numpy as np

from .bmv import bmv
from .hpsolb import hpsolb


def cauchy(
    n: int,
    x: np.ndarray,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    g: np.ndarray,
    iorder: np.ndarray,
    iwhere: np.ndarray,
    t: np.ndarray,
    d: np.ndarray,
    xcp: np.ndarray,
    m: int,
    wy: np.ndarray,
    ws: np.ndarray,
    sy: np.ndarray,
    wt: np.ndarray,
    theta: float,
    col: int,
    head: int,
    p: np.ndarray,
    c: np.ndarray,
    wbp: np.ndarray,
    v: np.ndarray,
    sbgnrm: float,
    epsmch: float,
    iprint: int = -1,
) -> tuple[int, int]:
    """Compute the Generalized Cauchy Point.

    Returns
    -------
    info : int
        0 success; nonzero if bmv reports singularity.
    nseg : int
        Number of piecewise-quadratic segments traversed.
    """
    info = 0

    # Early return for zero projected gradient.
    if sbgnrm <= 0.0:
        xcp[:] = x
        return info, 0

    bnded = True
    nfree = n + 1
    nbreak = 0
    ibkmin = 0
    bkmin = 0.0
    col2 = 2 * col
    f1 = 0.0

    # Initialize p = 0.
    for i in range(col2):
        p[i] = 0.0

    # Phase 1: bound status, Cauchy direction, breakpoints.
    for i in range(1, n + 1):                          # i = 1..n
        neggi = -g[i - 1]
        if iwhere[i - 1] != 3 and iwhere[i - 1] != -1:
            tl = 0.0
            tu = 0.0
            if nbd[i - 1] <= 2:
                tl = x[i - 1] - l[i - 1]
            if nbd[i - 1] >= 2:
                tu = u[i - 1] - x[i - 1]
            xlower = (nbd[i - 1] <= 2) and (tl <= 0.0)
            xupper = (nbd[i - 1] >= 2) and (tu <= 0.0)

            iwhere[i - 1] = 0
            if xlower:
                if neggi <= 0.0:
                    iwhere[i - 1] = 1
            elif xupper:
                if neggi >= 0.0:
                    iwhere[i - 1] = 2
            else:
                if abs(neggi) <= 0.0:
                    iwhere[i - 1] = -3

        pointr = head
        if iwhere[i - 1] != 0 and iwhere[i - 1] != -1:
            d[i - 1] = 0.0
        else:
            d[i - 1] = neggi
            f1 = f1 - neggi * neggi
            for j in range(1, col + 1):
                p[j - 1] += wy[i - 1, pointr - 1] * neggi
                p[col + j - 1] += ws[i - 1, pointr - 1] * neggi
                pointr = (pointr % m) + 1
            if (nbd[i - 1] <= 2) and (nbd[i - 1] != 0) and (neggi < 0.0):
                # Bounded toward lower.
                nbreak += 1
                tl = x[i - 1] - l[i - 1]
                iorder[nbreak - 1] = i
                t[nbreak - 1] = tl / (-neggi)
                if nbreak == 1 or t[nbreak - 1] < bkmin:
                    bkmin = t[nbreak - 1]
                    ibkmin = nbreak
            elif (nbd[i - 1] >= 2) and (neggi > 0.0):
                # Bounded toward upper.
                nbreak += 1
                tu = u[i - 1] - x[i - 1]
                iorder[nbreak - 1] = i
                t[nbreak - 1] = tu / neggi
                if nbreak == 1 or t[nbreak - 1] < bkmin:
                    bkmin = t[nbreak - 1]
                    ibkmin = nbreak
            else:
                # No breakpoint along d.
                nfree -= 1
                iorder[nfree - 1] = i
                if abs(neggi) > 0.0:
                    bnded = False

    # Adjust p for theta != 1.
    if theta != 1.0:
        for j in range(col):
            p[col + j] *= theta

    # Initialize xcp = x.
    xcp[:] = x

    # Trivial case: no breakpoints and all variables effectively fixed.
    if nbreak == 0 and nfree == n + 1:
        return info, 0

    # Phase 2: compute initial f2 and dtm.
    for j in range(col2):
        c[j] = 0.0

    f2 = -theta * f1
    f2_org = f2
    if col > 0:
        info = bmv(m, sy, wt, col, p[:col2].copy(), v[:col2])
        if info != 0:
            return info, 0
        f2 = f2 - float(np.dot(v[:col2], p[:col2]))
    dtm = -f1 / f2 if f2 != 0.0 else 0.0
    tsum = 0.0
    nseg = 1

    # Phase 3 / 4: piecewise quadratic walk if breakpoints exist.
    if nbreak == 0:
        # Skip walk; locate GCP in this (only) segment.
        pass
    else:
        nleft = nbreak
        iter_ = 1
        tj = 0.0
        in_walk = True
        while in_walk:
            tj0 = tj
            if iter_ == 1:
                tj = bkmin
                ibp = iorder[ibkmin - 1]
            else:
                if iter_ == 2:
                    if ibkmin != nbreak:
                        t[ibkmin - 1] = t[nbreak - 1]
                        iorder[ibkmin - 1] = iorder[nbreak - 1]
                hpsolb(nleft, t, iorder, iter_ - 2)
                tj = t[nleft - 1]
                ibp = iorder[nleft - 1]

            dt = tj - tj0

            if dtm < dt:
                # GCP in current segment.
                break

            tsum += dt
            nleft -= 1
            iter_ += 1
            dibp = d[ibp - 1]
            d[ibp - 1] = 0.0
            if dibp > 0.0:
                zibp = u[ibp - 1] - x[ibp - 1]
                xcp[ibp - 1] = u[ibp - 1]
                iwhere[ibp - 1] = 2
            else:
                zibp = l[ibp - 1] - x[ibp - 1]
                xcp[ibp - 1] = l[ibp - 1]
                iwhere[ibp - 1] = 1

            if nleft == 0 and nbreak == n:
                # All variables fixed; finish.
                dtm = dt
                in_walk = False
                # Final c update happens at very end (post-loop).
                # Skip the f1/f2 update + Phase 4 motion.
                tsum_skip_final_motion = True
                break

            nseg += 1
            dibp2 = dibp * dibp

            f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp
            f2 = f2 - theta * dibp2

            if col > 0:
                # c += dt * p
                for j in range(col2):
                    c[j] += dt * p[j]
                # wbp = W' e_{ibp} (with theta scaling on S half).
                pointr = head
                for j in range(1, col + 1):
                    wbp[j - 1] = wy[ibp - 1, pointr - 1]
                    wbp[col + j - 1] = theta * ws[ibp - 1, pointr - 1]
                    pointr = (pointr % m) + 1
                # v = M^{-1} wbp
                info = bmv(m, sy, wt, col, wbp[:col2].copy(), v[:col2])
                if info != 0:
                    return info, nseg
                wmc = float(np.dot(c[:col2], v[:col2]))
                wmp = float(np.dot(p[:col2], v[:col2]))
                wmw = float(np.dot(wbp[:col2], v[:col2]))
                # p -= dibp * wbp
                for j in range(col2):
                    p[j] -= dibp * wbp[j]
                f1 += dibp * wmc
                f2 += 2.0 * dibp * wmp - dibp2 * wmw

            f2 = max(epsmch * f2_org, f2)

            if nleft > 0:
                dtm = -f1 / f2
                continue                            # next iteration of walk
            elif bnded:
                f1 = 0.0
                f2 = 0.0
                dtm = 0.0
            else:
                dtm = -f1 / f2
            in_walk = False

        else:
            # `while` loop exited cleanly (no break) -- shouldn't reach here.
            pass

    # Phase 4: locate the GCP inside the current segment, then exit.
    if dtm <= 0.0:
        dtm = 0.0
    tsum += dtm

    # Move free variables and unfixed-bounded variables.
    # daxpy(n, tsum, d, 1, xcp, 1)
    for i in range(n):
        xcp[i] += tsum * d[i]

    # Update c = c + dtm * p (used for subsm).
    if col > 0:
        for j in range(col2):
            c[j] += dtm * p[j]

    return info, nseg
