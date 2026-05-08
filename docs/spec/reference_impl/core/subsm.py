"""subsm -- subspace minimization with bound projection and Morales/Nocedal 2011 safeguard.

Reference implementation following ``docs/spec/subroutines/subsm.md``.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular


def subsm(
    n: int,
    m: int,
    nsub: int,
    ind: np.ndarray,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    x: np.ndarray,
    d: np.ndarray,
    xp: np.ndarray,
    ws: np.ndarray,
    wy: np.ndarray,
    theta: float,
    xx: np.ndarray,
    gg: np.ndarray,
    col: int,
    head: int,
    wv: np.ndarray,
    wn: np.ndarray,
    iprint: int = -1,
) -> int:
    """Compute the subspace minimizer; return ``iword``.

    Mutates x (overwrites with new iterate), d (Newton direction), xp (saved x),
    wv (workspace).

    The F77 routine used to take an ``info`` output parameter for an
    "ill-conditioned K" status. Since K's conditioning is checked in
    ``formk`` / ``formt`` and the two ``dtrsm`` calls cannot fail on a
    non-singular factor, the parameter has been removed.
    """
    if nsub <= 0:
        return 0

    col2 = 2 * col

    # Phase 1a: wv = W' Z d.
    pointr = head
    for i in range(1, col + 1):
        temp1 = 0.0
        temp2 = 0.0
        for j in range(1, nsub + 1):
            k = ind[j - 1] - 1
            temp1 += wy[k, pointr - 1] * d[j - 1]
            temp2 += ws[k, pointr - 1] * d[j - 1]
        wv[i - 1] = temp1
        wv[col + i - 1] = theta * temp2
        pointr = (pointr % m) + 1

    # Phase 1b: wv := K^{-1} wv via wn factor (transpose then E-flip then no-transpose).
    if col2 > 0:
        wv[:col2] = solve_triangular(wn[:col2, :col2], wv[:col2], lower=False, trans="T")
        for i in range(col):
            wv[i] = -wv[i]
        wv[:col2] = solve_triangular(wn[:col2, :col2], wv[:col2], lower=False, trans="N")

    # Phase 1c: d = (1/theta) d + (1/theta^2) Z' W wv  (then a final division by theta).
    pointr = head
    for jy in range(1, col + 1):
        js = col + jy
        for i in range(1, nsub + 1):
            k = ind[i - 1] - 1
            d[i - 1] += (wy[k, pointr - 1] * wv[jy - 1] / theta
                          + ws[k, pointr - 1] * wv[js - 1])
        pointr = (pointr % m) + 1
    for i in range(nsub):
        d[i] = d[i] / theta

    # Phase 2: try projection; track if any bound is encountered.
    iword = 0
    xp[:] = x
    for i in range(1, nsub + 1):
        k = ind[i - 1] - 1
        dk = d[i - 1]
        xk = x[k]
        if nbd[k] != 0:
            if nbd[k] == 1:                            # lower only
                x[k] = max(l[k], xk + dk)
                if x[k] == l[k]:
                    iword = 1
            else:
                if nbd[k] == 2:                        # both bounds
                    xk2 = max(l[k], xk + dk)
                    x[k] = min(u[k], xk2)
                    if x[k] == l[k] or x[k] == u[k]:
                        iword = 1
                else:
                    if nbd[k] == 3:                    # upper only
                        x[k] = min(u[k], xk + dk)
                        if x[k] == u[k]:
                            iword = 1
        else:                                          # free
            x[k] = xk + dk

    # Phase 3: if no bound hit, done.
    if iword == 0:
        return iword

    # Phase 4: directional-derivative check (Morales/Nocedal 2011).
    dd_p = 0.0
    for i in range(n):
        dd_p += (x[i] - xx[i]) * gg[i]
    if dd_p > 0.0:
        # Positive dir derivative -> restore xp and fall through to backtracking.
        x[:] = xp
        # (Diagnostic; reference impl does not print.)
    else:
        return iword

    # Phase 5: backtracking with safeguarded step.
    alpha = 1.0
    temp1 = alpha
    ibd = 0
    for i in range(1, nsub + 1):
        k = ind[i - 1] - 1
        dk = d[i - 1]
        if nbd[k] != 0:
            if dk < 0.0 and nbd[k] <= 2:
                temp2 = l[k] - x[k]
                if temp2 >= 0.0:
                    temp1 = 0.0
                elif dk * alpha < temp2:
                    temp1 = temp2 / dk
            elif dk > 0.0 and nbd[k] >= 2:
                temp2 = u[k] - x[k]
                if temp2 <= 0.0:
                    temp1 = 0.0
                elif dk * alpha > temp2:
                    temp1 = temp2 / dk
            if temp1 < alpha:
                alpha = temp1
                ibd = i

    if alpha < 1.0:
        dk = d[ibd - 1]
        k = ind[ibd - 1] - 1
        if dk > 0.0:
            x[k] = u[k]
            d[ibd - 1] = 0.0
        elif dk < 0.0:
            x[k] = l[k]
            d[ibd - 1] = 0.0

    for i in range(1, nsub + 1):
        k = ind[i - 1] - 1
        x[k] = x[k] + alpha * d[i - 1]

    return iword
