"""cmprlb -- reduced gradient r at the Cauchy point.

Reference implementation following ``docs/spec/subroutines/cmprlb.md``.
Depends on ``core.bmv`` for the inverse-middle-matrix solve.
"""

from __future__ import annotations

import numpy as np

from .bmv import bmv


def cmprlb(
    n: int,
    m: int,
    x: np.ndarray,
    g: np.ndarray,
    ws: np.ndarray,
    wy: np.ndarray,
    sy: np.ndarray,
    wt: np.ndarray,
    z: np.ndarray,
    r: np.ndarray,
    wa: np.ndarray,
    index: np.ndarray,
    theta: float,
    col: int,
    head: int,
    nfree: int,
    cnstnd: bool,
) -> None:
    """Compute reduced gradient r at the Cauchy point.

    The F77 routine used to take an ``info`` output for forwarding bmv
    failures. Since LAPACK ``dtrsm`` cannot fail on a non-singular
    factor, that signaling is no longer needed.
    """
    # Path A: unconstrained with history -> r = -g.
    if (not cnstnd) and col > 0:
        for i in range(n):
            r[i] = -g[i]
        return

    # Path B: constrained or col == 0.
    for i in range(nfree):
        k = index[i] - 1                      # convert to 0-based
        r[i] = -theta * (z[k] - x[k]) - g[k]

    # Embedded bmv solve. Provides M^{-1} W'(z - x) in wa[1..2*col]
    # (Fortran 1-based) -> wa[0..2*col-1] in Python.
    if col > 0:
        v_in = wa[2 * m : 2 * m + 2 * col].copy()       # wa(2m+1..2m+2col)
        p_out = np.empty(2 * col, dtype=np.float64)
        bmv(m, sy, wt, col, v_in, p_out)
        wa[: 2 * col] = p_out

        pointr = head
        for j in range(1, col + 1):                       # 1-based j
            a1 = wa[j - 1]                                # wa(j)
            a2 = theta * wa[col + j - 1]                  # theta * wa(col + j)
            for i in range(nfree):
                k = index[i] - 1
                r[i] += wy[k, pointr - 1] * a1 + ws[k, pointr - 1] * a2
            pointr = (pointr % m) + 1
