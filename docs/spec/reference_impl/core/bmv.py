"""bmv -- compact L-BFGS middle-matrix solve M^{-1} v.

Reference implementation following ``docs/spec/subroutines/bmv.md``.
Uses ``scipy.linalg.solve_triangular`` for the two block solves;
the rest is explicit loops to preserve summation order.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular


def bmv(
    m: int,
    sy: np.ndarray,
    wt: np.ndarray,
    col: int,
    v: np.ndarray,
    p: np.ndarray,
) -> None:
    """Compute ``p = M^{-1} v`` using the compact L-BFGS factorization.

    Parameters
    ----------
    m : int
        Leading dimension of ``sy`` and ``wt``.
    sy : (m, m) float64
        Holds ``D`` on the diagonal and ``L`` (strict lower triangle).
    wt : (m, m) float64
        Upper Cholesky factor ``T`` of ``theta * S'S + L D^{-1} L'``.
    col : int
        Active size (0 <= col <= m).
    v : (2*col,) float64
        Right-hand side.
    p : (2*col,) float64
        Output (overwritten). On col==0, untouched.

    Notes
    -----
    The F77 routine used to take an ``info`` output parameter for the
    LINPACK ``dtrsl`` triangular-solve return status. Since LAPACK
    ``dtrsm`` cannot fail on a non-singular factor, the parameter is
    not present in this version.
    """
    if col == 0:
        return

    # Part I: solve [D^{1/2} 0; -L D^{-1/2}  J] [p1; p2] = [v1; v2]
    # Step 1: build (v2 + L D^{-1} v1) into p[col+1..2col]  (1-based).
    p[col] = v[col]                                   # p(col+1) = v(col+1)
    for i in range(2, col + 1):                       # i = 2..col, 1-based
        s = 0.0
        for k in range(1, i):                         # k = 1..i-1, 1-based
            s += sy[i - 1, k - 1] * v[k - 1] / sy[k - 1, k - 1]
        p[col + i - 1] = v[col + i - 1] + s

    # Step 2: solve J^T x = p[col+1..2col] in place. F77 dtrsm 'l','u','t','n'.
    p[col : col + col] = solve_triangular(
        wt[:col, :col], p[col : col + col], lower=False, trans="T"
    )

    # Step 3: p1 = v1 / sqrt(D).
    for i in range(1, col + 1):
        p[i - 1] = v[i - 1] / np.sqrt(sy[i - 1, i - 1])

    # Part II: solve [-D^{1/2} D^{-1/2} L'; 0 J'] [p1; p2] = [p1; p2]
    # Step 1: solve J x = p[col+1..2col] in place. F77 dtrsm 'l','u','n','n'.
    p[col : col + col] = solve_triangular(
        wt[:col, :col], p[col : col + col], lower=False, trans="N"
    )

    # Step 2: p1 = -D^{-1/2} p1_pre.
    for i in range(1, col + 1):
        p[i - 1] = -p[i - 1] / np.sqrt(sy[i - 1, i - 1])

    # Step 3: add the D^{-1} L' p2 correction.
    for i in range(1, col + 1):
        s = 0.0
        for k in range(i + 1, col + 1):
            s += sy[k - 1, i - 1] * p[col + k - 1] / sy[i - 1, i - 1]
        p[i - 1] = p[i - 1] + s
