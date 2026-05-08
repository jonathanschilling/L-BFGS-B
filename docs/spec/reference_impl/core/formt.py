"""formt -- form and Cholesky-factor T = theta*S'S + L D^{-1} L'.

Reference implementation following ``docs/spec/subroutines/formt.md``.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, LinAlgError


def formt(
    m: int,
    wt: np.ndarray,
    sy: np.ndarray,
    ss: np.ndarray,
    col: int,
    theta: float,
) -> int:
    """Form upper triangle of T and Cholesky-factor in place.

    On success, returns 0 with ``wt[:col, :col]`` holding the upper
    Cholesky factor. On failure, returns -3 (matches the F77 mapping
    of LAPACK ``dpotrf`` non-zero info).
    """
    # Phase 1: build T into upper triangle of wt[:col, :col].
    if col == 0:
        return 0

    for j in range(1, col + 1):
        wt[0, j - 1] = theta * ss[0, j - 1]

    for i in range(2, col + 1):
        for j in range(i, col + 1):
            s = 0.0
            for k in range(1, i):                # k = 1..i-1
                s += sy[i - 1, k - 1] * sy[j - 1, k - 1] / sy[k - 1, k - 1]
            wt[i - 1, j - 1] = s + theta * ss[i - 1, j - 1]

    # Phase 2: Cholesky factorization of the upper triangle.
    try:
        T_block = wt[:col, :col].copy()
        # scipy.linalg.cholesky needs the symmetric matrix; we have only the upper,
        # so mirror to lower for the call (or use lower=False which uses upper only).
        chol = cholesky(T_block, lower=False, check_finite=True)
    except LinAlgError:
        return -3

    # Write factor back into upper triangle. Lower stays whatever it was.
    for i in range(col):
        for j in range(i, col):
            wt[i, j] = chol[i, j]

    return 0
