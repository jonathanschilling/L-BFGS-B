"""formk -- form the LEL' factorization of the K matrix for subsm.

Reference implementation following ``docs/spec/subroutines/formk.md``.
Literal port of ``src/formk.f`` to preserve the cyclic-buffer indexing
and the active-set-drift modifications in their exact F77 order.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, solve_triangular, LinAlgError


def formk(
    n: int,
    nsub: int,
    ind: np.ndarray,
    nenter: int,
    ileave: int,
    indx2: np.ndarray,
    iupdat: int,
    updatd: bool,
    wn: np.ndarray,
    wn1: np.ndarray,
    m: int,
    ws: np.ndarray,
    wy: np.ndarray,
    sy: np.ndarray,
    theta: float,
    col: int,
    head: int,
) -> int:
    """Form K = [-D - Y'ZZ'Y/theta, L_a'-R_z'; L_a-R_z, theta*S'AA'S],
    then Cholesky-factor in two stages. Mutates wn (output) and wn1 (cached).
    Returns info = 0 / -1 / -2.
    """
    # Phase 1a: shift wn1 if iupdat>m and updatd.
    if updatd:
        if iupdat > m:
            for jy in range(1, m):                       # jy = 1..m-1
                js = m + jy
                # wn1[jy..m-1 (1-based)] gets shifted from wn1[jy+1..m]
                length = m - jy
                for k in range(length):
                    wn1[jy - 1 + k, jy - 1] = wn1[jy + k, jy]
                    wn1[js - 1 + k, js - 1] = wn1[js + k, js]
                # (2,1) block shift, length m-1
                for k in range(m - 1):
                    wn1[m + k, jy - 1] = wn1[m + 1 + k, jy]

        # Phase 1b: fill new row col of (1,1), (2,2), and entry (2,1) of new col.
        pbegin = 1
        pend = nsub
        dbegin = nsub + 1
        dend = n
        iy = col
        is_ = m + col
        ipntr = head + col - 1
        if ipntr > m:
            ipntr -= m
        jpntr = head
        for jy in range(1, col + 1):                    # jy = 1..col
            js = m + jy
            temp1 = 0.0
            temp2 = 0.0
            temp3 = 0.0
            for k in range(pbegin, pend + 1):
                k1 = ind[k - 1]
                temp1 += wy[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
            for k in range(dbegin, dend + 1):
                k1 = ind[k - 1]
                temp2 += ws[k1 - 1, ipntr - 1] * ws[k1 - 1, jpntr - 1]
                temp3 += ws[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
            wn1[iy - 1, jy - 1] = temp1
            wn1[is_ - 1, js - 1] = temp2
            wn1[is_ - 1, jy - 1] = temp3
            jpntr = (jpntr % m) + 1

        # New column col of (2,1), R_z part.
        jy = col
        jpntr = head + col - 1
        if jpntr > m:
            jpntr -= m
        ipntr = head
        for i in range(1, col + 1):
            is_i = m + i
            temp3 = 0.0
            for k in range(pbegin, pend + 1):
                k1 = ind[k - 1]
                temp3 += ws[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
            ipntr = (ipntr % m) + 1
            wn1[is_i - 1, jy - 1] = temp3
        upcl = col - 1
    else:
        upcl = col

    # Phase 2: modify old parts in (1,1) and (2,2) blocks for active-set drift.
    ipntr = head
    for iy in range(1, upcl + 1):
        is_iy = m + iy
        jpntr = head
        for jy in range(1, iy + 1):
            js = m + jy
            temp1 = 0.0
            temp2 = 0.0
            temp3 = 0.0
            temp4 = 0.0
            for k in range(1, nenter + 1):
                k1 = indx2[k - 1]
                temp1 += wy[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
                temp2 += ws[k1 - 1, ipntr - 1] * ws[k1 - 1, jpntr - 1]
            for k in range(ileave, n + 1):
                k1 = indx2[k - 1]
                temp3 += wy[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
                temp4 += ws[k1 - 1, ipntr - 1] * ws[k1 - 1, jpntr - 1]
            wn1[iy - 1, jy - 1] += temp1 - temp3
            wn1[is_iy - 1, js - 1] += -temp2 + temp4
            jpntr = (jpntr % m) + 1
        ipntr = (ipntr % m) + 1

    # Modify old parts in (2,1) block.
    ipntr = head
    for is_ in range(m + 1, m + upcl + 1):
        jpntr = head
        for jy in range(1, upcl + 1):
            temp1 = 0.0
            temp3 = 0.0
            for k in range(1, nenter + 1):
                k1 = indx2[k - 1]
                temp1 += ws[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
            for k in range(ileave, n + 1):
                k1 = indx2[k - 1]
                temp3 += ws[k1 - 1, ipntr - 1] * wy[k1 - 1, jpntr - 1]
            if is_ <= jy + m:
                wn1[is_ - 1, jy - 1] += temp1 - temp3
            else:
                wn1[is_ - 1, jy - 1] += -temp1 + temp3
            jpntr = (jpntr % m) + 1
        ipntr = (ipntr % m) + 1

    # Phase 3: form upper triangle of wn from wn1, sy, theta.
    m2 = 2 * m
    for iy in range(1, col + 1):
        is_iy = col + iy
        is1 = m + iy
        for jy in range(1, iy + 1):
            js = col + jy
            js1 = m + jy
            wn[jy - 1, iy - 1] = wn1[iy - 1, jy - 1] / theta
            wn[js - 1, is_iy - 1] = wn1[is1 - 1, js1 - 1] * theta
        for jy in range(1, iy):
            wn[jy - 1, is_iy - 1] = -wn1[is1 - 1, jy - 1]
        for jy in range(iy, col + 1):
            wn[jy - 1, is_iy - 1] = wn1[is1 - 1, jy - 1]
        wn[iy - 1, iy - 1] += sy[iy - 1, iy - 1]

    # First Cholesky on (1,1) block.
    try:
        chol_top = cholesky(wn[:col, :col].copy(), lower=False, check_finite=True)
    except LinAlgError:
        return -1
    for i in range(col):
        for j in range(i, col):
            wn[i, j] = chol_top[i, j]

    # Triangular solve: wn[1:col, col+1:2col] <- L^{-T} (current contents).
    for js in range(col + 1, 2 * col + 1):
        rhs = wn[:col, js - 1].copy()
        wn[:col, js - 1] = solve_triangular(
            wn[:col, :col], rhs, lower=False, trans="T"
        )

    # Form S'AA'S * theta + (L^{-1}(-L_a'+R_z'))' (L^{-1}(...)) in (2,2) block.
    for is_ in range(col + 1, 2 * col + 1):
        for js in range(is_, 2 * col + 1):
            wn[is_ - 1, js - 1] += float(np.dot(wn[:col, is_ - 1], wn[:col, js - 1]))

    # Second Cholesky on (2,2) block.
    try:
        block22 = wn[col : 2 * col, col : 2 * col].copy()
        chol_bot = cholesky(block22, lower=False, check_finite=True)
    except LinAlgError:
        return -2
    for i in range(col):
        for j in range(i, col):
            wn[col + i, col + j] = chol_bot[i, j]

    return 0
