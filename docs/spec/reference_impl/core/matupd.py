"""matupd -- update L-BFGS history (S, Y, sy, ss, theta).

Reference implementation following ``docs/spec/subroutines/matupd.md``.
Stateful: mutates ws, wy, sy, ss, and updates head, itail, col, theta
via mutable container.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MatupdState:
    """Mutable state passed in/out of matupd."""

    itail: int
    col: int
    head: int
    theta: float


def matupd(
    n: int,
    m: int,
    ws: np.ndarray,
    wy: np.ndarray,
    sy: np.ndarray,
    ss: np.ndarray,
    d: np.ndarray,
    r: np.ndarray,
    iupdat: int,
    state: MatupdState,
    rr: float,
    dr: float,
    stp: float,
    dtd: float,
) -> None:
    """Update the L-BFGS limited-memory history with a new (s, y) pair.

    See ``docs/spec/subroutines/matupd.md`` for the contract.
    """
    # Phase 1: slot pointers.
    if iupdat <= m:
        state.col = iupdat
        state.itail = ((state.head + iupdat - 2) % m) + 1
    else:
        state.itail = (state.itail % m) + 1
        state.head = (state.head % m) + 1

    # Phase 2: copy d -> ws[:, itail], r -> wy[:, itail]  (1-based itail).
    ws[:, state.itail - 1] = d
    wy[:, state.itail - 1] = r

    # Phase 3: theta = y'y / s'y.
    state.theta = rr / dr

    # Phase 4: if the buffer overflowed, shift sy and ss.
    if iupdat > m:
        for j in range(1, state.col):                # j = 1..col-1
            # ss[1..j, j] = ss[2..j+1, j+1]
            for k in range(1, j + 1):
                ss[k - 1, j - 1] = ss[k, j]
            # sy[j..col-1, j] = sy[j+1..col, j+1]
            for k in range(j, state.col):
                sy[k - 1, j - 1] = sy[k, j]

    # Phase 5: fill new row of sy and new column of ss.
    pointr = state.head
    for j in range(1, state.col):                    # j = 1..col-1
        # ddot(d, wy[:, pointr]) and ddot(ws[:, pointr], d)
        sy[state.col - 1, j - 1] = float(np.dot(d, wy[:, pointr - 1]))
        ss[j - 1, state.col - 1] = float(np.dot(ws[:, pointr - 1], d))
        pointr = (pointr % m) + 1

    # Phase 6: diagonal entries.
    if stp == 1.0:
        ss[state.col - 1, state.col - 1] = dtd
    else:
        ss[state.col - 1, state.col - 1] = stp * stp * dtd
    sy[state.col - 1, state.col - 1] = dr
