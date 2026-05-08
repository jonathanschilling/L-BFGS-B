"""projgr -- infinity norm of the projected gradient.

Reference implementation following ``docs/spec/subroutines/projgr.md``.
This module is part of the portability spec pack and is intended as a
language-agnostic reference, not as a production-quality optimizer.
"""

from __future__ import annotations

import numpy as np


def projgr(
    n: int,
    x: np.ndarray,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    g: np.ndarray,
) -> float:
    """Compute ``sbgnrm = ||g_proj||_inf`` (the projected gradient norm).

    Parameters
    ----------
    n : int
        Number of variables (must equal the length of ``x``, ``l``, ``u``,
        ``nbd``, ``g``).
    x : (n,) float64
        Current iterate. Must satisfy ``l[i] <= x[i] <= u[i]`` where the
        corresponding bound is active per ``nbd[i]``.
    l, u : (n,) float64
        Lower and upper bounds. Entries with no active bound (``nbd[i] = 0``,
        or ``nbd[i] = 3`` for ``l``, or ``nbd[i] = 1`` for ``u``) are
        ignored.
    nbd : (n,) int
        Bound types per variable: 0 = unbounded, 1 = lower only, 2 = both,
        3 = upper only.
    g : (n,) float64
        Gradient ``grad f(x)``.

    Returns
    -------
    float
        ``max_i |gi|`` where ``gi`` is the projected-gradient proxy
        defined in the spec.
    """
    sbgnrm = 0.0
    for i in range(n):
        gi = g[i]
        if nbd[i] != 0:
            if gi < 0.0:
                if nbd[i] >= 2:
                    gi = max(x[i] - u[i], gi)
            else:
                if nbd[i] <= 2:
                    gi = min(x[i] - l[i], gi)
        sbgnrm = max(sbgnrm, abs(gi))
    return sbgnrm
