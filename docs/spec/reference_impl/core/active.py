"""active -- initial projection and bound-status setup.

Reference implementation following ``docs/spec/subroutines/active.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ActiveResult:
    """Output flags from ``active``. ``x`` and ``iwhere`` are mutated in place."""

    prjctd: bool
    cnstnd: bool
    boxed: bool


def active(
    n: int,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    x: np.ndarray,
    iwhere: np.ndarray,
    iprint: int = -1,
) -> ActiveResult:
    """Project ``x`` to the feasible set and initialize ``iwhere``.

    Parameters
    ----------
    n : int
        Number of variables.
    l, u : (n,) float64
        Lower / upper bounds (consulted only where ``nbd`` enables them).
    nbd : (n,) int
        Bound types: 0, 1, 2, 3.
    x : (n,) float64
        Initial point. **Mutated in place** to the projected feasible point.
    iwhere : (n,) int
        Output buffer (any contents on entry are overwritten).
    iprint : int
        Diagnostic verbosity (informational; the reference impl does not
        emit diagnostics, since output format is port-defined).

    Returns
    -------
    ActiveResult
        Flags ``prjctd``, ``cnstnd``, ``boxed``.
    """
    # Phase 1: project x and count variables already at a bound.
    nbdd = 0
    prjctd = False
    for i in range(n):
        if nbd[i] > 0:
            if nbd[i] <= 2 and x[i] <= l[i]:
                if x[i] < l[i]:
                    prjctd = True
                    x[i] = l[i]
                nbdd += 1
            elif nbd[i] >= 2 and x[i] >= u[i]:
                if x[i] > u[i]:
                    prjctd = True
                    x[i] = u[i]
                nbdd += 1

    # Phase 2: classify each variable and summarise.
    cnstnd = False
    boxed = True
    for i in range(n):
        if nbd[i] != 2:
            boxed = False
        if nbd[i] == 0:
            iwhere[i] = -1
        else:
            cnstnd = True
            if nbd[i] == 2 and u[i] - l[i] <= 0.0:
                iwhere[i] = 3
            else:
                iwhere[i] = 0

    # Diagnostic output is port-defined; reference impl does not emit.
    _ = iprint
    _ = nbdd

    return ActiveResult(prjctd=prjctd, cnstnd=cnstnd, boxed=boxed)
