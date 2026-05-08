"""freev — free/active partition at the generalized Cauchy point.

Reference implementation following ``docs/spec/subroutines/freev.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FreevResult:
    """Outputs of ``freev`` (besides the in-place ``index``)."""

    nfree: int          # new number of free variables
    nenter: int         # number that entered the free set this iter
    ileave: int         # sentinel: indx2[ileave..n] are leaving variables
    wrk: bool           # rebuild-WN flag


def freev(
    n: int,
    nfree_in: int,
    index: np.ndarray,
    indx2: np.ndarray,
    iwhere: np.ndarray,
    updatd: bool,
    cnstnd: bool,
    iter: int,
    iprint: int = -1,
) -> FreevResult:
    """Build the free/active partition; detect entering/leaving variables.

    Parameters
    ----------
    n : int
    nfree_in : int
        Previous iteration's free-variable count (used only when
        ``iter > 0`` and ``cnstnd``).
    index : (n,) int
        On entry: previous partition (when ``iter > 0`` and ``cnstnd``).
        On exit: new partition built from ``iwhere``.
    indx2 : (n,) int
        Output buffer for the change record. Contents on entry are
        overwritten in the slots actually used; remaining slots may
        keep their entry values.
    iwhere : (n,) int
        Per-variable bound status from ``cauchy``: ``<= 0`` is free,
        ``> 0`` is at a bound.
    updatd, cnstnd : bool
    iter : int
    iprint : int
        Diagnostic verbosity (the reference impl does not emit prints).
    """
    # Phase 1: change detection (only when iter>0 and cnstnd).
    nenter = 0
    ileave = n + 1
    if iter > 0 and cnstnd:
        for i in range(nfree_in):
            k = index[i]
            if iwhere[k - 1] > 0:                  # k is 1-based
                ileave -= 1
                indx2[ileave - 1] = k
        for i in range(nfree_in, n):
            k = index[i]
            if iwhere[k - 1] <= 0:
                nenter += 1
                indx2[nenter - 1] = k

    # Phase 2: work flag.
    wrk = (ileave < n + 1) or (nenter > 0) or updatd

    # Phase 3: rebuild partition from iwhere.
    nfree = 0
    iact = n + 1
    for i in range(1, n + 1):                       # 1-based scan
        if iwhere[i - 1] <= 0:
            nfree += 1
            index[nfree - 1] = i
        else:
            iact -= 1
            index[iact - 1] = i

    _ = iprint  # diagnostics not emitted by reference impl

    return FreevResult(nfree=nfree, nenter=nenter, ileave=ileave, wrk=wrk)
