"""hpsolb -- partial heapsort: extract one minimum, leave heap.

Reference implementation following ``docs/spec/subroutines/hpsolb.md``.
Mirrors the F77 sift-up insertion order and sift-down tie-breaking so
that ``--strict`` conformance holds on tied minima.
"""

from __future__ import annotations

import numpy as np


def hpsolb(n: int, t: np.ndarray, iorder: np.ndarray, iheap: int) -> None:
    """Extract minimum to position ``n``; leave ``t[1..n-1]`` as a min-heap.

    Uses 1-based indexing in the algorithm description but Python
    conventions (0-based) in the array indexing. Two notations are
    related by ``t[i_one_based] == t_array[i_one_based - 1]``.
    """
    # Phase 1: build heap by sift-up insertion (skip if already a heap).
    if iheap == 0:
        for k in range(2, n + 1):                  # 1-based k = 2..n
            ddum = t[k - 1]
            indxin = iorder[k - 1]
            i = k
            while i > 1:
                j = i // 2                         # parent (integer divide truncates toward zero)
                if ddum < t[j - 1]:
                    t[i - 1] = t[j - 1]
                    iorder[i - 1] = iorder[j - 1]
                    i = j
                else:
                    break
            t[i - 1] = ddum
            iorder[i - 1] = indxin

    # Phase 2: extract root (only if n > 1).
    if n > 1:
        out = t[0]                                 # current minimum
        indxou = iorder[0]
        ddum = t[n - 1]                            # value to sift down
        indxin = iorder[n - 1]
        i = 1
        while True:
            j = 2 * i                              # left child
            if j > n - 1:
                break
            if t[(j + 1) - 1] < t[j - 1]:          # right child smaller -> use it
                j = j + 1
            if t[j - 1] < ddum:
                t[i - 1] = t[j - 1]
                iorder[i - 1] = iorder[j - 1]
                i = j
            else:
                break
        t[i - 1] = ddum
        iorder[i - 1] = indxin
        # Place extracted minimum at position n.
        t[n - 1] = out
        iorder[n - 1] = indxou
