"""errclb -- input validation.

Reference implementation following ``docs/spec/subroutines/errclb.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class ErrclbState:
    """Mutable state passed in/out of errclb (as in the F77 reverse-comm)."""

    task: str
    info: int
    k: int


def errclb(
    n: int,
    m: int,
    factr: float,
    l: np.ndarray,
    u: np.ndarray,
    nbd: np.ndarray,
    state: ErrclbState,
) -> None:
    """Validate optimizer inputs; mutate ``state`` on error.

    Follows the F77 sequential-overwrite semantics: later writes to
    ``state.task`` overwrite earlier ones, so multi-error inputs
    yield the *last* error's message in ``state.task``.

    Notes
    -----
    On success ``state`` is **not modified** (matching F77).
    """
    if n <= 0:
        state.task = "ERROR: N .LE. 0"
    if m <= 0:
        state.task = "ERROR: M .LE. 0"
    if factr < 0.0:
        state.task = "ERROR: FACTR .LT. 0"

    for i in range(n):              # 0-based here; i+1 is the F77 1-based index
        if nbd[i] < 0 or nbd[i] > 3:
            state.task = "ERROR: INVALID NBD"
            state.info = -6
            state.k = i + 1
        if nbd[i] == 2:
            if l[i] > u[i]:
                state.task = "ERROR: NO FEASIBLE SOLUTION"
                state.info = -7
                state.k = i + 1
