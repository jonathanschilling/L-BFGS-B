#!/usr/bin/env python3
"""Check a driver's JSON output against mathematically-derived bounds.

We do not record bit-level golden outputs because the L-BFGS-B
trajectory is sensitive to FP rounding (gfortran version, BLAS
implementation, thread count, CPU dispatch) and goldens captured on
one machine fail on another. Instead, each driver has a small spec
derived from its source code's stopping criteria and bound setup
(drivers/driver*.f / .f90). The test verifies the algorithm satisfied
its own contract:

    * final_task starts with the expected stopping reason
    * final_f and final_projg are within the bounds implied by the
      stopping criterion (with comfortable headroom)
    * final_x is feasible (every component within its declared bounds)

This catches:
    * algorithm fails to converge (wrong final_task, or values blow up)
    * algorithm converges to a wrong stationary point (final_f stays
      large; the unique critical point of bound-Rosenbrock is all-1's
      which gives f=0)
    * algorithm produces NaN / out-of-bounds x

It does not catch subtle convergence-quality regressions (e.g. an
extra iteration here or there). For this Fortran 77 library that's
the right trade.

Usage:  check_output.py <driver_name> <output.json>
        driver_name is one of {driver1, driver2, driver3}_{f77, f90}
"""

import json
import sys


# Bounds setup is identical across all three example drivers:
#   1-based odd i:  nbd=2, l=1,    u=100   (i.e. 0-based even index)
#   1-based even i: nbd=2, l=-100, u=100   (i.e. 0-based odd index)
def _x_bounds(i_zero_based):
    if i_zero_based % 2 == 0:
        return 1.0, 100.0
    return -100.0, 100.0


# Per-driver spec (the F77 and F90 versions of each driver call the
# library identically, so they share a spec).
#
# driver1: factr=1e7, pgtol=1e-5  -> stops via the L-BFGS-B internal
#   criterion (REL_REDUCTION_OF_F or NORM_OF_PROJECTED_GRADIENT).
#   Either way the resulting final_task starts with "CONVERGENCE: ".
#   With factr=1e7 the f-reduction criterion is the typical winner.
#
# driver2/driver3: factr=0, pgtol=0 (internal criteria suppressed) so
#   the only stopping condition that fires is the user's custom check
#   |proj g| <= 1e-10 * (1 + |f|), provided the test environment also
#   neutralises driver2's nfg cap (LBFGSB_NFG_LIMIT) and driver3's
#   wallclock cap (LBFGSB_TLIMIT). With those overrides, final_task is
#   exactly the corresponding driver-set string.
DRIVER_SPECS = {
    "driver1": {
        "n":           25,
        "task_prefix": "CONVERGENCE: ",
        "f_max":       1e-6,    # observed ~1e-9; ~1000x headroom
        "projg_max":   1e-2,    # not directly bounded by stopping; sanity
    },
    "driver2": {
        "n":           25,
        "task_prefix": "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL",
        "f_max":       1e-6,    # observed ~1e-15; algorithm well-converged
        "projg_max":   1e-9,    # criterion is < 1e-10*(1+|f|); 10x margin
    },
    "driver3": {
        "n":           1000,
        "task_prefix": "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL",
        "f_max":       1e-6,
        "projg_max":   1e-9,
    },
}


def canonical(driver_name):
    for suffix in ("_f77", "_f90"):
        if driver_name.endswith(suffix):
            return driver_name[: -len(suffix)]
    return driver_name


def main():
    if len(sys.argv) != 3:
        print("usage: check_output.py <driver_name> <output.json>", file=sys.stderr)
        return 2

    driver_name = sys.argv[1]
    spec_key = canonical(driver_name)
    if spec_key not in DRIVER_SPECS:
        print(f"FAIL: no spec for driver {driver_name!r}", file=sys.stderr)
        return 1
    spec = DRIVER_SPECS[spec_key]

    with open(sys.argv[2]) as f:
        out = json.load(f)

    failures = []

    if not out["final_task"].startswith(spec["task_prefix"]):
        failures.append(
            f"final_task: expected to start with {spec['task_prefix']!r}, "
            f"got {out['final_task']!r}"
        )

    f_val = out["final_f"]
    if not (0.0 <= f_val <= spec["f_max"]):
        failures.append(
            f"final_f: expected 0 <= f <= {spec['f_max']:g}, got {f_val!r}"
        )

    g_val = out["final_projg"]
    if not (0.0 <= g_val <= spec["projg_max"]):
        failures.append(
            f"final_projg: expected 0 <= |proj g| <= {spec['projg_max']:g}, "
            f"got {g_val!r}"
        )

    x = out["final_x"]
    n = spec["n"]
    if len(x) != n:
        failures.append(f"final_x: expected length {n}, got {len(x)}")
    else:
        oob = []
        for i, xi in enumerate(x):
            lo, hi = _x_bounds(i)
            if not (lo <= xi <= hi):
                oob.append((i, xi, lo, hi))
        if oob:
            lines = [f"final_x: {len(oob)} components outside their declared bounds"]
            for i, xi, lo, hi in oob[:5]:
                lines.append(f"  [{i}] expected [{lo}, {hi}], got {xi!r}")
            if len(oob) > 5:
                lines.append(f"  ... and {len(oob) - 5} more")
            failures.append("\n".join(lines))

    if failures:
        print(f"FAIL: {sys.argv[2]} ({driver_name})")
        for fline in failures:
            print(f"  {fline}")
        return 1

    max_dev = max(abs(xi - 1.0) for xi in x)
    print(
        f"OK ({driver_name}): "
        f"task={out['final_task'][:40]!r}, "
        f"f={f_val:.3e}, "
        f"projg={g_val:.3e}, "
        f"max|x-1|={max_dev:.3e}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
