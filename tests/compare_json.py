#!/usr/bin/env python3
"""Compare two L-BFGS-B aggregate JSON outputs.

Each JSON has three fields:
    final_task   (string,  exact match)
    final_f      (float,   IsCloseRelAbs metric)
    final_projg  (float,   IsCloseRelAbs metric)

The IsCloseRelAbs metric is |a-b| / (1 + |b|) < threshold. The "+1"
gracefully degrades to absolute tolerance when |b| is small (which is
typical for converged f or |proj g|, both near machine epsilon).

Default threshold is 1e-2 (1%): tight enough to catch real algorithmic
regressions, loose enough to absorb cross-toolchain FP-trajectory
divergence in iterative methods. Override with LBFGSB_TEST_TOLERANCE.

Exits 0 on match (and prints the values for visibility), 1 on mismatch.
"""
import json
import os
import sys

DEFAULT_THRESHOLD = 1e-2
THRESHOLD = float(os.environ.get("LBFGSB_TEST_TOLERANCE", DEFAULT_THRESHOLD))

EXPECTED_KEYS = {"final_task", "final_f", "final_projg"}


def is_close_rel_abs(a, b):
    return abs(float(a) - float(b)) / (1.0 + abs(float(b))) < THRESHOLD


def compare(ref, act):
    diffs = []
    if set(ref) != EXPECTED_KEYS or set(act) != EXPECTED_KEYS:
        diffs.append(
            f"structure mismatch: expected keys {sorted(EXPECTED_KEYS)}, "
            f"ref={sorted(ref)}, act={sorted(act)}"
        )
        return diffs

    if ref["final_task"] != act["final_task"]:
        diffs.append(
            f"final_task differs:\n"
            f"  ref: {ref['final_task']!r}\n"
            f"  act: {act['final_task']!r}"
        )

    for key in ("final_f", "final_projg"):
        if not is_close_rel_abs(act[key], ref[key]):
            diffs.append(
                f"{key} differs (threshold {THRESHOLD:g}):\n"
                f"  ref: {ref[key]!r}\n"
                f"  act: {act[key]!r}\n"
                f"  |a-b|/(1+|b|) = {abs(act[key]-ref[key])/(1+abs(ref[key])):g}"
            )

    return diffs


def main():
    if len(sys.argv) != 3:
        print("usage: compare_json.py <expected.json> <actual.json>", file=sys.stderr)
        return 2

    with open(sys.argv[1]) as f:
        ref = json.load(f)
    with open(sys.argv[2]) as f:
        act = json.load(f)

    diffs = compare(ref, act)
    if diffs:
        print(f"FAIL: {sys.argv[2]} differs from {sys.argv[1]}")
        for d in diffs:
            print(d)
        return 1

    print(
        f"OK (threshold {THRESHOLD:g}): "
        f"task={act['final_task']!r}, "
        f"f={act['final_f']:.4e}, "
        f"projg={act['final_projg']:.4e}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
