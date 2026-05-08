#!/usr/bin/env python3
"""Compare two L-BFGS-B aggregate JSON outputs.

Each JSON has four fields:
    final_task   (string,         exact match)
    final_f      (float,          IsCloseRelAbs metric)
    final_projg  (float,          IsCloseRelAbs metric)
    final_x      (list of float,  IsCloseRelAbs pointwise)

The IsCloseRelAbs metric is |a-b| / (1 + |b|) < threshold. The "+1"
gracefully degrades to absolute tolerance when |b| is small (which is
typical for converged f or |proj g|, both near machine epsilon).

Two thresholds, because final_x has very different cross-toolchain
behaviour than the scalar fields:

- THRESHOLD = 1e-9     for final_f, final_projg, final_task
- X_THRESHOLD = 2e-2   for final_x

The scalar fields are typically near machine epsilon at convergence, so
1e-9 leaves huge IsCloseRelAbs headroom. final_x is a different story:
the Rosenbrock optimum requires x(i+1) = x(i)^2 at every step, so tiny
FP rounding in early components is *squared* on each step along the
chain. The bulk of final_x is essentially bit-identical (median diff
~1e-14, p95 ~1e-12 in driver3) but the last few chain components blow
up to ~1% relative across toolchains. 2e-2 (2x margin over the worst
observed ~1e-2) tolerates this while still catching a regression that
puts the algorithm on the wrong stationary point - any other critical
point would differ by 10%+.

Override with LBFGSB_TEST_TOLERANCE / LBFGSB_TEST_X_TOLERANCE.

Exits 0 on match (and prints the values for visibility), 1 on mismatch.
"""
import json
import os
import sys

DEFAULT_THRESHOLD = 1e-9
DEFAULT_X_THRESHOLD = 2e-2
THRESHOLD = float(os.environ.get("LBFGSB_TEST_TOLERANCE", DEFAULT_THRESHOLD))
X_THRESHOLD = float(os.environ.get("LBFGSB_TEST_X_TOLERANCE", DEFAULT_X_THRESHOLD))

EXPECTED_KEYS = {"final_task", "final_f", "final_projg", "final_x"}
MAX_X_DIFFS_REPORTED = 10


def is_close_rel_abs(a, b, threshold):
    return abs(float(a) - float(b)) / (1.0 + abs(float(b))) < threshold


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
        if not is_close_rel_abs(act[key], ref[key], THRESHOLD):
            diffs.append(
                f"{key} differs (threshold {THRESHOLD:g}):\n"
                f"  ref: {ref[key]!r}\n"
                f"  act: {act[key]!r}\n"
                f"  |a-b|/(1+|b|) = {abs(act[key]-ref[key])/(1+abs(ref[key])):g}"
            )

    rx, ax = ref["final_x"], act["final_x"]
    if len(rx) != len(ax):
        diffs.append(f"final_x length differs: ref={len(rx)} act={len(ax)}")
    else:
        x_diffs = [
            (i, rx[i], ax[i], abs(ax[i] - rx[i]) / (1.0 + abs(rx[i])))
            for i in range(len(rx))
            if not is_close_rel_abs(ax[i], rx[i], X_THRESHOLD)
        ]
        if x_diffs:
            lines = [f"final_x: {len(x_diffs)}/{len(rx)} components differ "
                     f"(threshold {X_THRESHOLD:g})"]
            for i, r, a, m in x_diffs[:MAX_X_DIFFS_REPORTED]:
                lines.append(f"  [{i}] ref={r!r} act={a!r} |a-b|/(1+|b|)={m:g}")
            if len(x_diffs) > MAX_X_DIFFS_REPORTED:
                lines.append(f"  ... and {len(x_diffs) - MAX_X_DIFFS_REPORTED} more")
            diffs.append("\n".join(lines))

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
        f"OK (scalar tol {THRESHOLD:g}, x tol {X_THRESHOLD:g}): "
        f"task={act['final_task']!r}, "
        f"f={act['final_f']:.4e}, "
        f"projg={act['final_projg']:.4e}, "
        f"|x|={len(act['final_x'])}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
