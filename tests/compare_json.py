#!/usr/bin/env python3
"""Compare two L-BFGS-B JSON outputs using a relative-absolute closeness metric.

Metric: |a - b| / (1 + |b|) < THRESHOLD applied pointwise to every float.
Integers, strings and bools must match exactly.

Exits 0 on match, 1 on mismatch (with a diff report).
"""
import json
import sys

THRESHOLD = 1e-15
MAX_REPORTED_DIFFS = 30


def is_close_rel_abs(a, b):
    return abs(float(a) - float(b)) / (1.0 + abs(float(b))) < THRESHOLD


def walk(ref, act, path=""):
    if isinstance(ref, dict):
        if not isinstance(act, dict) or set(ref) != set(act):
            yield (path, "structure", ref, act)
            return
        for k in ref:
            yield from walk(ref[k], act[k], f"{path}.{k}")
    elif isinstance(ref, list):
        if not isinstance(act, list) or len(ref) != len(act):
            ref_len = len(ref) if isinstance(ref, list) else None
            act_len = len(act) if isinstance(act, list) else None
            yield (path, "length", ref_len, act_len)
            return
        for i, (r, a) in enumerate(zip(ref, act)):
            yield from walk(r, a, f"{path}[{i}]")
    elif isinstance(ref, bool):
        if act is not ref:
            yield (path, "exact", ref, act)
    elif isinstance(ref, float):
        if not isinstance(act, (int, float)) or not is_close_rel_abs(act, ref):
            yield (path, "numeric", ref, act)
    else:
        if act != ref:
            yield (path, "exact", ref, act)


def main():
    if len(sys.argv) != 3:
        print("usage: compare_json.py <expected.json> <actual.json>", file=sys.stderr)
        return 2

    with open(sys.argv[1]) as f:
        ref = json.load(f)
    with open(sys.argv[2]) as f:
        act = json.load(f)

    diffs = list(walk(ref, act))
    if not diffs:
        print(f"OK: {sys.argv[2]} matches {sys.argv[1]} within tolerance {THRESHOLD:g}")
        return 0

    print(f"FAIL: {len(diffs)} difference(s) between expected and actual:")
    for path, kind, r, a in diffs[:MAX_REPORTED_DIFFS]:
        print(f"  [{kind}] {path or '<root>'}: expected={r!r} actual={a!r}")
    if len(diffs) > MAX_REPORTED_DIFFS:
        print(f"  ... and {len(diffs) - MAX_REPORTED_DIFFS} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())
