#!/usr/bin/env python3
"""Branch coverage gate for L-BFGS-B src/.

Runs gcovr on the coverage data left behind by previously-executed test
binaries, parses the JSON, and reports any branches in src/*.f that are not
covered (or covered only in one direction).

Phase A behaviour (warn-only):
    Always exits 0. Prints a summary and lists uncovered branches so the
    maintainer can see the state of the suite.

Phase F behaviour (--fail):
    Exits 1 if any non-exempted branch is uncovered.

Files listed via --exclude-file are removed from the must-cover set entirely.
The plan excludes timer.f and prn{1,2,3}lb.f because they are not core
numerics.

Per-line/per-branch exemptions live in COVERAGE_EXCLUSIONS.md. Each entry
documents *why* a particular branch cannot or should not be covered. Format:

    ## src/<file>.f
    - line <N>, branch <0|1|both>: <one-paragraph justification>

This script parses that simple format. Anything not matching is treated as
prose and ignored.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


EXCLUSION_FILE_HEADER = re.compile(r"^##\s+src/(?P<file>[^\s]+\.f)\s*$")
EXCLUSION_LINE = re.compile(
    r"^-\s+line\s+(?P<line>\d+)"
    r"(?:,\s*branch\s+(?P<branch>0|1|both))?\s*[:.]"
)


def parse_exclusions(path: Path) -> dict[str, set[tuple[int, str]]]:
    """Return {filename: {(lineno, "0" | "1" | "both")}}.

    "both" means the line as a whole is exempted (both directions allowed
    uncovered). "0" / "1" exempt only that branch direction.
    """
    out: dict[str, set[tuple[int, str]]] = {}
    if not path.exists():
        return out
    current: str | None = None
    for raw in path.read_text().splitlines():
        m = EXCLUSION_FILE_HEADER.match(raw)
        if m:
            current = m.group("file")
            out.setdefault(current, set())
            continue
        if current is None:
            continue
        m = EXCLUSION_LINE.match(raw)
        if m:
            line = int(m.group("line"))
            branch = m.group("branch") or "both"
            out[current].add((line, branch))
    return out


def run_gcovr(root: Path, build: Path) -> dict | None:
    gcovr = shutil.which("gcovr")
    if gcovr is None:
        print(
            "gcovr not on PATH; install with `pip install --user gcovr` "
            "to enable the coverage gate.",
            file=sys.stderr,
        )
        return None
    out_json = build / "coverage" / "coverage.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gcovr,
        "--root", str(root),
        "--filter", str(root / "src") + "/",
        "--json", str(out_json),
        "--print-summary",
    ]
    try:
        subprocess.run(cmd, check=True, cwd=str(build))
    except subprocess.CalledProcessError as exc:
        print(f"gcovr failed: {exc}", file=sys.stderr)
        return None
    if not out_json.exists():
        return None
    with out_json.open() as f:
        return json.load(f)


def collect_uncovered(
    cov: dict,
    excluded_files: set[str],
    line_exclusions: dict[str, set[tuple[int, str]]],
) -> list[tuple[str, int, str]]:
    """Return list of (filename, lineno, branch_dir) for uncovered branches."""
    result: list[tuple[str, int, str]] = []
    for fileinfo in cov.get("files", []):
        path = fileinfo.get("file", "")
        name = Path(path).name
        if name in excluded_files:
            continue
        excl = line_exclusions.get(name, set())
        for line in fileinfo.get("lines", []):
            lineno = line.get("line_number")
            for branch in line.get("branches", []):
                count = branch.get("count", 0)
                if count > 0:
                    continue
                # Branches are listed as a flat per-line array; the index
                # within that array is the branch direction (0, 1, ...).
                # gcovr.json reports each branch separately with no explicit
                # direction label. Treat any uncovered entry as a gap.
                if (lineno, "both") in excl:
                    continue
                # Best-effort: mark direction by position; gcovr does not
                # always provide it. We default to "both" when ambiguous.
                direction = "both"
                if (lineno, direction) in excl:
                    continue
                result.append((name, lineno, direction))
    return result


def summary(cov: dict) -> str:
    total_b = total_b_covered = 0
    total_l = total_l_covered = 0
    for fileinfo in cov.get("files", []):
        for line in fileinfo.get("lines", []):
            total_l += 1
            if line.get("count", 0) > 0:
                total_l_covered += 1
            for branch in line.get("branches", []):
                total_b += 1
                if branch.get("count", 0) > 0:
                    total_b_covered += 1
    pct = lambda n, d: (100.0 * n / d) if d else 0.0
    return (
        f"Lines:    {total_l_covered}/{total_l} ({pct(total_l_covered, total_l):.1f}%)\n"
        f"Branches: {total_b_covered}/{total_b} ({pct(total_b_covered, total_b):.1f}%)"
    )


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True,
                   help="repo root (contains src/)")
    p.add_argument("--build", type=Path, required=True,
                   help="build directory (where .gcda files live)")
    p.add_argument("--exclusions", type=Path, required=True,
                   help="path to COVERAGE_EXCLUSIONS.md")
    p.add_argument("--exclude-file", action="append", default=[],
                   help="filename (basename) to exclude from must-cover set")
    p.add_argument("--warn-only", action="store_true",
                   help="exit 0 even if uncovered branches remain (default off)")
    p.add_argument("--fail", action="store_true",
                   help="exit non-zero on any uncovered, non-exempt branch")
    p.add_argument("--line-threshold", type=float, default=None,
                   help="minimum acceptable line-coverage fraction over the "
                        "in-scope files (e.g. 0.85). Used as a regression "
                        "detector: PR drops coverage below this -> fail.")
    args = p.parse_args(argv)

    if args.warn_only and args.fail:
        print("--warn-only and --fail are mutually exclusive", file=sys.stderr)
        return 2

    cov = run_gcovr(args.root, args.build)
    if cov is None:
        # No data — common in Phase A before any tests run.
        if args.fail:
            print("FAIL: no coverage data available")
            return 1
        print("WARN: no coverage data available; gate vacuously passes")
        return 0

    print(summary(cov))

    excluded_files = set(args.exclude_file)
    line_exclusions = parse_exclusions(args.exclusions)

    if args.line_threshold is not None:
        in_scope_total = in_scope_covered = 0
        for fileinfo in cov.get("files", []):
            name = Path(fileinfo.get("file", "")).name
            if name in excluded_files:
                continue
            for line in fileinfo.get("lines", []):
                in_scope_total += 1
                if line.get("count", 0) > 0:
                    in_scope_covered += 1
        frac = (in_scope_covered / in_scope_total) if in_scope_total else 1.0
        print(f"In-scope line coverage: "
              f"{in_scope_covered}/{in_scope_total} ({100*frac:.1f}%) "
              f"vs threshold {100*args.line_threshold:.1f}%")
        if frac < args.line_threshold:
            print(f"FAIL: line coverage {100*frac:.1f}% below "
                  f"threshold {100*args.line_threshold:.1f}%")
            return 1

    uncovered = collect_uncovered(cov, excluded_files, line_exclusions)

    if not uncovered:
        print("Branch coverage: 100% (net of exemptions). PASS.")
        return 0

    print(f"\n{len(uncovered)} uncovered branch(es) in src/ "
          f"(after excluding {sorted(excluded_files)}):")
    for fname, lineno, direction in uncovered[:50]:
        print(f"  src/{fname}:{lineno} [branch {direction}]")
    if len(uncovered) > 50:
        print(f"  ... and {len(uncovered) - 50} more")

    if args.warn_only:
        print("\n(warn-only mode: gate passes)")
        return 0
    if args.fail:
        return 1
    # Default behaviour: warn but pass, until Phase F flips it.
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
