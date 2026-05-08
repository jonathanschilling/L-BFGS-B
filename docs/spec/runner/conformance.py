#!/usr/bin/env python3
"""Conformance test runner for the L-BFGS-B portability spec pack.

Reads JSON test vectors from ``docs/spec/data/`` and runs each one
through the Python reference implementation under
``docs/spec/reference_impl/``. Reports pass / fail with diagnostics
on mismatch.

Modes
-----
``--strict`` (default)
    Compare output values bit-for-bit (IEEE-754 byte equality).
    Requires that the Python reference and the JSON expected values
    were produced under the same BLAS / LAPACK pin (see
    ``docs/spec/04_numerics.md``).

``--tolerance``
    Use per-subroutine numerical tolerances (see
    ``docs/spec/07_conformance.md``). Looser; for ports linked
    against a non-reference BLAS.

Usage
-----
    python3 docs/spec/runner/conformance.py [--strict|--tolerance]
                                            [--filter SUB] [--data DIR]

For ports in other languages: see ``docs/spec/runner/README.md`` for
the expected I/O protocol.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
from typing import Any, Callable

import numpy as np


# Make ``core.<sub>`` importable regardless of cwd.
ROOT = os.path.dirname(os.path.abspath(__file__))
SPEC_ROOT = os.path.dirname(ROOT)
REPO_ROOT = os.path.dirname(os.path.dirname(SPEC_ROOT))
sys.path.insert(0, SPEC_ROOT)
sys.path.insert(0, os.path.join(SPEC_ROOT, "reference_impl"))


# Default location of the F77 conformance driver (built by
# tests/conformance/CMakeLists.txt). Can be overridden via --driver-bin.
DEFAULT_FORTRAN_DRIVER = os.path.join(REPO_ROOT, "build", "tests",
                                      "conformance", "conformance_driver")


# ---------------------------------------------------------------------------
# Per-subroutine tolerances (used in --tolerance mode).
# These match docs/spec/07_conformance.md.
# ---------------------------------------------------------------------------
TOLERANCES = {
    "projgr":  {"abs": 0.0,    "rel": 0.0},          # exact arithmetic
    "errclb":  {"abs": 0.0,    "rel": 0.0},          # integer / string only
    "active":  {"abs": 0.0,    "rel": 0.0},          # exact arithmetic
    "freev":   {"abs": 0.0,    "rel": 0.0},          # integer only
    "hpsolb":  {"abs": 0.0,    "rel": 0.0},          # comparisons + assignments
    "cmprlb":  {"abs": 1e-14,  "rel": 1e-14},        # via bmv
    "bmv":     {"abs": 1e-14,  "rel": 1e-14},        # triangular solves
    "formt":   {"abs": 1e-14,  "rel": 1e-14},        # Cholesky
    "matupd":  {"abs": 0.0,    "rel": 0.0},          # exact (dot products inside BLAS)
    "dcstep":  {"abs": 1e-14,  "rel": 1e-14},        # cubic interpolation
    "formk":   {"abs": 1e-14,  "rel": 1e-14},        # Cholesky + dtrsm
    "dcsrch":  {"abs": 1e-14,  "rel": 1e-14},        # full line-search arithmetic
    "lnsrlb":  {"abs": 1e-14,  "rel": 1e-14},        # min/max / dot / projection
    "cauchy":  {"abs": 1e-13,  "rel": 1e-13},        # accumulated f1, f2 updates
    "subsm":   {"abs": 1e-13,  "rel": 1e-13},        # subspace solve
    "mainlb":  {"abs": 1e-12,  "rel": 1e-12},        # full algorithm loop
}


def _bit_equal(a: float, b: float) -> bool:
    """IEEE-754 byte equality (treats NaN as equal-to-itself)."""
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    return np.float64(a).tobytes() == np.float64(b).tobytes()


def _close(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    """Numerical closeness with absolute + relative tolerance."""
    return abs(a - b) <= abs_tol + rel_tol * max(abs(a), abs(b), 1.0)


def _compare(
    actual: Any,
    expected: Any,
    name: str,
    strict: bool,
    abs_tol: float,
    rel_tol: float,
    issues: list[str],
) -> bool:
    """Recursive structural compare. Appends issues for any mismatch."""
    if isinstance(expected, list):
        if not isinstance(actual, (list, np.ndarray)) or len(actual) != len(expected):
            issues.append(f"{name}: length mismatch ({len(actual)} != {len(expected)})")
            return False
        ok = True
        for i, (a, e) in enumerate(zip(actual, expected)):
            ok = _compare(a, e, f"{name}[{i}]", strict, abs_tol, rel_tol, issues) and ok
        return ok
    if isinstance(expected, bool):
        if bool(actual) != expected:
            issues.append(f"{name}: {actual!r} != {expected!r}")
            return False
        return True
    if isinstance(expected, int) and not isinstance(expected, bool):
        if int(actual) != expected:
            issues.append(f"{name}: {actual} != {expected}")
            return False
        return True
    if isinstance(expected, float):
        a = float(actual)
        e = float(expected)
        if strict:
            ok = _bit_equal(a, e)
        else:
            ok = _close(a, e, abs_tol, rel_tol)
        if not ok:
            issues.append(f"{name}: {a!r} != {e!r} ({'bit' if strict else 'tol'} diff)")
        return ok
    if isinstance(expected, str):
        # Strip trailing whitespace (F77 character buffers are space-padded).
        if str(actual).rstrip() != expected.rstrip():
            issues.append(f"{name}: {actual!r} != {expected!r}")
            return False
        return True
    issues.append(f"{name}: unsupported type {type(expected).__name__}")
    return False


# ---------------------------------------------------------------------------
# Per-subroutine handlers: convert JSON inputs -> Python ref call ->
# return actual outputs as a dict comparable to spec['expected'].
# ---------------------------------------------------------------------------
HANDLERS: dict[str, Callable[[dict], dict]] = {}


def handler(name: str):
    def deco(f):
        HANDLERS[name] = f
        return f
    return deco


@handler("projgr")
def _projgr(inp: dict) -> dict:
    from core.projgr import projgr
    n = inp["n"]
    sb = projgr(
        n,
        np.array(inp["x"], dtype=np.float64),
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        np.array(inp["g"], dtype=np.float64),
    )
    return {"sbgnrm": float(sb)}


@handler("errclb")
def _errclb(inp: dict) -> dict:
    from core.errclb import ErrclbState, errclb
    state = ErrclbState(task=inp["task_in"], info=inp["info_in"], k=inp["k_in"])
    errclb(
        inp["n"], inp["m"], inp["factr"],
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        state,
    )
    return {"task": state.task, "info": state.info, "k": state.k}


@handler("active")
def _active(inp: dict) -> dict:
    from core.active import active
    n = inp["n"]
    x = np.array(inp["x_in"], dtype=np.float64).copy()
    iwhere = np.zeros(n, dtype=np.int32)
    res = active(
        n,
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        x, iwhere, inp["iprint"],
    )
    return {
        "x": x.tolist(),
        "iwhere": iwhere.tolist(),
        "prjctd": res.prjctd, "cnstnd": res.cnstnd, "boxed": res.boxed,
    }


@handler("freev")
def _freev(inp: dict) -> dict:
    from core.freev import freev
    n = inp["n"]
    index = np.array(inp["index_in"], dtype=np.int32).copy()
    indx2 = np.array(inp["indx2_in"], dtype=np.int32).copy()
    res = freev(
        n, inp["nfree_in"], index, indx2,
        np.array(inp["iwhere"], dtype=np.int32),
        inp["updatd"], inp["cnstnd"], inp["iter"], inp["iprint"],
    )
    return {
        "nfree": res.nfree, "index": index.tolist(), "indx2": indx2.tolist(),
        "nenter": res.nenter, "ileave": res.ileave, "wrk": res.wrk,
    }


@handler("hpsolb")
def _hpsolb(inp: dict) -> dict:
    from core.hpsolb import hpsolb
    n = inp["n"]
    t = np.array(inp["t_in"], dtype=np.float64).copy()
    iorder = np.array(inp["iorder_in"], dtype=np.int32).copy()
    hpsolb(n, t, iorder, inp["iheap"])
    return {"t": t.tolist(), "iorder": iorder.tolist()}


@handler("cmprlb")
def _cmprlb(inp: dict) -> dict:
    from core.cmprlb import cmprlb
    n, m = inp["n"], inp["m"]
    r = np.array(inp["r_in"], dtype=np.float64).copy()
    info = cmprlb(
        n, m,
        np.array(inp["x"], dtype=np.float64),
        np.array(inp["g"], dtype=np.float64),
        np.array(inp["ws"], dtype=np.float64).reshape(n, m),
        np.array(inp["wy"], dtype=np.float64).reshape(n, m),
        np.array(inp["sy"], dtype=np.float64).reshape(m, m),
        np.array(inp["wt"], dtype=np.float64).reshape(m, m),
        np.array(inp["z"], dtype=np.float64),
        r,
        np.array(inp["wa_in"], dtype=np.float64).copy(),
        np.array(inp["index"], dtype=np.int32),
        inp["theta"], inp["col"], inp["head"], inp["nfree"], inp["cnstnd"],
    )
    return {"r": r.tolist(), "info": info}


@handler("bmv")
def _bmv(inp: dict) -> dict:
    from core.bmv import bmv
    p = np.array(inp["p_in"], dtype=np.float64).copy()
    info = bmv(
        inp["m"],
        np.array(inp["sy"], dtype=np.float64),
        np.array(inp["wt"], dtype=np.float64),
        inp["col"],
        np.array(inp["v"], dtype=np.float64),
        p,
    )
    return {"p": p.tolist(), "info": info}


@handler("formt")
def _formt(inp: dict) -> dict:
    from core.formt import formt
    m = inp["m"]
    wt = np.array(inp["wt_in"], dtype=np.float64).copy()
    info = formt(
        m, wt,
        np.array(inp["sy"], dtype=np.float64),
        np.array(inp["ss"], dtype=np.float64),
        inp["col"], inp["theta"],
    )
    if "wt_upper" in inp.get("__expected_keys__", []):
        return {"wt_upper": wt.tolist(), "info": info}
    # Default: return upper triangle padded with original entries on lower triangle.
    return {"wt_upper": wt.tolist(), "info": info}


@handler("matupd")
def _matupd(inp: dict) -> dict:
    from core.matupd import MatupdState, matupd
    n, m = inp["n"], inp["m"]
    ws = np.array(inp["ws_in"], dtype=np.float64).reshape(n, m).copy()
    wy = np.array(inp["wy_in"], dtype=np.float64).reshape(n, m).copy()
    sy = np.array(inp["sy_in"], dtype=np.float64).reshape(m, m).copy()
    ss = np.array(inp["ss_in"], dtype=np.float64).reshape(m, m).copy()
    state = MatupdState(itail=inp["itail_in"], col=inp["col_in"], head=inp["head_in"], theta=0.0)
    matupd(
        n, m, ws, wy, sy, ss,
        np.array(inp["d"], dtype=np.float64),
        np.array(inp["r"], dtype=np.float64),
        inp["iupdat"], state,
        inp["rr"], inp["dr"], inp["stp"], inp["dtd"],
    )
    return {
        "head": state.head, "itail": state.itail, "col": state.col, "theta": state.theta,
        "ws": ws.tolist(), "wy": wy.tolist(), "sy": sy.tolist(), "ss": ss.tolist(),
    }


@handler("dcstep")
def _dcstep(inp: dict) -> dict:
    from core.dcstep import DcstepState, dcstep
    state = DcstepState(
        stx=inp["stx_in"], fx=inp["fx_in"], dx=inp["dx_in"],
        sty=inp["sty_in"], fy=inp["fy_in"], dy=inp["dy_in"],
        stp=inp["stp_in"], brackt=inp["brackt_in"],
    )
    dcstep(state, inp["fp"], inp["dp"], inp["stpmin"], inp["stpmax"])
    return {
        "stx": state.stx, "fx": state.fx, "dx": state.dx,
        "sty": state.sty, "fy": state.fy, "dy": state.dy,
        "stp": state.stp, "brackt": state.brackt,
    }


@handler("formk")
def _formk(inp: dict) -> dict:
    from core.formk import formk
    n, m = inp["n"], inp["m"]
    ind = np.array(inp["ind"], dtype=np.int32)
    indx2 = np.array(inp["indx2"], dtype=np.int32)
    ws = np.array(inp["ws"], dtype=np.float64).reshape(n, m)
    wy = np.array(inp["wy"], dtype=np.float64).reshape(n, m)
    sy = np.array(inp["sy"], dtype=np.float64).reshape(m, m)
    wn = np.array(inp["wn_in"], dtype=np.float64).reshape(2 * m, 2 * m).copy()
    wn1 = np.array(inp["wn1_in"], dtype=np.float64).reshape(2 * m, 2 * m).copy()
    info = formk(
        n, inp["nsub"], ind, inp["nenter"], inp["ileave"], indx2,
        inp["iupdat"], inp["updatd"], wn, wn1, m, ws, wy, sy,
        inp["theta"], inp["col"], inp["head"],
    )
    return {"info": info, "wn": wn.tolist(), "wn1": wn1.tolist()}


@handler("dcsrch")
def _dcsrch(inp: dict) -> dict:
    from core.dcsrch import DcsrchState, dcsrch
    state = DcsrchState(task=inp["task_in"], isave=np.zeros(2, dtype=np.int32), dsave=np.zeros(13))
    f, g, stp = dcsrch(
        inp["f"], inp["g"], inp["stp"],
        inp["ftol"], inp["gtol"], inp["xtol"],
        inp["stpmin"], inp["stpmax"], state,
    )
    return {"task": state.task, "stp": stp, "f": f, "g": g}


@handler("lnsrlb")
def _lnsrlb(inp: dict) -> dict:
    from core.lnsrlb import LnsrlbState, lnsrlb
    n = inp["n"]
    state = LnsrlbState(
        fold=0.0, gd=0.0, gdold=0.0,
        stp=inp.get("stp_in", -42.0),
        dnorm=0.0, dtd=0.0, xstep=0.0, stpmx=0.0,
        ifun=inp["ifun_in"], iback=inp["iback_in"], nfgv=inp["nfgv_in"],
        info=inp["info_in"], task=inp["task_in"], csave="",
        isave=np.zeros(2, dtype=np.int32),
        dsave=np.zeros(13, dtype=np.float64),
    )
    x = np.array(inp["x"], dtype=np.float64).copy()
    f = lnsrlb(
        n,
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        x, inp["f_in"],
        np.array(inp["g"], dtype=np.float64).copy(),
        np.array(inp["d"], dtype=np.float64),
        np.array(inp["r_in"], dtype=np.float64).copy(),
        np.array(inp["t_in"], dtype=np.float64).copy(),
        np.array(inp["z"], dtype=np.float64),
        inp["iter"], inp["boxed"], inp["cnstnd"], state,
    )
    return {
        "stpmx": state.stpmx, "stp": state.stp, "dnorm": state.dnorm,
        "info": state.info, "task": state.task, "ifun": state.ifun,
        "x": x.tolist(),
    }


@handler("cauchy")
def _cauchy(inp: dict) -> dict:
    from core.cauchy import cauchy
    n, m, col = inp["n"], inp["m"], inp["col"]
    x = np.array(inp["x"], dtype=np.float64)
    iwhere = np.array(inp["iwhere_in"], dtype=np.int32).copy()
    iorder = np.zeros(n, dtype=np.int32)
    t = np.zeros(n)
    d = np.zeros(n)
    xcp = np.full(n, -42.0)
    ws = np.array(inp["ws"], dtype=np.float64).reshape(n, max(1, col))
    wy = np.array(inp["wy"], dtype=np.float64).reshape(n, max(1, col))
    sy = np.array(inp["sy"], dtype=np.float64).reshape(m, m)
    wt = np.array(inp["wt"], dtype=np.float64).reshape(m, m)
    p = np.zeros(max(2, 2 * col))
    c = np.zeros(max(2, 2 * col))
    wbp = np.zeros(max(2, 2 * col))
    v = np.zeros(max(2, 2 * col))
    info, nseg = cauchy(
        n, x,
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        np.array(inp["g"], dtype=np.float64),
        iorder, iwhere, t, d, xcp,
        m, wy, ws, sy, wt, inp["theta"], col, inp["head"],
        p, c, wbp, v, inp["sbgnrm"], inp["epsmch"], inp["iprint"],
    )
    return {"xcp": xcp.tolist(), "iwhere": iwhere.tolist(), "info": info, "nseg": nseg}


@handler("subsm")
def _subsm(inp: dict) -> dict:
    from core.subsm import subsm
    n, m = inp["n"], inp["m"]
    nsub = inp["nsub"]
    ind = np.array(inp["ind"], dtype=np.int32) if nsub > 0 else np.zeros(0, dtype=np.int32)
    x = np.array(inp["x_in"], dtype=np.float64).copy()
    d = np.array(inp["d_in"], dtype=np.float64).copy() if inp["d_in"] else np.zeros(max(1, nsub))
    xp = np.full(n, -42.0)
    ws = np.array(inp["ws"], dtype=np.float64).reshape(n, m)
    wy = np.array(inp["wy"], dtype=np.float64).reshape(n, m)
    wn = np.array(inp["wn_in"], dtype=np.float64).reshape(2 * m, 2 * m).copy()
    wv = np.zeros(2 * m)
    iword, info = subsm(
        n, m, nsub, ind,
        np.array(inp["l"], dtype=np.float64),
        np.array(inp["u"], dtype=np.float64),
        np.array(inp["nbd"], dtype=np.int32),
        x, d, xp, ws, wy, inp["theta"],
        np.array(inp["xx"], dtype=np.float64),
        np.array(inp["gg"], dtype=np.float64),
        inp["col"], inp["head"], wv, wn,
    )
    return {"x": x.tolist(), "iword": iword, "info": info}


# ---------------------------------------------------------------------------
# Fortran engine: subprocess-call ``conformance_driver``, marshal text I/O.
# ---------------------------------------------------------------------------

# Per-subroutine type schemas. For each input key, the entry is the
# conformance_driver type code (i, r, s, b, iv, rv, im, rm). Used when
# the JSON value is a list and we need to disambiguate vec vs mat, and
# when the JSON value is a Python int that should be marshalled as a
# real (rare; we don't have such cases yet).
SCHEMAS: dict[str, dict[str, str]] = {
    "projgr": {"n": "i", "x": "rv", "l": "rv", "u": "rv", "nbd": "iv", "g": "rv"},
    "errclb": {
        "n": "i", "m": "i", "factr": "r",
        "l": "rv", "u": "rv", "nbd": "iv",
        "task_in": "s", "info_in": "i", "k_in": "i",
    },
    "active": {
        "n": "i", "l": "rv", "u": "rv", "nbd": "iv",
        "x_in": "rv", "iprint": "i",
    },
    "freev": {
        "n": "i", "nfree_in": "i",
        "index_in": "iv", "indx2_in": "iv", "iwhere": "iv",
        "updatd": "b", "cnstnd": "b", "iter": "i", "iprint": "i",
    },
    "hpsolb": {"n": "i", "t_in": "rv", "iorder_in": "iv", "iheap": "i"},
    "cmprlb": {
        "n": "i", "m": "i", "col": "i", "head": "i", "nfree": "i",
        "theta": "r", "cnstnd": "b",
        "index": "iv", "x": "rv", "z": "rv", "g": "rv",
        "r_in": "rv", "wa_in": "rv",
        "ws": "rm", "wy": "rm", "sy": "rm", "wt": "rm",
    },
    "bmv": {
        "m": "i", "col": "i",
        "sy": "rm", "wt": "rm", "v": "rv", "p_in": "rv",
    },
    "formt": {
        "m": "i", "col": "i", "theta": "r",
        "sy": "rm", "ss": "rm", "wt_in": "rm",
    },
    "matupd": {
        "n": "i", "m": "i", "iupdat": "i",
        "head_in": "i", "itail_in": "i", "col_in": "i",
        "ws_in": "rm", "wy_in": "rm", "sy_in": "rm", "ss_in": "rm",
        "d": "rv", "r": "rv",
        "rr": "r", "dr": "r", "stp": "r", "dtd": "r",
    },
    "dcstep": {
        "stx_in": "r", "fx_in": "r", "dx_in": "r",
        "sty_in": "r", "fy_in": "r", "dy_in": "r",
        "stp_in": "r", "fp": "r", "dp": "r",
        "brackt_in": "b", "stpmin": "r", "stpmax": "r",
    },
    "formk": {
        "n": "i", "m": "i", "nsub": "i", "nenter": "i", "ileave": "i",
        "iupdat": "i", "col": "i", "head": "i",
        "updatd": "b", "theta": "r",
        "ind": "iv", "indx2": "iv",
        "ws": "rm", "wy": "rm", "sy": "rm",
        "wn_in": "rm", "wn1_in": "rm",
    },
    "dcsrch": {
        "f": "r", "g": "r", "stp": "r",
        "ftol": "r", "gtol": "r", "xtol": "r",
        "stpmin": "r", "stpmax": "r", "task_in": "s",
    },
    "lnsrlb": {
        "n": "i",
        "l": "rv", "u": "rv", "nbd": "iv",
        "x": "rv", "f_in": "r", "g": "rv", "d": "rv",
        "r_in": "rv", "t_in": "rv", "z": "rv",
        "stp_in": "r", "iter": "i",
        "ifun_in": "i", "iback_in": "i", "nfgv_in": "i", "info_in": "i",
        "task_in": "s", "boxed": "b", "cnstnd": "b",
    },
    "cauchy": {
        "n": "i", "m": "i", "col": "i", "head": "i", "iprint": "i",
        "theta": "r", "sbgnrm": "r", "epsmch": "r",
        "x": "rv", "l": "rv", "u": "rv", "nbd": "iv", "g": "rv",
        "iwhere_in": "iv",
        "ws": "rm", "wy": "rm", "sy": "rm", "wt": "rm",
    },
    "subsm": {
        "n": "i", "m": "i", "nsub": "i", "col": "i", "head": "i",
        "theta": "r",
        "ind": "iv", "l": "rv", "u": "rv", "nbd": "iv",
        "x_in": "rv", "d_in": "rv", "xx": "rv", "gg": "rv",
        "ws": "rm", "wy": "rm", "wn_in": "rm",
    },
}


def _format_real(x: float) -> str:
    """Format a float to a 17-digit IEEE-754 round-tripping representation."""
    return f"{float(x):.17g}"


def _write_input(sub: str, inputs: dict, path: str) -> None:
    """Write ``inputs`` in conformance_driver text format."""
    schema = SCHEMAS.get(sub)
    if schema is None:
        raise ValueError(f"no schema defined for subroutine '{sub}'")
    lines = [f"sub {sub}"]
    for key, val in inputs.items():
        typ = schema.get(key)
        if typ is None:
            # Skip extras not in the schema (e.g. JSON-only metadata).
            continue
        if typ == "i":
            lines.append(f"i {key} {int(val)}")
        elif typ == "r":
            lines.append(f"r {key} {_format_real(val)}")
        elif typ == "s":
            lines.append(f"s {key} {val}")
        elif typ == "b":
            lines.append(f"b {key} {'T' if val else 'F'}")
        elif typ == "iv":
            n = len(val)
            ints = " ".join(str(int(v)) for v in val)
            lines.append(f"iv {key} {n} {ints}")
        elif typ == "rv":
            n = len(val) if val else 0
            reals = " ".join(_format_real(v) for v in val) if n > 0 else ""
            lines.append(f"rv {key} {n} {reals}".rstrip())
        elif typ == "im":
            mat = np.array(val, dtype=np.int64)
            r, c = mat.shape
            flat = mat.T.ravel()              # column-major
            ints = " ".join(str(int(v)) for v in flat)
            lines.append(f"im {key} {r} {c} {ints}")
        elif typ == "rm":
            mat = np.array(val, dtype=np.float64)
            r, c = mat.shape
            flat = mat.T.ravel()              # column-major
            reals = " ".join(_format_real(v) for v in flat)
            lines.append(f"rm {key} {r} {c} {reals}")
        else:
            raise ValueError(f"unknown type code '{typ}' for {key}")
    lines.append("end")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _parse_output(path: str) -> dict:
    """Parse conformance_driver's text output into a dict.

    Output keys preserve their type-coded interpretation: int, real,
    string, bool, list-of-int, list-of-float, list-of-list (matrix
    de-flattened from column-major).
    """
    out: dict[str, Any] = {}
    with open(path) as fh:
        text = fh.read()
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("sub ") or line == "end":
            continue
        toks = line.split(None, 2)
        typ = toks[0]
        key = toks[1]
        rest = toks[2] if len(toks) > 2 else ""
        if typ == "i":
            out[key] = int(rest)
        elif typ == "r":
            out[key] = float(rest)
        elif typ == "s":
            out[key] = rest
        elif typ == "b":
            out[key] = (rest.strip() == "T")
        elif typ in ("iv", "rv"):
            parts = rest.split()
            n = int(parts[0])
            vals = parts[1:1 + n]
            if typ == "iv":
                out[key] = [int(v) for v in vals]
            else:
                out[key] = [float(v) for v in vals]
        elif typ in ("im", "rm"):
            parts = rest.split()
            r, c = int(parts[0]), int(parts[1])
            flat = parts[2:2 + r * c]
            cast = int if typ == "im" else float
            arr = np.array([cast(v) for v in flat]).reshape(c, r).T  # col-major -> row-major
            out[key] = arr.tolist()
        else:
            raise ValueError(f"unknown output type code '{typ}'")
    return out


def make_fortran_handler(driver_bin: str, sub: str) -> Callable[[dict], dict]:
    """Return a handler that runs the F77 conformance driver for ``sub``."""
    def _handler(inp: dict) -> dict:
        with tempfile.TemporaryDirectory() as td:
            in_path = os.path.join(td, "in.txt")
            out_path = os.path.join(td, "out.txt")
            _write_input(sub, inp, in_path)
            res = subprocess.run(
                [driver_bin, in_path, out_path],
                capture_output=True, text=True, timeout=30,
            )
            if res.returncode != 0:
                raise RuntimeError(
                    f"conformance_driver failed (code {res.returncode}):\n"
                    f"  stderr: {res.stderr}\n  stdout: {res.stdout}"
                )
            return _parse_output(out_path)
    return _handler


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run(data_dir: str, strict: bool, filter_sub: str | None,
        engine: str = "python", driver_bin: str = DEFAULT_FORTRAN_DRIVER) -> int:
    pattern = os.path.join(data_dir, "*.json")
    paths = sorted(glob.glob(pattern), key=lambda p: (
        re.sub(r"_case_\d+", "", p),
        int(re.search(r"_case_(\d+)", p).group(1)) if re.search(r"_case_(\d+)", p) else 0,
    ))

    if not paths:
        print(f"No JSON test vectors found under {data_dir}", file=sys.stderr)
        return 1

    if engine == "fortran":
        if not os.path.isfile(driver_bin) or not os.access(driver_bin, os.X_OK):
            print(f"Fortran conformance driver not found or not executable: {driver_bin}",
                  file=sys.stderr)
            print("Build it with: cmake --build build --target conformance_driver",
                  file=sys.stderr)
            return 2
        # Replace HANDLERS with subprocess-based handlers.
        active_handlers: dict[str, Callable[[dict], dict]] = {
            sub: make_fortran_handler(driver_bin, sub) for sub in SCHEMAS
        }
    else:
        active_handlers = HANDLERS

    n_pass = 0
    n_fail = 0
    n_skip = 0
    failures = []

    for p in paths:
        with open(p) as fh:
            spec = json.load(fh)
        name = spec["name"]
        sub = name.rsplit("_case_", 1)[0]
        if filter_sub and sub != filter_sub:
            continue
        if sub not in active_handlers:
            n_skip += 1
            continue
        try:
            actual = active_handlers[sub](spec["inputs"])
        except Exception as exc:
            n_fail += 1
            failures.append((name, [f"handler raised: {type(exc).__name__}: {exc}"]))
            continue

        tol = TOLERANCES.get(sub, {"abs": 1e-13, "rel": 1e-13})
        issues: list[str] = []
        ok = True
        for key, expected_val in spec["expected"].items():
            actual_val = actual.get(key)
            if actual_val is None:
                issues.append(f"{key}: handler did not return this output")
                ok = False
                continue
            ok = _compare(actual_val, expected_val, key, strict, tol["abs"], tol["rel"], issues) and ok
        if ok:
            n_pass += 1
        else:
            n_fail += 1
            failures.append((name, issues))

    mode = "strict" if strict else "tolerance"
    print(f"\n=== Conformance ({mode}, engine={engine}) ===")
    print(f"  Pass: {n_pass}")
    print(f"  Fail: {n_fail}")
    print(f"  Skip: {n_skip}")
    if failures:
        print("\nFailures:")
        for name, issues in failures:
            print(f"  {name}:")
            for issue in issues[:5]:
                print(f"    - {issue}")
            if len(issues) > 5:
                print(f"    ... and {len(issues) - 5} more")

    return 0 if n_fail == 0 else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--strict", action="store_true",
        help="Bit-equality (default).",
    )
    parser.add_argument(
        "--tolerance", action="store_true",
        help="Per-subroutine numerical tolerance (looser).",
    )
    parser.add_argument(
        "--data", default=os.path.join(SPEC_ROOT, "data"),
        help="Directory of JSON test vectors.",
    )
    parser.add_argument(
        "--filter", default=None,
        help="Restrict to a single subroutine (e.g. 'projgr').",
    )
    parser.add_argument(
        "--engine", choices=("python", "fortran"), default="python",
        help="Reference implementation to validate against. 'fortran' "
             "subprocess-calls tests/conformance/conformance_driver "
             "(must be built first with cmake).",
    )
    parser.add_argument(
        "--driver-bin", default=DEFAULT_FORTRAN_DRIVER,
        help="Path to the Fortran conformance driver binary.",
    )
    args = parser.parse_args()

    if args.tolerance and args.strict:
        print("--strict and --tolerance are mutually exclusive", file=sys.stderr)
        return 2

    strict = not args.tolerance                # default to strict
    return run(args.data, strict, args.filter, args.engine, args.driver_bin)


if __name__ == "__main__":
    sys.exit(main())
