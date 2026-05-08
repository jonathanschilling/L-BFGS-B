# projgr

## Purpose

Compute the infinity norm of the *projected gradient* at the current
iterate, for use in the convergence test `‖g_proj‖_∞ ≤ pgtol`.

The projected gradient measures how far the algorithm has converged
relative to the bounds: it is zero at any KKT point of the
bound-constrained problem and positive elsewhere.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Number of variables. |
| `x` | real vector, length `n` | Current iterate. Caller guarantees `l ≤ x ≤ u` componentwise (where bounds apply). |
| `l` | real vector, length `n` | Lower bounds. Entries with `nbd[i] ∈ {0, 3}` are not consulted. |
| `u` | real vector, length `n` | Upper bounds. Entries with `nbd[i] ∈ {0, 1}` are not consulted. |
| `nbd` | integer vector, length `n` | Bound types per variable: `0`, `1`, `2`, or `3` (see `02_api.md`). |
| `g` | real vector, length `n` | Gradient `∇f(x)` at the current iterate. |

### Logical output

| Name | Type | Description |
|------|------|-------------|
| `sbgnrm` | nonnegative real | Infinity norm of the projected gradient. |

### Preconditions

- `n ≥ 1`.
- `nbd[i] ∈ {0, 1, 2, 3}` for all `i` (caller validates via `errclb`).
- `x[i]` is feasible: `l[i] ≤ x[i]` if `nbd[i] ∈ {1, 2}` and
  `x[i] ≤ u[i]` if `nbd[i] ∈ {2, 3}`. (Caller ensures via initial
  projection in `active`.)

### Postconditions

- `sbgnrm ≥ 0`.
- Inputs `n`, `x`, `l`, `u`, `nbd`, `g` are unchanged.

## Algorithm

### Definition

For each variable `i`, compute a per-component projected-gradient
proxy `gi`:

```
gi = g[i]
if nbd[i] != 0:                            # has at least one bound
    if g[i] < 0:                           # gradient negative -> x wants to increase
        if nbd[i] >= 2:                    # has upper bound (nbd in {2, 3})
            gi = max(x[i] - u[i], g[i])
    else:                                  # gradient non-negative -> x wants to decrease
        if nbd[i] <= 2:                    # has lower bound (nbd in {1, 2})
            gi = min(x[i] - l[i], g[i])
sbgnrm = max(sbgnrm, |gi|)
```

The result is `sbgnrm = max_{i=1..n} |gi|`.

### Pseudocode (language-neutral)

```
input:  n, x, l, u, nbd, g
output: sbgnrm

sbgnrm = 0
for i = 1 to n:                            # ascending index order required
    gi = g[i]
    if nbd[i] != 0:
        if gi < 0:
            if nbd[i] >= 2:
                gi = max(x[i] - u[i], gi)
        else:                              # gi >= 0 (zero falls into the lower-bound branch)
            if nbd[i] <= 2:
                gi = min(x[i] - l[i], gi)
    sbgnrm = max(sbgnrm, abs(gi))
return sbgnrm
```

### Interpretation

The value `gi` is the negative of the standard "projected gradient
component" `(P[x - g, l, u] - x)[i]`, where `P[·, l, u]` is the
componentwise projection onto the box. To see this for a single
variable with both bounds (`nbd = 2`) and `g < 0`:

```
P[x - g, l, u] - x  =  P[x + |g|, l, u] - x      (since g < 0)
                   =  min(x + |g|, u) - x
                   =  min(|g|, u - x)
                   =  -max(g, x - u)             (negate)
                   =  -gi                         (matches the spec)
```

The infinity norm is invariant under sign flip, so `‖gi‖_∞` equals the
standard projected-gradient norm.

### Magic constants

None.

### Numerical safeguards

None required. The arithmetic is at most one subtraction and one
`min`/`max` per variable; no division, no square root.

### Order-of-operations dependencies

The reduction `sbgnrm = max(sbgnrm, |gi|)` is associative and
commutative for `max`, so the loop order does not affect the result.
However, **the F77 source iterates `i = 1..n` in ascending order**
and the `--strict` conformance mode requires the same ordering for
trace-level reproducibility (e.g., when comparing intermediate state
during debugging). Ports that only need the final scalar may iterate
in any order.

### Branching at `gi == 0`

The condition `gi < 0` distinguishes the negative-gradient and
non-negative-gradient branches. When `g[i] == 0` exactly (including
negative zero in IEEE-754), `g[i] < 0` evaluates to `false` and the
non-negative branch is taken. Negative zero (`-0.0`) compares equal to
positive zero in IEEE-754, so neither sign of zero changes the
result. The contribution of a zero gradient is always `0`.

### Branching at `nbd == 2`

A variable with `nbd[i] == 2` (both bounds) takes either the upper-
or lower-bound branch depending on the sign of `g[i]`, not both. The
F77 conditions `nbd >= 2` and `nbd <= 2` are mutually exclusive *given
the active branch* (both have `nbd != 0`).

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/projgr_case_1.json` | `nbd=0`, positive `g`: unbounded path |
| 2 | `data/projgr_case_2.json` | `nbd=0`, negative `g`: unbounded path with `|g|` |
| 3 | `data/projgr_case_3.json` | `nbd=2`, interior, `g > 0`: `min(x-l, g)` returns `g` |
| 4 | `data/projgr_case_4.json` | `nbd=2`, near lower bound, `g > 0`: `min(x-l, g)` returns `x-l` |
| 5 | `data/projgr_case_5.json` | `nbd=2`, near upper bound, `g < 0`: `max(x-u, g)` returns `x-u` |
| 6 | `data/projgr_case_6.json` | `nbd=1`, `g < 0`: lower-only, no upper-clip path |
| 7 | `data/projgr_case_7.json` | `nbd=3`, `g > 0`: upper-only, no lower-clip path |
| 8 | `data/projgr_case_8.json` | `n=4`, mixed `nbd`: `max` aggregation across components |

## Reference implementation

`reference_impl/core/projgr.py` — NumPy implementation that follows
this spec.

## Cross-references

- **Paper**: Conn/Gould/Toint 1988 §3 (active-set framework), the
  projected-gradient termination criterion is standard for box-
  constrained optimization. L-BFGS-B's specific formulation appears
  in `code.pdf` §3.
- **Related subroutines**: called by `mainlb` once per outer
  iteration. The result feeds the `pgtol` convergence test.
- **F77 source**: `src/projgr.f`.
- **Unit test**: `tests/unit/test_projgr.f90` (8 case_* blocks).
