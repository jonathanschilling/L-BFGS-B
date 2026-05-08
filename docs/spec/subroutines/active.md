# active

## Purpose

Initialize the bound-status state at the start of an optimization:

1. Project the user-supplied initial point `x0` onto the feasible
   box `l <= x <= u` if necessary.
2. Compute the per-variable bound-status array `iwhere` in its
   coarse-initial form (`-1` = unbounded, `0` = bounded but free,
   `3` = fixed by `l = u`).
3. Set the constraint summary flags `prjctd`, `cnstnd`, `boxed`.

Called once at the very beginning of `mainlb`, before the first
iteration.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Number of variables. |
| `l` | real vector, length `n` | Lower bounds (consulted only where `nbd[i] in {1, 2}`). |
| `u` | real vector, length `n` | Upper bounds (consulted only where `nbd[i] in {2, 3}`). |
| `nbd` | integer vector, length `n` | Bound-type codes per variable. Validated by `errclb` before this call. |
| `x` | real vector, length `n` | Initial iterate (may be infeasible on entry). |
| `iprint` | integer | Diagnostic verbosity flag (see `02_api.md`). |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `x` | real vector, length `n` | **In/out**: projected to feasible set componentwise. |
| `iwhere` | integer vector, length `n` | Initial bound status per variable: `-1`, `0`, or `3`. |
| `prjctd` | boolean | `true` if any `x[i]` was strictly outside its bounds and projected. |
| `cnstnd` | boolean | `true` if any variable has at least one bound (`nbd[i] != 0`). |
| `boxed` | boolean | `true` if every variable has both bounds (`nbd[i] = 2` for all `i`). |

### Preconditions

- Inputs validated by `errclb`: `n >= 1`, `nbd[i] in {0, 1, 2, 3}` for
  all `i`, and where `nbd[i] = 2`, `l[i] <= u[i]`.

### Postconditions

- `x` is feasible: `l[i] <= x[i]` if `nbd[i] in {1, 2}` and
  `x[i] <= u[i]` if `nbd[i] in {2, 3}`.
- `iwhere[i] in {-1, 0, 3}` (the coarse initial values; the finer
  values `1`, `2` are introduced later by `cauchy`).
- `l`, `u`, `nbd`, `n` are unchanged.

## Algorithm

### Phase 1: project `x` to feasible set

```
nbdd = 0                          # count of variables already at a bound
prjctd = false
for i = 1 to n:
    if nbd[i] > 0:                # has at least one bound
        if nbd[i] <= 2 and x[i] <= l[i]:    # has lower bound, and x at-or-below
            if x[i] < l[i]:
                prjctd = true
                x[i] = l[i]
            nbdd += 1
        else if nbd[i] >= 2 and x[i] >= u[i]:  # has upper, and x at-or-above
            if x[i] > u[i]:
                prjctd = true
                x[i] = u[i]
            nbdd += 1
```

The `else if` is critical: a variable cannot be both "at lower" and
"at upper" in the same iteration of the outer loop. The lower-bound
check fires first; if it fires, the upper-bound check is skipped for
that variable. (For `nbd = 2` with `l < u` and `x` interior, neither
`x <= l` nor `x >= u` is true, so the variable isn't counted in
`nbdd`.)

### Phase 2: initialize `iwhere`, `cnstnd`, `boxed`

```
cnstnd = false
boxed  = true
for i = 1 to n:
    if nbd[i] != 2:
        boxed = false
    if nbd[i] == 0:
        iwhere[i] = -1            # always free
    else:
        cnstnd = true
        if nbd[i] == 2 and u[i] - l[i] <= 0:
            iwhere[i] = 3         # fixed (l = u or worse)
        else:
            iwhere[i] = 0         # bounded but currently has slack
```

The `iwhere = 3` case captures both `l = u` (variable is pinned) and
the validation-edge case `l > u` (which `errclb` should have caught,
but `active` defends anyway).

### Pseudocode (combined)

```
input:  n, l, u, nbd, x (in/out), iprint
output: x (mutated), iwhere, prjctd, cnstnd, boxed

# Phase 1: project x
nbdd = 0
prjctd = false
for i = 1 to n:                                  # ascending order required
    if nbd[i] > 0:
        if nbd[i] <= 2 and x[i] <= l[i]:
            if x[i] < l[i]:
                prjctd = true
                x[i] = l[i]
            nbdd += 1
        else if nbd[i] >= 2 and x[i] >= u[i]:
            if x[i] > u[i]:
                prjctd = true
                x[i] = u[i]
            nbdd += 1

# Phase 2: classify and summarise
cnstnd = false
boxed  = true
for i = 1 to n:                                  # ascending order required
    if nbd[i] != 2:
        boxed = false
    if nbd[i] == 0:
        iwhere[i] = -1
    else:
        cnstnd = true
        if nbd[i] == 2 and u[i] - l[i] <= 0:
            iwhere[i] = 3
        else:
            iwhere[i] = 0

# Optional diagnostic prints (port-defined; see iprint policy in 02_api.md)
```

### Magic constants

None.

### Numerical safeguards

- The condition `u[i] - l[i] <= 0` uses subtraction rather than
  `u[i] <= l[i]` directly. For finite IEEE-754 doubles these are
  equivalent; for inputs where `u[i]` or `l[i]` is `+/-Inf`, `Inf - Inf
  = NaN`, but in this routine such pairs would never have `nbd = 2`
  (per `errclb` validation). No defensive handling required.
- The `nbdd` counter is informational (used only for diagnostic
  output) and does not affect the algorithm. Ports that don't print
  this diagnostic may omit it.

### Order-of-operations dependencies

The two phase loops both iterate `i = 1..n` in ascending order. The
order does not affect the final `x`, `iwhere` (each variable is
updated independently), `cnstnd`, or `boxed` (boolean OR / AND is
commutative). However, `nbdd` accumulates in order; a port that
reports `nbdd` in diagnostics must use ascending order to match the
F77 output.

`prjctd` is also set with OR semantics (any projection sets it
`true`); order doesn't matter for the final value.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/active_case_1.json` | All `nbd = 0`: cnstnd=false, boxed=false (n>0 with nbd!=2 means boxed=false), iwhere = -1 |
| 2 | `data/active_case_2.json` | All `nbd = 2` with `x` strictly interior: boxed=true, no projection |
| 3 | `data/active_case_3.json` | `nbd = 2`, `x = l` exactly: enters at-lower branch but no projection |
| 4 | `data/active_case_4.json` | `nbd = 2`, `x < l`: projects up |
| 5 | `data/active_case_5.json` | `nbd = 2`, `x = u` exactly: enters at-upper branch but no projection |
| 6 | `data/active_case_6.json` | `nbd = 2`, `x > u`: projects down |
| 7 | `data/active_case_7.json` | `nbd = 1`, `x < l`: projects up |
| 8 | `data/active_case_8.json` | `nbd = 3`, `x > u`: projects down |
| 9 | `data/active_case_9.json` | `nbd = 2` with `l = u`: variable fixed, iwhere = 3 |
| 10 | `data/active_case_10.json` | `iprint = 0` with all unbounded: triggers "unconstrained" diagnostic |
| 11 | `data/active_case_11.json` | `iprint > 0` with bounded var: triggers nbdd-count diagnostic |

## Reference implementation

`reference_impl/core/active.py` -- Python implementation.

## Cross-references

- **Paper**: `algorithm.pdf` sec.2 (active-set framework). The initial
  projection appears in `code.pdf` sec.3.
- **Related subroutines**: called once by `mainlb` at startup. Result
  feeds `cauchy` (refines `iwhere`), `freev`, and the projected-
  gradient test (`projgr`).
- **F77 source**: `src/active.f`.
- **Unit test**: `tests/unit/test_active.f90` (11 case_* blocks).
