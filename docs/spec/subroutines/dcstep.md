# dcstep

## Purpose

Compute the next safeguarded trial step in the More–Thuente line
search. Given the current bracket `[stx, sty]` (or just `stx` if no
bracket has been established yet), the function value/derivative at
the best step `stx`, the function value/derivative at the latest
trial `stp`, and the previous bracket endpoint info, produce the next
trial step `stp` clamped to `[stpmin, stpmax]`.

Also updates the bracket: which endpoint is replaced by the latest
trial, and whether the minimum has been definitively bracketed
(`brackt`).

This is the inner kernel of `dcsrch` (the line-search reverse-comm
driver). The algorithm follows More & Thuente 1994 §3 (case 1–4).

## Mathematical contract

### Logical inputs/outputs

All scalar (real or boolean), all in/out except where noted:

| Name | Type | Description |
|------|------|-------------|
| `stx` (in/out) | real | Best step so far (lowest `f`). Updated on success. |
| `fx`, `dx` (in/out) | real | `f` and `f'` at `stx`. Updated alongside `stx`. |
| `sty` (in/out) | real | Other bracket endpoint. May be updated. |
| `fy`, `dy` (in/out) | real | `f` and `f'` at `sty`. Updated alongside `sty`. |
| `stp` (in/out) | real | Latest trial step (in); next trial step (out). Always clamped to `[stpmin, stpmax]`. |
| `fp` (in) | real | `f` at the input `stp`. |
| `dp` (in) | real | `f'` at the input `stp`. |
| `brackt` (in/out) | boolean | True if a bracket has been established; the routine may set it true. |
| `stpmin`, `stpmax` (in) | real | Step-length bounds (clamps the new `stp`). |

### Preconditions

- If `brackt = true`, then `min(stx, sty) < stp < max(stx, sty)`.
- The derivative at `stx` is in the descent direction:
  `dx · (stp - stx) < 0`.

### Postconditions

- `stp` is the next trial step, clamped to `[stpmin, stpmax]` (and to
  the safeguarded interval when `brackt = true`).
- `stx`, `fx`, `dx`, `sty`, `fy`, `dy` are updated according to
  More-Thuente case logic (see below).
- `brackt` may transition from `false` to `true` (never the reverse).

## Algorithm

The algorithm dispatches on three boolean conditions:

```
case_higher_f      = (fp > fx)
case_opposite_sign = (sgnd < 0)               where sgnd = dp * sign(dx)
case_dec_grad      = (|dp| < |dx|)
```

into one of four cases:

- **Case 1 — Higher function value**: `fp > fx`. The minimum is
  bracketed by `[stx, stp]`. Compute cubic and quadratic step
  candidates; if cubic is closer to `stx`, take it; else take the
  average. **Set `brackt = true`.**
- **Case 2 — Lower `f`, opposite-sign derivatives**: `fp ≤ fx and sgnd < 0`.
  The minimum is bracketed by `[stx, stp]` (different bracket from
  case 1). Compute cubic and secant; pick whichever is farther from
  `stp`. **Set `brackt = true`.**
- **Case 3 — Lower `f`, same sign, magnitude of `f'` decreases**:
  `fp ≤ fx and sgnd ≥ 0 and |dp| < |dx|`. Cubic step is computed
  only if the cubic tends to ∞ in the search direction or its
  minimum lies beyond `stp`. Otherwise the cubic is set to the secant.
  Bracket-aware safeguarding clamps the result.
- **Case 4 — Lower `f`, same sign, magnitude of `f'` does not decrease**:
  `fp ≤ fx and sgnd ≥ 0 and |dp| ≥ |dx|`. If bracketed, cubic step
  using `sty` data instead of `stx`. Else: clamp to `stpmin` or
  `stpmax` depending on the side of the trial.

After computing `stpf`, update the bracket data:

```
if fp > fx:
    sty = stp; fy = fp; dy = dp                 # update sty endpoint
else:
    if sgnd < 0:
        sty = stx; fy = fx; dy = dx             # swap sty to old stx
    stx = stp; fx = fp; dx = dp                 # update stx to new best
```

Finally `stp = stpf`.

### Cubic-step computation (used in Cases 1, 2, 3)

```
theta = 3 * (fx - fp) / (stp - stx) + dx + dp
s = max(|theta|, |dx|, |dp|)
gamma = s * sqrt((theta/s)^2 - (dx/s) * (dp/s))    # may be 0 if discriminant <= 0
```

The `s`-normalization avoids overflow for large derivatives. In Case
3 the discriminant is clamped to non-negative via
`max(0, (theta/s)^2 - ...)` to handle numerical noise.

The sign of `gamma` is flipped depending on the side of `stp`
relative to `stx` (see source for exact conditions per case).

The cubic step is then `stx + (gamma - dx + theta) / ((gamma - dx) + gamma + dp) * (stp - stx)` (Case 1) or analogous.

### Secant / quadratic step

- **Case 1 quadratic**: `stpq = stx + (dx / ((fx - fp)/(stp - stx) + dx)) / 2 * (stp - stx)`.
- **Case 2 secant**: `stpq = stp + (dp / (dp - dx)) * (stx - stp)`.
- **Case 3 secant**: same form as Case 2.
- **Case 4 cubic** uses `sty`, `fy`, `dy` instead of `stx`, `fx`, `dx`.

### Bracket safeguarding (Case 3 only when `brackt`)

After computing `stpf` from the cubic-or-secant, clamp it to a
safeguarded interval:

```
if stp > stx:
    stpf = min(stp + 0.66 * (sty - stp), stpf)
else:
    stpf = max(stp + 0.66 * (sty - stp), stpf)
```

The factor `0.66` (variable `p66` in the F77 source) prevents the
new step from getting too close to the bracket boundary too quickly.

### Magic constants

| Constant | F77 | Value | Meaning |
|----------|-----|-------|---------|
| `p66` | `p66` | `0.66d0` | Bracket-safeguard multiplier. From More-Thuente 1994 §3.5. |
| `two` | `two` | `2.0d0` | Used in Case 1 quadratic step. |
| `three` | `three` | `3.0d0` | Cubic-interpolation coefficient. |

The value `0.66` is in the **magic-constants table** in `04_numerics.md`.

### Numerical safeguards

- The `s` normalization in cubic computation prevents overflow for
  large derivatives.
- In Case 3, the cubic discriminant is clamped to `≥ 0` via `max(0, ...)`.
- The new `stpf` is always clamped to `[stpmin, stpmax]` either
  explicitly (Case 3 unbracketed, Case 4 unbracketed) or implicitly
  via the bracket interval safeguarding.
- In Case 3 with `r ≥ 0` and `gamma = 0`, the cubic is set to
  `stpmax` or `stpmin` depending on the side; this avoids a
  divide-by-zero.

### Order-of-operations dependencies

Every arithmetic operation in `dcstep` is single-pass; no reductions
involved. **Order is fully fixed by the F77 source** — the cubic
computation, quadratic computation, choice between them, and
clamping are all sequential. Ports must follow the F77 sequence to
match bit-for-bit.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/dcstep_case_1.json` | Case 1: `fp > fx`, `stp > stx`. Sets `brackt = true`. |
| 2 | `data/dcstep_case_2.json` | Case 1: `fp > fx`, `stp < stx` — gamma sign flip. |
| 3 | `data/dcstep_case_3.json` | Case 2: opposite-sign derivatives. Bracketed. |
| 4 | `data/dcstep_case_4.json` | Case 3: decreasing-grad with `brackt = true`. |
| 5 | `data/dcstep_case_5.json` | Case 3: decreasing-grad with `brackt = false`. |
| 6 | `data/dcstep_case_6.json` | Case 3: `r ≥ 0` cubic-direction-clamp branch. |
| 7 | `data/dcstep_case_7.json` | Case 4: `brackt = true`, cubic from `sty`. |
| 8 | `data/dcstep_case_8.json` | Case 4: `brackt = false`, `stp > stx`: clamp to `stpmax`. |
| 9 | `data/dcstep_case_9.json` | Case 4: `brackt = false`, `stp < stx`: clamp to `stpmin`. |
| 10 | `data/dcstep_case_10.json` | Case 3: `stp < stx`, `r ≥ 0` ⟹ `stpmin` branch. |
| 11 | `data/dcstep_case_11.json` | Case 3: `brackt = true`, `stp < stx`, max safeguard. |

## Reference implementation

`reference_impl/core/dcstep.py`.

## Cross-references

- **Paper**: More & Thuente 1994 §3 (the algorithm). MINPACK-1 (1983)
  is the original implementation.
- **Related subroutines**: called by `dcsrch` to generate each new
  trial step.
- **F77 source**: `src/dcstep.f`.
- **Unit test**: `tests/unit/test_dcstep.f90` (11 case_* blocks).
