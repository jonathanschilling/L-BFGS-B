# errclb

## Purpose

Validate user-supplied inputs to the optimizer. Detected errors are
reported through `task` (a human-readable error message), `info`
(a numeric error code), and `k` (the index of the offending variable
for array-valued errors).

Called once at startup, before any iteration.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | integer | Claimed problem dimension. May be invalid (`≤ 0`). |
| `m` | integer | Claimed memory parameter. May be invalid (`≤ 0`). |
| `factr` | real | Function-decrease tolerance. May be invalid (`< 0`). |
| `l` | real vector, length `n` if `n > 0` | Lower bounds. |
| `u` | real vector, length `n` if `n > 0` | Upper bounds. |
| `nbd` | integer vector, length `n` if `n > 0` | Bound-type codes per variable. May be invalid (outside `{0, 1, 2, 3}`). |

### Logical outputs

| Name | Type | Description |
|------|------|-------------|
| `task` | string | Set to an `'ERROR: ...'` string on validation failure; **untouched** on success. |
| `info` | integer | `-6` if any `nbd[i] ∉ {0, 1, 2, 3}`; `-7` if any `nbd[i] = 2` with `l[i] > u[i]`; **untouched** otherwise. |
| `k` | integer | Index of the *last* offending variable (1-based); **untouched** if no array-valued error. |

### Preconditions

None. The whole point of this routine is to validate, so it must be
robust to garbage inputs (subject to the array length being at least
`n`, which is the caller's responsibility).

### Postconditions

- `task` may be unchanged (success), or contain an `'ERROR: ...'`
  message.
- If `task` is unchanged, `info` and `k` are also unchanged.

## Algorithm

### Decision flow

The validation runs sequentially and is **not short-circuiting**:
multiple checks may execute even after the first failure, and later
writes to `task` overwrite earlier ones. The order of checks matters
for the final `task` value when multiple errors are present.

```
1. if n <= 0:        task = 'ERROR: N .LE. 0'
2. if m <= 0:        task = 'ERROR: M .LE. 0'
3. if factr < 0:     task = 'ERROR: FACTR .LT. 0'
4. for i = 1 to n:
       if nbd[i] not in {0, 1, 2, 3}:
           task = 'ERROR: INVALID NBD'
           info = -6
           k = i
       if nbd[i] = 2 and l[i] > u[i]:
           task = 'ERROR: NO FEASIBLE SOLUTION'
           info = -7
           k = i
```

### "Last error wins" semantics

When multiple errors are present:

- If `n ≤ 0`, `m ≤ 0`, and `factr < 0` are all true:
  `task = 'ERROR: FACTR .LT. 0'` (check 3 wins).
- If both `n ≤ 0` and `nbd[1] = -1`: the nbd loop runs (because the
  loop bound `n` is an integer; if `n ≤ 0` the loop does nothing in
  F77 — see below); the final `task` depends on which checks fire.
- Within the nbd loop: if both an invalid `nbd` and an infeasible
  `nbd = 2` occur on different indices, the *last* such index wins
  for `k`, and the last-fired check determines `task`.

**This is the F77 behavior; ports must reproduce it for `--strict`
conformance.** A port that returns early on first error will pass the
"is the input valid?" test but produce different `task`/`info`/`k`
values on multi-error inputs.

### Behavior when `n ≤ 0`

The F77 loop `do i = 1, n` with `n ≤ 0` executes zero iterations
(see `06_portability_notes.md`). So the nbd-array checks are skipped
when `n ≤ 0`, even though `task` has been set by check 1. This is
**defensive** behavior: the array `nbd` may be empty or malformed
when `n ≤ 0`, and skipping the loop avoids out-of-bounds reads.

### Pseudocode

```
input:  n, m, factr, l, u, nbd, task (in/out), info (in/out), k (in/out)
output: task, info, k mutated as needed

if n <= 0:
    task = 'ERROR: N .LE. 0'
if m <= 0:
    task = 'ERROR: M .LE. 0'
if factr < 0:
    task = 'ERROR: FACTR .LT. 0'

for i = 1 to n:                            # ascending order required
    if nbd[i] < 0 or nbd[i] > 3:
        task = 'ERROR: INVALID NBD'
        info = -6
        k = i
    if nbd[i] == 2:
        if l[i] > u[i]:
            task = 'ERROR: NO FEASIBLE SOLUTION'
            info = -7
            k = i
```

### Magic constants

None.

### Numerical safeguards

None. All comparisons are integer or simple real comparisons.

### Order-of-operations dependencies

- The three top-level checks (`n`, `m`, `factr`) are in fixed order;
  `factr < 0` overwrites a prior `n` or `m` error.
- The `nbd` loop iterates `i = 1..n` in ascending order; for
  multi-error inputs, the last error's index is recorded in `k`.

Both orderings are required for `--strict` conformance.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/errclb_case_1.json` | All inputs valid (no error path) |
| 2 | `data/errclb_case_2.json` | `n = 0` triggers the n-check |
| 3 | `data/errclb_case_3.json` | `m = 0` triggers the m-check |
| 4 | `data/errclb_case_4.json` | `factr < 0` triggers the factr-check |
| 5 | `data/errclb_case_5.json` | `nbd[i] < 0` triggers invalid-nbd; `info=-6` |
| 6 | `data/errclb_case_6.json` | `nbd[i] > 3` triggers invalid-nbd from above |
| 7 | `data/errclb_case_7.json` | `nbd=2` with `l > u` triggers infeasible; `info=-7` |
| 8 | `data/errclb_case_8.json` | `nbd=2` with `l ≤ u`: feasible, no error |

## Reference implementation

`reference_impl/core/errclb.py` — Python implementation that follows
this spec.

## Mapping to other languages

A natural port idiom is to *replace* `errclb` with input validation
that raises an exception (Python: `ValueError`; Java: `IllegalArgumentException`;
Rust: `Err`; Go: `error`). The exception type can encode `info`; the
exception message can be the `task` string; the offending index can
be a field on the exception.

For `--strict` conformance, the port must:
- Detect all the same error conditions.
- Produce the same `task` string (or equivalent enum value), `info`,
  and `k` *as if* it had run the F77 sequential overwrite logic.

Ports that prefer fail-fast semantics (return on first error) lose
bit-conformance on multi-error inputs but preserve the "is the input
valid?" property. Whether this matters depends on the port's audience.

## Cross-references

- **Paper**: not algorithmic; this is implementation-level input
  validation.
- **Related subroutines**: called by `mainlb` once at startup before
  `active`. If `errclb` sets `task = 'ERROR: ...'`, `mainlb` returns
  immediately without further computation.
- **F77 source**: `src/errclb.f`.
- **Unit test**: `tests/unit/test_errclb.f90` (8 case_* blocks).
