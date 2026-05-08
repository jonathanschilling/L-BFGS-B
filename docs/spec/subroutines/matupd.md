# matupd

## Purpose

Append the latest step pair `(s, y)` to the L-BFGS history and
update the cached compact-representation matrices `sy`, `ss`. Also
recompute `theta = y'y / s'y`.

When the buffer is at capacity (`col = m` already), the *oldest* pair
is dropped to make room -- the storage is cyclic, with `head`/`itail`
indices marking the oldest and newest columns of `S` / `Y`.

## Mathematical contract

### Logical inputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Number of variables. |
| `m` | positive integer | Memory parameter. |
| `iupdat` | positive integer | Total number of L-BFGS updates *including this one* (incremented by `mainlb` before the call). |
| `head` (in/out) | integer, `1 <= head <= m` | Cyclic index of oldest column. Updated when `iupdat > m`. |
| `itail` (in/out) | integer, `1 <= itail <= m` | Cyclic index of newest column. Updated unconditionally. |
| `col` (in/out) | integer, `0 <= col <= m` | Number of stored pairs. Updated to `min(iupdat, m)`. |
| `d` | real vector, length `n` | Search direction (the new `s`-vector is `stp * d`; F77 stores `d` directly and folds `stp` into `ss[col, col]`). |
| `r` | real vector, length `n` | The new `y`-vector `g_new - g`. |
| `stp` | positive real | Line-search step length used to compute the new `s = stp * d`. |
| `dr` | real | Curvature `d' r = s' y / stp` (caller validated `dr > 0`). |
| `rr` | nonneg real | `||r||^2 = y' y`. |
| `dtd` | nonneg real | `||d||^2 = d' d`. |
| `ws`, `wy` (in/out) | real matrices `n * m` | History columns. The column at `itail` is overwritten. |
| `sy`, `ss` (in/out) | real matrices `m * m` | Cached `S'Y` / `S'S` (lower / upper packing -- see `01_glossary.md`). |
| `theta` (out) | real | Reset to `rr / dr`. |

### Logical outputs

All in-place mutations of the in/out arguments above. There is no
separate return value.

### Preconditions

- `dr > 0`. The caller (`mainlb`) checks the curvature condition
  `s' y > eps * ||y||^2` *before* calling `matupd`; `dr` here is
  `s' y / stp`, and the caller's check ensures the corresponding
  `s' y` is positive.
- `(s = stp * d, y = r)` are the step and gradient-difference at the
  iteration that just completed.

### Postconditions

- `col = min(iupdat, m)`.
- The newest column of `ws`, `wy` (at index `itail`) holds `d`, `r`.
- `theta = rr / dr`.
- `sy`, `ss` upper-and-lower packed matrices reflect the new history.

## Algorithm

### Phase 1: update slot pointers and `col`

```
if iupdat <= m:                                # buffer not yet full
    col = iupdat
    itail = ((head + iupdat - 2) mod m) + 1    # append slot
else:                                          # buffer already full -> cyclic shift
    itail = (itail mod m) + 1                  # advance newest pointer
    head  = (head  mod m) + 1                  # advance oldest pointer (drop oldest)
```

When `iupdat = 1`, the formula gives `itail = ((head - 1) mod m) + 1 = head`,
which is correct (single column at slot 1 if `head = 1`).

### Phase 2: copy `d`, `r` into `ws`, `wy` at column `itail`

```
ws[:, itail] = d
wy[:, itail] = r
```

(In F77 these are `dcopy` calls; ports use whatever idiom their
language provides.)

### Phase 3: update `theta`

```
theta = rr / dr           # = y'y / s'y, the standard scaling for B_0
```

### Phase 4: if buffer overflowed, shift `sy` and `ss`

When `iupdat > m`, the stored `(s, y)` pairs have all advanced one
slot; the cached `sy` and `ss` matrices must be left-and-up-shifted
by one to reflect the new ordering (oldest pair dropped).

```
if iupdat > m:
    for j = 1 to col - 1:
        # ss column j gets ss column j+1 with the top entry dropped
        for k = 1 to j:
            ss[k, j] = ss[k + 1, j + 1]
        # sy column j (subdiagonal entries) gets sy column j+1 with one top entry dropped
        for k = j to col - 1:
            sy[k, j] = sy[k + 1, j + 1]
```

In F77 these are `dcopy` calls with stride `1`; the indexing is
clearer in the explicit form above.

### Phase 5: fill the new last row of `sy` and last column of `ss`

```
pointr = head
for j = 1 to col - 1:
    sy[col, j] = d * wy[:, pointr]              # s_new ' y_j
    ss[j, col] = ws[:, pointr] * d              # s_j ' s_new
    pointr = (pointr mod m) + 1
```

The chronological walk through `pointr = head, head+1, ..., wrapping`
hits the columns in the order they were inserted; the inner products
fill the `(col, j)` and `(j, col)` slots for all the *previous* pairs.

### Phase 6: diagonal entries

```
if stp == 1.0:
    ss[col, col] = dtd                          # ||s||^2 = ||d||^2
else:
    ss[col, col] = stp^2 * dtd                   # ||s||^2 = stp^2 * ||d||^2
sy[col, col] = dr                               # s'y for the new pair
```

The branch on `stp == 1.0` exact equality is intentional: it skips
the multiplication when not needed, which is bit-equivalent to
multiplying by 1.0 anyway. **For `--strict` conformance ports must
preserve this branch** (the dispatch on exact `1.0` is observable in
diagnostic traces, though the numerical output is identical for IEEE
multiplication-by-one).

### Pseudocode (combined)

See above sections; combined into one routine in source order.

### Magic constants

The exact-equality test `stp == 1.0`. Documented above.

### Numerical safeguards

- Caller checks `s' y > eps * ||y||^2` and skips `matupd` if violated;
  `matupd` itself does not validate.
- `dr` should be positive on entry (otherwise `theta = rr/dr` is
  garbage). The unit tests use `dr > 0`.

### Order-of-operations dependencies

- The shift loop (Phase 4) must execute before the new-row loop
  (Phase 5), because Phase 5 indexes into the shifted `ws`, `wy`.
- The new-row loop walks `j = 1..col-1` in ascending order, with
  `pointr` advancing cyclically. `--strict` requires this order.
- The dot products (`ddot`) follow BLAS reference conventions.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/matupd_case_1.json` | First update (`iupdat = 1`), `stp = 0.5`: append, `ss[1,1] = stp^2*dtd` |
| 2 | `data/matupd_case_2.json` | Second update (`iupdat = 2`), still appending |
| 3 | `data/matupd_case_3.json` | Third update with `iupdat = 3 > m = 2`: cyclic overflow |
| 4 | `data/matupd_case_4.json` | First update with `stp = 1.0`: skips the multiply branch |

(The F77 unit test `test_matupd.f90` has 2 case_* blocks but each
exercises multiple sub-cases sequentially; we split into 4 JSON
vectors to give one logical operation per case.)

## Reference implementation

`reference_impl/core/matupd.py`.

## Cross-references

- **Paper**: `code.pdf` sec.2.2 (compact representation update),
  Byrd/Nocedal/Schnabel 1994 sec.3.
- **Related subroutines**: called by `mainlb` after each successful
  iteration where the curvature condition holds. Updated state is
  consumed by `formt` (recompute `T`), `bmv`, `cauchy`, `subsm`.
- **F77 source**: `src/matupd.f`.
- **Unit test**: `tests/unit/test_matupd.f90` (2 case_* blocks, 4
  sub-cases).
