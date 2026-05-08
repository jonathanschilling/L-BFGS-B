# hpsolb

## Purpose

Extract the **minimum** of array `t` and place it at position `n`,
leaving the remaining `n-1` elements as a min-heap rooted at position
`1`. Tracks original indices in a parallel array `iorder`.

Used by `cauchy` to process breakpoints in order of increasing
magnitude: the first call with `iheap = 0` builds the heap and
extracts the smallest breakpoint; subsequent calls with `iheap = 1`
extract the next-smallest in `O(log n)` each.

The algorithm is HEAPSORT (Williams 1964, CACM Algorithm 232).

## Mathematical contract

### Logical inputs/outputs

| Name | Type | Description |
|------|------|-------------|
| `n` | positive integer | Array length. |
| `t` | real vector, length `n` | **In/out**: on entry, the values to process; on exit, `t[n]` is the minimum, and `t[1..n-1]` is a min-heap of the remaining values. |
| `iorder` | integer vector, length `n` | **In/out**: on entry, the indices `iorder[i]` corresponding to `t[i]` (typically `iorder[i] = i` on the first call); on exit, permuted such that `(t[i], iorder[i])` pairs are preserved as a multiset, and `iorder[n]` is the index of the minimum. |
| `iheap` | integer | `0` if `t` is not yet a heap (this routine builds the heap first); non-zero if `t` is already a min-heap. |

### Preconditions

- `n >= 1`.
- `iorder[i]` are arbitrary integers (typically distinct, but the
  routine treats them as opaque tags).

### Postconditions

- `t[n]` is the minimum of the original `t[1..n]`.
- `iorder[n]` is the index (from the original parallel `iorder`) of
  that minimum element.
- `t[1..n-1]` is a valid min-heap.
- The multiset `{(t[i], iorder[i]) : 1 <= i <= n}` after the call
  equals the multiset before the call (only the order changes).

## Algorithm

### Phase 1: build heap (only if `iheap = 0`)

Repeated sift-up insertion: for each `k = 2, 3, ..., n`, insert the
pair `(t[k], iorder[k])` into the existing heap of size `k-1` rooted
at index 1.

```
if iheap == 0:
    for k = 2 to n:
        ddum   = t[k]
        indxin = iorder[k]
        i = k
        while i > 1:
            j = i / 2                 # parent
            if ddum < t[j]:           # strictly less -> swap up
                t[i] = t[j]
                iorder[i] = iorder[j]
                i = j
                # continue
            else:
                break                 # ddum is in correct place
        t[i] = ddum
        iorder[i] = indxin
```

The inner sift-up uses **strict less** (`ddum < t[j]`); equal-valued
elements stay in their original relative order ("stable for ties").

### Phase 2: extract minimum (only if `n > 1`)

```
if n > 1:
    out    = t[1]                     # the current minimum
    indxou = iorder[1]
    ddum   = t[n]                     # element to sift down from root
    indxin = iorder[n]
    i = 1
    while True:
        j = 2 * i                     # left child
        if j > n - 1:
            break                     # no children in the remaining heap
        if t[j+1] < t[j]:             # right child smaller than left
            j = j + 1
        if t[j] < ddum:
            t[i] = t[j]
            iorder[i] = iorder[j]
            i = j
        else:
            break
    t[i] = ddum
    iorder[i] = indxin
    t[n] = out
    iorder[n] = indxou
```

The sift-down compares against `n - 1` (not `n`), because position
`n` is reserved for the extracted minimum. Ties: prefer the left
child (`<` is strict).

### Pseudocode (combined)

See above.

### Magic constants

None.

### Numerical safeguards

- All comparisons are real-vs-real `<`. Equal values are treated as
  not-less, so the algorithm is stable for ties.
- No arithmetic on `t` values other than comparison (no addition,
  multiplication). No floating-point safeguards needed.

### Order-of-operations dependencies

- The sift-up sequence (insert order `k = 2, 3, ..., n`) is *not* an
  algorithmic invariant -- any valid heap-build (e.g., bottom-up
  Floyd's heapify) produces a valid heap. **However**, the F77
  source uses sift-up insertion in this specific order, and ties in
  `t` are broken by insertion order. For `--strict` conformance the
  port must match this insertion order.
- The sift-down ties also matter: F77 prefers the left child on
  equal values. Ports must match.

If a port uses a standard-library priority queue, the resulting
`iorder[n]` may differ from the F77 result on tied minima. In
practice, breakpoints in `cauchy` are extremely unlikely to be
exactly equal (they are computed from differences of distinct
floating-point values); but for `--strict` conformance, ports must
implement the F77 sift-up/sift-down explicitly.

## Test vectors

| Case | File | Branch exercised |
|------|------|------------------|
| 1 | `data/hpsolb_case_1.json` | `n = 1`, `iheap = 0`: heap-build empty, no extract |
| 2 | `data/hpsolb_case_2.json` | `n = 1`, `iheap = 1`: skip build, no extract |
| 3 | `data/hpsolb_case_3.json` | Unsorted `n = 4`, `iheap = 0`: full build + extract |
| 4 | `data/hpsolb_case_4.json` | Pre-built heap `n = 4`, `iheap = 1`: right child smaller branch |
| 5 | `data/hpsolb_case_5.json` | Sift-down stops early (root is already small) |
| 6 | `data/hpsolb_case_6.json` | `iheap = 0` with min already at position 1: build is no-op |

## Reference implementation

`reference_impl/core/hpsolb.py`.

## Cross-references

- **Paper**: Williams 1964, CACM Algorithm 232 (HEAPSORT). The cited
  reference in the F77 source.
- **Related subroutines**: called repeatedly by `cauchy` to extract
  breakpoints in order.
- **F77 source**: `src/hpsolb.f`.
- **Unit test**: `tests/unit/test_hpsolb.f90` (6 case_* blocks).
