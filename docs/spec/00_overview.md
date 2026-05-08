# Algorithm overview

L-BFGS-B minimizes a smooth function `f : ℝⁿ → ℝ` subject to box bounds
`l ≤ x ≤ u`, using a limited-memory quasi-Newton approximation to the
Hessian and a projected gradient framework for the bounds. Per
iteration:

1. **Identify the active bound set.** Variables that the search
   direction would push into a bound are fixed at that bound; the rest
   are *free*.
2. **Compute the generalized Cauchy point** `x_c`: the minimizer of a
   quadratic model of `f` along the projected gradient direction. This
   determines which bounds become active in this iteration.
3. **Subspace minimization.** Approximately solve the box-constrained
   quadratic subproblem on the free variables, using the L-BFGS Hessian
   approximation. This yields a search direction `d`.
4. **Line search** along `d` (More–Thuente Wolfe-condition line
   search), projected onto the box. The accepted step `α*d` produces
   the next iterate `x_{k+1}`.
5. **Update the L-BFGS history**: append `(s_k, y_k) = (x_{k+1} - x_k,
   g_{k+1} - g_k)` to the limited-memory pair list, dropping the
   oldest pair if the memory bound `m` is reached.
6. **Test convergence** on `‖projected gradient‖_∞ < pgtol` and on
   relative function decrease `(f_k - f_{k+1}) / max(|f_k|, |f_{k+1}|, 1)
   < factr * epsmch`.

The algorithm runs until convergence, an iteration/evaluation limit
is hit, or the line search fails (which can indicate the function is
unbounded below, the gradient is inconsistent with the function, or
the problem is too ill-conditioned).

The "limited memory" part: instead of storing the `n×n` Hessian
approximation, we store the last `m` pairs `(s_i, y_i)`, with `m`
typically 5–17. The compact representation
(Byrd/Nocedal/Schnabel 1994) lets us compute Hessian-vector products
in `O(mn)` time without ever materializing the dense Hessian.

## Paper map

L-BFGS-B is described across three papers. **`code.pdf` is the
authoritative reference** — it documents the algorithm as actually
implemented in this codebase. The others give derivations or fixes.

### `algorithm.pdf` — Byrd/Lu/Nocedal/Zhu, "A Limited Memory Algorithm for Bound Constrained Optimization" (SIAM J. Sci. Comput. 1995)

The original algorithm derivation. Read this for:
- The mathematical justification for the active-set strategy.
- The generalized Cauchy point definition (eq 4.6).
- The compact L-BFGS representation
  (also covered in Byrd/Nocedal/Schnabel 1994).
- The original subspace minimization scheme.

**Caveat:** the implemented code follows the variant in `code.pdf`,
not this paper. The paper states this explicitly: *"the algorithm we
present below differs in a few important details from the one
originally described."* See `05_deviations.md`.

### `code.pdf` — Zhu/Byrd/Nocedal, "L-BFGS-B: Algorithm 778" (ACM TOMS 1997)

The algorithm as implemented. **This is the reference for the spec.**
Covers:
- The iteration as actually coded, including the line search wrapper
  (`lnsrlb`), the curvature test, and the refresh logic.
- The reverse-comm interface (`setulb`).
- Termination tests (`factr`, `pgtol`).
- Workspace-array layout (relevant only for the F77 implementation;
  ports do not need to follow this).

### `acm-remark.pdf` — Morales/Nocedal, "Remark on Algorithm 778" (ACM TOMS 2011)

A correction to the subspace minimization in `subsm` (the
"Cauchy-direction" variant). The original 1997 code had a bug in how
it handled bounds during the subspace solve; this paper documents the
fix. The fix is included in `src/subsm.f`.

## Reading order if you are porting

1. Read `algorithm.pdf` §1–§5 for intuition (1–2 hours).
2. Skim `code.pdf` §1–§3 for the coded algorithm and the termination
   tests (30 minutes).
3. Skim `acm-remark.pdf` (15 minutes) — only matters for `subsm`.
4. Then read this pack starting from `01_glossary.md`.

You do not need to read `code.pdf` §4 (Fortran implementation details)
unless you want to mimic the F77 reverse-comm interface, in which case
also read this pack's `08_legacy_reverse_comm.md`.

## Algorithm complexity per iteration

For dimension `n` and memory `m` (typically `n ≫ m`):

| Step | Work |
|------|------|
| Identify active set (`active`/`cauchy`) | `O(n + m²)` |
| Cauchy point breakpoint sort (`hpsolb`) | `O(n log n)` worst case |
| Subspace minimization (`subsm`) | `O(m² n)` |
| Line search (`lnsrlb`) | `O(n)` per evaluation, plus user f/g cost |
| L-BFGS update (`matupd`, `formt`) | `O(m² + mn)` |

Memory: `O(mn + m²)`. The F77 `wa` workspace is sized `2mn + 5n + 11m² + 8m`
double-precision values; the `iwa` integer workspace is `3n` integers.
Ports may use whatever data structures fit their language as long as
the per-iteration math matches.
