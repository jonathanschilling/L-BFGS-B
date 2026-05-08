# Glossary: notation across paper, code, and spec

This document maps the notation used in the three reference papers
(`algorithm.pdf`, `code.pdf`, `acm-remark.pdf`) to the F77 identifiers
in `src/` and the spec notation used throughout this pack.

When papers disagree on notation, the **Spec** column gives the
canonical name used here; per-subroutine specs in `subroutines/` use
this notation.

## Scalars and counters

| Spec | Paper | F77 | Description |
|------|-------|-----|-------------|
| `n` | `n` | `n` | Number of variables (problem dimension). |
| `m` | `m` | `m` | Maximum number of L-BFGS pairs to retain. Typical 5-17. |
| `col` | `k` (alg) / `m_k` (code) | `col` | Current number of stored L-BFGS pairs (`0 <= col <= m`). Increases each iteration until reaching `m`, then stays at `m`. |
| `iter` | `k` | `iter` | Iteration counter (0-based; `iter=0` is the initial point). |
| `theta` | `theta` | `theta` | Hessian scaling parameter; `theta = (y_k^T y_k) / (s_k^T y_k)` after each successful update. Initialized to 1. |
| `f` | `f(x)` | `f` | Function value at current iterate. |
| `factr` | `factr` | `factr` | Function-decrease tolerance (relative). Convergence when `(f_old - f_new) / max(\|f_old\|, \|f_new\|, 1) <= factr * eps`. Typical `10^7` (low accuracy), `10^1` (machine accuracy). |
| `pgtol` | `pgtol` | `pgtol` | Projected-gradient tolerance. Convergence when `||g_proj||_Inf <= pgtol`. |
| `eps` / `epsmch` | `eps_M` | `epsmch` | Machine epsilon (~= `2.22 * 10^{-16}` for IEEE-754 double). |

## Vectors (current state)

| Spec | Paper | F77 | Shape | Description |
|------|-------|-----|-------|-------------|
| `x` | `x_k` | `x` | `n` | Current iterate. |
| `g` | `g_k = gradf(x_k)` | `g` | `n` | Gradient at current iterate. |
| `g_proj` | `\|P[x - g, l, u] - x\|_Inf` | (computed in `projgr`) | scalar | Infinity norm of projected gradient. |
| `s` | `s_k = x_{k+1} - x_k` | (incremental) | `n` | Step from previous iterate. |
| `y` | `y_k = g_{k+1} - g_k` | (incremental) | `n` | Gradient difference. |
| `d` | `d_k` | `d` | `n` | Search direction (subspace minimizer). |
| `xc` | `x^c` | `z` | `n` | Generalized Cauchy point. |
| `xp` | `x_{k-1}` | `xp` | `n` | Previous iterate (kept for refresh logic). |
| `r` | (varies) | `r` | `n` | Reduced gradient at the Cauchy point (used by `subsm`). |
| `t` | -- | `t` | `n` | Scratch vector inside `mainlb`. |

## Bounds

| Spec | F77 | Description |
|------|-----|-------------|
| `l` | `l` | Lower-bound vector (length `n`). Entries irrelevant where `nbd(i) in {0, 3}`. |
| `u` | `u` | Upper-bound vector (length `n`). Entries irrelevant where `nbd(i) in {0, 1}`. |
| `nbd[i]` | `nbd(i)` | Bound code per variable: `0` = unbounded, `1` = lower only, `2` = both, `3` = upper only. |

## Bound-status flags (per variable)

After `active` runs once, then refined by `cauchy` each iteration.

| Value | Meaning |
|-------|---------|
| `iwhere[i] = -1` | Variable has no bounds (`nbd[i] = 0`); always free. |
| `iwhere[i] = 0` | Has bounds, currently free (interior). |
| `iwhere[i] = 1` | Currently at the lower bound. |
| `iwhere[i] = 2` | Currently at the upper bound. |
| `iwhere[i] = 3` | Variable is fixed (`nbd[i] = 2` and `l[i] = u[i]`). |

## L-BFGS history

The algorithm retains the last `col` step pairs `(s_i, y_i)` where
indexing is by *insertion order*. The most recent pair is at column
`itail` of `S` and `Y`; the oldest is at column `head`.

| Spec | Paper | F77 | Shape | Description |
|------|-------|-----|-------|-------------|
| `S` | `S_k = [s_{k-col}, ..., s_{k-1}]` | `ws` | `n * col` | Matrix of `s`-vectors. Columns ordered by insertion; F77 uses `head`/`itail` to track cyclic position when `col = m`. |
| `Y` | `Y_k` | `wy` | `n * col` | Matrix of `y`-vectors, same ordering as `S`. |
| `D` | `D_k = diag(s_i^T y_i)` | (diagonal of `sy`) | `col * col` | Diagonal of dot products. **Stored as the diagonal of `sy`**. |
| `L` | `L_k` (strictly lower triangular of `S_k^T Y_k`) | (lower triangle of `sy`) | `col * col` | Strict lower triangle: `L_{i,j} = s_i^T y_j` for `i > j`. **Stored as the strict lower triangle of `sy`**. |
| `R` | `R_k` (upper triangular of `S_k^T Y_k`) | not stored | `col * col` | Used in derivations; not materialized in the F77 code (the relevant info is in `D` + `L`). |
| `S^TS` | `S_k^T S_k` | `ss` | `col * col` | Symmetric Gram matrix. F77 stores the upper triangle in `ss`. |

### F77 packing of `sy` and `ss`

The F77 implementation packs both `D` and `L` into a single `m * m`
array `sy`:
- `sy(i, i) = D[i] = s_i^T y_i` (diagonal)
- `sy(i, j) = L[i,j] = s_i^T y_j` for `i > j` (strict lower triangle)
- `sy(i, j)` for `i < j` is **unused** (held but not read).

Similarly `ss(i, j)` for `i <= j` holds `s_i^T s_j` (upper triangle of
`S^TS`); the strict lower part is unused.

**Ports may use any storage they prefer** -- separate `D` (vector) and
`L` (lower-triangular matrix) variables, a full symmetric `S^TS` matrix,
sparse formats, etc. -- as long as the entries above are accessible.
The packing above is an F77 implementation detail.

## Compact representation

Following Byrd/Nocedal/Schnabel 1994 sec.3 and `code.pdf` sec.2.

| Spec | Paper | F77 | Shape | Description |
|------|-------|-----|-------|-------------|
| `W` | `W_k = [Y_k, theta*S_k]` | (formed implicitly) | `n * 2col` | Wide matrix. Not stored explicitly in F77; columns are taken on the fly from `wy` and `ws`. |
| `M_inv` | `M_k^{-1} = [[-D, L^T], [L, theta*S^TS]]` | (formed implicitly via `sy`, `ss`, `theta`) | `2col * 2col` | Middle matrix inverse. F77 never materializes this; `bmv` solves `M_k v = result` using the components. |
| `T` | `J_k`: upper Cholesky factor of `theta*S^TS + L D^{-1} L^T` | `wt` | `col * col` | Upper-triangular Cholesky factor. Computed in `formt` from `sy`, `ss`, `theta`. Used by `bmv` for the matrix-vector solve. |
| `B` | `B_k = thetaI - W M^{-1} W^T` | (not stored) | `n * n` | Hessian approximation. Never materialized; only Hessian-vector products are computed. |
| `K` | `K = [[-D - Y^T_a Z B Z^T Y_a, ...], ...]` (subspace minimization) | `wn` | `2col * 2col` | Reduced Hessian for the subspace minimization. Built by `formk`, used by `subsm`. |
| `K_chol` | upper Cholesky of `K` | `snd` | `2col * 2col` | Cholesky factor of `K` for back-substitution in `subsm`. |

## F77 workspace arrays (mapped to logical objects)

Ports do **not** need to use a single contiguous workspace -- these
mappings exist so that readers of the F77 source can locate the
logical objects within `wa` / `iwa`.

`wa` is the double-precision workspace, sized `2mn + 5n + 11m^2 + 8m`.
The F77 driver layout (see `setulb.f`):

| F77 slice | Spec object | Shape |
|-----------|-------------|-------|
| `wa[1:mn]` | `S` flattened (column-major) | `n * m` |
| `wa[mn+1:2mn]` | `Y` flattened | `n * m` |
| `wa[2mn+1:2mn+m^2]` | `sy` packed | `m * m` |
| `wa[2mn+m^2+1:2mn+2m^2]` | `ss` packed | `m * m` |
| `wa[2mn+2m^2+1:2mn+3m^2]` | `T` (Cholesky factor) | `m * m` |
| `wa[2mn+3m^2+1:2mn+7m^2]` | `K` workspace | `2m * 2m` |
| `wa[2mn+7m^2+1:2mn+11m^2]` | `K_chol` workspace | `2m * 2m` |
| `wa[2mn+11m^2+1:2mn+11m^2+5n]` | working vectors `z, r, d, t, xp` | `5 * n` |
| `wa[remaining]` | scratch for `cauchy`, `bmv`, etc. | `8m` |

`iwa` is the integer workspace, sized `3n`:

| F77 slice | Spec object | Shape |
|-----------|-------------|-------|
| `iwa[1:n]` | `index` (free-variable permutation) | `n` |
| `iwa[n+1:2n]` | `iwhere` after `active`/`cauchy` | `n` |
| `iwa[2n+1:3n]` | `indx2` (entering/leaving) | `n` |

## Reverse-comm state arrays *(legacy interface only)*

These are how the F77 `setulb` interface preserves state across calls.
They are **not part of the canonical algorithm spec**; ports using a
callback-based interface have no analogue.

| F77 | Type | Description |
|-----|------|-------------|
| `task` | `character(60)` | Reverse-comm state: `&#39;START&#39;`, `&#39;FG&#39;`, `&#39;NEW_X&#39;`, `&#39;CONV...&#39;`, `&#39;ABNORMAL&#39;`, `&#39;ERROR&#39;`, `&#39;STOP&#39;`. |
| `csave` | `character(60)` | Inner reverse-comm state for `dcsrch` line search. |
| `lsave[1:4]` | `logical(4)` | `prjctd`, `cnstnd`, `boxed`, `updatd` flags preserved across calls. |
| `isave[1:44]` | `integer(44)` | Saved integer state (workspace offsets, iter counters, indices). |
| `dsave[1:29]` | `real(29)` | Saved real state (line-search bracket, theta, sbgnrm, etc.). |

For details on how these reconstruct the algorithm state, see
`08_legacy_reverse_comm.md`.

## Termination codes

The F77 `task` string on exit:

| `task` | Meaning |
|--------|---------|
| `CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL` | `||g_proj||_Inf <= pgtol`. |
| `CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH` | Function decrease below `factr * eps`. |
| `ABNORMAL TERMINATION IN LNSRLB` | Line search failed (more than `maxls` evaluations, or step too small/large). |
| `ERROR: ...` | Input validation failed (see `errclb.md`). |
| `STOP: ...` | User-requested halt via `task = &#39;STOP&#39;`. |

Ports may use enum values, exceptions, or other idiomatic mechanisms;
the strings above are the F77 representation.

## Index conventions in this pack

Spec pseudocode is **1-based** to match Fortran source for ease of
cross-reference. Ports in 0-based languages must adjust accordingly.
Where the distinction matters numerically (e.g., a loop bound on a
`break` condition), the spec calls it out explicitly.
