# Porting from Fortran 77: practical gotchas

The L-BFGS-B reference implementation in `src/` is Fortran 77 (with
Doxygen comment markup). Port authors reading the F77 source
alongside the spec — to verify behavior or to confirm operation order
— need to be aware of F77 conventions that differ from C-family
languages.

This document is **not** a Fortran tutorial. It catalogues the
specific traps that have bitten ports of similar codes. Each section
flags whether the convention matters for *spec conformance* (i.e., a
port must reproduce it) or only for *reading* the F77 source
correctly.

## Array memory layout: column-major

Fortran arrays are **column-major** (Fortran-natural / "first index
varies fastest"). C, C++, Java, Rust, Go, NumPy default-row-major.

For a 2D array `A(i, j)` declared `double precision A(M, N)`:
- F77 in-memory order: `A(1,1), A(2,1), ..., A(M,1), A(1,2), ...`
- C-equivalent declaration `double A[N][M]`: in-memory order
  `A[0][0], A[0][1], ..., A[0][M-1], A[1][0], ...` — **different**.

For the C-equivalent layout (so the ports can index naturally), the
declaration order is *reversed*: F77 `A(M, N)` ↔ C `A[N][M]`.

**Affects spec conformance**:
- Operation order in inner loops sometimes depends on memory layout
  (e.g., a column-by-column accumulation reduces to `daxpy` calls in
  F77). Per-subroutine specs document the order independent of layout
  choice.
- BLAS/LAPACK calls assume column-major. Ports using a row-major BLAS
  binding must transpose accordingly. The standard CBLAS bindings
  expose a layout flag (`CblasColMajor`/`CblasRowMajor`) — ports must
  pick `CblasColMajor` to match the F77 reference.

**For reading F77 source**: when you see `A(i, j)`, mentally translate
to "row `i`, column `j`, with column-`j` contiguous in memory". Inner
loops over `i` are the cache-friendly direction in F77.

## 1-based indexing

F77 arrays default to 1-based indices: `A(1)` is the first element,
`A(N)` is the last. Loops written as `do i = 1, N` are inclusive on
both ends.

**Spec convention**: this pack uses 1-based indexing in pseudocode
to match the F77 source for ease of cross-reference. Ports in 0-based
languages (C, Python, Rust, ...) translate to 0-based in the obvious
way. Where the index value matters numerically (rare — usually only
for diagnostic output), the spec calls it out.

## Character strings: `character*60` and trailing-space padding

F77 strings are fixed-length character arrays. `character*60 task`
declares a 60-character buffer; assignments shorter than 60 characters
are **right-padded with spaces**. Comparisons (`task .eq. 'FG'`) match
on the full padded buffer, so `'FG' || 58 spaces == 'FG'` (with
`||` denoting concatenation).

The L-BFGS-B convention: the first 1–8 non-space characters are the
state code; the remainder is human-readable diagnostic text:

```
'CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL          '
 ^^^^^^^^^^^                                                 (first 11 = state code)
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   (rest = description)
```

In F77 callers, the convention is to test only a prefix:
`task(1:2) .eq. 'FG'`, `task(1:5) .eq. 'NEW_X'`, etc.

**For ports**:
- Use the language's native string type. Right-padding is a F77
  artifact, not part of the algorithm.
- The state codes themselves matter only for the legacy reverse-comm
  interface (see `08_legacy_reverse_comm.md`). Callback-based ports
  use `info` codes instead.

## Logical type: `.true.` / `.false.`

F77 logicals (`logical :: prjctd`) are not integers, but in-memory
representations vary by compiler (gfortran uses `1` for true, `0` for
false, in a 4-byte slot). Comparisons use `.eqv.` / `.neqv.`, not
`==`.

**For ports**: use the language's native boolean. No conformance
implications.

## `save` semantics

F77 local variables can be declared `save`, which gives them static
storage with a value persisting across calls. L-BFGS-B uses `save` in
`dcsrch` and a few utility routines.

**For ports**: the reverse-comm `setulb` interface relies on `save`
implicitly via the `isave`/`dsave` arrays passed in by the caller.
Callback-based ports do not need static storage; the algorithm's
state lives in normal local variables of the running optimization
function.

In the F77 source, `save` variables are sometimes initialized only
once (via `data` statements) and then mutated; ports reading the
source should track this carefully. Per-subroutine specs flag any
relevant `save` use.

## `implicit none` *(or absence thereof)*

F77 has implicit typing: a variable starting with `i`, `j`, `k`, `l`,
`m`, or `n` is integer by default; otherwise real. The L-BFGS-B
sources **do not** declare `implicit none` (though we now compile with
`-fimplicit-none`, see `CMakeLists.txt`).

This means a typo in a F77 variable name silently introduces a new
variable rather than producing a compiler error. Two such bugs have
been caught in this codebase historically (the `bmv one` bug, the
`test_cmprlb theta` bug). Reading the F77 source, double-check that
every variable you see is declared somewhere.

**For ports**: declare all variables. This is automatic in any modern
language.

## Integer divide

F77 integer division truncates toward zero (`5 / 2 = 2`, `-5 / 2 = -2`).
Most modern languages match this; Python's `/` on integers used to do
true division (Python 2) and now does true division everywhere
(Python 3 — use `//` for truncating).

**For ports**: use the language's truncating-divide operator (or floor
divide for non-negative operands; the difference only matters for
negative operands, which L-BFGS-B does not produce in indexing
contexts).

## Logical operators: no short-circuit evaluation

F77's `.and.` and `.or.` are **not short-circuit**. Both operands are
evaluated.

```fortran
if (i .ge. 1 .and. a(i) .gt. 0) ...   ! evaluates a(i) even if i < 1!
```

This rarely produces visible differences in L-BFGS-B because the F77
source is structured to avoid the issue (with explicit nested `if`),
but porters reading the source should not assume short-circuit
behavior even where intermediate inputs would be invalid.

**For ports**: use the language's short-circuit operators (`&&`, `||`,
`and`, `or`); this is safer and produces equivalent behavior in
L-BFGS-B's actual code paths (which are structured for non-SC).

## `do` loops with start > end

F77 `do i = start, end, step` where `step > 0` and `start > end`
**executes zero iterations** (compatible with C-family `for` loops on
this point).

```fortran
do i = 5, 4
   ...           ! loop body NOT executed
end do
```

This is the same as `for(i=5; i<=4; i++)` in C — body never runs. No
conformance gotcha; mentioned only because some ancient F77 dialects
behaved differently.

## `goto` and arithmetic-`if`

F77 makes liberal use of `goto`. L-BFGS-B has many: `cauchy.f` has
~10 labels, `mainlb.f` has ~20. Most are forward gotos that simulate
`break` / `continue` of inner loops.

The Doxygen-doc'd `c>` comments and structured indentation usually
clarify intent. When in doubt, the per-subroutine specs in
`subroutines/` describe the control flow in language-neutral terms.

There are no `assigned goto` or `computed goto` constructs in
L-BFGS-B, only labeled gotos.

## `common` blocks

L-BFGS-B has **no `common` blocks**; all state is passed explicitly
through arguments. Ports do not need to worry about this F77 idiom.

(Some older numerical Fortran code used `common` for shared storage;
L-BFGS-B avoids this, which made the F77-to-callback translation in
`mainlb` straightforward.)

## `equivalence` declarations

L-BFGS-B has **no `equivalence`** statements (which let two F77
variables share the same memory location). Ports do not need to model
aliasing.

## `data` statement initialization

F77 `data` initializes variables once before execution. Used in
L-BFGS-B for constants like `parameter (zero=0.0d0)`. Any modern
language's `const` / `static` initialization is equivalent.

## Reverse communication (covered separately)

The most language-specific F77 idiom in L-BFGS-B is the
reverse-communication protocol of `setulb`. This is captured in
`08_legacy_reverse_comm.md` as an optional appendix. It is **not part
of the canonical algorithm spec**; ports targeting modern languages
should use callbacks instead (see `02_api.md`).

## Source-file format: fixed-form columns

The F77 sources in `src/*.f` use **fixed-form Fortran**: columns 1–5
are reserved for labels, column 6 for continuation marker, columns
7–72 for code, columns beyond 72 ignored. The newer drivers in
`drivers/*.f90` use free-form.

Port authors do not need to worry about column conventions, but
should be aware that:
- A line continuation in F77 fixed form is a non-blank, non-zero
  character in column 6 (the `+` characters in `mainlb.f`).
- Long F77 expressions span multiple lines via this mechanism;
  reading them requires recognizing that each `+` line continues the
  prior.

## Doxygen `c>` markup

The L-BFGS-B sources use Doxygen-aware F77 comments: lines starting
with `c>` (or `!>` in F90 sources) are treated as Doxygen
documentation. Most subroutines have an `@param` block at the top
documenting each argument, and a brief `@brief` description. These
are the canonical per-parameter docs; the spec's per-subroutine files
in `subroutines/` are the canonical *behavioral* docs.

## Compiler flags used by this repo

The project's `CMakeLists.txt` sets, for `gfortran`:

```
-fimplicit-none -Wall -Wextra -Wno-compare-reals
```

Plus, when `LBFGSB_COVERAGE=ON`:

```
--coverage -O0 -g
```

The `-Wno-compare-reals` is intentional: L-BFGS-B legitimately uses
exact floating-point comparison (`==` / `/=`) on real values to detect
algorithmic states (e.g., variable exactly at a bound, exact unit
step). This is checked-and-correct usage; the compiler warning is a
stylistic check that does not apply here. Per-subroutine specs flag
each such comparison and explain why it's correct.

For bit-for-bit conformance, additionally:

```
-O2 -fno-fast-math -fno-associative-math -frounding-math
```

(see `04_numerics.md`).
