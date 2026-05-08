# L-BFGS-B Portability Specification Pack

This directory contains a language-neutral specification of the
L-BFGS-B algorithm sufficient to produce a faithful port to any
language with a BLAS/LAPACK binding (C, C++, Java, Rust, Go, Python,
...). The pack is the contract between the canonical Fortran 77
implementation in `src/` and any future port.

A port that passes the conformance suite in this pack reproduces the
F77 reference's numerical behavior **bit-for-bit** when both are built
against the same reference BLAS/LAPACK (see `04_numerics.md` for the
exact pinning policy). Ports using a different BLAS can validate
against the looser `--tolerance` mode of the conformance runner.

## Audience

This pack is written for someone who:
- Wants to implement L-BFGS-B in another language without reading the
  Fortran source.
- Has access to a BLAS/LAPACK binding for their target language.
- Is comfortable with double-precision floating-point arithmetic and
  the basics of nonlinear optimization.

It is *not* a tutorial on the algorithm. For background, read the
papers in `docs/`:
- `algorithm.pdf` (Byrd/Lu/Nocedal/Zhu 1995): the original derivation.
- `code.pdf` (Zhu/Byrd/Nocedal 1997): ACM TOMS 778 — the algorithm as
  implemented, including line search and termination details.
- `acm-remark.pdf` (Morales/Nocedal 2011): the subsequent fix to the
  subspace minimization.

`00_overview.md` summarizes which paper covers what.

## How to read this pack

Recommended order:

1. **`00_overview.md`** — algorithm in one page, paper map.
2. **`01_glossary.md`** — notation table (paper ↔ Fortran ↔ spec).
3. **`02_api.md`** — the canonical algorithm interface (callback-based).
4. **`03_state.md`** — logical state model: what the algorithm holds in
   memory between iterations, decoupled from the F77 1D workspace
   packing.
5. **`04_numerics.md`** — tolerances, magic constants, BLAS policy.
6. **`06_portability_notes.md`** — F77 → other-language gotchas.
7. **Per-subroutine specs** in `subroutines/` — one file per
   computational primitive; read in dependency order
   (`projgr` → `bmv` → ... → `mainlb`).
8. **`07_conformance.md`** — what "passes" means per subroutine.
9. **`05_deviations.md`** — places the F77 code deviates from the
   published algorithm (and why).

If your port needs to expose a Fortran-77-compatible reverse-comm
interface (the legacy `setulb` API with `task` strings and
`isave`/`dsave` packing), see the optional **`08_legacy_reverse_comm.md`**.
Most modern ports do not need this and can skip it.

## Pack layout

```
docs/spec/
├── README.md                   # This file
├── 00_overview.md              # Algorithm summary, paper map
├── 01_glossary.md              # Paper ↔ code ↔ spec notation
├── 02_api.md                   # Abstract callback-based interface
├── 03_state.md                 # Logical state model
├── 04_numerics.md              # Tolerances, magic constants, BLAS policy
├── 05_deviations.md            # Code vs papers
├── 06_portability_notes.md     # F77 → other-language gotchas
├── 07_conformance.md           # What "passes" per subroutine
├── 08_legacy_reverse_comm.md   # OPTIONAL: F77 setulb compatibility
├── subroutines/                # 16 per-subroutine specs
├── data/                       # Language-neutral test vectors (JSON)
├── runner/                     # Reference conformance test runner
└── reference_impl/             # Python+NumPy reference implementation
```

## Conformance

A port is **conformant** if it passes the conformance suite in `runner/`
in `--strict` mode against the reference BLAS pin. There is also a
looser `--tolerance` mode for ports that must use a non-reference BLAS;
see `07_conformance.md` for per-subroutine tolerance values.

The conformance suite consumes the JSON test vectors in `data/`; each
vector exercises one branch of one subroutine, matching one `case_*`
block in `tests/unit/test_*.f90`.

## Scope

**In scope** (16 subroutines, all numerical):

| Tier | Subroutines |
|------|-------------|
| 1 — utilities | `projgr`, `errclb`, `active`, `freev`, `hpsolb`, `cmprlb` |
| 2 — math primitives | `bmv`, `formt`, `matupd`, `dcstep`, `formk` |
| 3 — line search | `dcsrch`, `lnsrlb` |
| 4 — algorithm core | `cauchy`, `subsm` |
| 5 — algorithm loop | `mainlb` |

**Out of scope** (port-defined behavior):

- `timer.f` — wallclock helper. Ports may use any monotonic clock.
- `prn1lb.f`, `prn2lb.f`, `prn3lb.f` — diagnostic printing. Ports may
  emit any human-readable diagnostics format.
- The reverse-comm wrapper `setulb.f` — captured optionally in
  `08_legacy_reverse_comm.md` for F77 ABI compatibility only.

## Getting help

If something in this pack is ambiguous, contradicts the F77 source, or
fails to specify behavior you encounter — **the spec is wrong**. Open
an issue or PR against this repository; ambiguity in the spec is a
defect.
