# Appendix: F77 reverse-communication interface

> **This appendix is OPTIONAL.** It is needed only by ports that want
> to expose the F77 `setulb` ABI for compatibility with existing
> callers. Ports targeting modern languages should use the
> callback-based interface in `02_api.md` instead.

The F77 reference implementation exposes the optimizer through a
*reverse-communication* protocol: instead of the user passing function
and gradient callbacks, the user calls `setulb` repeatedly; each call
returns control with a `task` string indicating what the user should
do next (typically: "evaluate `f` and `g` at the current `x`, then
call back").

Reverse communication exists because Fortran 77 lacks first-class
function values. Modern languages (including Fortran 2003 and later)
have closures, so most ports do not need this protocol. Use this
appendix only if:

- Existing F77 / scientific code calls `setulb` and you want a
  drop-in replacement, or
- Your language has restricted callback support (some embedded /
  RT-Java environments, certain FPGA toolchains).

## Protocol overview

### State diagram

```
                                       +------+
                                       v      |
   START --> (compute f,g) --> FG -->  user  |
              [task='FG'              evaluates
               first time]            f, g    |
                                       |      |
                                       v      |
                                    NEW_X ----+   (outer iter complete;
                                       |           user may inspect, then
                                       |           re-call setulb)
                                       |
                                       v
                  CONV* / ABNORMAL* / ERROR* / STOP*  (terminal)
```

### Task strings

`task` is a `character*60` buffer. The first 1-8 non-space characters
are the state code; the rest is human-readable description. Test
state with prefix matching: `task(1:2) == 'FG'`, `task(1:5) == 'NEW_X'`,
etc.

| `task` prefix | When set | Meaning |
|---------------|----------|---------|
| `START` | By user before first call | Initialize a new optimization. |
| `FG` | By `setulb` on return | User must compute `f` and `g` at the current `x` and call `setulb` again. |
| `NEW_X` | By `setulb` on return | A new iterate is available in `x`; the user may inspect, log, etc., then call `setulb` again to continue iterating. |
| `CONV...` | By `setulb` on return | Convergence reached. Specific suffix indicates which test fired. Terminal. |
| `ABNORMAL...` | By `setulb` on return | Line search failed or other abnormal termination. Terminal. |
| `ERROR...` | By `setulb` on return | Input validation failed. Terminal. |
| `STOP` | By user any time | Request immediate termination. `setulb` will tidy up and return with a terminal `task`. |

Specific terminal task strings:

```
'CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL'
'CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH'
'ABNORMAL TERMINATION IN LNSRLB'
'ERROR: N <= 0'
'ERROR: M <= 0'
'ERROR: FACTR < 0'
'ERROR: INVALID NBD'
'ERROR: NO FEASIBLE SOLUTION'
'STOP: CPU EXCEEDING THE TIME LIMIT.'
```

(See `subroutines/errclb.md` for the full input-validation case list.)

## Calling sequence (worked example)

```fortran
character*60 :: task, csave
logical      :: lsave(4)
integer      :: isave(44), iwa(3*n)
real(8)      :: dsave(29), wa(2*m*n + 5*n + 11*m*m + 8*m)
real(8)      :: x(n), l(n), u(n), g(n), f, factr, pgtol
integer      :: nbd(n), iprint

! Initial setup
x = ...        ! initial point
l = ...; u = ...; nbd = ...
factr = 1.0d7
pgtol = 1.0d-5
iprint = -1
task = 'START'

do
  call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, &
              wa, iwa, task, iprint, csave, lsave, isave, dsave)
  if (task(1:2) == 'FG') then
    f = ...      ! evaluate f(x)
    g = ...      ! evaluate gradient at x
    cycle        ! call setulb again with the new f, g
  else if (task(1:5) == 'NEW_X') then
    ! a new iterate is in x; optionally inspect, then continue
    cycle
  else
    exit         ! terminal task: CONV / ABNORMAL / ERROR / STOP
  end if
end do

! task contains the termination message; isave/dsave have stats
```

## Workspace arrays

`setulb` requires the user to provide all storage; nothing is
allocated internally. Sizes:

| Array | Type | Size | Purpose |
|-------|------|------|---------|
| `wa` | `real(8)` | `2*m*n + 5*n + 11*m^2 + 8*m` | All double-precision working storage. |
| `iwa` | `integer` | `3*n` | All integer working storage (`index`, `iwhere`, `indx2`). |
| `task` | `character*60` | -- | Reverse-comm state code. |
| `csave` | `character*60` | -- | Inner reverse-comm state for `dcsrch`. |
| `lsave` | `logical` | `4` | `prjctd`, `cnstnd`, `boxed`, `updatd`. |
| `isave` | `integer` | `44` | Saved integer state across calls. |
| `dsave` | `real(8)` | `29` | Saved real state across calls. |

The user must not modify any of these between calls except `f` and
`g` when `task = 'FG'`, and `task` itself for `'START'` or `'STOP'`.

## `wa` partitioning

Set during the `'START'` call by `setulb`; remains stable across all
subsequent calls. The first 16 slots of `isave` hold the offsets:

| `isave[k]` | Stores | Logical object | F77 name |
|-----------|--------|----------------|----------|
| `1` | `m*n` | (size constant) | -- |
| `2` | `m^2` | (size constant) | -- |
| `3` | `4*m^2` | (size constant) | -- |
| `4` | offset | `S` matrix (history of `s`-vectors) | `ws` |
| `5` | offset | `Y` matrix (history of `y`-vectors) | `wy` |
| `6` | offset | `sy` (`S^TY` packed: `D` + `L`) | `sy` |
| `7` | offset | `ss` (`S^TS` upper triangle) | `ss` |
| `8` | offset | `T` (Cholesky factor) | `wt` |
| `9` | offset | `K` (`2m * 2m` reduced Hessian) | `wn` |
| `10` | offset | `K_chol` (Cholesky of `K`) | `snd` |
| `11` | offset | `xc` (Cauchy point) | `z` |
| `12` | offset | `r` (reduced gradient) | `r` |
| `13` | offset | `d` (search direction) | `d` |
| `14` | offset | `t` (scratch) | `t` |
| `15` | offset | `xp` (previous iterate) | `xp` |
| `16` | offset | scratch (8*m) | `wa` (inner) |

All offsets are 1-based positions within the user's `wa` array.

## `iwa` partitioning

Set deterministically by `setulb`'s call to `mainlb`:

| `iwa[k]` | Logical object |
|---------|-----------------|
| `iwa[1..n]` | `index` (free-variable permutation) |
| `iwa[n+1..2n]` | `iwhere` (per-variable bound status) |
| `iwa[2n+1..3n]` | `indx2` (entering/leaving lists) |

## `lsave` semantics

Available on exit with `task = 'NEW_X'`:

| `lsave[k]` | Meaning |
|-----------|---------|
| `1` | `prjctd`: initial `x0` was infeasible and projected. |
| `2` | `cnstnd`: at least one variable has a bound. |
| `3` | `boxed`: every variable has both lower and upper bounds. |
| `4` | `updatd`: L-BFGS history was updated in this iteration. |

## `isave` saved state

Beyond the workspace offsets in `1..16`, `isave` carries diagnostic
counters and saved internal state. The slots documented for users
(available on exit with `task = 'NEW_X'`):

| `isave[k]` | Meaning |
|-----------|---------|
| `22` | Total intervals explored across all Cauchy point searches. |
| `26` | Total skipped BFGS updates before the current iteration. |
| `30` | Iteration counter (`iter`). |
| `31` | Total BFGS updates *prior to* the current iteration. |
| `33` | Cauchy intervals explored in the current iteration. |
| `34` | Total `f` and `g` evaluations (cumulative). |
| `36` | `f`/`g` evaluations in the current iteration's line search. |
| `37` | Subspace-argmin location flag: `0` if within box, `1` if beyond box (clipped). |
| `38` | Number of free variables in the current iteration. |
| `39` | Number of active constraints in the current iteration. |
| `40` | Sentinel `n + 1 - isave[40]` = number of variables that *left* the active set this iteration. |
| `41` | Number of variables that *entered* the active set this iteration. |

Slots 17-21, 23-25, 27-29, 32, 35, 42-44 are private to `mainlb` and
its callees (saving control-flow state for resume-after-`'FG'`); the
user must not read or modify them.

## `dsave` saved state

Available on exit with `task = 'NEW_X'`:

| `dsave[k]` | Meaning |
|-----------|---------|
| `1` | Current `theta`. |
| `2` | `f` at the previous iterate. |
| `3` | `factr * epsmch`. |
| `4` | 2-norm of the search direction `d`. |
| `5` | `epsmch` as computed by `mainlb`. |
| `7` | Accumulated wallclock spent on Cauchy point. |
| `8` | Accumulated wallclock on subspace minimization. |
| `9` | Accumulated wallclock on line search. |
| `11` | Slope of `f` along `d` at the current line-search point. |
| `12` | Maximum relative step length imposed by the box geometry. |
| `13` | `||g_proj||_Inf`. |
| `14` | Relative step length accepted in the line search. |
| `15` | Slope of `f` along `d` at the start of the line search. |
| `16` | Square of the 2-norm of `d`. |

Slots 6, 10, 17-29 are private to internal routines.

## `csave` saved state

Holds the inner reverse-comm state for the More-Thuente line search
(`dcsrch`). The user must not modify it. When the line search needs
another `f`/`g` evaluation, the outer state machine in `mainlb` sets
`task = 'FG'` and the user provides values; `dcsrch` resumes from the
state captured in `csave`.

## Implementation notes for ports

A port that exposes this protocol must:

1. **On `task = 'START'`**: validate inputs (per `errclb.md`), partition
   `wa`/`iwa` (or whatever storage equivalent the port uses), do the
   initial projection, and yield with `task = 'FG'` requesting the
   first `f`/`g`.
2. **On subsequent calls**: dispatch on the previous `task` value (or
   on a hidden program-counter saved in `isave`) to resume the
   correct internal location. The F77 source uses a `goto`-based
   state machine; ports without `goto` typically use a switch on a
   saved enum.
3. **Preserve all working state** in `isave`/`dsave`/`lsave`/`csave`
   between calls. The user is allowed to inspect the documented
   slots; the port must not move them.
4. **On `task = 'STOP'`**: the user can set this any time. The port
   tidies up (no work in progress to abort cleanly) and returns with
   `task = 'STOP: <reason>'` as a terminal state.

## Mapping back to the callback interface

A port that already implements the canonical callback API in
`02_api.md` can wrap it in a reverse-comm shim by:

1. Running the callback-based optimizer **on a separate stack /
   coroutine / thread** that yields when the callback is invoked.
2. The shim's `setulb`-equivalent function calls `resume()` on the
   coroutine; when the coroutine yields with the trial `x`, the shim
   sets `task = 'FG'` and returns to the user.
3. On the next call, the shim writes the user-supplied `f` and `g`
   into the coroutine's communication slot and `resume()`s again.

Most languages support this shape via:
- Python: generators (`yield`) or `asyncio` coroutines.
- C/C++: `setjmp`/`longjmp`, fibers, or `std::generator` (C++23).
- Rust: `Generator` trait or hand-rolled state machine.
- Go: goroutines + channels.
- Java: `Thread.suspend`/`resume` (deprecated; use a worker thread
  with synchronization), or virtual threads (Loom).

Ports that need the legacy ABI but can't use coroutines must
implement the protocol natively as the F77 source does -- saving
program-counter and local state in `isave`/`dsave`. This is more
work than the coroutine wrapper but has no extra runtime cost.

## Testing

The conformance suite includes JSON test vectors that exercise the
reverse-comm protocol end-to-end (driving the optimizer through full
runs of the integration drivers). These vectors test the legacy
interface only; ports that don't implement this appendix can skip
them.
