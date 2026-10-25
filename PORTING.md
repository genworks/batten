# Porting NLib to Batten

Batten is a routine-by-routine translation of NLib, from NVIDIA's open-source
SMLib (usd-brep `source/SMLib`, Apache-2.0).  This document is the contract
every ported routine follows.  The bar is simple to state: **a ported routine
returns the same bits as the C routine it came from, on every input, or fails
on exactly the inputs where the C fails.**  The test suite checks that against
NLib itself.

## What to port from

- Only usd-brep's Apache-2.0 source, `source/SMLib/src/NL_*.cpp` and its
  headers -- including routines NVIDIA fences in `#if NLIB_UNUSED`, which are
  part of that source.  Never port from any other copy of NLib or SMLib.
- `tools/manifest.sexp` lists every routine to port: C name, Lisp name, source
  file and line, wave, status (`:done`, `:core`, `:port`, `:vanish`) and callees.

## Where things go

- One Lisp file per NLib source file: `NL_CrvFit.cpp` -> `src/nl/crv-fit.lisp`.
- Its tests: `tests/nl/crv-fit.lisp`.
- These files already exist; add to them.  Do not edit `batten.asd`,
  `src/package.lisp`, `src/types.lisp`, `src/constants.lisp`, `tests/harness.lisp`,
  `tests/oracle*.lisp`, `tests/marshal.lisp` or another module's files.  If a
  shared helper is missing, say so in your report instead.
- Everything lives in package `BATTEN`, unexported.  `src/api.lisp` holds the
  small exported API; ported routines do not need exporting.

## Names

The Lisp name is mechanical (`tools/nlib_names.py`; the manifest has it):
drop `N_`, hyphenate at case changes, lower-case.

    N_CrvEval          -> crv-eval
    N_Combine2CPts     -> combine-2-c-pts
    N_XYZTo2dCPtArray  -> xyz-to-2d-c-pt-array

File-private (`static`) C functions get their file stem: `ST_findStep` in
`NL_CrvNurb.cpp` -> `crv-nurb--st-find-step`.  File-private variables
(`NL_PRIVATE NL_REAL KN_TOL = 1.0e-6;`) become constants or special variables
named the same way: `+crv-adv--kn-tol+`, or `*crv-adv--...*` if written to.
Your own helpers, which are not NLib routines, start with `%`.

## Types

| C | Batten |
|---|---|
| `NL_REAL`, `NL_PARAMETER` | `double-float` |
| `NL_INDEX`, `NL_INTEGER`, `NL_DEGREE` | `fixnum` |
| `NL_FLAG`, `NL_BOOLEAN` arguments | `fixnum`, using the constants in `src/constants.lisp` (`+left+`, `+true+`, `+udir+` ...) |
| `NL_POINT`, `NL_VECTOR` | `point`, a `(simple-array double-float (3))` |
| `NL_CPOINT` | `cpoint`, `(simple-array double-float (4))` -- `(wx wy wz w)` |
| `NL_POINT *` array, `NL_EPOLYGON` | `point-array`, `(* 3)` |
| `NL_CPOINT *` array, `NL_CPOLYGON` | `cpoint-array`, `(* 4)` |
| `NL_POINT **`, `NL_ENET` | `point-net`, `(* * 3)` |
| `NL_CPOINT **`, `NL_CNET` | `cpoint-net`, `(* * 4)` |
| `NL_REAL *` | `real-vector` |
| `NL_REAL **` | `real-matrix` (2D) |
| `NL_INDEX *` | `index-vector` |
| `NL_KNOTVECTOR` | `real-vector` of the knots |
| `NL_CURVE` | `curve` (degree, knots, cpoints) |
| `NL_SURFACE` | `surface` (degree-u, degree-v, knots-u, knots-v, cpoints) |
| `NL_CFUN`, `NL_SFUN` | `cfun`, `sfun` (degree(s), knots, `vals`) |
| `NL_INTERVAL`, `NL_RECTANGLE`, `NL_MINMAXBOX` | `real-vector` of 2, 4, 6 in C field order |
| `NL_LINESEG`, `NL_PLANE` | `lineseg`, `plane` structures |
| `NL_RMATRIX` | `rmatrix` structure (n, m, rm, mtp, bw); `(make-rmatrix)` is the NULL matrix; `+mt-full+` ... `+mt-banded+` |
| `NL_REAL ***` | `real-3d-array`, `(* * *)` |
| `NL_INTEGER **` | `index-matrix`, a fixnum `(* *)` |

`NL_NOZ` and `NL_NOW` stay inside the coordinate arrays as `+no-z+` and
`+no-w+`, exactly where the C has them.  NLib's globals are `*mtol*`, `*ptol*`,
`*ltol*` ... (`src/types.lisp`); `NL_PascalTri` is `+pascal-tri+`.  NLib's
macros: `nl-min`, `nl-max`, `nl-min3`, `nl-max3`, `nl-sign`, `nl-siab` -- use
these, not `min`/`max`, where the C uses the macros.

## Signatures

Arguments in C order, with these changes:

- **Errors.**  An `NL_FLAG` return that reports failure becomes a condition:
  where the C does `NL_ERROR(NL_GEO_ERR)` or `N_ErrSet(...)` and returns
  `NL_YES`, signal `(nurbs-error :geo-err 'routine-name)` -- the keyword is the
  NLib error constant's name.  Callers that `NL_OUT` on a callee's error just
  let the condition propagate.  **Where a C caller ignores a callee's error, or
  inspects it and carries on, the Lisp caller catches `nurbs-error` and carries
  on the same way** (see `srf-eval-pt-normal-at-pole`).  Where the C leaves
  through `NL_QUIT` with no error code set, signal `(nurbs-error :quit 'routine)`.
- **`NL_STACKS *`** arguments disappear; allocate Lisp arrays instead.
- **Scalar outputs** (`NL_REAL *`, `NL_INDEX *`, `NL_FLAG *` holding one value)
  are returned as values, after the C return value if that is not an error
  flag.
- **Point outputs** (`NL_POINT *`, `NL_VECTOR *`, `NL_CPOINT *` for one point)
  and **arrays the caller allocates** (the C documents "MEMORY MUST BE
  ALLOCATED IN THE CALLING ROUTINE", or sizes them by NL_MAXDEG) are passed in
  by the caller and written.
- **Objects the C callee allocates** (an `NL_CURVE *` it sizes with
  `N_CrvSizeArrays`, an `NL_REAL **` it allocates, curve arrays) are returned
  as fresh Lisp objects, in C argument order among the returned values.
- **In-out arguments** (modified in place in C) are modified in place.
- A routine with no scalar outputs returns its first output argument, or
  nothing in particular.
- `NL_BOOLEAN` *return values* of predicates become Lisp booleans.
- `NL_STRING rname` arguments, where kept, are the routine name for errors.

C passes `NL_POINT` and `NL_CPOINT` by value; Batten passes arrays by
reference.  **Read every input coordinate before writing any output**, so that
an output may be the same array as an input (`with-point` / `with-cpoint` in
`src/nl/math.lisp` do this).

## Arithmetic

Same bits means same operations:

- Keep the C's grouping exactly: `a*b + c*d + e*f` is
  `(+ (+ (* a b) (* c d)) (* e f))`.  Operand order within one `+` or `*` does
  not matter (IEEE addition and multiplication commute exactly); grouping does.
- Keep divisions as divisions and reciprocal multiplications as such
  (`x * (1/w)` and `x / w` round differently).
- Same comparisons (`<` vs `<=`), same literal tolerances (`1e-24`), same
  early exits.
- Integer-to-double conversions are explicit: `(float i 1d0)`.  C integer
  division truncates: `(truncate a b)`.
- `sqrt` of a double known non-negative: `(sqrt (the (double-float 0d0) x))`.
  `fabs` is `abs`, `floor`/`ceil` of doubles are `ffloor`/`fceiling`.
  `pow(x, y)` is `(expt x y)` with Y a double, even when it holds a whole
  number.  `sin`, `cos`, `atan2` (`(atan y x)`), `acos`, `exp`, `log` call the
  same libm on Linux.
- Do not "simplify" or merge computations, even when they look redundant; two
  NLib routines computing the same quantity differently must stay different.
- Floating-point traps are masked, as in C (`with-c-float-environment`): a
  division by zero gives an infinity, an invalid operation a NaN.  Port the C's
  own guards (`float-op-is-bad` ...) and nothing more; the test harness runs
  both sides in this environment.

## Shared helpers

Use these rather than writing private copies (`src/nl/math.lisp`):

- Rows and net entries: C hands `Pw[i]` or `Pw[i][j]` straight to
  `N_CopyCPt`, `N_Combine2CPts` and the like.  `%row-c-pt`, `%c-pt-row`,
  `%row-pt`, `%pt-row`, `%net-c-pt`, `%c-pt-net`, `%net-pt`, `%pt-net` move one
  row or entry to or from a scratch point; `%copy-c-pt-row`, `%copy-pt-row`,
  `%copy-c-pt-net`, `%combine-2-c-pt-rows` and `%vector-blend-c-pt-row` do the
  NLib operation on rows directly, with NLib's arithmetic.
- `%check-mem-real-matrix`: `N_CheckMemRealMatrix` for an optional
  caller-allocated `NL_RMATRIX` output.
- `knot-vector-is-end-param`, `knot-vector-is-param-out-of-bounds`
  (`basis-basic.lisp`); `crv-decompose-bez` (`crv-tool.lisp`);
  `%reverse-knots`, NLib's inline `NI_ReverseKnots` (`crv-nurb.lisp`).
- A routine that was first written as a private `%helper` and is now in
  the manifest under its own name: port it under that name (rename the
  helper and its callers in your own files) and give it a test.

Look in the lower-level files for what is already ported before writing a
helper: `git grep -n "(defun NAME"` in your worktree.

## Style

- No `LOOP`: `dolist`, `dotimes`, `do`/`do*`, mapping functions, `reduce`.
- Declare types; put hot loops in `(fast ...)` (speed 3, safety 0); scratch
  arrays `dynamic-extent`, sized as the C sizes them where that is a constant.
- Every ported routine has a docstring whose last line names its source:
  `NLib: N_CrvEval.`
- Each source file starts with the SPDX header both NVIDIA and Genworks hold
  copyright under; the existing stubs have it.

## Routines that vanish

Allocators, frees, `N_InitNurbs`/`N_EndNurbs`, accessors that only unpack a
struct (`N_CrvGetCPtsAndKnots` ...), printers and size checks have no Lisp
counterpart (manifest status `:vanish`).  Where a routine you port calls one,
use the Lisp object directly.

## Tests

For each non-static routine you port, register a differential test **named
exactly its Lisp name** in `tests/nl/<stem>.lisp`.  The core tests
(`tests/nl/math.lisp`, `srf-nurb.lisp` ...) are worked examples.

    (define-oracle-test crv-eval (:file "crv-nurb" :trials 4000)
      (let* ((c (random-curve)) (u (random-parameter (curve-knots c))) (flg (random-side)))
        (agree (list :curve c :u u :flg flg)
               (list (batten::crv-eval c u flg (make-point)))
               (with-foreign ((cur c :curve))
                 (with-foreign-out ((pt :point))
                   (values (c-flag-error-p (c-crv-eval cur u flg pt))
                           (list (from-foreign pt :point))))))))

- `c-<lisp-name>` calls the C routine (bindings generated in
  `tests/oracle-bindings.lisp`, each with its C signature as a comment).
  Structs C passes by value are passed by pointer.
- `with-foreign`, `with-foreign-out`, `from-foreign` move objects across
  (`tests/marshal.lisp` lists the types, `:rmatrix`, `:cnet`, `:index-matrix`
  and `:real-3d` among them); `with-scalars` gives cells for scalar outputs,
  and `:pointer` cells for pointer outputs (`NL_REAL **`, `NL_CURVE ***` ...);
  `with-stacks` gives an `NL_STACKS` for C routines that take one -- read
  results back inside it.
- NLib's constants (`+udir+`, `+epoint+`, `+mt-banded+` ...) are visible in
  tests unqualified; other Batten internals are `batten::name`.
- Inputs: `random-curve`, `random-surface`, `random-cfun`, `random-sfun`,
  `random-points`, `random-parameter` ... (`tests/random.lisp`).  Cover the
  branches the C has: degrees, rational and not, 2D, parameters at knots and
  ends, degenerate cases the C handles.  Keep inputs inside the C's domain --
  a segmentation fault in the C kills the test process.
- Static routines are tested through their callers.

Run them:

    tools/run-tests.sh crv-fit          # one file's tests, plus its coverage
    tools/run-tests.sh crv-fit 200      # fewer trials while working
    tools/run-tests.sh all              # everything

Each run is a fresh SBCL in a container (about three seconds), so several can
run at once.  A routine is done when its test shows `0 MISMATCH` and is not
`CRASHED`.  "both-fail" is agreement: both sides failed on that input.

## Bit-identical, or explain

If a routine cannot be made bit-identical -- a libm function that differs, an
uninitialized read in the C, a C bug that depends on memory layout -- do not
paper over it with a tolerance.  Leave the test reporting the mismatch and
describe the cause in your report.
