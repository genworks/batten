# Batten

NURBS curves and surfaces in portable Common Lisp: a port of NLib, the NURBS
library inside NVIDIA's open-source SMLib
([usd-brep](https://github.com/NVIDIA-Omniverse/usd-brep)), under the same
Apache License 2.0.

A batten is the long flexible strip that shipwrights and loftsmen bent
through weights to draw fair curves: the physical ancestor of the spline.

**Status: release 0.1.0 (October 2026).**  Every NLib routine the Gendl/GDL
stack reaches is ported, and NLib's tessellator -- 473 routines across 29 NLib source
files: basis functions, curve and surface evaluation, knot insertion and
removal, Bezier decomposition, interpolation and approximation, offsets,
conics, matrices, triangulation of trimmed surfaces and the rest
(`tools/manifest.sexp` lists them).  504 differential tests cover them and
the API (the file-private statics through their callers), and all agree with
NLib.

Two packages: `batten` exports the readable API (`API.md`, `src/api/`) --
construction, queries and operations on curves and surfaces, geometry
utilities, tessellation -- and `batten.nlib` exports the port itself, NLib's
own routines under their Lisp names (N_CrvEval is `crv-eval`), for anyone
who wants NLib's interface.

## Getting started

Batten is plain portable Common Lisp with no dependencies.  Clone it where
ASDF or Quicklisp finds local projects (`~/common-lisp/`,
`~/quicklisp/local-projects/`), then:

```lisp
(asdf:load-system "batten")          ; or (ql:quickload "batten")

(defpackage :batten-example (:use :common-lisp :batten))
(in-package :batten-example)

;; a cubic through four points, chord-length parameterized
(defparameter *c* (interpolate-curve '((0 0 0) (1 2 0) (3 3 0) (4 0 0))))
(curve-point *c* 0.5)                 ; => #(2.417893218813452d0 2.892766952966369d0 0.0d0)
(curve-length *c*)                    ; => 8.265825079520468d0
(curve-closest-point *c* '(2 1 0))    ; => parameter, point, distance

;; a ruled surface between two lines, and its triangles
(defparameter *s* (ruled-surface (line-segment '(0 0 0) '(4 0 0))
                                 (line-segment '(0 3 1) '(4 3 1))))
(surface-point *s* 0.5 0.5)           ; => #(2.0d0 1.5d0 0.5d0)
(tessellate-surface *s*)              ; => points, triangles, parameters
```

`API.md` lists the whole API; every function's docstring is its reference.

## How it is held to NLib

Each routine is a translation of one NLib routine, operation for operation,
so that its results agree with NLib's C to the bit -- not merely within a
tolerance.  The test suite checks this directly: it loads NLib from an SMLib
build through CFFI and runs every ported routine against the C routine it
came from, on random curves and parameters (`batten/tests`, which needs CFFI
and the library).  Each docstring names its NLib source, e.g.
`NLib: N_CrvDerivs`.

To run the tests, build the oracle from a checkout of
[usd-brep](https://github.com/NVIDIA-Omniverse/usd-brep) (a C++17 compiler is
needed; it compiles SMLib once, a few minutes), then run them in a
throwaway container of the public Gendl SBCL image:

```sh
tools/build-oracle.sh /path/to/usd-brep/source/SMLib oracle
tools/run-tests.sh all          # or one NLib file's tests: tools/run-tests.sh crv-nurb
```

Because every routine performs the same IEEE double operations in the same
order, the results are also identical across Common Lisp implementations.
`batten-random:evaluation-digest` hashes the bits of 50,000 random
evaluations; seven implementations produce the same digest -- SBCL, Clozure
CL, Allegro CL, ECL, CLISP, dotcl and EGCL.

## Conventions

They follow NLib, which follows *The NURBS Book* (Piegl and Tiller):

- A curve has a degree, a knot vector `U[0..m]` and homogeneous control points
  `Pw[0..n]`, stored flat as `(wx wy wz w)` per point.
- `w` equal to `+no-w+` marks a non-rational curve; `z` equal to `+no-z+`
  marks a 2D one.
- At a knot, `:side :left` takes the span ending there and `:right` the span
  beginning there.
- Routines ending in `!` write into an array the caller supplies and cons
  nothing for the result; their plain counterparts return fresh arrays.
- Failures signal `nurbs-error`, whose code names the NLib error.

## Benchmarks

`batten/bench` times the routines on the same seeded curves on any
implementation; `batten/bench-c` sets them beside NLib's C.

```lisp
(asdf:load-system "batten/bench")
(batten-bench:run-bench)
```

## License

Apache License 2.0; see `LICENSE` and `NOTICE`.  Batten is a modified version
of NVIDIA's NLib and is not affiliated with or endorsed by NVIDIA.
