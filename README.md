# Batten

NURBS curves and surfaces in portable Common Lisp: a port of NLib, the NURBS
library inside NVIDIA's open-source SMLib
([usd-brep](https://github.com/NVIDIA-Omniverse/usd-brep)), under the same
Apache License 2.0.

A batten is the long flexible strip that shipwrights and loftsmen bent
through weights to draw fair curves: the physical ancestor of the spline.

**Status: not yet released.**  Every NLib routine the Gendl/GDL stack reaches
is ported -- 406 routines across 28 NLib source files: basis functions, curve
and surface evaluation, knot insertion and removal, Bezier decomposition,
interpolation and approximation, offsets, conics, matrices and the rest
(`tools/manifest.sexp` lists them).  377 differential tests cover them (the
file-private statics through their callers), and all agree with NLib.  What remains
before a release is integration into GDL and a performance pass; the readable
exported API (`src/api.lisp`) covers only evaluation so far, and the ported
routines are reached as `batten::name`.

## How it is held to NLib

Each routine is a translation of one NLib routine, operation for operation,
so that its results agree with NLib's C to the bit -- not merely within a
tolerance.  The test suite checks this directly: it loads NLib from an SMLib
build through CFFI and runs every ported routine against the C routine it
came from, on random curves and parameters (`batten/tests`, which needs CFFI
and the library).  Each docstring names its NLib source, e.g.
`NLib: N_CrvDerivs`.

Because every routine performs the same IEEE double operations in the same
order, the results are also identical across Common Lisp implementations.
`batten-random:evaluation-digest` hashes the bits of 50,000 random
evaluations; SBCL, Clozure CL and Allegro CL produce the same digest.

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
