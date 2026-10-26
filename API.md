# Batten's readable API

Batten has two layers.  The **port** is every NLib routine under its
mechanical name (`crv-refine`, `fit-crv-interp` ...), with NLib's argument
order and conventions, held bit-identical to the C; package `batten.nlib`
exports it.  The **API**, exported from `batten`, is what a Lisp programmer
should reach for: readable names, keyword arguments, fresh results,
conditions.  Each API function is a thin wrapper over one or a few port
routines, so it inherits their exactness; each has a differential test
against the C routines it wraps (`tests/api/`).

The API lives in `src/api/`, a file per section below.  The tables give the
names as built; each function's docstring is its reference (lambda list,
keywords and defaults, values), ending with the port routines it calls.
Where the built API departs from the first sketch, the section's "Changed"
note says how.

## Principles

- **Objects are the port's structures**: `curve`, `surface`, `cfun`, `sfun`,
  points as `(simple-array double-float (3))`.  No CLOS layer over them; a
  few generic functions only where the same operation applies to curves and
  surfaces alike (`domain`, `bounding-box`, `copy`, `transform` ...).
- **Pure by default.**  Functions return fresh objects and leave their
  arguments alone.  A `!` variant writes into a result the caller supplies
  and conses nothing (`curve-point!`), for inner loops.
- **Any real, any point sequence**: parameters and lengths may be any real,
  points and vectors any sequence of three reals.
- **Keywords for everything optional**: `:side :left/:right` at a knot,
  `:tolerance`, `:degree`, `:parameterization :chord-length/:centripetal/
  :uniform`.  Defaults follow NLib's own (`*mtol*`, `*ptol*` ...) where NLib
  has one.
- **Errors are conditions.**  `nurbs-error` with NLib's code (`:par-err`,
  `:geo-err` ...) and the routine; no status returns.
- **The C float environment is established at the boundary**: every API
  function runs its port routines under `with-c-float-environment`, so
  callers never see a stray `division-by-zero` from NLib's own
  compute-then-test idiom (Allegro excepted; see `types.lisp`).
- **Names say what, not how**: `curve-closest-point`, not
  `crv-closest-pt`; `interpolate-curve`, not `fit-crv-interp`.  No API name
  is a port routine's or a COMMON-LISP symbol.

## Shared by curves and surfaces

Generic functions (`src/api/generics.lisp`), with methods for `curve` and
`surface`: `domain`, `bounding-box`, `copy`, `transform`, `translate`,
`insert-knots`, `elevate-degree`, `reparameterize`, `split`.  The surface
methods of `insert-knots`, `elevate-degree` and `split` take
`:direction :u` or `:v`.

## Evaluation

`curve-point`, `curve-point!`, `curve-derivatives`, `curve-derivatives!`,
`surface-point`, `surface-point!`, `surface-derivatives`,
`surface-derivatives!`, `find-span`, `curve-rational-p`, `surface-rational-p`
(`crv-eval`, `crv-derivs`, `srf-eval-pt`, `srf-derivs`, `basis-find-span`,
`is-crv-rat`, `is-srf-rat`).

## Construction

| API | Port routines |
|---|---|
| `(make-curve degree knots points &key weights two-d)` | — |
| `(make-surface degree-u degree-v knots-u knots-v net &key weights)` | `%make-surface` |
| `(line-segment start end)` | `crv-line-from-pt-and-vector` |
| `(circular-arc center radius &key x-axis y-axis start-angle end-angle degree angle-tolerance)` | `create-circ-arc-tol` |
| `(elliptical-arc center major-radius minor-radius &key ...)` | `create-elliptical-arc-ang-type-tol` |
| `(conic-arc start start-tangent end end-tangent through)` | `create-conic-arc` |
| `(interpolate-curve points &key degree parameterization parameters start-tangent end-tangent start-derivative end-derivative)` | `fit-crv-interp`, `fit-crv-knots-and-tangents` |
| `(approximate-curve points &key degree tolerance ... control-points)` | `fit-crv-approx`, `fit-crv-approx-knots` |
| `(fit-arc points &key tolerance ...)`, `(fit-line points)` | `fit-arc-to-pts`, `line-fit-pts` |
| `(interpolate-surface net &key degree-u degree-v ...)` | `fit-srf-to-pts`, `fit-srf-to-pts-knots` |
| `(approximate-surface net &key tolerance ...)` | `fit-srf-approx-tol`, `fit-srf-lst-sq-approx` |
| `(ruled-surface curve-1 curve-2 &key across)` | `create-ruled-srf-from-boundary-crvs` |
| `(extruded-surface curve direction length &key across)` | `create-srf-extrude-crv` |
| `(revolved-surface curve axis-point axis-direction &key angle degree)` | `create-revolved-srf` |
| `(skinned-surface curves &key degree ...)` | `create-skin-srf`, `create-skin-srf-params` |
| `(gordon-surface u-curves vs v-curves us &key ...)` | `create-gordon-srf` |
| `(coons-surface boundary derivatives &key check-corners)` | `create-coons-boundary-crvs` |
| `(sphere center radius &key start-angle end-angle sweep-angle ...)` | `create-sphere` |

Changed: `line-curve` became `line-segment` (a test helper of that name
exists), and it ends exactly at END where NLib's START + (END − START) can
miss by a bit.  Angles are in degrees, as in NLib.  `circular-arc` takes
NLib's center, radius and axes rather than a normal and a start point (which
would add arithmetic NLib does not do); Gordon takes the curves' parameters
and Coons the four cross-boundary derivative curves, as NLib needs them.

## Curve queries

| API | Port routines |
|---|---|
| `(curve-tangent curve u &key side)` → `(values tangent point)` | `crv-eval-tangent` |
| `(curve-curvature curve u &key side)` | `crv-eval-curvature` |
| `(curve-frenet-frame curve u &key side)` → `(values point tangent normal binormal)` | `crv-eval-frenet-frame` |
| `(curve-length curve &key from to tolerance relative)` | `crv-arc-length` |
| `(curve-point-at-length curve fraction &key tolerance relative length)` → `(values point parameter length)` | `crv-percentage-along-pt` |
| `(curve-closest-point curve point &key tolerance cosine-tolerance start)` → `(values u point distance)` | `crv-closest-pt` |
| `(curve-project-points curve points &key tolerance)` → `(values parameters projections)` | `crv-project-pts` |
| `(curve-polyline curve tolerance &key relative)` → `(values points parameters)` | `approx-crv-with-polyline` |
| `curve-closed-p`, `curve-planar-p`, `curve-line-p`, `curve-degenerate-p`, `curves-coincident-p` | `crv-is-closed`, `crv-is-planar`, `crv-is-line`, `crv-is-degen`, `crvs-are-coincident` |
| `(curve-type curve &key tolerance)` → `:line`, `:circle` or `:freeform`, with its data | `crv-get-type` |

Changed: `curve-type` tells lines and circles from free-form curves only
(NLib's `crv-get-type` does; nothing ported detects a conic).  Zero or
negative tolerances are refused where NLib would subdivide for ever.
`curve-closest-point` can signal `:con-err` where NLib's Newton iteration
does not converge (at corners, at a closed curve's seam).

## Curve operations

| API | Port routines |
|---|---|
| `insert-knots`, `elevate-degree`, `split`, `reparameterize`, `translate`, `transform` (methods) | `crv-refine`, `crv-elevate-degree`, `crv-extract-crv-seg`, `crv-reparam`, `crv-translate`, `crv-transform` |
| `(remove-knots curve &key tolerance)` | `crv-remove-knots` |
| `(subcurve curve u0 u1)` | `crv-extract-crv-seg` |
| `(bezier-segments curve)` | `crv-decompose-bez` |
| `(reverse-curve curve)` | `crv-reverse` |
| `(extend-curve curve distance &key end continuity relative)` | `crv-extend-by-dist` |
| `(offset-curve curve distance &key normal tolerance degree parameterization keep-end-tangents)` | `crv-offset-pt-sampling` |
| `(make-compatible curves &key knot-tolerance)` | `crvs-make-compatible`, `crvs-make-compatible-knot-tol` |
| `curve-to-2d`, `curve-to-3d`, `curve-to-rational`, `(curve-to-polynomial curve &key tolerance ...)` | `crv-3d-to-2d`, `crv-2d-to-3d`, `crv-non-rat-to-rat`, `approx-nurbs-with-non-rat-crv` |

Changed: `rational` and `polynomial` became `curve-to-rational` and
`curve-to-polynomial` (`rational` is a COMMON-LISP symbol).  A curve already
3D, rational or polynomial comes back from those as a copy (the NLib routines
would alter it).

## Surface queries and operations

| API | Port routines |
|---|---|
| `domain`, `bounding-box`, `copy`, `transform`, `translate`, `insert-knots`, `elevate-degree`, `reparameterize`, `split` (methods) | `srf-get-parameter-bounds`, `srf-get-b-box`, `srf-copy`, `transform-c-pt`, `srf-translate`, `srf-insert-knots`, `srf-elevate-degree`, `srf-reparam`, `srf-split` |
| `(surface-normal s u v &key u-side v-side)` → `(values point normal u-tangent v-tangent)`; poles handled | `srf-eval-pt-pt-deriv-normal` |
| `(surface-grid s us vs &key u-side v-side)` | `srf-eval-pt-grid` |
| `(surface-closest-point s point &key u v tolerance cosine-tolerance error-p)` → `(values u v point distance converged-p)` | `get-closest-pt-on-srf` |
| `(isocurve s &key u v)` | `srf-extract-iso-crv` |
| `(subsurface s u0 u1 v0 v1)` | `srf-extract-patch` |
| `(bezier-patches s)` | `srf-decompose-to-bez` |
| `(offset-surface s distance &key tolerance ...)` | `srf-offset` |
| `(surface-swap-uv s)` | `swap-uv` |
| `surface-closed-p`, `surface-flat-p` | `srf-is-closed`, `srf-is-flat-cheap` |
| `(surface-type s &key tolerance)` → `(values type properties)` | `srf-type` |

Changed: `swap-uv` became `surface-swap-uv` (`swap-uv` is the port's).
`surface-closest-point` gained `:error-p` and a fifth value, whether NLib's
iteration converged; with `:error-p nil` it returns NLib's best point.

## Geometry utilities and scalar B-spline functions

| API | Port routines |
|---|---|
| `distance`, `dot`, `cross`, `(normalize v)` → `(values unit magnitude)`, `angle-between` (degrees) | `dist-pt-pt`, `vector-dot`, `vector-cross`, `vector-normalize`, `vectors-angle` |
| `project-point-to-line`, `project-point-to-plane` | `project-pt-line`, `project-pt-plane` |
| `intersect-lines`, `intersect-line-plane`, `intersect-planes` | `isect-line-line`, `isect-line-plane`, `isect-plane-plane` |
| `make-cfun`, `cfun-value`, `cfun-derivatives`, `(cfun-derivative f &key order)`, `cfun*`, `cfun+`, `cfun-` | `c-func-from-knots`, `c-func-eval`, `c-func-derivs`, `crv-func-deriv`, `crv-func-multiply-crv-func`, `crv-func-sum-diff-crv-func` |
| `make-sfun` | `s-func-from-knot-vectors` |

A line is a point and a direction, a plane a point and a normal.  Changed:
`cfun-` added; for `sfun` only the constructor, as nothing ported evaluates
one.

## Tessellation

`(tessellate-surface s &key tolerance trim-tolerance parameter-tolerance
outer-loops holes normals)` → `(values points triangles parameters
[normals])`: NLib's tessellator (`tess-trimmed-srf`, `tess-get-triangles`),
untrimmed or trimmed by 2D loops in parameter space, holes included.

## Decisions on the sketch's open questions

1. **Package layout**: both.  `batten` exports the API; `batten.nlib`
   (`src/nlib.lisp`, generated by `tools/gen-nlib-package.lisp`) re-exports
   the port -- every ported routine that is not file-private, and NLib's
   constants -- as a documented surface for anyone who wants NLib's own
   interface.
2. **Generic functions**: yes, for the operations above; the `!` functions
   and the port are there for inner loops.
3. **Gendl integration**: done -- Gendl's `batten-kernel` implements the
   NURBS curve and surface part of surf's geometry-kernel protocol on this
   API and the port, tessellation included; breps and solids stay with a
   kernel that has them.
