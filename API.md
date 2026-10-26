# Batten's readable API — a sketch

Batten has two layers.  The **port** is every NLib routine under its
mechanical name (`batten::crv-refine`, `batten::fit-crv-interp` ...), with
NLib's argument order and conventions, held bit-identical to the C.  The
**API**, exported from `batten`, is what a Lisp programmer should reach for:
readable names, keyword arguments, fresh results, conditions.  Each API
function is a thin wrapper over one or a few port routines, so it inherits
their exactness.  This file proposes that API; only evaluation exists today
(`src/api.lisp`).

## Principles

- **Objects are the port's structures**: `curve`, `surface`, `cfun`, `sfun`,
  points as `(simple-array double-float (3))`.  No CLOS layer over them; a
  few generic functions only where the same operation applies to curves and
  surfaces alike (`domain`, `bounding-box`, `copy`, `transform`).
- **Pure by default.**  Functions return fresh objects and leave their
  arguments alone.  A `!` variant writes into a result the caller supplies
  and conses nothing (`curve-point!`), for inner loops.
- **Keywords for everything optional**: `:side :left/:right` at a knot,
  `:tolerance`, `:degree`, `:parameterization :chord-length/:centripetal/
  :uniform`.  Defaults follow NLib's own (`*mtol*`, `*ptol*` ...).
- **Errors are conditions.**  `nurbs-error` with NLib's code (`:par-err`,
  `:geo-err` ...) and the routine; no status returns.
- **The C float environment is established at the boundary**: every API
  function runs its port routines under `with-c-float-environment`, so
  callers never see a stray `division-by-zero` from NLib's own
  compute-then-test idiom (Allegro excepted; see `types.lisp`).
- **Names say what, not how**: `curve-closest-point`, not
  `crv-closest-pt`; `interpolate-curve`, not `fit-crv-interp`.

## Construction

| API | Port routines |
|---|---|
| `(make-curve degree knots points &key weights two-d)` | exists |
| `(make-surface degree-u degree-v knots-u knots-v net &key weights)` | `%make-surface` |
| `(line-curve p q)` | `crv-line-from-pt-and-vector` |
| `(circular-arc center normal start sweep &key radius)` | `create-circ-arc`, `create-circ-arc-tol` |
| `(elliptical-arc ...)`, `(conic-arc p0 t0 p2 t2 p)` | `create-elliptical-arc`, `create-conic-arc` |
| `(interpolate-curve points &key degree parameterization tangents)` | `fit-crv-interp`, `fit-crv-knots-and-tangents`, `fit-crv-derivs` |
| `(approximate-curve points &key degree tolerance)` | `fit-crv-approx`, `fit-crv-approx-knots`, `fit-remove-knots` |
| `(fit-arc points)`, `(fit-line points)` | `fit-arc-to-pts`, `line-fit-pts` |
| `(interpolate-surface net &key degree-u degree-v)` | `fit-srf-to-pts`, `fit-srf-to-pts-knots` |
| `(approximate-surface net &key tolerance)` | `fit-srf-approx-tol`, `fit-srf-lst-sq-approx` |
| `(ruled-surface c1 c2)` | `create-ruled-srf-from-boundary-crvs` |
| `(extruded-surface curve direction length)` | `create-srf-extrude-crv` |
| `(revolved-surface curve axis-point axis-direction &key angle)` | `create-revolved-srf` |
| `(skinned-surface curves &key degree)` | `create-skin-srf`, `create-skin-srf-params` |
| `(gordon-surface u-curves v-curves)` | `create-gordon-srf` |
| `(coons-surface boundary)` | `create-coons-boundary-crvs` |
| `(sphere center radius)` | `create-sphere` |

## Curve queries

| API | Port routines |
|---|---|
| `(domain curve)` → `(values u0 u1)` | `crv-get-param-bounds` |
| `(curve-point curve u &key side)`, `curve-point!` | exists (`crv-eval`) |
| `(curve-derivatives curve u n &key side)`, `curve-derivatives!` | exists (`crv-derivs`) |
| `(curve-tangent curve u)` | `crv-eval-tangent` |
| `(curve-curvature curve u)` | `crv-eval-curvature` |
| `(curve-frenet-frame curve u)` → `(values point tangent normal binormal)` | `crv-eval-frenet-frame` |
| `(curve-length curve &key from to tolerance)` | `crv-arc-length` |
| `(curve-point-at-length curve fraction)` | `crv-percentage-along-pt` |
| `(curve-closest-point curve point &key tolerance)` → `(values u point distance)` | `crv-closest-pt` |
| `(curve-project-points curve points)` | `crv-project-pts` |
| `(bounding-box curve)` | `crv-get-b-box` |
| `(curve-polyline curve tolerance)` | `approx-crv-with-polyline` |
| predicates: `curve-closed-p`, `curve-planar-p`, `curve-line-p`, `curve-degenerate-p`, `curve-rational-p`, `curves-coincident-p` | `crv-is-closed`, `crv-is-planar`, `crv-is-line`, `crv-is-degen`, `is-crv-rat`, `crvs-are-coincident` |
| `(curve-type curve)` → `:line`, `:circle`, `:conic`, `:freeform` | `crv-get-type`, `crv-get-conic-data` |

## Curve operations

| API | Port routines |
|---|---|
| `(insert-knots curve knots)` | `crv-refine`, `crv-insert-knot` |
| `(remove-knots curve &key tolerance)` | `crv-remove-knots`, `crv-remove-all-knots`, `crv-remove-knot-max-err` |
| `(elevate-degree curve by)` | `crv-elevate-degree` |
| `(split-curve curve u)` → two curves; `(subcurve curve u0 u1)` | `crv-extract-crv-seg` |
| `(bezier-segments curve)` | `crv-decompose-bez` |
| `(reverse-curve curve)` | `crv-reverse` |
| `(reparameterize curve u0 u1)` | `crv-reparam-to-interval`, `crv-reparam` |
| `(extend-curve curve distance &key end continuity)` | `crv-extend-by-dist` |
| `(offset-curve curve distance &key normal tolerance)` | `crv-offset-pt-sampling` |
| `(make-compatible curves)` | `crvs-make-compatible`, `crvs-make-compatible-knot-tol` |
| `(translate curve vector)`, `(transform curve matrix)` | `crv-translate`, `crv-transform` |
| `(curve-to-2d curve)`, `(curve-to-3d curve)`, `(rational curve)`, `(polynomial curve &key tolerance)` | `crv-3d-to-2d`, `crv-2d-to-3d`, `crv-non-rat-to-rat`, `approx-nurbs-with-non-rat-crv` |

## Surface queries and operations

| API | Port routines |
|---|---|
| `(domain surface)` → `(values u0 u1 v0 v1)` | `srf-get-parameter-bounds` |
| `(surface-point s u v &key u-side v-side)`, `surface-point!` | exists (`srf-eval-pt`) |
| `(surface-derivatives s u v nu nv)`, `surface-derivatives!` | exists (`srf-derivs`) |
| `(surface-normal s u v)` → `(values point normal)`; poles handled | `srf-eval-pt-pt-deriv-normal`, `srf-eval-pt-normal-at-pole` |
| `(surface-grid s us vs)` | `srf-eval-pt-grid` |
| `(surface-closest-point s point)` → `(values u v point distance)` | `get-closest-pt-on-srf` |
| `(bounding-box surface)` | `srf-get-b-box` |
| `(isocurve s :u value)` / `(isocurve s :v value)` | `srf-extract-iso-crv` |
| `(subsurface s u0 u1 v0 v1)`, `(split-surface s :u value)` | `srf-extract-patch`, `srf-split` |
| `(insert-knots s knots :direction :u)`, `(elevate-degree s by :direction :v)` | `srf-insert-knots`, `srf-elevate-degree` |
| `(bezier-patches s)` | `srf-decompose-to-bez` |
| `(offset-surface s distance &key tolerance)` | `srf-offset` |
| `(swap-uv s)`, `(reparameterize s rectangle)` | `swap-uv`, `srf-reparam-to-interval` |
| predicates: `surface-closed-p`, `surface-flat-p`, `surface-rational-p` | `srf-is-closed`, `srf-is-flat-cheap`, `is-srf-rat` |
| `(surface-type s)` → `:plane`, `:cylinder`, `:cone`, `:sphere`, `:revolved` ... | `srf-type` |

## Geometry utilities

Points and vectors are plain double arrays; a short set is worth exporting
under readable names: `distance`, `dot`, `cross`, `normalize`,
`project-point-to-line`, `project-point-to-plane`, `intersect-lines`,
`intersect-line-plane`, `intersect-planes`, `angle-between`
(`dist-pt-pt`, `vector-dot`, `vector-cross`, `vector-normalize`,
`project-pt-line`, `project-pt-plane`, `isect-line-line`,
`isect-line-plane`, `isect-plane-plane`, `vectors-angle`).

## Scalar B-spline functions

`(make-cfun degree knots values)`, `(cfun-value f u)`,
`(cfun-derivatives f u n)`, `(cfun-derivative f)`, `(cfun* f g)`,
`(cfun+ f g)` (`c-func-eval`, `c-func-derivs`, `crv-func-deriv`,
`crv-func-multiply-crv-func`, `crv-func-sum-diff-crv-func`); the same
shape for `sfun`.

## Open questions

1. **Package layout**: one `batten` package, or `batten` (API) plus
   `batten.nlib` (the port, exported under its mechanical names for anyone
   who wants NLib's exact interface)?  The second makes the port a
   documented, supported surface instead of `batten::` internals.
2. **Generic functions** (`domain`, `bounding-box`, `copy`, `transform`,
   `insert-knots`, `elevate-degree`, `reparameterize` across curves and
   surfaces) cost a dispatch per call; worth it for readability at this
   level, with the `!` and port layers for speed.
3. **Gendl integration**: Gendl's `surf` objects (`b-spline-curve`,
   `b-spline-surface`, `fitted-curve` ...) call a geometry-kernel protocol
   -- `(make-b-spline-curve *geometry-kernel* ...)`,
   `(get-nurbs-from-curve *geometry-kernel* ...)` -- whose only
   implementation today is GDL's SMLib frosting.  A Batten kernel would
   implement the NURBS curve and surface part of that protocol in Lisp,
   giving Gendl working curves and surfaces on every implementation; breps
   and solids stay SMLib's.  The API above is what those methods would
   call.
