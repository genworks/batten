;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: Asking and changing a surface: domain, normals (poles
;;;; included), grids, closest points, bounding box, iso-curves, patches and
;;;; splits, knots and degree, Bezier patches, offsets, its type and its
;;;; predicates.

(in-package :batten)

;;; ------------------------------------------------------------------ helpers

(defun %surface-direction-flag (direction)
  "NLib's flag for a parameter direction, :u or :v."
  (ecase direction (:u +udir+) (:v +vdir+)))

(defun %surface-real-vector (reals)
  "REALS, any sequence of reals, as a fresh real-vector."
  (map 'real-vector #'%real reals))

(defun %surface-matrix (matrix)
  "MATRIX, a 4x4 array of reals or an rmatrix, as the rmatrix N_TransformCPt reads."
  (if (rmatrix-p matrix)
      matrix
      (let ((rm (make-real-matrix 4 4)))
        (dotimes (i 4)
          (dotimes (j 4)
            (setf (aref rm i j) (%real (aref matrix i j)))))
        (make-rmatrix 3 3 rm +mt-full+ 0))))

(defun %surface-sample-parameters (lo hi count)
  "COUNT+1 evenly spaced parameters from LO to HI, the last exactly HI."
  (declare (type double-float lo hi) (type fixnum count))
  (let ((v (make-real-vector (1+ count)))
        (step (/ (- hi lo) (float count 1d0))))
    (declare (type real-vector v) (type double-float step))
    (dotimes (i count)
      (setf (aref v i) (+ lo (* (float i 1d0) step))))
    (setf (aref v count) hi)
    v))

(defun %surface-closest-start (surface point)
  "The parameters (values u v) of the point nearest POINT among a grid of
2(n+1)+1 by 2(m+1)+1 evenly spaced samples of SURFACE: the start of the
Newton iteration of surface-closest-point.  Must run in the C float
environment."
  (declare (type surface surface) (type point point))
  (multiple-value-bind (u0 u1 v0 v1) (srf-get-parameter-bounds surface)
    (let* ((k (* 2 (1+ (surface-n surface))))
           (l (* 2 (1+ (surface-m surface))))
           (us (%surface-sample-parameters u0 u1 k))
           (vs (%surface-sample-parameters v0 v1 l))
           (grid (srf-eval-pt-grid surface us vs k l +left+ +left+ (make-point-net (1+ k) (1+ l))))
           (q (make-point))
           (best +bigd+) (bu u0) (bv v0))
      (declare (type fixnum k l) (type real-vector us vs) (type point-net grid)
               (type double-float best bu bv))
      (dotimes (i (1+ k))
        (dotimes (j (1+ l))
          (set-point q (aref grid i j 0) (aref grid i j 1) (aref grid i j 2))
          (let ((d (dist-pt-pt q point)))
            (declare (type double-float d))
            (when (< d best)
              (setf best d bu (aref us i) bv (aref vs j))))))
      (values bu bv))))

(defun %surface-type-properties (stp dir aa bb cc dd xx yy r1 r2 ac al h cur-a cur-b)
  "The outputs of N_SrfType (srf-type) as surface-type returns them:
(values type properties)."
  (let ((direction (cond ((= dir +udir+) :u) ((= dir +vdir+) :v) (t nil))))
    (cond
      ((= stp +nplane+)
       (values :plane (list :corners (list aa bb cc dd))))
      ((or (= stp +ncylinder+) (= stp +ncone+))
       (values (if (= stp +ncylinder+) :cylinder :cone)
               (list :base-center aa :axis bb :x-axis xx :y-axis yy
                     :radius r1 :end-radius r2 :height h :angle al
                     :profile cur-a :direction direction)))
      ((or (= stp +nsphere+) (= stp +ntorus+))
       (values (if (= stp +nsphere+) :sphere :torus)
               (list :axis-point aa :axis bb :angle al
                     :center dd :normal cc :x-axis xx :y-axis yy
                     :radius r1 :axis-distance r2 :profile-sweep ac
                     :profile cur-a :direction direction)))
      ((= stp +nrevolution+)
       (values :revolved
               (list :axis-point aa :axis bb :angle al :profile-normal cc
                     :profile cur-a :direction direction)))
      ((= stp +nextrusion+)
       (values :extruded
               (list :vector aa :length r1 :profile cur-a :direction direction)))
      ((= stp +nruled+)
       (values :ruled (list :profile cur-a :opposite cur-b :direction direction)))
      (t (values :freeform nil)))))

;;; ------------------------------------------------------- the generic functions

(defmethod domain ((surface surface))
  "The parameter domain of SURFACE: (values u0 u1 v0 v1).
Port: srf-get-parameter-bounds."
  (with-c-float-environment
    (srf-get-parameter-bounds surface)))

(defmethod bounding-box ((surface surface))
  "The axis-aligned box around the control points of SURFACE (which contains
the surface): (values min-point max-point), fresh.
Port: srf-get-b-box."
  (let ((box (make-real-vector 6)))
    (declare (type real-vector box))
    (with-c-float-environment
      (srf-get-b-box surface box))
    (values (make-point (aref box 0) (aref box 2) (aref box 4))
            (make-point (aref box 1) (aref box 3) (aref box 5)))))

(defmethod copy ((surface surface))
  "A copy of SURFACE sharing nothing with it.
Port: srf-copy."
  (with-c-float-environment
    (srf-copy surface)))

(defmethod transform ((surface surface) matrix)
  "A fresh copy of SURFACE with every control point transformed by MATRIX, a
4x4 array of reals (rotation and scale upper left, translation in the last
column, projection in the last row) or an rmatrix.  A 2D surface stays 2D.
Port: srf-copy, transform-c-pt."
  (let ((rma (%surface-matrix matrix)))
    (with-c-float-environment
      (let* ((result (srf-copy surface))
             (pw (surface-cpoints result))
             (cw (make-cpoint)))
        (declare (type cpoint-net pw) (dynamic-extent cw))
        (dotimes (i (array-dimension pw 0))
          (dotimes (j (array-dimension pw 1))
            (%net-c-pt pw i j cw)
            (transform-c-pt cw rma cw)
            (%c-pt-net cw pw i j)))
        result))))

(defmethod translate ((surface surface) vector)
  "A fresh copy of SURFACE moved by VECTOR (any three reals).
Port: srf-copy, srf-translate."
  (let ((tv (%point vector)))
    (with-c-float-environment
      (srf-translate (srf-copy surface) tv))))

(defmethod insert-knots ((surface surface) knots &key direction)
  "A fresh copy of SURFACE with KNOTS (a sequence of reals strictly inside
the domain, in any order, repeats allowed up to the degree) inserted as
DIRECTION (:u or :v) knots; the shape is unchanged.  No knots gives a copy.
Port: srf-insert-knots (srf-copy for no knots)."
  (let ((dir (%surface-direction-flag direction))
        (knx (sort (%surface-real-vector knots) #'<)))
    (declare (type real-vector knx))
    (with-c-float-environment
      (if (zerop (length knx))
          (srf-copy surface)
          (srf-insert-knots surface knx dir)))))

(defmethod elevate-degree ((surface surface) by &key direction)
  "A fresh copy of SURFACE with its degree in DIRECTION (:u or :v) raised by
BY, a non-negative integer; the shape is unchanged.
Port: srf-elevate-degree."
  (let ((dir (%surface-direction-flag direction)))
    (with-c-float-environment
      (srf-elevate-degree surface by dir))))

(defmethod reparameterize ((surface surface) &rest bounds)
  "A fresh copy of SURFACE with its domain mapped linearly onto BOUNDS, the
four reals u0 u1 v0 v1; a direction whose ends already match is left alone.
Port: srf-copy, srf-reparam (srf-reparam-to-interval)."
  (destructuring-bind (u0 u1 v0 v1) bounds
    (let ((u0 (%real u0)) (u1 (%real u1)) (v0 (%real v0)) (v1 (%real v1)))
      (with-c-float-environment
        (srf-reparam (srf-copy surface) u0 u1 v0 v1)))))

(defmethod split ((surface surface) at &key direction)
  "SURFACE cut at AT, a u value (DIRECTION :u) or a v value (:v) strictly
inside the domain, into two fresh surfaces: (values before after), the part
below AT and the part above.  Where AT is a knot of multiplicity above the
degree (an invalid knot vector) NLib writes nothing and this returns
(values nil nil).
Port: srf-split."
  (let ((dir (%surface-direction-flag direction))
        (at (%real at)))
    (with-c-float-environment
      (srf-split surface at dir))))

;;; ------------------------------------------------------------------ queries

(defun surface-normal (surface u v &key (u-side :left) (v-side :left))
  "The point of SURFACE at (U, V) and the unit normal there: (values point
normal u-tangent v-tangent), fresh, the last two the unit first derivatives.
At a pole, where a first derivative vanishes, the normal comes from the
mixed derivative (and the vanished tangent is zero); signals nurbs-error
:geo-err where no normal exists.  U-SIDE and V-SIDE (:left or :right) pick
the side at a knot.
Port: srf-eval-pt-pt-deriv-normal (srf-eval-pt-normal-at-pole at a pole)."
  (let ((u (%real u)) (v (%real v))
        (ufl (side-flag u-side)) (vfl (side-flag v-side))
        (p (make-point)) (su (make-point)) (sv (make-point)) (n (make-point)))
    (with-c-float-environment
      (srf-eval-pt-pt-deriv-normal surface u v ufl vfl p su sv n))
    (values p n su sv)))

(defun surface-grid (surface us vs &key (u-side :left) (v-side :left))
  "The points of SURFACE at every pair of parameters from US and VS
(ascending sequences of reals in the domain), a fresh point-net: entry
[i][j] is the point at (US[i], VS[j]).  U-SIDE and V-SIDE (:left or :right)
pick the side at a knot.
Port: srf-eval-pt-grid."
  (let* ((u (%surface-real-vector us))
         (v (%surface-real-vector vs))
         (nn (1- (length u))) (mm (1- (length v))))
    (declare (type real-vector u v) (type fixnum nn mm))
    (with-c-float-environment
      (srf-eval-pt-grid surface u v nn mm (side-flag u-side) (side-flag v-side)
                        (make-point-net (max 1 (1+ nn)) (max 1 (1+ mm)))))))

(defun surface-closest-point (surface point &key u v (tolerance *mtol*) (cosine-tolerance 1d-5)
                                                 (error-p t))
  "The point of SURFACE closest to POINT (any three reals): (values u v
closest-point distance converged-p), the point fresh.  Newton iteration from
U and V when given, otherwise from the nearest of a grid of 2(n+1)+1 by
2(m+1)+1 evenly spaced samples of the surface (n+1, m+1 its control point
counts).  TOLERANCE is the point coincidence tolerance and COSINE-TOLERANCE
the zero cosine tolerance of the iteration.  When the iteration does not
converge, signals nurbs-error :con-err, or with ERROR-P nil returns the best
point it found and nil for CONVERGED-P.
Port: get-closest-pt-on-srf (srf-eval-pt-grid for the start)."
  (let ((p (%point point))
        (top (%real tolerance)) (toc (%real cosine-tolerance))
        (q (make-point))
        (uv (make-real-vector 2)))
    (with-c-float-environment
      (multiple-value-bind (su sv) (if (and u v)
                                       (values (%real u) (%real v))
                                       (%surface-closest-start surface p))
        (let ((u0 (if u (%real u) su))
              (v0 (if v (%real v) sv)))
          (setf (aref uv 0) u0 (aref uv 1) v0)
          (handler-bind ((nurbs-error
                           (lambda (e)
                             (when (and (not error-p) (eq (nurbs-error-code e) :con-err))
                               (return-from surface-closest-point
                                 (values (aref uv 0) (aref uv 1) q (dist-pt-pt q p) nil))))))
            (multiple-value-bind (cu cv)
                (get-closest-pt-on-srf surface p u0 v0 top toc q (make-point-net 3 3) uv)
              (values cu cv q (dist-pt-pt q p) t))))))))

(defun isocurve (surface &key u v)
  "The iso-parametric curve of SURFACE at U = U (give :u), running in v, or
at V = V (give :v), running in u, as a fresh curve.  Give exactly one.
Port: srf-extract-iso-crv."
  (unless (if u (null v) v)
    (error "ISOCURVE takes exactly one of :U and :V."))
  (with-c-float-environment
    (if u
        (srf-extract-iso-crv surface (%real u) +vdir+)
        (srf-extract-iso-crv surface (%real v) +udir+))))

(defun subsurface (surface u0 u1 v0 v1)
  "The part of SURFACE over [U0, U1] x [V0, V1], as a fresh surface.  A
corner parameter within 1e-8 of a knot is snapped to it.
Port: srf-extract-patch."
  (let ((u0 (%real u0)) (u1 (%real u1)) (v0 (%real v0)) (v1 (%real v1)))
    (with-c-float-environment
      (srf-extract-patch surface u0 u1 v0 v1))))

(defun bezier-patches (surface)
  "The Bezier patches of SURFACE, as a fresh 2D array of surfaces: element
[i][j] is the i-th patch in u and j-th in v, each of SURFACE's degrees with
clamped knots over its own span.
Port: srf-decompose-to-bez."
  (with-c-float-environment
    (values (srf-decompose-to-bez surface))))

(defun offset-surface (surface distance &key (tolerance *mtol*) (type-tolerance *mtol*)
                                            (degree-u 3) (degree-v 3)
                                            (parameterization :inherited))
  "An approximation, to within TOLERANCE, of the surface at DISTANCE from
SURFACE along its normal, as a fresh surface.  SURFACE is first classified to
within TYPE-TOLERANCE (see surface-type): a plane is offset exactly; a
surface of revolution or an extrusion of a planar profile offsets its
profile; anything else is sampled and the offset points interpolated with
degrees DEGREE-U, DEGREE-V at SURFACE's own parameters (PARAMETERIZATION
:inherited), or approximated at :uniform, :chord-length or :centripetal
ones, then knots are removed within TOLERANCE.
Port: srf-offset."
  (let ((d (%real distance)) (tol (%real tolerance)) (eps (%real type-tolerance))
        (par (ecase parameterization
               (:inherited +inherited+) (:uniform +uniform+)
               (:chord-length +chordlength+) (:centripetal +centripetal+))))
    (with-c-float-environment
      (srf-offset surface d degree-u degree-v par eps tol))))

(defun surface-swap-uv (surface)
  "A fresh copy of SURFACE with the roles of u and v exchanged: the result at
(v, u) is SURFACE at (u, v).
Port: srf-copy, swap-uv."
  (with-c-float-environment
    (swap-uv (srf-copy surface))))

(defun surface-closed-p (surface &key direction)
  "True when SURFACE is closed in DIRECTION: :u (its first and last rows of
control points coincide, to *mtol*), :v (its first and last columns), or nil
for either.
Port: srf-is-closed."
  (with-c-float-environment
    (if direction
        (srf-is-closed surface (%surface-direction-flag direction))
        (or (srf-is-closed surface +udir+) (srf-is-closed surface +vdir+)))))

(defun surface-flat-p (surface &key (tolerance *mtol*))
  "True when every control point of SURFACE lies within TOLERANCE of the
plane through the centroid of its corners (a surface closed in u or v, or
one whose corners span no plane, is not flat).
Port: srf-is-flat-cheap."
  (let ((eps (%real tolerance)))
    (with-c-float-environment
      (= (srf-is-flat-cheap surface eps) +yes+))))

(defun surface-type (surface &key (tolerance *mtol*))
  "What SURFACE is, to within TOLERANCE: (values type properties).  TYPE is
:plane, :cylinder, :cone, :sphere, :torus, :revolved, :extruded, :ruled or
:freeform, tried in that order (a plane first, then a surface of revolution,
then a ruled one); PROPERTIES a plist of fresh objects describing it:
  :plane      :corners, the corner points (0,0) (n,0) (n,m) (0,m)
  :cylinder, :cone  :base-center, :axis, :x-axis, :y-axis, :radius,
              :end-radius, :height, :angle (of revolution, degrees), :profile
  :sphere, :torus   :axis-point, :axis, :angle, :center, :normal, :x-axis,
              :y-axis (of the profile circle), :radius (its radius),
              :axis-distance (of its center), :profile-sweep (degrees), :profile
  :revolved   :axis-point, :axis, :angle, :profile-normal (zero unless the
              profile is planar), :profile
  :extruded   :vector, :length, :profile
  :ruled      :profile, :opposite (the two boundary curves)
and for all but a plane :direction, NLib's construction direction (:u or :v).
Port: srf-type."
  (let ((tol (%real tolerance))
        (aa (make-point)) (bb (make-point)) (cc (make-point))
        (dd (make-point)) (xx (make-point)) (yy (make-point)))
    (with-c-float-environment
      (multiple-value-bind (cur-a cur-b r1 r2 ac al h stp dir)
          (srf-type surface tol aa bb cc dd xx yy)
        (%surface-type-properties stp dir aa bb cc dd xx yy r1 r2 ac al h cur-a cur-b)))))
