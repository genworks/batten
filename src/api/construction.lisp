;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: Making curves and surfaces: from points and knots,
;;;; lines and arcs and conics, interpolation and approximation, and the classic
;;;; constructions (ruled, extruded, revolved, skinned, Gordon, Coons).

(in-package :batten)

;;; ------------------------------------------------------------ input coercions

;;; What this section accepts beyond %real, %point and %point-array: point
;;; arrays as well as sequences of points, nets as 3D arrays or sequences of
;;; rows, and NLib's flags under keyword names.  Each returns a fresh object,
;;; so the port routines may write into it.

(defun %construction-points (points)
  "POINTS -- a sequence of points (as %point takes them), or a 2D array of
rows of three reals -- as a fresh point-array."
  (if (and (arrayp points) (= (array-rank points) 2))
      (let* ((n (array-dimension points 0))
             (array (make-point-array n)))
        (unless (= (array-dimension points 1) 3)
          (nurbs-error :inp-err '%construction-points))
        (dotimes (i n array)
          (dotimes (k 3)
            (setf (aref array i k) (float (aref points i k) 1d0)))))
      (%point-array points)))

(defun %construction-net (net routine)
  "NET -- a sequence of rows, each a sequence of points, or a 3D array of
rows of points -- as a fresh point-net: entry [i][j] is row I's point J.
Signals :inp-err (for ROUTINE) when the rows differ in length."
  (if (and (arrayp net) (= (array-rank net) 3))
      (let* ((rows (array-dimension net 0))
             (cols (array-dimension net 1))
             (out (make-point-net rows cols)))
        (unless (= (array-dimension net 2) 3)
          (nurbs-error :inp-err routine))
        (dotimes (i rows out)
          (dotimes (j cols)
            (dotimes (k 3)
              (setf (aref out i j k) (float (aref net i j k) 1d0))))))
      (let* ((rows (length net))
             (cols (if (zerop rows) 0 (length (elt net 0))))
             (out (make-point-net rows cols))
             (i 0))
        (declare (type fixnum rows cols i))
        (map nil #'(lambda (row)
                     (unless (= (length row) cols)
                       (nurbs-error :inp-err routine))
                     (let ((j 0))
                       (declare (type fixnum j))
                       (map nil #'(lambda (p)
                                    (let ((p (%point p)))
                                      (dotimes (k 3)
                                        (setf (aref out i j k) (aref p k))))
                                    (incf j))
                            row))
                     (incf i))
             net)
        out)))

(defun %construction-reals (reals)
  "REALS, a sequence of reals, as a fresh real-vector."
  (map 'real-vector #'%real reals))

(defun %construction-curves (curves)
  "Fresh copies of the curves in the sequence CURVES, as a simple-vector."
  (map 'simple-vector #'crv-copy curves))

(defun %construction-parameterization (parameterization)
  (ecase parameterization
    (:chord-length +chordlength+)
    (:centripetal +centripetal+)
    (:uniform +uniform+)))

(defun %construction-direction (direction)
  (ecase direction (:u +udir+) (:v +vdir+)))

(defun %construction-circle-kind (degree)
  "NLib's flag for a rational circle of DEGREE: +quadratic+, +quartic+ or
+quintic+ (2, 4 or 5)."
  (ecase degree (2 +quadratic+) (4 +quartic+) (5 +quintic+)))

(defun %construction-knot-kind (knots)
  (ecase knots (:single +single+) (:multiple +multiple+)))

(defun %construction-flag (generalized-boolean)
  (if generalized-boolean +yes+ +no+))

;;; ---------------------------------------------------------------- from parts

(defun make-surface (degree-u degree-v knots-u knots-v net &key weights)
  "A NURBS surface of degrees DEGREE-U and DEGREE-V on the knot vectors
KNOTS-U and KNOTS-V, with the control net NET: a sequence of rows (one per
control point in u), each a sequence of points (one per control point in v),
or a 3D array of them.  With WEIGHTS, a sequence of rows of reals shaped like
NET, the surface is rational: control points (w*x w*y w*z w); without, w is
+no-w+.  Signals :inp-err when the counts of knots, control points and
weights do not agree.
Port: %make-surface."
  (let* ((net (%construction-net net 'make-surface))
         (n+1 (array-dimension net 0))
         (m+1 (array-dimension net 1))
         (knots-u (%construction-reals knots-u))
         (knots-v (%construction-reals knots-v))
         (cpoints (make-cpoint-net n+1 m+1)))
    (declare (type point-net net) (type fixnum n+1 m+1) (type cpoint-net cpoints))
    (unless (and (= (length knots-u) (+ n+1 degree-u 1))
                 (= (length knots-v) (+ m+1 degree-v 1)))
      (nurbs-error :inp-err 'make-surface))
    (let ((w (when weights
               (unless (= (length weights) n+1)
                 (nurbs-error :inp-err 'make-surface))
               (let ((w (make-real-matrix n+1 m+1)) (i 0))
                 (declare (type fixnum i))
                 (map nil #'(lambda (row)
                              (unless (= (length row) m+1)
                                (nurbs-error :inp-err 'make-surface))
                              (let ((j 0))
                                (declare (type fixnum j))
                                (map nil #'(lambda (x) (setf (aref w i j) (%real x)) (incf j)) row))
                              (incf i))
                      weights)
                 w))))
      (dotimes (i n+1)
        (dotimes (j m+1)
          (if w
              (let ((wij (aref w i j)))
                (setf (aref cpoints i j 0) (* wij (aref net i j 0))
                      (aref cpoints i j 1) (* wij (aref net i j 1))
                      (aref cpoints i j 2) (* wij (aref net i j 2))
                      (aref cpoints i j 3) wij))
              (setf (aref cpoints i j 0) (aref net i j 0)
                    (aref cpoints i j 1) (aref net i j 1)
                    (aref cpoints i j 2) (aref net i j 2)
                    (aref cpoints i j 3) +no-w+)))))
    (%make-surface degree-u degree-v knots-u knots-v cpoints)))

;;; ------------------------------------------------------- lines, arcs, conics

(defun line-segment (start end)
  "The straight line from START to END, a degree-1 curve on [0, 1].  Its end
control point is END exactly (NLib builds it as START + (END - START)).
Port: crv-line-from-pt-and-vector."
  (let ((p (%point start))
        (q (%point end)))
    (declare (type point p q))
    (with-c-float-environment
      (let* ((v (make-point (- (aref q 0) (aref p 0)) (- (aref q 1) (aref p 1)) (- (aref q 2) (aref p 2))))
             (curve (crv-line-from-pt-and-vector p v))
             (pw (curve-cpoints curve)))
        (declare (type point v) (type cpoint-array pw))
        (setf (aref pw 1 0) (aref q 0) (aref pw 1 1) (aref q 1) (aref pw 1 2) (aref q 2))
        curve))))

(defun circular-arc (center radius &key (x-axis '(1 0 0)) (y-axis '(0 1 0))
                                        (start-angle 0) (end-angle 360) (degree 2)
                                        (angle-tolerance 0.01d0))
  "The circle, or circular arc, of CENTER and RADIUS in the plane of X-AXIS
and Y-AXIS (orthogonal; their lengths do not matter), from START-ANGLE to
END-ANGLE in DEGREES, measured from X-AXIS towards Y-AXIS.  DEGREE is 2 (a
quadratic rational curve of up to four arcs), 4 (quartic, its parameter
closer to arc length) or 5 (quintic, no interior knots).  A sweep within
ANGLE-TOLERANCE degrees of 360 makes a closed circle.
Port: create-circ-arc-tol."
  (let ((c (%point center))
        (x (copy-seq (%point x-axis)))
        (y (copy-seq (%point y-axis))))
    (with-c-float-environment
      (create-circ-arc-tol c x y (%real radius) (%real start-angle) (%real end-angle)
                           (%real angle-tolerance) (%construction-circle-kind degree)))))

(defun elliptical-arc (center major-radius minor-radius
                       &key (x-axis '(1 0 0)) (y-axis '(0 1 0))
                            (start-angle 0) (end-angle 360) (angles-on :circle) (degree 2)
                            (angle-tolerance 0.01d0))
  "The ellipse, or elliptical arc, of CENTER with MAJOR-RADIUS along X-AXIS
and MINOR-RADIUS along Y-AXIS, from START-ANGLE to END-ANGLE in DEGREES.
ANGLES-ON :circle takes the angles on the circle of the major radius (the
eccentric anomaly), :ellipse on the ellipse itself (the polar angle).  The
circle made by circular-arc (DEGREE 2, 4 or 5, ANGLE-TOLERANCE) is scaled
onto the ellipse.
Port: create-elliptical-arc-ang-type-tol."
  (let ((c (%point center))
        (x (copy-seq (%point x-axis)))
        (y (copy-seq (%point y-axis))))
    (with-c-float-environment
      (create-elliptical-arc-ang-type-tol c x y (%real major-radius) (%real minor-radius)
                                          (%real start-angle) (%real end-angle)
                                          (%real angle-tolerance)
                                          (ecase angles-on (:circle 0) (:ellipse 1))
                                          (%construction-circle-kind degree)))))

(defun conic-arc (start start-tangent end end-tangent through)
  "The conic arc from START, tangent to START-TANGENT, to END, tangent to
END-TANGENT, passing through the point THROUGH: a rational quadratic curve of
one, two or four segments.
Port: create-conic-arc."
  (let ((p0 (%point start)) (t0 (%point start-tangent))
        (p2 (%point end)) (t2 (%point end-tangent))
        (p (%point through)))
    (with-c-float-environment
      (create-conic-arc p0 t0 p2 t2 p))))

;;; ------------------------------------------------------------- curve fitting

(defun interpolate-curve (points &key (degree 3) (parameterization :chord-length) parameters
                                      start-tangent end-tangent start-derivative end-derivative)
  "The non-rational curve of DEGREE through POINTS (a sequence of points or a
point array), at parameters by PARAMETERIZATION (:chord-length, :centripetal
or :uniform, on [0, 1]) or at the given PARAMETERS (a sequence of reals, one
per point), with knots averaged from the parameters.  At either end the
curve may be given a direction (START-TANGENT, END-TANGENT: scaled by NLib to
the length of the data polygon) or a derivative used as it is
(START-DERIVATIVE, END-DERIVATIVE), not both kinds.  Signals :inp-err for a
degree above the number of points less one, or below 2 with end conditions
or PARAMETERS.
Port: fit-crv-interp, fit-crv-knots-and-tangents."
  (when (and (or start-tangent end-tangent) (or start-derivative end-derivative))
    (nurbs-error :inp-err 'interpolate-curve))
  (let* ((pts (%construction-points points))
         (n (1- (array-dimension pts 0)))
         (par (%construction-parameterization parameterization))
         (ts (or start-tangent start-derivative))
         (te (or end-tangent end-derivative)))
    (declare (type point-array pts) (type fixnum n par))
    (with-c-float-environment
      (if (or ts te parameters)
          (let ((u (if parameters
                       (%construction-reals parameters)
                       (fit-calc-crv-param-values pts n +epoint+ par (make-real-vector (1+ n))))))
            (unless (= (length u) (1+ n))
              (nurbs-error :inp-err 'interpolate-curve))
            (values (fit-crv-knots-and-tangents pts n u
                                                (and ts (copy-seq (%point ts)))
                                                (and te (copy-seq (%point te)))
                                                (if (or start-derivative end-derivative)
                                                    +derivative+
                                                    +tangent+)
                                                nil 0d0 degree)))
          (fit-crv-interp pts n degree par)))))

(defun approximate-curve (points &key (degree 3) (tolerance *mtol*) (start-degree 1) closed
                                      (knots :multiple) control-points
                                      (parameterization :chord-length))
  "A non-rational curve of DEGREE near POINTS (a sequence of points or a point
array), in one of two ways.  Without CONTROL-POINTS: within TOLERANCE of every
point, with as few control points as NLib's fit finds -- the interpolant of
START-DEGREE (chord length parameters) has knots removed, and is refitted one
degree higher on single or multiple knots (KNOTS :single or :multiple) until
DEGREE; with CLOSED and the first point equal to the last, the curve closes
with matching end derivatives (START-DEGREE at least 2).  With
CONTROL-POINTS, a count: the least-squares curve of that many control points
through the first and last points, at parameters by PARAMETERIZATION
(:chord-length, :centripetal or :uniform), its knots set by the parameter
clusters.
Port: fit-crv-approx; fit-calc-crv-param-values, fit-calc-knot-vector-crv-approx, fit-crv-approx-knots."
  (let* ((pts (%construction-points points))
         (m (1- (array-dimension pts 0))))
    (declare (type point-array pts) (type fixnum m))
    (with-c-float-environment
      (if control-points
          (let* ((n (1- control-points))
                 (u (fit-calc-crv-param-values pts m +epoint+
                                               (%construction-parameterization parameterization)
                                               (make-real-vector (1+ m))))
                 (knt (fit-calc-knot-vector-crv-approx u m n degree (make-real-vector (+ n degree 2)))))
            (fit-crv-approx-knots pts m +epoint+ u knt n degree))
          (fit-crv-approx pts m (%construction-flag closed) start-degree degree (%real tolerance)
                          (%construction-knot-kind knots))))))

(defun fit-arc (points &key (tolerance *mtol*) (error-measure :maximum) (arc t) (degree 2))
  "The least-squares circle or circular arc to POINTS (a sequence of points
or a point array) in the plane z = 0 -- or the least-squares line segment
when the points are collinear (their :maximum or :average distance from the
line, ERROR-MEASURE, below TOLERANCE) or the line fits better on average.
ARC true gives the arc spanning the points (leaving out the largest gap in
their angles about the center), false the full circle; DEGREE is the
circle's, as circular-arc takes it.  Returns (values curve kind
average-error maximum-error), KIND :line or :circle.
Port: fit-arc-to-pts."
  (let* ((pts (%construction-points points))
         (n (1- (array-dimension pts 0))))
    (declare (type point-array pts) (type fixnum n))
    (multiple-value-bind (curve kind average maximum)
        (with-c-float-environment
          (fit-arc-to-pts pts n (%real tolerance)
                          (ecase error-measure (:maximum +maximum+) (:average +average+))
                          (%construction-flag arc) (%construction-circle-kind degree)))
      (values curve (if (= kind +nline+) :line :circle) average maximum))))

(defun fit-line (points)
  "The best fitting line to POINTS (a sequence of points or a point array):
(values point direction average-distance maximum-distance), the point the
centroid, the direction a unit vector, and the distances the points'.
Port: line-fit-pts."
  (let* ((pts (%construction-points points))
         (n (1- (array-dimension pts 0)))
         (c (make-point))
         (v (make-point)))
    (declare (type point-array pts) (type fixnum n))
    (multiple-value-bind (average maximum) (with-c-float-environment (line-fit-pts pts n c v))
      (values c v average maximum))))

;;; ----------------------------------------------------------- surface fitting

(defun interpolate-surface (net &key (degree-u 3) (degree-v 3) (parameterization :chord-length)
                                     parameters-u parameters-v knots-u knots-v)
  "The non-rational surface of DEGREE-U and DEGREE-V through the grid of
points NET (rows in u of points in v: a sequence of sequences or a 3D array),
each row and column on an isoparameter line.  The parameters are by
PARAMETERIZATION (:chord-length, :centripetal or :uniform) unless given
(PARAMETERS-U, one per row, PARAMETERS-V, one per column), the knots averaged
from them unless given (KNOTS-U, KNOTS-V).  Signals :inp-err for a degree
above the number of rows (or columns) less one.
Port: fit-srf-to-pts; fit-calc-srf-param-values, fit-crv-calc-knot-vector, fit-srf-to-pts-knots."
  (let* ((pnet (%construction-net net 'interpolate-surface))
         (n (1- (array-dimension pnet 0)))
         (m (1- (array-dimension pnet 1)))
         (par (%construction-parameterization parameterization)))
    (declare (type point-net pnet) (type fixnum n m par))
    (with-c-float-environment
      (if (or parameters-u parameters-v knots-u knots-v)
          (let ((u (make-real-vector (1+ n)))
                (v (make-real-vector (1+ m))))
            (unless (and parameters-u parameters-v)
              (fit-calc-srf-param-values pnet n m +epoint+ par u v))
            (when parameters-u (setf u (%construction-reals parameters-u)))
            (when parameters-v (setf v (%construction-reals parameters-v)))
            (unless (and (= (length u) (1+ n)) (= (length v) (1+ m)))
              (nurbs-error :inp-err 'interpolate-surface))
            (let ((knu (if knots-u
                           (%construction-reals knots-u)
                           (fit-crv-calc-knot-vector u n degree-u (make-real-vector (+ n degree-u 2)))))
                  (knv (if knots-v
                           (%construction-reals knots-v)
                           (fit-crv-calc-knot-vector v m degree-v (make-real-vector (+ m degree-v 2))))))
              (fit-srf-to-pts-knots pnet n m +epoint+ u v knu knv degree-u degree-v)))
          (fit-srf-to-pts pnet n m degree-u degree-v par)))))

(defun approximate-surface (net &key (degree-u 3) (degree-v 3) (tolerance *mtol*)
                                     (start-degree-u 1) (start-degree-v 1)
                                     (knots :multiple) (knots-from-removal t)
                                     control-points-u control-points-v
                                     (parameterization :chord-length))
  "A non-rational surface of DEGREE-U and DEGREE-V near the grid of points
NET (rows in u of points in v: a sequence of sequences or a 3D array), each
row and column of points on an isoparameter line, in one of two ways.
Without control point counts: within TOLERANCE of every point -- the
interpolant of START-DEGREE-U, START-DEGREE-V (chord length parameters) has
knots removed and is refitted a degree higher, on single or multiple knots
(KNOTS :single or :multiple) taken from the removal (KNOTS-FROM-REMOVAL) or
from the parameters, until the degrees are reached.  With
CONTROL-POINTS-U and/or CONTROL-POINTS-V (counts; the other defaults to the
number of rows or columns): the least-squares surface of that many control
points, at parameters by PARAMETERIZATION (:chord-length, :centripetal or
:uniform).
Port: fit-srf-approx-tol, fit-srf-lst-sq-approx."
  (let* ((pnet (%construction-net net 'approximate-surface))
         (r (1- (array-dimension pnet 0)))
         (s (1- (array-dimension pnet 1))))
    (declare (type point-net pnet) (type fixnum r s))
    (with-c-float-environment
      (if (or control-points-u control-points-v)
          (fit-srf-lst-sq-approx pnet r s
                                 (if control-points-u (1- control-points-u) r)
                                 (if control-points-v (1- control-points-v) s)
                                 degree-u degree-v
                                 (%construction-parameterization parameterization))
          (fit-srf-approx-tol pnet r s start-degree-u degree-u start-degree-v degree-v
                              (%real tolerance) (%construction-knot-kind knots)
                              (%construction-flag knots-from-removal))))))

;;; ---------------------------------------------------- classic constructions

(defun ruled-surface (curve-1 curve-2 &key (across :v))
  "The ruled surface between CURVE-1 and CURVE-2, linear in the parameter
direction ACROSS (:u or :v), where CURVE-1 is at 0 and CURVE-2 at 1 (knots
0 0 1 1); the curves run in the other direction, made compatible (on copies:
the arguments are left alone).
Port: create-ruled-srf-from-boundary-crvs."
  (let ((c1 (crv-copy curve-1))
        (c2 (crv-copy curve-2)))
    (with-c-float-environment
      (create-ruled-srf-from-boundary-crvs c1 c2 (%construction-direction across)))))

(defun extruded-surface (curve direction length &key (across :v))
  "CURVE swept the distance LENGTH along DIRECTION (a vector; its length does
not matter): a surface linear in the parameter direction ACROSS (:u or :v),
knots 0 0 1 1, the curve at 0; the curve runs in the other direction.
Port: create-srf-extrude-crv."
  (let ((c (crv-copy curve))
        (w (copy-seq (%point direction))))
    (with-c-float-environment
      (create-srf-extrude-crv c w (%real length) (%construction-direction across)))))

(defun revolved-surface (curve axis-point axis-direction &key (angle 360) (degree 2))
  "CURVE revolved about the axis through AXIS-POINT along AXIS-DIRECTION by
ANGLE degrees (clockwise looking along AXIS-DIRECTION): the circles run in
u, rational of DEGREE (2 or 5; 4 below 241 degrees), the curve in v.
Port: create-revolved-srf."
  (let ((c (crv-copy curve))
        (s (copy-seq (%point axis-point)))
        (tv (copy-seq (%point axis-direction))))
    (with-c-float-environment
      (create-revolved-srf c s tv (%real angle) (%construction-circle-kind degree)))))

(defun skinned-surface (curves &key degree (across :v) closed (make-compatible t)
                                    (knot-tolerance 0d0) parameters knots)
  "The surface through the cross-sectional CURVES (a sequence of two or
more), skinned in the parameter direction ACROSS (:u or :v) with DEGREE
(default 3, or one less than the number of curves if that is less); the
curves run in the other direction.  With MAKE-COMPATIBLE they are made
compatible first (on copies), knots closer than KNOT-TOLERANCE merged (0: a
thousandth of the domain); otherwise they must be compatible.  The curves
sit at chord-length parameters unless PARAMETERS are given (one per curve;
KNOTS, the knot vector across, are then averaged from them unless given).
With CLOSED and the first and last curve the same, the surface closes
smoothly (then give KNOTS of K+DEGREE+4 values with PARAMETERS).
Port: create-skin-srf, create-skin-srf-params."
  (let* ((cs (%construction-curves curves))
         (k (1- (length cs)))
         (deg (or degree (min 3 k)))
         (dir (%construction-direction across)))
    (declare (type simple-vector cs) (type fixnum k deg dir))
    (with-c-float-environment
      (if (or parameters knots)
          (let* ((tt (if parameters
                         (%construction-reals parameters)
                         (nurbs-error :inp-err 'skinned-surface)))
                 (knt (if knots
                          (%construction-reals knots)
                          (fit-crv-calc-knot-vector tt k deg (make-real-vector (+ k deg 2))))))
            (unless (= (length tt) (1+ k))
              (nurbs-error :inp-err 'skinned-surface))
            (create-skin-srf-params cs k (%construction-flag make-compatible)
                                    (%construction-flag closed) deg dir tt knt))
          (create-skin-srf cs k (%construction-flag make-compatible) (%construction-flag closed)
                           deg dir (%real knot-tolerance))))))

(defun gordon-surface (u-curves vs v-curves us &key degree-u degree-v tensor-degree-u tensor-degree-v
                                                     (preserve-parameters t))
  "The Gordon surface through a network of non-rational curves: U-CURVES
(each set compatible), U-CURVES[i] lying at v = VS[i], and V-CURVES,
V-CURVES[j] at u = US[j].  The v-curves are lofted in u with DEGREE-U
(default 3, or one less than their number if that is less), the u-curves in
v with DEGREE-V, and the interpolant of the network's points with
TENSOR-DEGREE-U, TENSOR-DEGREE-V (defaults: the lofting degrees) is
subtracted.  With PRESERVE-PARAMETERS the lofts and the interpolant use US
and VS; otherwise chord-length parameters.
Port: create-gordon-srf."
  (let* ((cu (%construction-curves u-curves))
         (cv (%construction-curves v-curves))
         (k (1- (length cu)))
         (l (1- (length cv)))
         (v-prms (%construction-reals vs))
         (u-prms (%construction-reals us))
         (dlu (or degree-u (min 3 l)))
         (dlv (or degree-v (min 3 k))))
    (declare (type simple-vector cu cv) (type fixnum k l dlu dlv))
    (unless (and (= (length v-prms) (1+ k)) (= (length u-prms) (1+ l)))
      (nurbs-error :inp-err 'gordon-surface))
    (with-c-float-environment
      (create-gordon-srf cu v-prms k cv u-prms l dlu dlv
                         (or tensor-degree-u dlu) (or tensor-degree-v dlv)
                         (%construction-flag preserve-parameters)))))

(defun coons-surface (boundary derivatives &key check-corners)
  "The bicubically blended Coons surface of four non-rational BOUNDARY
curves, a sequence (left right bottom top) -- left and right on one
v-domain, at u = start and end; bottom and top on one u-domain -- and their
cross-boundary derivative curves DERIVATIVES, in the same order.  With
CHECK-CORNERS the derivatives must match the boundary tangents at the
corners, and the twists from either side agree.
Port: create-coons-boundary-crvs."
  (let ((b (%construction-curves boundary))
        (d (%construction-curves derivatives)))
    (declare (type simple-vector b d))
    (unless (and (= (length b) 4) (= (length d) 4))
      (nurbs-error :inp-err 'coons-surface))
    (with-c-float-environment
      (create-coons-boundary-crvs (svref b 0) (svref d 0) (svref b 1) (svref d 1)
                                  (svref b 2) (svref d 2) (svref b 3) (svref d 3)
                                  (%construction-flag check-corners)))))

(defun sphere (center radius &key (start-angle 0) (end-angle 180) (sweep-angle 360)
                                  (degree-u 2) (degree-v 2))
  "The sphere of CENTER and RADIUS, or a patch of it: the meridian arc from
START-ANGLE to END-ANGLE (degrees from the +z pole, 0 <= start < end <=
180), rational of DEGREE-V (2, 4 or 5), revolved about the z axis through
CENTER by SWEEP-ANGLE degrees with circles of DEGREE-U (2 or 5; 4 below 241
degrees).
Port: create-sphere."
  (let ((c (copy-seq (%point center))))
    (with-c-float-environment
      (create-sphere c (%real radius) (%real start-angle) (%real end-angle) (%real sweep-angle)
                     (%construction-circle-kind degree-v) (%construction-circle-kind degree-u)))))
