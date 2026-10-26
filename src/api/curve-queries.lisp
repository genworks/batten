;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: Asking a curve: domain, tangents, curvature and the
;;;; Frenet frame, length, closest points, bounding box, polylines, its type and
;;;; its predicates.

(in-package :batten)

;;; ------------------------------------------------------- the shared generics

(defmethod domain ((curve curve))
  "The parameter domain of CURVE, (values u0 u1).
Port: crv-get-param-bounds."
  (with-c-float-environment
    (crv-get-param-bounds curve)))

(defmethod bounding-box ((curve curve))
  "The axis-aligned box around CURVE's control points (which contains the
curve), (values min-point max-point), both fresh.
Port: crv-get-b-box."
  (let ((box (make-real-vector 6)))
    (declare (type real-vector box))
    (with-c-float-environment
      (crv-get-b-box curve box))
    ;; NLib's box is (xl xr yb yt zn zf)
    (values (make-point (aref box 0) (aref box 2) (aref box 4))
            (make-point (aref box 1) (aref box 3) (aref box 5)))))

(defmethod copy ((curve curve))
  "A copy of CURVE sharing nothing with it.
Port: crv-copy."
  (with-c-float-environment
    (crv-copy curve)))

;;; ------------------------------------------------- tangent, curvature, frame

(defun curve-tangent (curve u &key (side :left))
  "The unit tangent of CURVE at parameter U, with the point there:
(values tangent point), both fresh.  SIDE (:left or :right) picks the span
at a knot.  Signals a nurbs-error (:geo-err) where the first derivative is
zero.
Port: crv-eval-tangent."
  (declare (type curve curve))
  (let ((p (make-point)) (tv (make-point)))
    (with-c-float-environment
      (crv-eval-tangent curve (%real u) (side-flag side) p tv))
    (values tv p)))

(defun curve-curvature (curve u &key (side :left))
  "The curvature of CURVE at parameter U: |C' x C''| / |C'|^3, the reciprocal
of the osculating circle's radius (0 on a straight stretch).  SIDE (:left or
:right) picks the span at a knot.  Signals a nurbs-error where the first
derivative is zero.
Port: crv-eval-curvature."
  (declare (type curve curve))
  (let ((p (make-point)) (tv (make-point)) (c (make-point)) (nv (make-point)))
    (declare (dynamic-extent p tv c nv))
    (with-c-float-environment
      (crv-eval-curvature curve (%real u) (side-flag side) p tv +no+ c nv))))

(defun curve-frenet-frame (curve u &key (side :left))
  "The Frenet frame of CURVE at parameter U: (values point tangent normal
binormal), all fresh, the three vectors of unit length.  SIDE (:left or
:right) picks the span at a knot.  Signals a nurbs-error for a curve of
degree 1 (:deg-err), and where the first derivative or the curvature is zero
(:geo-err).
Port: crv-eval-frenet-frame."
  (declare (type curve curve))
  (let ((p (make-point)) (tv (make-point)) (nv (make-point)) (bv (make-point)))
    (with-c-float-environment
      (crv-eval-frenet-frame curve (%real u) (side-flag side) p tv nv bv))
    (values p tv nv bv)))

;;; ------------------------------------------------------------------- length

(defun %positive-tolerance (tolerance routine)
  "TOLERANCE as a double-float, or a nurbs-error (:inp-err) for ROUTINE
unless it is positive: NLib subdivides for ever on a zero tolerance."
  (let ((tol (%real tolerance)))
    (declare (type double-float tol))
    (unless (> tol 0d0)
      (nurbs-error :inp-err routine))
    tol))

(defun curve-length (curve &key from to (tolerance 1d-6) (relative t))
  "The arc length of CURVE between parameters FROM and TO (by default its
whole domain).  TOLERANCE bounds the error: a fraction of the measured
stretch's size (the diagonal of its control points' box) when RELATIVE, the
default, else a distance.  A line or circular arc is measured exactly.
Signals a nurbs-error (:inp-err) when TO is before FROM or TOLERANCE is not
positive.
Port: crv-arc-length."
  (declare (type curve curve))
  (let ((tol (%positive-tolerance tolerance 'curve-length)))
    (with-c-float-environment
      (multiple-value-bind (u0 u1) (crv-get-param-bounds curve)
        (crv-arc-length curve (if from (%real from) u0) (if to (%real to) u1)
                        tol (if relative +relative+ +absolute+))))))

(defun curve-point-at-length (curve fraction &key (tolerance 1d-6) (relative t) length)
  "The point of CURVE at FRACTION (0 to 1) of its arc length from the start:
(values point parameter length), the point fresh, LENGTH the curve's whole
length.  TOLERANCE bounds the error in length: a fraction of the curve's size
(the diagonal of its control points' box) when RELATIVE, the default, else a
distance.  Pass the curve's LENGTH when it is known, to save measuring it.
Signals a nurbs-error (:inp-err) for a FRACTION outside [0, 1] or a
TOLERANCE that is not positive.
Port: crv-percentage-along-pt."
  (declare (type curve curve))
  (let ((tol (%positive-tolerance tolerance 'curve-point-at-length))
        (p (make-point)))
    (multiple-value-bind (len up)
        (with-c-float-environment
          (crv-percentage-along-pt curve (%real fraction) tol *mtol*
                                   (if relative +relative+ +absolute+)
                                   (if length (%real length) -1d0) p))
      (values p up len))))

;;; ----------------------------------------------------------- closest points

(defun %polyline-nearest-parameter (points parameters p)
  "The parameter of the point of the polyline POINTS (a point-array of its
vertices, PARAMETERS their parameters) nearest point P, interpolated
linearly along its leg: the first nearest, the first vertex for a polyline
of one."
  (declare (type point-array points) (type real-vector parameters) (type point p))
  (let ((a (make-point)) (b (make-point)) (q (make-point))
        (best-u (aref parameters 0))
        (best-d 0d0))
    (declare (dynamic-extent a b q) (type double-float best-u best-d))
    (set-point a (aref points 0 0) (aref points 0 1) (aref points 0 2))
    (setf best-d (dist-pt-pt p a))
    (dotimes (j (1- (array-dimension points 0)))
      (set-point a (aref points j 0) (aref points j 1) (aref points j 2))
      (set-point b (aref points (1+ j) 0) (aref points (1+ j) 1) (aref points (1+ j) 2))
      (let ((alf (project-pt-line-param p a b q)))
        (declare (type double-float alf))
        ;; the nearest point of the segment, not of its line
        (cond ((< alf 0d0) (setf alf 0d0) (copy-pt a q))
              ((> alf 1d0) (setf alf 1d0) (copy-pt b q)))
        (let ((d (dist-pt-pt p q)))
          (declare (type double-float d))
          (when (< d best-d)
            (setf best-d d
                  best-u (+ (aref parameters j)
                            (* alf (- (aref parameters (1+ j)) (aref parameters j)))))))))
    best-u))

(defun %closest-point-start (curve p)
  "A starting parameter for Newton's closest-point iteration from point P to
CURVE: that of the nearest point of a polyline within a thousandth of the
curve's size (the diagonal of its control points' box; never under *mtol*)."
  (declare (type curve curve) (type point p))
  (let ((tol (nl-max (* 1d-3 (the double-float (crv-get-b-box-max-diag-dist curve)))
                     (the double-float *mtol*))))
    (multiple-value-bind (points parameters) (approx-crv-with-polyline curve tol +absolute+ +both+)
      (%polyline-nearest-parameter points parameters p))))

(defun curve-closest-point (curve point &key (tolerance *mtol*) (cosine-tolerance tolerance) start)
  "The point of CURVE closest to POINT: (values u closest distance), CLOSEST
fresh.  Newton's iteration runs from START, a parameter, or by default from
the point nearest POINT of a polyline through the curve (within a thousandth
of the curve's size), until the curve point is within TOLERANCE of POINT or
the cosine between the tangent and the offset to POINT is within
COSINE-TOLERANCE of zero; it wraps around a closed curve.  Signals a
nurbs-error (:con-err) when the iteration does not converge, as it can where
the nearest point is a corner (a kink, or the seam of a closed curve whose
tangents do not meet).
Port: crv-closest-pt, started by approx-crv-with-polyline (sized by
crv-get-b-box-max-diag-dist) and project-pt-line-param."
  (declare (type curve curve))
  (let ((p (%point point))
        (q (make-point)))
    (with-c-float-environment
      (let* ((u0 (if start (%real start) (%closest-point-start curve p)))
             (u (crv-closest-pt curve p u0 (%real tolerance) (%real cosine-tolerance) q)))
        (values u q (dist-pt-pt p q))))))

(defun curve-project-points (curve points &key (tolerance 1d-3))
  "Project POINTS (a point-array, or a sequence of points, lying along CURVE
in its direction) onto a polyline through CURVE within TOLERANCE times the curve's length:
(values parameters projections), a fresh real-vector and point-array, each
point projected onto the first polyline leg it falls on from the previous
point's leg onward, its parameter interpolated along the leg.  Signals a
nurbs-error (:num-err) when a point falls on no leg, and (:inp-err) for a
TOLERANCE that is not positive.
Port: crv-project-pts."
  (declare (type curve curve))
  (let ((tol (%positive-tolerance tolerance 'curve-project-points))
        ;; a point-array is passed through: crv-project-pts only reads it
        (pts (if (typep points 'point-array) points (%point-array points))))
    (multiple-value-bind (projections parameters)
        (with-c-float-environment
          (crv-project-pts curve pts (1- (array-dimension pts 0)) tol))
      (values parameters projections))))

(defun curve-polyline (curve tolerance &key relative)
  "A polyline through points of CURVE deviating from it by at most TOLERANCE,
a distance (or with RELATIVE, a fraction of the curve's length): (values
points parameters), a fresh point-array and real-vector of the vertices and
their parameters on the curve.  A curve of degree 1 gives its control points.
Signals a nurbs-error (:inp-err) for a TOLERANCE that is not positive.
Port: approx-crv-with-polyline."
  (declare (type curve curve))
  (let ((tol (%positive-tolerance tolerance 'curve-polyline)))
    (multiple-value-bind (points parameters)
        (with-c-float-environment
          (approx-crv-with-polyline curve tol (if relative +relative+ +absolute+) +both+))
      (values points parameters))))

;;; --------------------------------------------------------------- predicates

(defun curve-closed-p (curve)
  "True when CURVE's first and last control points coincide (to within
*mtol*).
Port: crv-is-closed."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-is-closed curve)))

(defun curve-planar-p (curve &key (tolerance *mtol*))
  "When every control point of CURVE lies within TOLERANCE of one plane,
(values t point normal): a point on the plane (the control points' centroid)
and its unit normal, both fresh; else nil.  As in NLib, a curve that is a
point or a straight line, or that has fewer than three control points, is
not planar (its plane is not determined).
Port: crv-is-planar."
  (declare (type curve curve))
  (let ((pp (make-point)) (nn (make-point)))
    (if (= (the fixnum (with-c-float-environment
                         (crv-is-planar curve (%real tolerance) pp nn)))
           +true+)
        (values t pp nn)
        nil)))

(defun curve-line-p (curve &key (tolerance *mtol*))
  "True when CURVE is a straight line: its control points lie within
TOLERANCE of the segment between its end points.
Port: crv-is-line."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-is-line curve (%real tolerance))))

(defun curve-degenerate-p (curve)
  "True when CURVE is a single point: all its control points lie within
*mtol* of their centroid.
Port: crv-is-degen."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-is-degen curve)))

(defun curves-coincident-p (curve-1 curve-2 &key (tolerance *mtol*))
  "True when CURVE-1 and CURVE-2 are the same curve: made compatible (one
degree, one knot vector, the domains mapped onto one), their control points
pair up within TOLERANCE (in homogeneous space).  False as well when they
cannot be made compatible.  The curves themselves are left alone.
Port: crvs-are-coincident (on copies made by crv-copy)."
  (declare (type curve curve-1 curve-2))
  (with-c-float-environment
    (crvs-are-coincident (crv-copy curve-1) (crv-copy curve-2) (%real tolerance))))

;;; --------------------------------------------------------------------- type

(defun curve-type (curve &key (tolerance *mtol*))
  "What CURVE is, to within TOLERANCE, with what describes it:
  (values :line start direction) -- DIRECTION the vector to the end point;
  (values :circle center normal radius sweep x-axis y-axis) -- SWEEP in
    degrees, the unit X-AXIS and Y-AXIS spanning the circle's plane;
  :freeform -- anything else.
The points and vectors are fresh.
Port: crv-get-type."
  (declare (type curve curve))
  (let ((aa (make-point)) (bb (make-point)) (xx (make-point)) (yy (make-point)))
    (multiple-value-bind (radius sweep type)
        (with-c-float-environment
          (crv-get-type curve (%real tolerance) aa bb xx yy))
      (declare (type double-float radius sweep) (type fixnum type))
      (cond ((= type +nline+) (values :line aa bb))
            ((= type +ncircle+) (values :circle aa bb radius sweep xx yy))
            (t :freeform)))))
