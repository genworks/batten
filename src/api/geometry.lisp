;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: Points and vectors under readable names, and the
;;;; scalar B-spline functions (cfun, sfun).

(in-package :batten)

;;; ------------------------------------------------------- points and vectors

;;; A point or a vector is any sequence of three reals (a Batten point is used
;;; as it is).  A line is a point and a direction, L(t) = point + t*direction;
;;; a bounded line is the segment 0 <= t <= 1.  A plane is a point and a
;;; normal.  Results are fresh points.

(defun %geometry-line (point direction bounded)
  "The port's line through POINT along DIRECTION, bounded (0 <= t <= 1) when
BOUNDED.  The port only reads it, so a Batten point is shared, not copied."
  (%make-lineseg (%point point) (%point direction) (if bounded +bounded+ +unbounded+)))

(defun %geometry-plane (point normal)
  "The port's plane through POINT with NORMAL (read only, like %GEOMETRY-LINE)."
  (%make-plane (%point point) (%point normal)))

(defun distance (p q)
  "The Euclidean distance between points P and Q.
Port: dist-pt-pt."
  (let ((p (%point p)) (q (%point q)))
    (declare (type point p q))
    (with-c-float-environment (dist-pt-pt p q))))

(defun dot (a b)
  "The dot product of vectors A and B.
Port: vector-dot."
  (let ((a (%point a)) (b (%point b)))
    (declare (type point a b))
    (with-c-float-environment (vector-dot a b))))

(defun cross (a b)
  "The cross product A x B, a fresh vector.
Port: vector-cross."
  (let ((a (%point a)) (b (%point b)))
    (declare (type point a b))
    (with-c-float-environment (vector-cross a b (make-point)))))

(defun normalize (v)
  "V scaled to unit length, a fresh vector, and its length as a second value:
(values unit length).  Signals nurbs-error :geo-err for a zero vector.
Port: vector-normalize."
  (let ((v (%point v)) (unit (make-point)))
    (declare (type point v unit))
    (let ((magnitude (with-c-float-environment (vector-normalize v unit))))
      (values unit magnitude))))

(defun angle-between (a b)
  "The angle between vectors A and B in degrees, 0 to 180.  Signals
nurbs-error :geo-err when either is a zero vector.
Port: vectors-angle."
  (let ((a (%point a)) (b (%point b)))
    (declare (type point a b))
    (with-c-float-environment (vectors-angle a b))))

(defun project-point-to-line (point line-point line-direction &key bounded)
  "The foot of the perpendicular from POINT to the line through LINE-POINT
along LINE-DIRECTION: (values foot parameter), FOOT a fresh point and
PARAMETER its t on the line, FOOT = LINE-POINT + t*LINE-DIRECTION.  With
BOUNDED the line is the segment 0 <= t <= 1, and a foot outside it gives
NIL.  Signals nurbs-error :geo-err for a zero direction.
Port: project-pt-line."
  (let ((line (%geometry-line line-point line-direction bounded))
        (point (%point point))
        (foot (make-point)))
    (declare (type point point foot))
    (multiple-value-bind (found parameter) (with-c-float-environment (project-pt-line line point foot))
      (if (= found +true+) (values foot parameter) nil))))

(defun project-point-to-plane (point plane-point plane-normal)
  "The projection of POINT onto the plane through PLANE-POINT with normal
PLANE-NORMAL (along the normal), a fresh point.  Signals nurbs-error
:geo-err for a zero normal.
Port: project-pt-plane."
  (let ((plane (%geometry-plane plane-point plane-normal))
        (point (%point point)))
    (declare (type point point))
    (with-c-float-environment (project-pt-plane plane point (make-point)))))

(defun intersect-lines (point1 direction1 point2 direction2 &key bounded1 bounded2)
  "Where the line through POINT1 along DIRECTION1 meets the line through
POINT2 along DIRECTION2, in the least-squares sense: (values point t1 t2),
POINT (fresh) the midpoint of the two lines' closest points, at T1 on the
first line and T2 on the second (Li(t) = POINTi + t*DIRECTIONi).  Skew lines
are thus \"met\" at their closest approach.  NIL when the lines are parallel,
or when a closest point falls outside a line made a segment (0 <= t <= 1) by
BOUNDED1 or BOUNDED2.  Signals nurbs-error :geo-err for a zero direction.
Port: isect-line-line."
  (let ((line1 (%geometry-line point1 direction1 bounded1))
        (line2 (%geometry-line point2 direction2 bounded2))
        (point (make-point)))
    (declare (type point point))
    (multiple-value-bind (found t1 t2) (with-c-float-environment (isect-line-line line1 line2 point))
      (if (= found +true+) (values point t1 t2) nil))))

(defun intersect-line-plane (line-point line-direction plane-point plane-normal &key bounded)
  "Where the line through LINE-POINT along LINE-DIRECTION meets the plane
through PLANE-POINT with normal PLANE-NORMAL: (values point parameter), POINT
fresh and PARAMETER its t on the line, POINT = LINE-POINT +
t*LINE-DIRECTION.  NIL when the line is parallel to the plane or, with
BOUNDED (the segment 0 <= t <= 1), when the segment does not reach it.
Signals nurbs-error :geo-err for a zero direction or normal.
Port: isect-line-plane."
  (let ((line (%geometry-line line-point line-direction bounded))
        (plane (%geometry-plane plane-point plane-normal))
        (point (make-point)))
    (declare (type point point))
    (multiple-value-bind (found parameter) (with-c-float-environment (isect-line-plane line plane point))
      (if (= found +true+) (values point parameter) nil))))

(defun intersect-planes (point1 normal1 point2 normal2)
  "The line where the plane through POINT1 with normal NORMAL1 meets the
plane through POINT2 with normal NORMAL2: (values point direction), two fresh
vectors, DIRECTION the cross product of the unit normals (not itself unit
length).  NIL when the planes are parallel.  Signals nurbs-error :geo-err for
a zero normal.
Port: isect-plane-plane."
  (let ((plane1 (%geometry-plane point1 normal1))
        (plane2 (%geometry-plane point2 normal2))
        (point (make-point))
        (direction (make-point)))
    (declare (type point point direction))
    (if (= (with-c-float-environment (isect-plane-plane plane1 plane2 point direction)) +true+)
        (values point direction)
        nil)))

;;; ------------------------------------------------------ scalar B-spline functions

;;; A cfun is a scalar B-spline function of one parameter, f(u) =
;;; Sum_i N_i,p(u) * value_i; an sfun one of two, f(u, v).  Their knots and
;;; values are real-vectors (a real-matrix of values for an sfun).

(defun %scalar-real-vector (reals)
  "REALS, any sequence of reals, as a fresh real-vector."
  (map 'real-vector #'%real reals))

(defun make-cfun (degree knots values)
  "A scalar B-spline function of degree DEGREE on KNOTS (a sequence of reals,
clamped: DEGREE+1 equal knots at each end) with control VALUES (a sequence of
reals); KNOTS has (length VALUES) + DEGREE + 1 elements.  The arrays are
fresh.  Signals nurbs-error :inp-err for a wrong count and :knt-err for an
invalid knot vector.
Port: knot-vector-is-valid, c-func-from-knots."
  (check-type degree (integer 0 #.+maxdeg+))
  (let ((knots (%scalar-real-vector knots))
        (vals (%scalar-real-vector values)))
    (declare (type real-vector knots vals))
    (unless (= (length knots) (+ (length vals) degree 1))
      (nurbs-error :inp-err 'make-cfun))
    (with-c-float-environment
      (knot-vector-is-valid knots degree 'make-cfun)
      (c-func-from-knots (%make-cfun 0 (make-real-vector 0) (make-real-vector 0))
                         vals (1- (length vals)) degree knots (1- (length knots))))))

(defun cfun-value (f u &key (side :left))
  "The value of scalar function F at U.  At a knot, SIDE (:left or :right)
picks the span on that side.  Signals nurbs-error :par-err for U outside the
domain.
Port: c-func-eval."
  (with-c-float-environment (c-func-eval f (%real u) (side-flag side))))

(defun cfun-derivatives (f u n &key (side :left))
  "The value of scalar function F at U and its derivatives to order N, a
fresh real-vector of N+1: element k is the kth derivative.  At a knot, SIDE
(:left or :right) picks the span on that side.  Signals nurbs-error :par-err
for U outside the domain.
Port: c-func-derivs."
  (check-type n (integer 0 #.most-positive-fixnum))
  (with-c-float-environment (c-func-derivs f (%real u) (side-flag side) n (make-real-vector (1+ n)))))

(defun cfun-derivative (f &key (order 1))
  "The derivative of scalar function F of order ORDER (0 to F's degree), as a
fresh scalar function of degree p - ORDER.  Signals nurbs-error :inp-err for
an ORDER outside 0..p and :der-err when F's knots make the derivative
undefined.
Port: crv-func-deriv."
  (check-type order fixnum)
  (with-c-float-environment (crv-func-deriv f order)))

(defun cfun* (f g)
  "The product of scalar functions F and G, a fresh function of degree p + q.
When their domains differ, the one with the shorter domain is rescaled onto
the other's (in a copy: F and G are untouched).  Signals nurbs-error
:deg-err when p + q exceeds the maximum degree and :inp-err when a degree-0
factor, or the partner of one, has more than one span.
Port: crv-func-copy, crv-func-multiply-crv-func."
  (with-c-float-environment
    (crv-func-multiply-crv-func (crv-func-copy f) (crv-func-copy g))))

(defun cfun+ (f g)
  "The sum of scalar functions F and G, a fresh function: both are raised to
the higher degree and refined to common knots first (in copies), the one with
the shorter domain rescaled onto the other's.
Port: crv-func-sum-diff-crv-func."
  (with-c-float-environment (crv-func-sum-diff-crv-func f g +plus+ +no+)))

(defun cfun- (f g)
  "The difference F - G of scalar functions F and G, a fresh function, made
as CFUN+ makes the sum.
Port: crv-func-sum-diff-crv-func."
  (with-c-float-environment (crv-func-sum-diff-crv-func f g +minus+ +no+)))

(defun %scalar-real-matrix (rows)
  "ROWS, a 2D array of reals or a sequence of equally long sequences of reals,
as a fresh real-matrix."
  (if (and (arrayp rows) (= (array-rank rows) 2))
      (let ((matrix (make-real-matrix (array-dimension rows 0) (array-dimension rows 1))))
        (dotimes (k (array-total-size rows) matrix)
          (setf (row-major-aref matrix k) (%real (row-major-aref rows k)))))
      (let* ((n (length rows))
             (m (if (zerop n) 0 (length (elt rows 0))))
             (matrix (make-real-matrix n m))
             (i 0))
        (map nil #'(lambda (row)
                     (unless (= (length row) m)
                       (nurbs-error :inp-err 'make-sfun))
                     (let ((j 0))
                       (map nil #'(lambda (x) (setf (aref matrix i j) (%real x)) (incf j)) row))
                     (incf i))
             rows)
        matrix)))

(defun make-sfun (degree-u degree-v knots-u knots-v values)
  "A scalar tensor-product B-spline function of degrees DEGREE-U and DEGREE-V
on KNOTS-U and KNOTS-V (sequences of reals, clamped) with control VALUES, a
2D array of reals or a sequence of rows: VALUES[i][j] goes with the ith
u-basis function and the jth v-basis function, so KNOTS-U has (rows) +
DEGREE-U + 1 elements and KNOTS-V (columns) + DEGREE-V + 1.  The arrays are
fresh.  Signals nurbs-error :inp-err for wrong counts and :knt-err for an
invalid knot vector.
Port: knot-vector-is-valid, s-func-from-knot-vectors."
  (check-type degree-u (integer 0 #.+maxdeg+))
  (check-type degree-v (integer 0 #.+maxdeg+))
  (let ((knots-u (%scalar-real-vector knots-u))
        (knots-v (%scalar-real-vector knots-v))
        (vals (%scalar-real-matrix values)))
    (declare (type real-vector knots-u knots-v) (type real-matrix vals))
    (unless (and (= (length knots-u) (+ (array-dimension vals 0) degree-u 1))
                 (= (length knots-v) (+ (array-dimension vals 1) degree-v 1)))
      (nurbs-error :inp-err 'make-sfun))
    (with-c-float-environment
      (knot-vector-is-valid knots-u degree-u 'make-sfun)
      (knot-vector-is-valid knots-v degree-v 'make-sfun)
      (s-func-from-knot-vectors (%make-sfun 0 0 (make-real-vector 0) (make-real-vector 0)
                                            (make-real-matrix 0 0))
                                vals degree-u degree-v knots-u knots-v))))
