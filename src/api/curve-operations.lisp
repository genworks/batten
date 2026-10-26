;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: New curves from old: knots inserted and removed,
;;;; degree raised, split and cut, reversed, reparameterized, extended, offset,
;;;; made compatible, moved, and changed between 2D and 3D and rational and
;;;; polynomial.

(in-package :batten)

;;; Every function here copies before calling a port routine that works in
;;; place (crv-reparam, crv-translate, crv-transform, crv-3d-to-2d ...), so the
;;; caller's curve is never changed and every result is fresh.

;;; ------------------------------------------------------------ keywords

(defun %crv-ops-end (end)
  (ecase end (:start +start+) (:end +end+)))

(defun %crv-ops-continuity (continuity)
  (ecase continuity (:g1 +g-1+) (:g2 +g-2+) (:extrapolate +cmax+)))

(defun %crv-ops-parameterization (parameterization)
  (ecase parameterization
    (:inherited +inherited+)
    (:uniform +uniform+)
    (:chord-length +chordlength+)
    (:centripetal +centripetal+)))

(defun %crv-ops-matrix (matrix)
  "MATRIX, a 4x4 array of reals or a sequence of four rows of four reals, as
a fresh full rmatrix."
  (let ((rm (make-real-matrix 4 4)))
    (declare (type real-matrix rm))
    (if (arrayp matrix)
        (progn
          (unless (equal (array-dimensions matrix) '(4 4))
            (nurbs-error :inp-err 'transform))
          (dotimes (i 4)
            (dotimes (j 4)
              (setf (aref rm i j) (%real (aref matrix i j))))))
        (progn
          (unless (and (= (length matrix) 4) (every (lambda (row) (= (length row) 4)) matrix))
            (nurbs-error :inp-err 'transform))
          (dotimes (i 4)
            (dotimes (j 4)
              (setf (aref rm i j) (%real (elt (elt matrix i) j)))))))
    (make-rmatrix 3 3 rm)))

;;; ------------------------------------------------------------ knots

(defmethod insert-knots ((curve curve) knots &key direction)
  "A fresh copy of CURVE with KNOTS inserted: a sequence of reals strictly
inside the domain, in any order, repeats allowed (a knot's multiplicity
should stay at most the degree).  The shape is unchanged.  No knots gives a
copy.  DIRECTION is for surfaces and ignored here.  Signals :par-err when a
knot is not inside the domain.
Port: crv-refine."
  (declare (ignore direction))
  (let ((x (sort (map 'real-vector #'%real knots) #'<)))
    (declare (type real-vector x))
    (with-c-float-environment
      (if (zerop (length x))
          (crv-copy curve)
          (crv-refine curve x)))))

(defun remove-knots (curve &key (tolerance *mtol*))
  "A fresh copy of CURVE with every knot removed that can go without moving
the curve by more than TOLERANCE (the error bounds accumulate over the
removals; for a rational curve they are scaled by the weights, and a removal
that would take a weight outside [*wmin*, *wmax*] is refused).  Knots go
cheapest first.  A copy when no knot can go.
Port: crv-remove-knots."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-remove-knots curve (%real tolerance))))

;;; ------------------------------------------------------------ degree

(defmethod elevate-degree ((curve curve) by &key direction)
  "A fresh copy of CURVE with its degree raised BY (0 gives a copy); the
shape is unchanged.  DIRECTION is for surfaces and ignored here.  Signals
:deg-err when BY is negative or the new degree exceeds *dmax*, and :knt-err
when CURVE's knot vector is not valid.
Port: crv-elevate-degree."
  (declare (ignore direction))
  (with-c-float-environment
    (crv-elevate-degree curve by)))

;;; ------------------------------------------------------------ pieces

(defun subcurve (curve u0 u1)
  "The piece of CURVE between parameters U0 and U1 (U0 < U1), a fresh curve
on [U0, U1] with clamped ends.  Signals :inp-err when U1 <= U0.
Port: crv-extract-crv-seg."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-extract-crv-seg curve (%real u0) (%real u1))))

(defmethod split ((curve curve) at &key direction)
  "CURVE cut at parameter AT, strictly inside its domain, into two fresh
curves: (values before after), on [u0, AT] and [AT, u1].  DIRECTION is for
surfaces and ignored here.  Signals :inp-err when AT is not inside the
domain.
Port: crv-extract-crv-seg."
  (declare (ignore direction))
  (let* ((knots (curve-knots curve))
         (u0 (aref knots 0))
         (u1 (aref knots (1- (length knots))))
         (at (%real at)))
    (declare (type real-vector knots) (type double-float u0 u1 at))
    (with-c-float-environment
      (values (crv-extract-crv-seg curve u0 at)
              (crv-extract-crv-seg curve at u1)))))

(defun bezier-segments (curve)
  "The Bezier segments of CURVE, in order: a fresh list of degree-p curves,
one for each span of non-zero length, each with knots (a ... a b ... b).
Port: crv-decompose-bez."
  (declare (type curve curve))
  (with-c-float-environment
    (coerce (values (crv-decompose-bez curve)) 'list)))

;;; ------------------------------------------------------------ direction and domain

(defun reverse-curve (curve)
  "CURVE traced the other way, fresh: the same points, its parameter running
from end to start over the same domain.
Port: crv-reverse."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-reverse curve nil)))

(defmethod reparameterize ((curve curve) &rest bounds)
  "A fresh copy of CURVE with its domain mapped linearly onto [U0, U1], BOUNDS
being U0 U1: the same points at the new parameters.
Port: crv-copy, crv-reparam."
  (destructuring-bind (u0 u1) bounds
    (let ((u0 (%real u0)) (u1 (%real u1)))
      (with-c-float-environment
        (crv-reparam (crv-copy curve) u0 u1)))))

;;; ------------------------------------------------------------ extension and offset

(defun extend-curve (curve distance &key (end :end) (continuity :g1) relative)
  "CURVE extended at END (:start or :end) by about the arc length DISTANCE,
or by DISTANCE times CURVE's arc length when RELATIVE, as a fresh curve.
CONTINUITY says how the extension joins: :g1 a straight segment along the end
tangent, :g2 the end segment reflected, :extrapolate the end polynomial piece
carried on (no knot at the join).
Port: crv-extend-by-dist."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-extend-by-dist curve (%real distance) (if relative +relative+ +absolute+)
                        (%crv-ops-end end) (%crv-ops-continuity continuity))))

(defun offset-curve (curve distance &key (normal '(0 0 1)) (tolerance *mtol*) (degree 3)
                                         (parameterization :inherited)
                                         (keep-end-tangents (/= degree 1)))
  "The offset of the planar, tangent-continuous CURVE by DISTANCE, within
TOLERANCE: a fresh curve, each point moved DISTANCE along NORMAL x T (NORMAL
the normal of CURVE's plane, T the unit tangent; a negative DISTANCE goes the
other way).  A line or a circle is offset exactly; otherwise the offset is
sampled, interpolated by a curve of DEGREE at CURVE's own parameters
(PARAMETERIZATION :inherited) or ones computed from the samples (:uniform,
:chord-length, :centripetal), keeping CURVE's end tangent directions when
KEEP-END-TANGENTS, and its knots removed within TOLERANCE.  A 2D CURVE gives
a 2D offset.
Port: crv-offset-pt-sampling."
  (declare (type curve curve))
  (let ((normal (copy-seq (%point normal))))
    (with-c-float-environment
      (crv-offset-pt-sampling curve normal (%real distance) degree
                              (if keep-end-tangents +yes+ +no+)
                              (%crv-ops-parameterization parameterization)
                              (%real tolerance)))))

;;; ------------------------------------------------------------ compatibility

(defun make-compatible (curves &key knot-tolerance)
  "Fresh copies of CURVES (a sequence), made compatible, as a list: all of
the same rationality and dimension, of the highest degree among them, on one
knot vector (each on [0, 1] when their domains differ), every curve's shape
unchanged.  Knots closer than KNOT-TOLERANCE count as one, at their average;
without it, a thousandth of the domain.
Port: crv-copy, crvs-make-compatible, crvs-make-compatible-knot-tol."
  (let ((copies (map 'simple-vector #'crv-copy curves)))
    (declare (type simple-vector copies))
    (when (plusp (length copies))
      (with-c-float-environment
        (if knot-tolerance
            (crvs-make-compatible-knot-tol copies (1- (length copies)) (%real knot-tolerance))
            (crvs-make-compatible copies (1- (length copies))))))
    (coerce copies 'list)))

;;; ------------------------------------------------------------ moving

(defmethod translate ((curve curve) vector)
  "A fresh copy of CURVE moved by VECTOR (any three reals).  A 2D curve
comes back 3D, in the plane z = 0 before the move.
Port: crv-copy, crv-translate."
  (let ((vector (%point vector)))
    (with-c-float-environment
      (crv-translate (crv-copy curve) vector))))

(defmethod transform ((curve curve) matrix)
  "A fresh copy of CURVE with its control points mapped by MATRIX, a 4x4
matrix (a 2D array, or four rows of four reals) applied to homogeneous
points: rotation and scale upper left, translation in the last column,
projection in the last row.  A 2D curve comes back 3D.
Port: crv-copy, crv-transform."
  (let ((rma (%crv-ops-matrix matrix)))
    (with-c-float-environment
      (crv-transform (crv-copy curve) rma))))

;;; ------------------------------------------------------------ 2D, 3D, rational, polynomial

(defun curve-to-2d (curve)
  "A fresh 2D copy of CURVE: its projection onto the xy plane.
Port: crv-copy, crv-3d-to-2d."
  (declare (type curve curve))
  (with-c-float-environment
    (crv-3d-to-2d (crv-copy curve))))

(defun curve-to-3d (curve)
  "A fresh 3D copy of CURVE: a 2D curve placed in the plane z = 0, a 3D one
copied as it is.
Port: crv-copy, crv-2d-to-3d."
  (declare (type curve curve))
  (with-c-float-environment
    (let ((copy (crv-copy curve)))
      (if (crv-is-3d curve)
          copy
          (crv-2d-to-3d copy)))))

(defun curve-to-rational (curve)
  "A fresh rational copy of CURVE: a polynomial curve given weights 1, a
rational one copied as it is.
Port: crv-copy, crv-non-rat-to-rat."
  (declare (type curve curve))
  (with-c-float-environment
    (let ((copy (crv-copy curve)))
      (if (is-crv-rat curve)
          copy
          (crv-non-rat-to-rat copy)))))

(defun curve-to-polynomial (curve &key (tolerance *mtol*) (degree 3)
                                       (parameterization :inherited)
                                       (keep-end-tangents (/= degree 1)))
  "A fresh polynomial (non-rational) curve within TOLERANCE of CURVE, on
CURVE's domain; a polynomial CURVE is copied as it is.  A rational CURVE is
sampled span by span, as densely as TOLERANCE and its curvature need, the
samples interpolated by a curve of DEGREE at CURVE's own parameters
(PARAMETERIZATION :inherited) or ones computed from the samples (:uniform,
:chord-length, :centripetal), keeping CURVE's end tangents when
KEEP-END-TANGENTS, and then as many knots removed as TOLERANCE allows.
TOLERANCE must be positive.
Port: crv-copy, approx-nurbs-with-non-rat-crv."
  (declare (type curve curve))
  (with-c-float-environment
    (if (is-crv-rat curve)
        (approx-nurbs-with-non-rat-crv curve (%real tolerance) degree
                                       (%crv-ops-parameterization parameterization)
                                       (if keep-end-tangents +yes+ +no+))
        (crv-copy curve))))
