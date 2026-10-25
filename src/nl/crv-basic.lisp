;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(declaim (inline is-crv-rat))
(defun is-crv-rat (cur)
  "True when curve CUR is rational (its first control point has a weight).
NLib: N_IsCrvRat."
  (declare (type curve cur))
  (/= (aref (curve-cpoints cur) 0 3) +no-w+))

(defun crv-get-param-bounds (cur)
  "The parameter domain of CUR, (values U[0] U[m]).  NLib: N_CrvGetParamBounds."
  (declare (type curve cur))
  (let ((k (curve-knots cur)))
    (values (aref k 0) (aref k (1- (length k))))))

(defun crv-is-3d (cur)
  "True when curve CUR is three dimensional: its first control point's wz is
not +no-z+.
NLib: N_CrvIs3d."
  (declare (type curve cur))
  (let ((wz (aref (curve-cpoints cur) 0 2)))
    (declare (type double-float wz))
    (/= wz +no-z+)))

(defun crv-from-c-polygon-and-knot-vector (cur pol p knt)
  "Make curve CUR (modified in place) the curve of control polygon POL (a
cpoint-array), degree P and knot vector KNT (a real-vector), sharing the
arrays as the C shares the pointers.  Return CUR.
NLib: N_CrvFromCPolygonAndKnotVector."
  (declare (type curve cur) (type cpoint-array pol) (type fixnum p) (type real-vector knt))
  (setf (curve-cpoints cur) pol
        (curve-degree cur) p
        (curve-knots cur) knt)
  cur)

(defun crv-get-conic-data (ctp as ae)
  "For a conic of type CTP (+quadratic+, +quartic+ or +quintic+) from start
angle AS to end angle AE (degrees), return (values ae n p m): AE adjusted
for the angle tolerance and to lie in (as, as+360], the highest control
point index N, the degree P and the highest knot index M.
NLib: N_CrvGetConicData."
  (declare (type fixnum ctp) (type double-float as ae))
  (let ((angle-tol-deg 1.0d-6)
        (n 0) (p 0))
    (declare (type double-float angle-tol-deg) (type fixnum n p))
    (do () ((not (< ae (+ as angle-tol-deg))))
      (setf ae (+ ae 360d0)))
    (do () ((not (> (- ae as) (- 360d0 angle-tol-deg))))
      (setf ae (- ae 360d0)))
    (when (< (- ae as) angle-tol-deg)
      (setf ae (+ as 360d0)))
    (let ((theta (- ae as)))
      (declare (type double-float theta))
      (cond ((= ctp +quadratic+)
             (setf p 2
                   n (cond ((< theta 91d0) 2)
                           ((< theta 181d0) 4)
                           ((< theta 271d0) 6)
                           (t 8))))
            ((= ctp +quartic+)
             (setf p 4
                   n (cond ((< theta 121d0) 4)
                           ((< theta 241d0) 8)
                           (t 12))))
            ((= ctp +quintic+)
             (setf p 5 n 5))))
    (values ae n p (+ n p 1))))

(defun crv-set-c-pts-and-knots (cur pw u)
  "Set the control points of curve CUR to PW (a cpoint-array) and its knots
to U (a real-vector), sharing the arrays as the C shares the pointers.
Return CUR.
NLib: N_CrvSetCPtsAndKnots."
  (declare (type curve cur) (type cpoint-array pw) (type real-vector u))
  (setf (curve-cpoints cur) pw
        (curve-knots cur) u)
  cur)

(defun crv-is-closed (cur)
  "True when the first and last control points of curve CUR coincide in
Euclidean space to within *mtol*.
NLib: N_CrvIsClosed."
  (declare (type curve cur))
  (let ((pw (curve-cpoints cur))
        (n (curve-n cur))
        (p0 (make-point))
        (pn (make-point)))
    (declare (type cpoint-array pw) (type fixnum n) (dynamic-extent p0 pn))
    (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw 0)
      (set-point p0 x y z))
    (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw n)
      (set-point pn x y z))
    (< (dist-pt-pt p0 pn) (the double-float *mtol*))))

(defun crv-is-degen (cur)
  "True when curve CUR degenerates to a point: every control point lies
within *mtol* of the control points' Euclidean centroid.
NLib: N_CrvIsDegen."
  (declare (type curve cur))
  (let* ((pw (curve-cpoints cur))
         (n (curve-n cur))
         (mtol *mtol*)
         (fac (/ 1d0 (+ (float n 1d0) 1d0)))
         (q (make-point))
         (m (make-point)))
    (declare (type cpoint-array pw) (type fixnum n) (type double-float mtol fac)
             (dynamic-extent q m))
    (fast
      (dotimes (i (1+ n))
        (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw i)
          (set-point q x y z))
        (sum-2-pts m q m))
      (scale-pt fac m m)
      (dotimes (i (1+ n) t)
        (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw i)
          (set-point q x y z))
        (when (> (dist-pt-pt q m) mtol)
          (return nil))))))

(defun crv-copy (curp)
  "A fresh copy of curve CURP: its degree, control points and knots.
NLib: N_CrvCopy."
  (declare (type curve curp))
  (let* ((pw (curve-cpoints curp))
         (up (curve-knots curp))
         (n+1 (array-dimension pw 0))
         (m+1 (length up))
         (qw (make-cpoint-array n+1))
         (uq (make-real-vector m+1)))
    (declare (type cpoint-array pw qw) (type real-vector up uq) (type fixnum n+1 m+1))
    (fast
      (dotimes (i n+1)
        (setf (aref qw i 0) (aref pw i 0)
              (aref qw i 1) (aref pw i 1)
              (aref qw i 2) (aref pw i 2)
              (aref qw i 3) (aref pw i 3)))
      (dotimes (i m+1)
        (setf (aref uq i) (aref up i))))
    (%make-curve (curve-degree curp) uq qw)))

(defun create-crv-from-num-and-denom (num den)
  "A fresh NURBS curve whose homogeneous coordinates (wx wy wz) are the
control points of curve NUM and whose weights are the values of the cfun
DEN, or non-rational when DEN is nil.  NUM and DEN must have the same
degree, control point count and knot count, else :inp-err.
NLib: N_CreateCrvFromNumAndDenom."
  (declare (type curve num) (type (or null cfun) den))
  (let* ((nw (curve-cpoints num))
         (un (curve-knots num))
         (nn (curve-n num))
         (mn (curve-m num))
         (pn (curve-degree num)))
    (declare (type cpoint-array nw) (type real-vector un) (type fixnum nn mn pn))
    (when den
      (let ((nd (1- (length (cfun-vals den))))
            (md (1- (length (cfun-knots den))))
            (pd (cfun-degree den)))
        (declare (type fixnum nd md pd))
        (when (or (/= nn nd) (/= mn md) (/= pn pd))
          (nurbs-error :inp-err 'create-crv-from-num-and-denom))))
    (let ((pw (make-cpoint-array (1+ nn)))
          (uc (make-real-vector (1+ mn)))
          (fu (if den (cfun-vals den) nil)))
      (declare (type cpoint-array pw) (type real-vector uc) (type (or null real-vector) fu))
      (dotimes (i (1+ nn))
        (setf (aref pw i 0) (aref nw i 0)
              (aref pw i 1) (aref nw i 1)
              (aref pw i 2) (aref nw i 2)
              (aref pw i 3) (if fu (aref fu i) +no-w+)))
      (fast
        (dotimes (i (1+ mn))
          (setf (aref uc i) (aref un i))))
      (%make-curve pn uc pw))))
