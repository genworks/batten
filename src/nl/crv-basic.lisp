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

;;; ------------------------------------------------------------ level 1, group 1

(defun crv-from-c-pts-and-knots (cur pw n p u m)
  "Make curve CUR (modified in place) the curve of control points PW[0..N],
degree P and knots U[0..M], sharing the arrays where they are exactly that
long, as the C shares the pointers.  Return CUR.
NLib: N_CrvFromCPtsAndKnots."
  (declare (type curve cur) (type cpoint-array pw) (type fixnum n p m) (type real-vector u))
  (let ((pol (c-polygon-from-c-pts pw n))
        (knt (knot-vector-from-real-array u m)))
    (declare (type cpoint-array pol) (type real-vector knt))
    (crv-from-c-polygon-and-knot-vector cur pol p knt)))

(defun crv-compress (cur)
  "Give curve CUR (modified in place) fresh control point and knot arrays
holding copies of its control points and knots.  Return CUR.
NLib: N_CrvCompress."
  (declare (type curve cur))
  (let* ((pw (curve-cpoints cur))
         (u (curve-knots cur))
         (n (curve-n cur))
         (m (curve-m cur))
         (qw (make-cpoint-array (1+ n)))
         (v (make-real-vector (1+ m))))
    (declare (type cpoint-array pw qw) (type real-vector u v) (type fixnum n m))
    (fast
      (dotimes (i (1+ n))
        (%copy-c-pt-row pw i qw i))
      (dotimes (i (1+ m))
        (setf (aref v i) (aref u i))))
    (crv-set-c-pts-and-knots cur qw v)))

(defun crv-is-line (cur tol)
  "True when curve CUR is a straight line to within TOL: every interior
control point projects inside the segment between the end points (extended
by TOL at each end) and lies within TOL of it.  When the end points are
within TOL of each other, every control point must also lie within TOL of
the first one.
NLib: N_CrvIsLine."
  (declare (type curve cur) (type double-float tol))
  (let ((pw (curve-cpoints cur))
        (n (curve-n cur))
        (q (make-point)) (a (make-point)) (b (make-point))
        (line (make-point)) (aprime (make-point)) (bprime (make-point)))
    (declare (type cpoint-array pw) (type fixnum n)
             (dynamic-extent q a b line aprime bprime))
    (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw 0)
      (set-point a x y z))
    (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw n)
      (set-point b x y z))
    (when (<= (dist-pt-pt a b) tol)
      (do ((i 1 (1+ i))) ((>= i n))
        (declare (type fixnum i))
        (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw i)
          (set-point q x y z))
        (when (> (dist-pt-pt a q) tol)
          (return-from crv-is-line nil))))
    (vector-dir a b line)
    (let ((length (vector-magnitude line)))
      (declare (type double-float length))
      ;; the C ignores N_VectorNormalizeRef's error, leaving LINE unscaled
      (handler-case (vector-normalize-ref line)
        (nurbs-error () nil))
      (vector-pt-along-vector a (+ length tol) line bprime)
      (vector-reverse-in-place line)
      (vector-pt-along-vector b (+ length tol) line aprime))
    (do ((i 1 (1+ i))) ((>= i n) t)
      (declare (type fixnum i))
      (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw i)
        (set-point q x y z))
      (let ((d (handler-case (dist-perp-pt-line-seg q aprime bprime)
                 (nurbs-error () nil))))
        (when (null d)
          (return nil))
        (locally (declare (type double-float d))
          (when (or (< d 0d0) (> d tol))
            (return nil)))))))

(defun crvs-are-combatible (cur k)
  "True when the curves CUR[0..K] (a simple-vector of curves) are compatible:
all rational or all non-rational, all 2D or all 3D, with the same highest
control point and knot indexes, and knots equal to within *ptol*.
NLib: N_CrvsAreCombatible."
  (declare (type simple-vector cur) (type fixnum k))
  (let ((allra t) (allnra t) (all2d t) (all3d t))
    (dotimes (i (1+ k))
      (unless (is-crv-rat (svref cur i))
        (setf allra nil)
        (return)))
    (dotimes (i (1+ k))
      (when (is-crv-rat (svref cur i))
        (setf allnra nil)
        (return)))
    (when (and (not allra) (not allnra))
      (return-from crvs-are-combatible nil))
    (dotimes (i (1+ k))
      (when (crv-is-3d (svref cur i))
        (setf all2d nil)
        (return)))
    (dotimes (i (1+ k))
      (unless (crv-is-3d (svref cur i))
        (setf all3d nil)
        (return)))
    (when (and (not all2d) (not all3d))
      (return-from crvs-are-combatible nil))
    (let ((n0 (curve-n (svref cur 0)))
          (m0 (curve-m (svref cur 0))))
      (declare (type fixnum n0 m0))
      (do ((i 1 (1+ i))) ((> i k))
        (declare (type fixnum i))
        (unless (and (= (curve-n (svref cur i)) n0) (= (curve-m (svref cur i)) m0))
          (return-from crvs-are-combatible nil))))
    (let ((u0 (curve-knots (svref cur 0)))
          (ptol *ptol*))
      (declare (type real-vector u0) (type double-float ptol))
      (do ((i 1 (1+ i))) ((> i k) t)
        (declare (type fixnum i))
        (let* ((u (curve-knots (svref cur i)))
               (m (curve-m (svref cur i))))
          (declare (type real-vector u) (type fixnum m))
          (do ((j 0 (1+ j))) ((> j m))
            (declare (type fixnum j))
            (when (> (abs (- (aref u j) (aref u0 j))) ptol)
              (return-from crvs-are-combatible nil))))))))

(defun crv-get-num-and-denom (cur)
  "The numerator and denominator of curve CUR, as two values: a fresh
non-rational curve whose control points are CUR's homogeneous (wx wy wz),
and, when CUR is rational, a fresh cfun whose values are CUR's weights (nil
when CUR is non-rational).
NLib: N_CrvGetNumAndDenom."
  (declare (type curve cur))
  (let* ((pw (curve-cpoints cur))
         (uc (curve-knots cur))
         (n (curve-n cur))
         (m (curve-m cur))
         (p (curve-degree cur))
         (rat (is-crv-rat cur))
         (nw (make-cpoint-array (1+ n)))
         (un (make-real-vector (1+ m)))
         (fu (if rat (make-real-vector (1+ n)) nil))
         (ud (if rat (make-real-vector (1+ m)) nil))
         (pt (make-cpoint)))
    (declare (type cpoint-array pw nw) (type real-vector uc un) (type fixnum n m p)
             (type (or null real-vector) fu ud) (dynamic-extent pt))
    (dotimes (i (1+ n))
      (%row-c-pt pw i pt)
      (multiple-value-bind (xw yw zw w) (c-pt-to-wx-wy-wz pt)
        (c-pt-from-wx-wy-wz xw yw zw +no-w+ pt)
        (%c-pt-row pt nw i)
        (when rat
          (setf (aref fu i) w))))
    (dotimes (i (1+ m))
      (setf (aref un i) (aref uc i))
      (when rat
        (setf (aref ud i) (aref uc i))))
    (values (%make-curve p un nw)
            (if rat (%make-cfun p ud fu) nil))))

(defun crv-get-e-polygon (cur k l)
  "The Euclidean polygon of control points PW[K..L] of curve CUR, a fresh
point-array of L-K+1 vertices.  Signals :ind-err unless 0 <= K <= L <= n.
NLib: N_CrvGetEPolygon."
  (declare (type curve cur) (type fixnum k l))
  (let ((pw (curve-cpoints cur))
        (n (curve-n cur)))
    (declare (type cpoint-array pw) (type fixnum n))
    (when (or (> k l) (< k 0) (> l n))
      (nurbs-error :ind-err 'crv-get-e-polygon))
    (let ((p (make-point-array (1+ (- l k)))))
      (declare (type point-array p))
      (do ((i k (1+ i))) ((> i l))
        (declare (type fixnum i))
        (multiple-value-bind (x y z) (%cpoint-row-to-euclid pw i)
          (setf (aref p (- i k) 0) x
                (aref p (- i k) 1) y
                (aref p (- i k) 2) z)))
      (e-polygon-from-pts (- l k) p))))

(defun crv-3d-to-2d (cur)
  "Make curve CUR 2D, in place, by setting every control point's wz to
+no-z+.  Return CUR.
NLib: N_Crv3dTo2d."
  (declare (type curve cur))
  (let ((pw (curve-cpoints cur)))
    (declare (type cpoint-array pw))
    (dotimes (i (array-dimension pw 0))
      (setf (aref pw i 2) +no-z+)))
  cur)

(defun crv-2d-to-3d (cur)
  "Make curve CUR 3D, in place, by setting every control point's wz to 0.
Return CUR.
NLib: N_Crv2dTo3d."
  (declare (type curve cur))
  (let ((pw (curve-cpoints cur)))
    (declare (type cpoint-array pw))
    (dotimes (i (array-dimension pw 0))
      (setf (aref pw i 2) 0d0)))
  cur)

(defun crv-rat-to-non-rat (cur)
  "Make curve CUR non-rational, in place, by setting every control point's
weight to +no-w+.  Return CUR.
NLib: N_CrvRatToNonRat."
  (declare (type curve cur))
  (let ((pw (curve-cpoints cur)))
    (declare (type cpoint-array pw))
    (dotimes (i (array-dimension pw 0))
      (setf (aref pw i 3) +no-w+)))
  cur)

(defun crv-non-rat-to-rat (cur)
  "Make curve CUR rational, in place, by setting every control point's weight
to 1.  Return CUR.
NLib: N_CrvNonRatToRat."
  (declare (type curve cur))
  (let ((pw (curve-cpoints cur)))
    (declare (type cpoint-array pw))
    (dotimes (i (array-dimension pw 0))
      (setf (aref pw i 3) 1d0)))
  cur)

(defun crv-reparam-to-interval (cur i)
  "Map the knots of curve CUR, in place, linearly onto interval I (a
real-vector (a b)): the first p+1 knots become a, the last p+1 become b.
Nothing changes when the domain already is [a, b].  Return CUR.
NLib: N_CrvReparamToInterval."
  (declare (type curve cur) (type real-vector i))
  (let* ((p (curve-degree cur))
         (u (curve-knots cur))
         (m (curve-m cur)))
    (declare (type fixnum p m) (type real-vector u))
    (multiple-value-bind (a b) (interval-get-data i)
      (declare (type double-float a b))
      (when (or (/= a (aref u 0)) (/= b (aref u m)))
        (let* ((u0 (aref u 0))
               (fac (/ (- b a) (- (aref u m) (aref u 0)))))
          (declare (type double-float u0 fac))
          (fast
            (do ((j 0 (1+ j))) ((> j p))
              (declare (type fixnum j))
              (setf (aref u j) a))
            (do ((j (1+ p) (1+ j))) ((> j (- m p 1)))
              (declare (type fixnum j))
              (setf (aref u j) (+ (* fac (- (aref u j) u0)) a)))
            (do ((j (- m p) (1+ j))) ((> j m))
              (declare (type fixnum j))
              (setf (aref u j) b)))))))
  cur)
