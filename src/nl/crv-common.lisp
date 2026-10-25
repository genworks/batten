;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvCommon.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------ level 2, group 5

(declaim (inline %crv-common-sqrt))
(defun %crv-common-sqrt (x)
  "C sqrt of X, known non-negative or NaN (a NaN goes through as the C's)."
  (declare (type double-float x))
  (if (>= x 0d0)
      (sqrt (the (double-float 0d0) x))
      (- x x)))

(defun %crv-common-row-weight (p w pw i)
  "N_Weight( P, w, &Pw[i] ): write point P with weight W into row I of
cpoint-array PW."
  (declare (type point p) (type double-float w) (type cpoint-array pw) (type fixnum i))
  (let ((cw (make-cpoint)))
    (declare (dynamic-extent cw))
    (weight p w cw)
    (%c-pt-row cw pw i)))

(defun create-quadratic-arc-tol (c x y r as ae at)
  "A fresh quadratic rational NURBS circle or circular arc of center C, axes
X and Y (orthogonal; normalized here, on copies) and radius R, from angle AS
to angle AE in degrees, made of up to four Bezier arcs.  When the sweep
exceeds 360 - AT degrees the last control point is made the first.  An arc
too flat for its end tangents to meet gets the midpoint of its ends as the
middle control point.  Signals :geo-err for a zero axis, and :inp-err when
the end tangents of an arc cannot be intersected.
The C leaves the intersection flag as it was when N_IsectLineLine fails;
here it starts out true (uninitialized in the C, whose first arc then reads
garbage), and a failure the C carries past (flag false) is signaled at the
end when it came from the last arc, as the C returns the last error.
NLib: N_CreateQuadraticArcTol."
  (declare (type point c x y) (type double-float r as ae at))
  (multiple-value-bind (ae n p m) (crv-get-conic-data +quadratic+ as ae)
    (declare (type double-float ae) (type fixnum n p m) (ignore p))
    (let* ((narcs (truncate n 2))
           (theta (- ae as))
           (dtheta (/ theta (float narcs 1d0)))
           (angle-tol-deg at)
           (pw (make-cpoint-array (1+ n)))
           (u (make-real-vector (1+ m)))
           (xx (make-point)) (yy (make-point))
           (p0 (make-point)) (p1 (make-point)) (p2 (make-point))
           (t0 (make-point)) (t2 (make-point))
           (l0 (%make-lineseg (make-point) (make-point) +unbounded+))
           (l2 (%make-lineseg (make-point) (make-point) +unbounded+))
           (ints +true+)
           (pending nil)
           (angle 0d0) (dangle 0d0) (w1 0d0) (cosa 0d0) (sina 0d0) (alf 0d0) (bet 0d0))
      (declare (type fixnum narcs ints) (type double-float theta dtheta angle-tol-deg
                                                         angle dangle w1 cosa sina alf bet)
               (type cpoint-array pw) (type real-vector u)
               (dynamic-extent xx yy p0 p1 p2 t0 t2 l0 l2))
      (copy-pt x xx)
      (copy-pt y yy)
      (vector-normalize-ref xx)
      (vector-normalize-ref yy)
      (setf angle (/ (* +pi+ dtheta) 360d0)
            w1 (cos angle))
      (setf angle (/ (* +pi+ as) 180d0)
            cosa (cos angle)
            sina (sin angle)
            alf (* r cosa)
            bet (* r sina))
      (translate-sum-2-pts c alf xx bet yy p0)
      (combine-2-pts (- sina) xx cosa yy t0)
      (%crv-common-row-weight p0 1d0 pw 0)
      (setf dangle (/ (* +pi+ dtheta) 180d0))
      (let ((j 0))
        (declare (type fixnum j))
        (do ((i 1 (1+ i))) ((> i narcs))
          (declare (type fixnum i))
          (setf angle (+ angle dangle)
                cosa (cos angle)
                sina (sin angle)
                alf (* r cosa)
                bet (* r sina))
          (translate-sum-2-pts c alf xx bet yy p2)
          (combine-2-pts (- sina) xx cosa yy t2)
          (create-line-start-dir-vector l0 p0 t0 +unbounded+)
          (create-line-start-dir-vector l2 p2 t2 +unbounded+)
          (setf pending nil)
          (handler-case (setf ints (values (isect-line-line l0 l2 p1)))
            (nurbs-error (e) (setf pending e)))
          (cond ((= ints +false+)
                 ;; (the C also takes N_DistPtPt( P0, P2 ) here, unused)
                 (combine-2-pts 0.5d0 p0 0.5d0 p2 p1))
                (pending
                 (nurbs-error :inp-err 'create-quadratic-arc-tol)))
          (%crv-common-row-weight p1 w1 pw (+ j 1))
          (%crv-common-row-weight p2 1d0 pw (+ j 2))
          (when (< i narcs)
            (copy-pt p2 p0)
            (copy-pt t2 t0))
          (incf j 2)))
      (when (> theta (- 360d0 angle-tol-deg))
        (%copy-c-pt-row pw 0 pw n))
      (dotimes (i 3)
        (setf (aref u i) 0d0))
      (case narcs
        (2 (setf (aref u 3) 0.5d0 (aref u 4) 0.5d0))
        (3 (setf (aref u 3) (/ 1d0 3d0) (aref u 4) (/ 1d0 3d0)
                 (aref u 5) (/ 2d0 3d0) (aref u 6) (/ 2d0 3d0)))
        (4 (setf (aref u 3) 0.25d0 (aref u 4) 0.25d0
                 (aref u 5) 0.50d0 (aref u 6) 0.50d0
                 (aref u 7) 0.75d0 (aref u 8) 0.75d0)))
      (dotimes (i 3)
        (setf (aref u (+ n i 1)) 1d0))
      (when pending
        (error pending))
      (%make-curve 2 u pw))))

(defun create-quartic-arc-tol (c x y r as ae at)
  "A fresh quartic rational NURBS circle or circular arc of center C, axes X
and Y (orthogonal; normalized here, on copies) and radius R, from angle AS
to angle AE in degrees: up to three quadratic Bezier arcs, each
reparameterized to a quartic so that the parameter is closer to arc length.
When the sweep exceeds 360 - AT degrees the last control point is made the
first.  Signals :geo-err for a zero axis, :num-err when the
reparameterization cannot be solved, and :inp-err when the end tangents of
an arc are parallel.
NLib: N_CreateQuarticArcTol."
  (declare (type point c x y) (type double-float r as ae at))
  (multiple-value-bind (ae n p m) (crv-get-conic-data +quartic+ as ae)
    (declare (type double-float ae) (type fixnum n p m) (ignore p))
    (let* ((narcs (truncate n 4))
           (theta (- ae as))
           (dtheta (/ theta (float narcs 1d0)))
           (angle-tol-deg at)
           (pw (make-cpoint-array (1+ n)))
           (u (make-real-vector (1+ m)))
           (xx (make-point)) (yy (make-point))
           (p0 (make-point)) (p1 (make-point)) (p2 (make-point))
           (t0 (make-point)) (t2 (make-point))
           (l0 (%make-lineseg (make-point) (make-point) +unbounded+))
           (l2 (%make-lineseg (make-point) (make-point) +unbounded+))
           (qw (make-cpoint)) (aw (make-cpoint)) (bw (make-cpoint)) (cw (make-cpoint))
           (angle 0d0) (dangle 0d0) (alf 0d0) (bet 0d0) (disc 0d0) (wf 0d0) (fw 0d0)
           (cosa 0d0) (sina 0d0) (w1 0d0) (x2 0d0) (a 0d0) (b 0d0) (cc 0d0) (d 0d0)
           (ang 0d0) (cs 0d0) (u1 0d0) (u2 0d0) (v1 0d0) (v2 0d0))
      (declare (type fixnum narcs)
               (type double-float theta dtheta angle-tol-deg angle dangle alf bet disc wf fw
                     cosa sina w1 x2 a b cc d ang cs u1 u2 v1 v2)
               (type cpoint-array pw) (type real-vector u)
               (dynamic-extent xx yy p0 p1 p2 t0 t2 l0 l2 qw aw bw cw))
      (copy-pt x xx)
      (copy-pt y yy)
      (vector-normalize-ref xx)
      (vector-normalize-ref yy)
      (setf dtheta (/ theta (float narcs 1d0))
            angle (/ (* +pi+ dtheta) 180d0)
            x2 (cos angle)
            angle (* 0.5d0 angle)
            w1 (cos angle))
      (setf angle (/ (* +pi+ as) 180d0)
            cosa (cos angle)
            sina (sin angle)
            alf (* r cosa)
            bet (* r sina))
      (translate-sum-2-pts c alf xx bet yy p0)
      (combine-2-pts (- sina) xx cosa yy t0)
      (%crv-common-row-weight p0 1d0 pw 0)
      ;; the parameters u1, u2 of the quadratic Bezier arc at 1/4 and 3/4 of
      ;; the sweep
      (dotimes (i 2)
        (setf ang (if (= i 0)
                      (/ (* (* 0.25d0 +pi+) dtheta) 180d0)
                      (/ (* (* 0.75d0 +pi+) dtheta) 180d0))
              cs (cos ang))
        (setf a (- (+ (- 1d0 (* 2d0 w1)) x2) (* (* 2d0 cs) (- 1d0 w1)))
              b (+ (* 2d0 (- w1 1d0)) (* (* 2d0 cs) (- 1d0 w1)))
              cc (- 1d0 cs))
        (setf disc (- (* b b) (* (* 4d0 a) cc)))
        (when (< disc 0d0)
          (nurbs-error :num-err 'create-quartic-arc-tol))
        ;; a NaN sweep makes DISC NaN: the C goes on to its :num-err below
        (setf disc (%crv-common-sqrt disc)
              v1 (/ (+ (- b) disc) (* 2d0 a))
              v2 (/ (- (- b) disc) (* 2d0 a)))
        (cond ((and (>= v1 0d0) (<= v1 1d0)) (setf a v1))
              ((and (>= v2 0d0) (<= v2 1d0)) (setf a v2))
              (t (nurbs-error :num-err 'create-quartic-arc-tol)))
        (if (= i 0) (setf u1 a) (setf u2 a)))
      ;; the rational quadratic Bezier reparameterization
      (setf a (- u1 u2))
      (setf wf (/ (- (+ 8d0 (* 10d0 a))) (* 6d0 a))
            fw (/ (- (* u1 (+ 10d0 (* 6d0 wf))) 1d0) 6d0))
      (setf d (/ (* (* 2d0 fw) fw) 3d0)
            cc (/ (- (+ (- 1d0 (* 6d0 fw)) (* (* 4d0 wf) fw)) (* (* 4d0 fw) fw)) 3d0)
            b (/ (+ (+ (- (+ (* -6d0 wf) (* (* 2d0 wf) wf)) (* (* 4d0 wf) fw)) (* 6d0 fw))
                    (* (* 2d0 fw) fw))
                 3d0)
            a (- wf fw))
      (setf dangle (/ (* +pi+ dtheta) 180d0))
      (let ((j 0))
        (declare (type fixnum j))
        (do ((i 1 (1+ i))) ((> i narcs))
          (declare (type fixnum i))
          (setf angle (+ angle dangle)
                cosa (cos angle)
                sina (sin angle)
                alf (* r cosa)
                bet (* r sina))
          (translate-sum-2-pts c alf xx bet yy p2)
          (combine-2-pts (- sina) xx cosa yy t2)
          (create-line-start-dir-vector l0 p0 t0 +unbounded+)
          (create-line-start-dir-vector l2 p2 t2 +unbounded+)
          (when (= (the fixnum (values (isect-line-line l0 l2 p1))) +false+)
            (nurbs-error :inp-err 'create-quartic-arc-tol))
          (weight p1 w1 qw)
          (%crv-common-row-weight p2 1d0 pw (+ j 4))
          ;; the quadratic segment is Pw[j], Qw, Pw[j+4]
          (%row-c-pt pw j aw)
          (%row-c-pt pw (+ j 4) bw)
          (combine-2-c-pts a aw fw qw cw)
          (%c-pt-row cw pw (+ j 1))
          (combine-4-c-pts b aw cc qw d bw 2d0 cw cw)
          (%c-pt-row cw pw (+ j 2))
          (combine-2-c-pts a qw fw bw cw)
          (%c-pt-row cw pw (+ j 3))
          (when (< i narcs)
            (copy-pt p2 p0)
            (copy-pt t2 t0))
          (incf j 4)))
      (when (> theta (- 360d0 angle-tol-deg))
        (%copy-c-pt-row pw 0 pw n))
      (dotimes (i 5)
        (setf (aref u i) 0d0))
      (case narcs
        (2 (do ((i 5 (1+ i))) ((> i 8))
             (declare (type fixnum i))
             (setf (aref u i) 0.5d0)))
        (3 (do ((i 5 (1+ i))) ((> i 8))
             (declare (type fixnum i))
             (setf (aref u i) (/ 1d0 3d0)))
           (do ((i 9 (1+ i))) ((> i 12))
             (declare (type fixnum i))
             (setf (aref u i) (/ 2d0 3d0)))))
      (dotimes (i 5)
        (setf (aref u (+ n i 1)) 1d0))
      (%make-curve 4 u pw))))

;;; ------------------------------------------------------------ level 3, group 5

(defun create-quintic-arc-tol (c x y r as ae at)
  "A fresh quintic rational NURBS circle or circular arc of center C, axes X
and Y (normalized here, on copies) and radius R, from angle AS to angle AE in
degrees, with no interior knots: a quintic full circle starting at the point
opposite the arc's middle is built, and when the sweep is at most 360 - AT
degrees the arc is extracted from it between the parameters of its end
points.  Signals :geo-err for a zero axis, and what crv-closest-pt and
crv-extract-crv-seg signal.
NLib: N_CreateQuinticArcTol."
  (declare (type point c x y) (type double-float r as ae at))
  (multiple-value-bind (ae n p m) (crv-get-conic-data +quintic+ as ae)
    (declare (type double-float ae) (type fixnum n p m))
    (let* ((pw (make-cpoint-array (1+ n)))
           (u (make-real-vector (1+ m)))
           (angle-tol-deg at)
           (xx (make-point)) (yy (make-point))
           (ps (make-point)) (pe (make-point)) (pm (make-point)) (a (make-point))
           (tv (make-point)) (sv (make-point))
           (angle 0d0) (dangle 0d0) (alf 0d0) (bet 0d0) (cosa 0d0) (sina 0d0)
           (am 0d0) (us 0d0) (ue 0d0) (u0 0d0))
      (declare (type cpoint-array pw) (type real-vector u)
               (type double-float angle-tol-deg angle dangle alf bet cosa sina am us ue u0)
               (dynamic-extent xx yy ps pe pm a tv sv))
      ;; start, end and mid points
      (setf am (* 0.5d0 (+ as ae)))
      (setf am (+ am 180d0))
      (copy-pt x xx)
      (copy-pt y yy)
      (vector-normalize-ref xx)
      (vector-normalize-ref yy)
      (setf angle (/ (* +pi+ as) 180d0)
            alf (* r (cos angle))
            bet (* r (sin angle)))
      (translate-sum-2-pts c alf xx bet yy ps)
      (setf angle (/ (* +pi+ ae) 180d0)
            alf (* r (cos angle))
            bet (* r (sin angle)))
      (translate-sum-2-pts c alf xx bet yy pe)
      (setf angle (/ (* +pi+ am) 180d0)
            cosa (cos angle)
            sina (sin angle)
            alf (* r cosa)
            bet (* r sina))
      (translate-sum-2-pts c alf xx bet yy pm)
      (combine-2-pts (- sina) xx cosa yy tv)
      (combine-2-pts (- cosa) xx (- sina) yy sv)
      (vector-normalize-ref tv)
      (vector-normalize-ref sv)
      ;; quintic full circle
      (%crv-common-row-weight pm 5d0 pw 0)
      (%crv-common-row-weight pm 5d0 pw 5)
      (combine-2-pts 1d0 pm (* 4d0 r) tv a)
      (%crv-common-row-weight a 1d0 pw 1)
      (combine-2-pts 1d0 pm (* -4d0 r) tv a)
      (%crv-common-row-weight a 1d0 pw 4)
      (translate-sum-2-pts c (* 3d0 r) sv (* 2d0 r) tv a)
      (%crv-common-row-weight a 1d0 pw 2)
      (translate-sum-2-pts c (* 3d0 r) sv (* -2d0 r) tv a)
      (%crv-common-row-weight a 1d0 pw 3)
      (dotimes (i 6)
        (setf (aref u i) 0d0
              (aref u (+ i 6)) 1d0))
      (let ((cur (%make-curve p u pw)))
        ;; full circle
        (when (> (- ae as) (- 360d0 angle-tol-deg))
          (%copy-c-pt-row pw 0 pw 5)
          (return-from create-quintic-arc-tol cur))
        ;; extract the arc
        (setf dangle (- as am))
        (setf dangle (+ 360d0 dangle))
        (setf u0 (/ dangle 360d0))
        (setf us (crv-closest-pt cur ps u0 *mtol* *mtol* a))
        (setf dangle (- ae am))
        (when (< dangle 0d0)
          (setf dangle (+ 360d0 dangle)))
        (setf u0 (/ dangle 360d0))
        (setf ue (crv-closest-pt cur pe u0 *mtol* *mtol* a))
        (crv-extract-crv-seg cur us ue)))))

(defun create-conic-arc (p0 t0 p2 t2 p)
  "A fresh rational quadratic NURBS conic arc from P0 (tangent T0) to P2
(tangent T2) through point P: one Bezier arc, or the arc split at its
shoulder point into two or four, by its middle weight and the angle between
the tangents.  When the middle point cannot be found
(bez-conic-arc-from-pts-and-tangents fails), a straight 'conic' with the
middle point halfway between P0 and P2 and weight 1 is made instead.  Signals
:wei-err when the middle weight is -1 or less.  A NaN classification (zero
tangents) leaves no arc chosen: the C then makes a curve with one control
point, P2, and that is returned too.
NLib: N_CreateConicArc."
  (declare (type point p0 t0 p2 t2 p))
  (let ((p1 (make-point))
        (wp 0d0)
        (one-arc nil) (two-arcs nil) (four-arcs nil)
        (n 0))
    (declare (type double-float wp) (type fixnum n))
    ;; one arc first
    (handler-case (setf wp (bez-conic-arc-from-pts-and-tangents p0 t0 p2 t2 p p1))
      (nurbs-error ()
        ;; a conic 'line'
        (set-point p1
                   (/ (+ (aref p0 0) (aref p2 0)) 2d0)
                   (/ (+ (aref p0 1) (aref p2 1)) 2d0)
                   (/ (+ (aref p0 2) (aref p2 2)) 2d0))
        (setf wp 1d0)))
    (when (<= wp -1d0)
      (nurbs-error :wei-err 'create-conic-arc))
    ;; classify by the number of arcs
    (if (>= wp 1d0)
        (setf one-arc t)
        (let* ((dot (vector-dot t0 t2))
               (m1 (vector-magnitude t0))
               (m2 (vector-magnitude t2))
               (cosa (/ (- dot) (* m1 m2))))
          (declare (type double-float dot m1 m2 cosa))
          (cond ((and (>= wp 0d0) (<= cosa 0.5d0)) (setf one-arc t))
                ((and (= wp 0d0) (> cosa (- 1d0 (the double-float *ptol*)))) (setf two-arcs t))
                ((and (>= wp 0d0) (> cosa 0.5d0)) (setf two-arcs t))
                ((and (< wp 0d0) (>= cosa 0d0)) (setf two-arcs t))
                ((and (< wp 0d0) (< cosa 0d0)) (setf four-arcs t)))))
    (cond (one-arc (setf n 2))
          (two-arcs (setf n 4))
          (four-arcs (setf n 8)))
    (let* ((m (+ n 3))
           (pw (make-cpoint-array (1+ n)))
           (u (make-real-vector (1+ m))))
      (declare (type fixnum m) (type cpoint-array pw) (type real-vector u))
      ;; end points and end knots
      (%crv-common-row-weight p0 1d0 pw 0)
      (%crv-common-row-weight p2 1d0 pw n)
      (dotimes (i 3)
        (setf (aref u i) 0d0
              (aref u (+ n i 1)) 1d0))
      (cond
        (one-arc
         (%crv-common-row-weight p1 wp pw 1))
        (t
         (let ((pp (make-point-array 3)) (q (make-point-array 3)) (r (make-point-array 3))
               (pt (make-point)))
           (declare (dynamic-extent pp q r pt))
           (%pt-row p0 pp 0)
           (%pt-row p1 pp 1)
           (%pt-row p2 pp 2)
           (multiple-value-bind (wq wr) (bez-conic-split pp wp q r)
             (declare (type double-float wq wr))
             (flet ((put (arr i w k)
                      (%row-pt arr i pt)
                      (%crv-common-row-weight pt w pw k)))
               (cond
                 (two-arcs
                  (put q 1 wq 1)
                  (put q 2 1d0 2)
                  (put r 1 wr 3)
                  (setf (aref u 3) 0.5d0
                        (aref u 4) 0.5d0))
                 (four-arcs
                  (let ((q1 (make-point-array 3)) (q2 (make-point-array 3))
                        (r1 (make-point-array 3)) (r2 (make-point-array 3)))
                    (declare (dynamic-extent q1 q2 r1 r2))
                    (multiple-value-bind (wq1 wq2) (bez-conic-split q wq q1 q2)
                      (declare (type double-float wq1 wq2))
                      (multiple-value-bind (wr1 wr2) (bez-conic-split r wr r1 r2)
                        (declare (type double-float wr1 wr2))
                        (put q1 1 wq1 1)
                        (put q1 2 1d0 2)
                        (put q2 1 wq2 3)
                        (put q2 2 1d0 4)
                        (put r1 1 wr1 5)
                        (put r1 2 1d0 6)
                        (put r2 1 wr2 7)
                        (setf (aref u 3) 0.25d0 (aref u 4) 0.25d0
                              (aref u 5) 0.50d0 (aref u 6) 0.50d0
                              (aref u 7) 0.75d0 (aref u 8) 0.75d0)))))))))))
      (%make-curve 2 u pw))))

;;; ------------------------------------------------------------- L4-g0

(defun create-circ-arc-tol (c x y r as ae at ctp)
  "A fresh NURBS circle or circular arc of center C, axes X and Y, radius R,
from angle AS to angle AE in degrees with angle tolerance AT, of the kind CTP
names: +quadratic+ (create-quadratic-arc-tol), +quartic+
(create-quartic-arc-tol) or +quintic+ (create-quintic-arc-tol).  Signals
:cal-err for any other CTP.
NLib: N_CreateCircArcTol."
  (declare (type point c x y) (type double-float r as ae at) (type fixnum ctp))
  (cond ((= ctp +quadratic+) (create-quadratic-arc-tol c x y r as ae at))
        ((= ctp +quartic+) (create-quartic-arc-tol c x y r as ae at))
        ((= ctp +quintic+) (create-quintic-arc-tol c x y r as ae at))
        (t (nurbs-error :cal-err 'create-circ-arc-tol))))

;;; ------------------------------------------------------------- L5-g0

(defun create-circ-arc (c x y r as ae ctp)
  "A fresh NURBS circle or circular arc of center C, axes X and Y, radius R,
from angle AS to angle AE in degrees, of the kind CTP names: create-circ-arc-tol
with the angle tolerance 0.01 degrees.
NLib: N_CreateCircArc."
  (declare (type point c x y) (type double-float r as ae) (type fixnum ctp))
  (create-circ-arc-tol c x y r as ae 0.01d0 ctp))

(defun create-elliptical-arc-ang-type-tol (c x y r1 r2 as ae at atp ctp)
  "A fresh NURBS ellipse or elliptical arc of center C, axes X and Y
(normalized here, on copies), major radius R1 along X and minor radius R2
along Y, from angle AS to angle AE in degrees (ATP 0: angles of the point on
the circle of radius R1, projected to the ellipse; otherwise angles of the
point on the ellipse), with angle tolerance AT and conic kind CTP: the circle
of radius R1 is made by create-circ-arc-tol and its control points are scaled
toward the X axis by R2/R1.  Signals :geo-err for a zero axis, :inp-err when a
control point cannot be projected onto the X axis, and what
create-circ-arc-tol signals.
NLib: N_CreateEllipticalArcAngTypeTol."
  (declare (type point c x y) (type double-float r1 r2 as ae at) (type fixnum atp ctp))
  (let ((xx (make-point)) (yy (make-point))
        (lx (%make-lineseg (make-point) (make-point) +unbounded+))
        (p (make-point)) (q (make-point))
        (apnt (make-real-vector 2))
        (fac 0d0) (omf 0d0) (angle 0d0) (alf 0d0) (bet 0d0) (e 0d0)
        (tmp-angle 0d0) (ap 0d0))
    (declare (dynamic-extent xx yy lx p q apnt)
             (type double-float fac omf angle alf bet e tmp-angle ap))
    (copy-pt x xx)
    (copy-pt y yy)
    (vector-normalize-ref xx)
    (vector-normalize-ref yy)
    (create-line-start-dir-vector lx c xx +unbounded+)
    ;; (the C also sets fac = r1/r2, omf = 1 - fac here, unused)
    ;; start and end angles
    (do () ((not (<= ae as)))
      (setf ae (+ 360d0 ae)))
    (do () ((not (>= (- ae as) 360.1d0)))
      (setf ae (- ae 360d0)))
    (when (> (- ae as) 360d0)
      (setf ae (+ as 360d0)))
    (dotimes (i 2)
      (setf angle (/ (* +pi+ (if (= i 0) as ae)) 180d0))
      (cond ((= atp 0)
             ;; point on circle at angle, projected to the ellipse
             (setf alf (* r1 (cos angle))
                   bet (* r2 (sin angle))))
            (t
             ;; point on the ellipse at angle
             (setf e (tan angle))
             (setf alf (/ (* r1 r2) (%crv-common-sqrt (+ (* r2 r2) (* (* (* r1 r1) e) e)))))
             (when (< (cos angle) 0d0)
               (setf alf (- alf)))
             (setf bet (- (* r1 r1) (* alf alf)))
             (if (< bet 1d-12)
                 (setf bet 0d0)
                 (setf bet (* (/ r2 r1) (%crv-common-sqrt bet))))
             (when (< (sin angle) 0d0)
               (setf bet (- bet)))))
      ;; polar angle of the vector
      (setf tmp-angle (atan (/ bet r2) (/ alf r1)))
      (setf ap (/ (* 180d0 tmp-angle) +pi+))
      ;; angle snapping up to NL_PTOL degrees
      (when (<= (abs ap) (the double-float *ptol*))
        (setf ap 0d0))
      (when (< ap 0d0)
        (setf ap (+ 360d0 ap)))
      (setf (aref apnt i) ap))
    ;; the circle or circular arc
    (let* ((cur (create-circ-arc-tol c xx yy r1 (aref apnt 0) (aref apnt 1) at ctp))
           (pw (curve-cpoints cur))
           (cw (make-cpoint))
           (w 0d0))
      (declare (dynamic-extent cw) (type cpoint-array pw) (type double-float w))
      ;; the affine map to the ellipse's control points
      (setf fac (/ r2 r1)
            omf (- 1d0 fac))
      (dotimes (i (array-dimension pw 0))
        (%row-c-pt pw i cw)
        (c-pt-to-pt-euclid cw p)
        (setf w (c-pt-get-w cw))
        (when (= (the fixnum (values (project-pt-line lx p q))) +false+)
          (nurbs-error :inp-err 'create-elliptical-arc-ang-type-tol))
        (combine-2-pts fac p omf q p)
        (weight p w cw)
        (%c-pt-row cw pw i))
      cur)))

;;; ------------------------------------------------------------- L6-g1

(defun calc-circ-weights (as ae ctp w)
  "Write into W (caller-allocated, 9 entries suffice) the weights of the unit
circle or circular arc from angle AS to angle AE in degrees of conic kind CTP
(+quadratic+ or +quintic+; whatever create-circ-arc accepts), and return the
highest index written.  Signals what create-circ-arc signals.
NLib: N_CalcCircWeights."
  (declare (type double-float as ae) (type fixnum ctp) (type real-vector w))
  (let* ((cura (create-circ-arc (make-point 0d0 0d0 0d0) (make-point 1d0 0d0 0d0)
                                (make-point 0d0 1d0 0d0) 1d0 as ae ctp))
         (aw (curve-cpoints cura))
         (n (1- (array-dimension aw 0))))
    (declare (type cpoint-array aw) (type fixnum n))
    (do ((i 0 (1+ i))) ((> i n))
      (declare (type fixnum i))
      (setf (aref w i) (aref aw i 3)))
    n))

(defun create-elliptical-arc (c x y r1 r2 as ae ctp)
  "A fresh NURBS ellipse or elliptical arc: create-elliptical-arc-ang-type-tol
with angle tolerance 0.01 degrees and ATP 0 (AS and AE are angles of the point
on the circle of radius R1).
NLib: N_CreateEllipticalArc."
  (declare (type point c x y) (type double-float r1 r2 as ae) (type fixnum ctp))
  (create-elliptical-arc-ang-type-tol c x y r1 r2 as ae 0.01d0 0 ctp))
