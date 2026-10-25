;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfCommon.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun create-srf-corner-pts (p00 p10 p01 p11)
  "A fresh bilinear, non-rational surface through corner points P00, P10, P01
and P11 (Pij is the corner at u = i, v = j), on [0,1] x [0,1].
NLib: N_CreateSrfCornerPts."
  (declare (type point p00 p10 p01 p11))
  (let ((pw (make-cpoint-net 2 2))
        (u (make-real-vector 4))
        (v (make-real-vector 4))
        (cp (make-cpoint)))
    (declare (type cpoint-net pw) (type real-vector u v) (type cpoint cp))
    (flet ((put (p i j)
             (pt-to-c-pt p cp)
             (dotimes (k 4) (setf (aref pw i j k) (aref cp k)))))
      (put p00 0 0)
      (put p10 1 0)
      (put p01 0 1)
      (put p11 1 1))
    (setf (aref u 0) 0d0 (aref u 1) 0d0 (aref u 2) 1d0 (aref u 3) 1d0
          (aref v 0) 0d0 (aref v 1) 0d0 (aref v 2) 1d0 (aref v 3) 1d0)
    (%make-surface 1 1 u v pw)))

;;; ------------------------------------------------------------------- L2-g3

(defun create-srf-extrude-crv (cur w d dir)
  "A fresh generalized cylinder: curve CUR extruded the distance D along the
direction W (normalized; W itself is left alone, as the C takes it by value).
DIR +udir+ makes the surface linear in u (control points Qw[0][j] the curve's,
Qw[1][j] moved), +vdir+ linear in v (Qw[i][0] and Qw[i][1]); the linear
direction has knots (0 0 1 1).  A 2D CUR is made 3D in place first.  Signals
:cal-err for another DIR and :geo-err for a zero W (CUR already made 3D).
NLib: N_CreateSrfExtrudeCrv."
  (declare (type curve cur) (type point w) (type double-float d) (type fixnum dir))
  (let* ((pw (curve-cpoints cur))
         (n (1- (array-dimension pw 0)))
         (p (curve-degree cur))
         (u (curve-knots cur))
         (m (1- (length u)))
         (ni 0) (mi 0) (pi* 0) (qi 0) (ri 0) (si 0))
    (declare (type cpoint-array pw) (type real-vector u)
             (type fixnum n p m ni mi pi* qi ri si))
    ;; highest indexes
    (cond ((= dir +udir+)
           (setf ni 1 mi n pi* 1 qi p ri 3 si m))
          ((= dir +vdir+)
           (setf ni n mi 1 pi* p qi 1 ri m si 3))
          (t (nurbs-error :cal-err 'create-srf-extrude-crv)))
    (let ((qw (make-cpoint-net (1+ ni) (1+ mi)))
          (uq (make-real-vector (1+ ri)))
          (vq (make-real-vector (1+ si)))
          (wv (make-point (aref w 0) (aref w 1) (aref w 2)))
          (cp (make-cpoint))
          (a (make-point)))
      (declare (type cpoint-net qw) (type real-vector uq vq) (type point wv a) (type cpoint cp))
      ;; prepare for extrusion
      (unless (crv-is-3d cur)
        (crv-2d-to-3d cur))
      (vector-normalize-ref wv)
      (flet ((extrude (k i0 j0 i1 j1)
               ;; Qw[i0][j0] = Pw[k]; Qw[i1][j1] = Pw[k] moved along W
               (%row-c-pt pw k cp)
               (%c-pt-net cp qw i0 j0)
               (let ((wt (c-pt-to-pt-and-w cp a)))
                 (declare (type double-float wt))
                 (vector-pt-along-vector a d wv a)
                 (if (is-crv-rat cur)
                     (weight a wt cp)
                     (pt-to-c-pt a cp))
                 (%c-pt-net cp qw i1 j1))))
        (cond
          ((= dir +udir+)
           (dotimes (j (1+ n))
             (extrude j 0 j 1 j))
           (dotimes (i 2)
             (setf (aref uq i) 0d0
                   (aref uq (+ 2 i)) 1d0))
           (dotimes (j (1+ m))
             (setf (aref vq j) (aref u j))))
          (t
           (dotimes (i (1+ n))
             (extrude i i 0 i 1))
           (dotimes (i (1+ m))
             (setf (aref uq i) (aref u i)))
           (dotimes (j 2)
             (setf (aref vq j) 0d0
                   (aref vq (+ 2 j)) 1d0)))))
      (%make-surface pi* qi uq vq qw))))

;;; ------------------------------------------------------------- L4-g1

(defun create-ruled-srf-from-boundary-crvs (cur-p cur-q dir)
  "A fresh ruled surface between curves CUR-P (at u = 0, or v = 0) and CUR-Q
(at u = 1, or v = 1): S = (1-u) Np/Dp + u Nq/Dq.  DIR +udir+ rules in u (the
surface is linear in u, knots (0 0 1 1)), +vdir+ in v.  When either curve is
rational both are brought over the common denominator Dp*Dq (a non-rational
one is made rational in place for the while, and made non-rational again);
otherwise copies of the two are degree elevated and refined to compatible
knots, the shorter domain rescaled onto the longer.  Signals :inp-err for
another DIR, and what the curve routines signal.
NLib: N_CreateRuledSrfFromBoundaryCrvs."
  (declare (type curve cur-p cur-q) (type fixnum dir))
  (let ((rat-p (is-crv-rat cur-p))
        (rat-q (is-crv-rat cur-q))
        (cur-pw nil) (cur-qw nil)
        (n 0) (p 0) (m 0) (q 0)
        (up (make-real-vector 0)))
    (declare (type fixnum n p m q) (type real-vector up))
    (flet ((locals-p ()
             (setf n (curve-n cur-pw) p (curve-degree cur-pw)
                   m (curve-m cur-pw) up (curve-knots cur-pw)))
           (locals-q ()
             (setf n (curve-n cur-qw) q (curve-degree cur-qw) m (curve-m cur-qw))))
      ;; make the boundary curves compatible
      (cond
        ((or rat-p rat-q)
         ;; make a non-rational curve rational
         (unless rat-p (crv-non-rat-to-rat cur-p))
         (unless rat-q (crv-non-rat-to-rat cur-q))
         ;; numerators and denominators
         (multiple-value-bind (num-p den-p) (crv-get-num-and-denom cur-p)
           (multiple-value-bind (num-q den-q) (crv-get-num-and-denom cur-q)
             ;; product of the denominators, and the numerator terms
             (let* ((den (crv-func-multiply-crv-func den-p den-q))
                    (num-pq (crv-func-multiply-crv den-q num-p))
                    (num-qp (crv-func-multiply-crv den-p num-q)))
               (setf cur-pw (create-crv-from-num-and-denom num-pq den))
               (locals-p)
               (setf cur-qw (create-crv-from-num-and-denom num-qp den))
               ;; the C reads Q's indexes and degree into n, p and m
               (setf n (curve-n cur-qw) p (curve-degree cur-qw) m (curve-m cur-qw)))))
         (unless rat-p (crv-rat-to-non-rat cur-p))
         (unless rat-q (crv-rat-to-non-rat cur-q)))
        (t
         ;; working curves
         (setf cur-pw (crv-copy cur-p))
         (setf cur-qw (crv-copy cur-q))
         (locals-p)
         (locals-q)
         ;; degree elevate
         (when (/= p q)
           (if (< p q)
               (progn (setf cur-pw (crv-elevate-degree cur-pw (- q p)))
                      (locals-p))
               (progn (setf cur-qw (crv-elevate-degree cur-qw (- p q)))
                      (locals-q))))
         ;; refine
         (multiple-value-bind (kxp kxq)
             (get-compatible-knot-array-mult (curve-knots cur-pw) (curve-knots cur-qw) p)
           (declare (type real-vector kxp kxq))
           (when (>= (1- (length kxp)) 0)
             (setf cur-pw (crv-refine cur-pw kxp))
             (locals-p))
           (when (>= (1- (length kxq)) 0)
             (setf cur-qw (crv-refine cur-qw kxq))
             (locals-q))))))
    ;; surface indexes
    (let ((ns 0) (ms 0) (ps 0) (qs 0) (rs 0) (ss 0))
      (declare (type fixnum ns ms ps qs rs ss))
      (cond ((= dir +udir+) (setf ns 1 ps 1 rs 3 ms n qs p ss m))
            ((= dir +vdir+) (setf ns n ps p rs m ms 1 qs 1 ss 3))
            (t (nurbs-error :inp-err 'create-ruled-srf-from-boundary-crvs)))
      (let ((sw (make-cpoint-net (1+ ns) (1+ ms)))
            (u (make-real-vector (1+ rs)))
            (v (make-real-vector (1+ ss)))
            (pw (curve-cpoints cur-pw))
            (qw (curve-cpoints cur-qw))
            (cp (make-cpoint)))
        (declare (type cpoint-net sw) (type real-vector u v) (type cpoint-array pw qw)
                 (dynamic-extent cp))
        (cond
          ((= dir +udir+)
           (dotimes (j (1+ n))
             (%c-pt-net (%row-c-pt pw j cp) sw 0 j)
             (%c-pt-net (%row-c-pt qw j cp) sw 1 j))
           (dotimes (i 2)
             (setf (aref u i) 0d0
                   (aref u (+ 2 i)) 1d0))
           (dotimes (j (1+ m))
             (setf (aref v j) (aref up j))))
          (t
           (dotimes (i (1+ n))
             (%c-pt-net (%row-c-pt pw i cp) sw i 0)
             (%c-pt-net (%row-c-pt qw i cp) sw i 1))
           (dotimes (i (1+ m))
             (setf (aref u i) (aref up i)))
           (dotimes (j 2)
             (setf (aref v j) 0d0
                   (aref v (+ 2 j)) 1d0))))
        (%make-surface ps qs u v sw)))))

;;; ------------------------------------------------------------- L7-g1

(defun create-revolved-srf (cur s tv al ctp)
  "A fresh surface of revolution: curve CUR revolved by the angle AL (degrees,
clockwise looking along TV) about the axis through S along TV.  The u
direction is the circle, of the kind CTP names (+quadratic+ or +quintic+; the
C's weight array holds 9, so +quartic+ only below 241 degrees), the v
direction CUR's.  Each control point is revolved about its projection on the
axis, weighted with its own weight; one on the axis (closer than *mtol*) stays
there, with the circle's weights.  CUR is made rational and 3D in place first,
as the C does.  Signals :inp-err when a projection fails, :cal-err when every
control point is on the axis or for another CTP, and what the projection and
the circle signal.
NLib: N_CreateRevolvedSrf."
  (declare (type curve cur) (type point s tv) (type double-float al) (type fixnum ctp))
  (multiple-value-bind (al n p r) (crv-get-conic-data ctp 0d0 al)
    (declare (type double-float al) (type fixnum n p r))
    (let ((pw (curve-cpoints cur))
          (m (curve-n cur))
          (q (curve-degree cur))
          (v (curve-knots cur))
          (ss (curve-m cur)))
      (declare (type cpoint-array pw) (type fixnum m q ss) (type real-vector v))
      ;; make the profile curve 3D rational
      (unless (is-crv-rat cur)
        (crv-non-rat-to-rat cur))
      (unless (crv-is-3d cur)
        (crv-2d-to-3d cur))
      (let ((qw (make-cpoint-net (1+ n) (1+ m)))
            (uq (make-real-vector (1+ r)))
            (vq (make-real-vector (1+ ss)))
            (wi (make-real-vector 9))
            (lst (%make-lineseg (make-point) (make-point) +unbounded+))
            (allrzero t)
            (cura nil)
            (a (make-point)) (c (make-point)) (x (make-point)) (y (make-point))
            (cp (make-cpoint)))
        (declare (type cpoint-net qw) (type real-vector uq vq wi) (type lineseg lst)
                 (type point a c x y) (type cpoint cp))
        ;; the circle weights
        (calc-circ-weights 0d0 al ctp wi)
        ;; the surface control points
        (create-line-start-dir-vector lst s tv +unbounded+)
        (dotimes (j (1+ m))
          (c-pt-to-pt-euclid (%row-c-pt pw j cp) a)
          (let ((wp (aref pw j 3)))
            (declare (type double-float wp))
            (when (= (the fixnum (values (project-pt-line lst a c))) +false+)
              (nurbs-error :inp-err 'create-revolved-srf))
            (vector-diff a c x)
            (let ((rad (vector-magnitude x)))
              (declare (type double-float rad))
              (if (< rad (the double-float *mtol*))
                  (dotimes (i (1+ n))
                    (%c-pt-net (weight c (* (aref wi i) wp) cp) qw i j))
                  (progn
                    (setf allrzero nil)
                    (vector-cross tv x y)
                    (setf cura (create-circ-arc c x y rad 0d0 al ctp))
                    (let ((aw (curve-cpoints cura)))
                      (declare (type cpoint-array aw))
                      (dotimes (i (1+ n))
                        (%c-pt-net (scale-c-pt wp (%row-c-pt aw i cp) cp) qw i j))))))))
        (when allrzero
          (nurbs-error :cal-err 'create-revolved-srf))
        ;; the knots
        (let ((u (curve-knots cura)))
          (declare (type real-vector u))
          (dotimes (i (1+ r))
            (setf (aref uq i) (aref u i))))
        (dotimes (j (1+ ss))
          (setf (aref vq j) (aref v j)))
        (%make-surface p q uq vq qw)))))
