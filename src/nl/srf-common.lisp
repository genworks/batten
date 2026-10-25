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
