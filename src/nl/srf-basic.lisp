;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(declaim (inline is-srf-rat))
(defun is-srf-rat (sur)
  "True when surface SUR is rational.  NLib: N_IsSrfRat."
  (declare (type surface sur))
  (/= (aref (surface-cpoints sur) 0 0 3) +no-w+))

(defun srf-get-parameter-bounds (sur)
  "The parameter domain of SUR, (values ul ur vb vt).  NLib: N_SrfGetParameterBounds."
  (declare (type surface sur))
  (let ((ku (surface-knots-u sur)) (kv (surface-knots-v sur)))
    (values (aref ku 0) (aref ku (1- (length ku)))
            (aref kv 0) (aref kv (1- (length kv))))))

(defun srf-weights-are-valid (sur rname)
  "Signal :wei-err in routine RNAME when a weight of rational surface SUR lies
outside [*wmin*, *wmax*].  Returns SUR otherwise.  NLib: N_SrfWeightsAreValid."
  (declare (type surface sur))
  (when (is-srf-rat sur)
    (let ((pw (surface-cpoints sur))
          (wmin *wmin*) (wmax *wmax*))
      (declare (type cpoint-net pw) (type double-float wmin wmax))
      (dotimes (i (1+ (surface-n sur)))
        (dotimes (j (1+ (surface-m sur)))
          (let ((w (aref pw i j 3)))
            (declare (type double-float w))
            (when (or (< w wmin) (> w wmax))
              (nurbs-error :wei-err rname)))))))
  sur)

(defun srf-from-c-net-and-knot-vectors (sur net p q knu knv)
  "Make SUR the surface of control net NET (a cpoint-net), degrees P and Q
and knot vectors KNU and KNV, sharing them as the C shares its pointers.
Returns SUR.  NLib: N_SrfFromCNetAndKnotVectors."
  (declare (type surface sur) (type cpoint-net net) (type fixnum p q)
           (type real-vector knu knv))
  (setf (surface-cpoints sur) net
        (surface-degree-u sur) p
        (surface-degree-v sur) q
        (surface-knots-u sur) knu
        (surface-knots-v sur) knv)
  sur)

(defun srf-get-net-and-knot-vectors (sur)
  "The control net, degrees and knot vectors of SUR, the surface's own
objects: (values net p q knu knv).  NLib: N_SrfGetNetAndKnotVectors."
  (declare (type surface sur))
  (values (surface-cpoints sur) (surface-degree-u sur) (surface-degree-v sur)
          (surface-knots-u sur) (surface-knots-v sur)))

(defun srf-get-c-pts-knot-vector-and-knots (sur)
  "The control points, knot vectors and knot arrays of SUR, the surface's own
objects: (values pw knu knv u v).  A Batten knot vector is its knot array, so
KNU is U and KNV is V.  NLib: N_SrfGetCPtsKnotVectorAndKnots."
  (declare (type surface sur))
  (values (surface-cpoints sur) (surface-knots-u sur) (surface-knots-v sur)
          (surface-knots-u sur) (surface-knots-v sur)))

(defun srf-set-c-pts-and-knots (sur pw u v)
  "Give SUR the control points PW and the knot arrays U and V, sharing them as
the C shares its pointers.  Returns SUR.  NLib: N_SrfSetCPtsAndKnots."
  (declare (type surface sur) (type cpoint-net pw) (type real-vector u v))
  (setf (surface-cpoints sur) pw
        (surface-knots-u sur) u
        (surface-knots-v sur) v)
  sur)

(defun srf-is-closed (sur dir)
  "True when SUR is closed in direction DIR (+udir+ or +vdir+): the Euclidean
control points of its first and last rows (+udir+) or columns (+vdir+) lie
within *mtol* of each other.  False for any other DIR.  NLib: N_SrfIsClosed."
  (declare (type surface sur) (type fixnum dir))
  (let* ((pw (surface-cpoints sur))
         (n (surface-n sur)) (m (surface-m sur))
         (mtol *mtol*)
         (cw (make-cpoint))
         (ps (make-point)) (pe (make-point)))
    (declare (type cpoint-net pw) (type fixnum n m) (type double-float mtol)
             (dynamic-extent cw ps pe))
    (flet ((euclid (i j p)
             (setf (aref cw 0) (aref pw i j 0) (aref cw 1) (aref pw i j 1)
                   (aref cw 2) (aref pw i j 2) (aref cw 3) (aref pw i j 3))
             (c-pt-to-pt-euclid cw p)))
      (declare (inline euclid))
      (cond
        ((= dir +udir+)
         (dotimes (j (1+ m) t)
           (euclid 0 j ps)
           (euclid n j pe)
           (when (> (dist-pt-pt ps pe) mtol)
             (return nil))))
        ((= dir +vdir+)
         (dotimes (i (1+ n) t)
           (euclid i 0 ps)
           (euclid i m pe)
           (when (> (dist-pt-pt ps pe) mtol)
             (return nil))))
        (t nil)))))

(defun 1d-pt-set-is-closed (a m ptp tol)
  "True when the first and last (index M) points of A lie within TOL of each
other: A is a point-array for PTP +epoint+ (Euclidean distance) or a
cpoint-array for +hpoint+ (homogeneous distance).  False for any other PTP.
NLib: N_1dPtSetIsClosed."
  (declare (type fixnum m ptp) (type double-float tol))
  (cond
    ((= ptp +epoint+)
     (let ((a a) (p0 (make-point)) (pm (make-point)))
       (declare (type point-array a) (dynamic-extent p0 pm))
       (set-point p0 (aref a 0 0) (aref a 0 1) (aref a 0 2))
       (set-point pm (aref a m 0) (aref a m 1) (aref a m 2))
       (not (> (dist-pt-pt p0 pm) tol))))
    ((= ptp +hpoint+)
     (let ((a a) (p0 (make-cpoint)) (pm (make-cpoint)))
       (declare (type cpoint-array a) (dynamic-extent p0 pm))
       (set-cpoint p0 (aref a 0 0) (aref a 0 1) (aref a 0 2) (aref a 0 3))
       (set-cpoint pm (aref a m 0) (aref a m 1) (aref a m 2) (aref a m 3))
       (not (> (dist-cpt-cpt-homo p0 pm) tol))))
    (t nil)))

(defun srf-copy (sur-p)
  "A fresh copy of surface SUR-P.  NLib: N_SrfCopy."
  (declare (type surface sur-p))
  (let* ((n (surface-n sur-p)) (m (surface-m sur-p))
         (pw (surface-cpoints sur-p))
         (up (surface-knots-u sur-p)) (vp (surface-knots-v sur-p))
         (qw (make-cpoint-net (1+ n) (1+ m)))
         (uq (make-real-vector (length up)))
         (vq (make-real-vector (length vp))))
    (declare (type fixnum n m) (type cpoint-net pw qw) (type real-vector up vp uq vq))
    (dotimes (i (1+ n))
      (dotimes (j (1+ m))
        (dotimes (c 4) (setf (aref qw i j c) (aref pw i j c)))))
    (dotimes (i (length up)) (setf (aref uq i) (aref up i)))
    (dotimes (j (length vp)) (setf (aref vq j) (aref vp j)))
    (%make-surface (surface-degree-u sur-p) (surface-degree-v sur-p) uq vq qw)))
