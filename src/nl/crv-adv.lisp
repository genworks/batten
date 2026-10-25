;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvAdv.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; Control point rows: %copy-c-pt-row, %combine-2-c-pt-rows and
;;; %vector-blend-c-pt-row in math.lisp.  Bezier decomposition:
;;; crv-decompose-bez in crv-tool.lisp.

;;; ------------------------------------------------------------------ routines

(defun crv-derivs-at-knot (cur u flg der dw)
  "Write the derivatives of curve CUR at U in homogeneous space, orders 0
through DER, into DW, a cpoint-array of at least (1+ DER) rows: DW[k] is the
kth derivative.  FLG +left+ means u is in [u_j, u_j+1) (right derivatives),
+right+ u in (u_j, u_j+1] (left derivatives).  Return DW.
NLib: N_CrvDerivsAtKnot."
  (declare (type curve cur) (type double-float u) (type fixnum flg der) (type cpoint-array dw))
  (let ((knots (curve-knots cur))
        (pw (curve-cpoints cur))
        (p (curve-degree cur)))
    (declare (type real-vector knots) (type cpoint-array pw) (type fixnum p))
    (knot-vector-is-param-out-of-bounds knots u 'crv-derivs-at-knot)
    (let* ((bd (make-real-matrix (1+ der) (1+ p)))
           (spn (basis-derivs knots p u flg der bd)))
      (declare (type real-matrix bd) (type fixnum spn))
      (fast
        (dotimes (k (1+ der))
          (setf (aref dw k 0) 0d0 (aref dw k 1) 0d0 (aref dw k 2) 0d0 (aref dw k 3) 0d0)
          (dotimes (i (1+ p))
            (%vector-blend-c-pt-row (aref bd k i) pw (+ (- spn p) i) dw k))))
      dw)))

(defun crv-get-max-2nd-deriv (cur)
  "An upper bound on the magnitude of the second derivative of curve CUR,
from sampling each Bezier segment at its nodes.  NLib: N_CrvGetMax2ndDeriv."
  (declare (type curve cur))
  (let* ((p (curve-degree cur))
         (bez (crv-decompose-bez cur))
         (k (1- (length bez)))
         (cd (make-point-array 3))
         (c2 (make-point))
         (u (make-real-vector (1+ p)))
         (duu 0d0))
    (declare (type fixnum p k) (type simple-vector bez) (type point-array cd)
             (type point c2) (type real-vector u) (type double-float duu))
    (do ((i 0 (1+ i))) ((> i k))
      (declare (type fixnum i))
      (let* ((al (if (= i 0) 0 1))
             (b (svref bez i))
             (uk (curve-knots b))
             (r (1- (length uk)))
             (ui (/ (- (aref uk r) (aref uk 0)) (float p 1d0))))
        (declare (type fixnum al r) (type real-vector uk) (type double-float ui))
        (setf (aref u 0) (aref uk 0)
              (aref u p) (aref uk r))
        (do ((a 1 (1+ a))) ((> a (1- p)))
          (declare (type fixnum a))
          (setf (aref u a) (+ (aref uk 0) (* (float a 1d0) ui))))
        (do ((a al (1+ a))) ((> a p))
          (declare (type fixnum a))
          (crv-derivs b (aref u a) +left+ 2 cd)
          (setf (aref c2 0) (aref cd 2 0) (aref c2 1) (aref cd 2 1) (aref c2 2) (aref cd 2 2))
          (let ((muu (pt-magnitude c2)))
            (declare (type double-float muu))
            (when (> muu duu)
              (setf duu muu))))))
    duu))

(defun crv-non-rat-get-first-deriv (cur d)
  "The Dth derivative curve of NON-RATIONAL curve CUR, 0 <= D <= degree, as a
fresh curve of degree p-D.  Signals :inp-err for a D out of range or a rational
curve, :der-err for a zero-length knot span in the difference, :num-err when a
coefficient would overflow.  NLib: N_CrvNonRatGetFirstDeriv."
  (declare (type curve cur) (type fixnum d))
  (let* ((n (curve-n cur))
         (pw (curve-cpoints cur))
         (p (curve-degree cur))
         (m (curve-m cur))
         (u (curve-knots cur)))
    (declare (type fixnum n p m) (type cpoint-array pw) (type real-vector u))
    (when (or (< d 0) (> d p))
      (nurbs-error :inp-err 'crv-non-rat-get-first-deriv))
    (when (is-crv-rat cur)
      (nurbs-error :inp-err 'crv-non-rat-get-first-deriv))
    (let* ((nd (- n d))
           (pd (- p d))
           (md (+ nd pd 1))
           (dw (make-cpoint-array (1+ nd)))
           (v (make-real-vector (1+ md)))
           (tw (make-cpoint-array (1+ n))))
      (declare (type fixnum nd pd md) (type cpoint-array dw tw) (type real-vector v))
      (fast
        (dotimes (i (1+ n))
          (%copy-c-pt-row pw i tw i))
        ;; control points of the derivative curve
        (do ((k 1 (1+ k))) ((> k d))
          (declare (type fixnum k))
          (do ((i 0 (1+ i))) ((> i (- n k)))
            (declare (type fixnum i))
            (when (= (aref u (+ i p 1)) (aref u (+ i k)))
              (nurbs-error :der-err 'crv-non-rat-get-first-deriv))
            (when (float-op-is-bad (+ (- (float p 1d0) (float k 1d0)) 1d0)
                                   (- (aref u (+ i p 1)) (aref u (+ i k)))
                                   +division+)
              (nurbs-error :num-err 'crv-non-rat-get-first-deriv))
            (let ((alf (/ (+ (- (float p 1d0) (float k 1d0)) 1d0)
                          (- (aref u (+ i p 1)) (aref u (+ i k))))))
              (declare (type double-float alf))
              (%combine-2-c-pt-rows alf tw (1+ i) (- alf) tw i tw i))))
        (do ((i 0 (1+ i))) ((> i nd))
          (declare (type fixnum i))
          (%copy-c-pt-row tw i dw i))
        ;; the knot vector
        (let ((j -1))
          (declare (type fixnum j))
          (do ((i 0 (1+ i))) ((> i pd))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u 0)))
          (do ((i (1+ p) (1+ i))) ((> i n))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u i)))
          (do ((i 0 (1+ i))) ((> i pd))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u m)))))
      (%make-curve pd v dw))))

(defun crv-planar-offset-get-deriv (cur u d nv ufl od)
  "Write the first and second derivatives at U of the offset by distance D of
planar curve CUR into rows 1 and 2 of OD, a point-array of at least 3 rows
allocated by the caller (row 0 is left alone).  NV is the normal to the
curve's plane: the curve is offset in the direction NV x T.  UFL +left+ means
u is in [u_j, u_j+1) (right derivatives), +right+ u in (u_j, u_j+1] (left
derivatives).  Signals as crv-derivs and crv-get-curvature-deriv do, and
:num-err when the curvature's division would overflow.  Returns OD.
NLib: N_CrvPlanarOffsetGetDeriv."
  (declare (type curve cur) (type double-float u d) (type point nv) (type fixnum ufl)
           (type point-array od))
  (let ((dd (make-point-array 3))
        (d1 (make-point)) (d2 (make-point)) (b (make-point)) (o (make-point))
        (num 0d0) (den 0d0) (kap 0d0) (kpp 0d0) (f1 0d0) (f2 0d0) (sgn 0d0))
    (declare (type point-array dd) (type point d1 d2 b o)
             (type double-float num den kap kpp f1 f2 sgn)
             (dynamic-extent dd d1 d2 b o))
    ;; derivatives, the curvature and its derivative
    (crv-derivs cur u ufl 2 dd)
    (%point-array-row dd 1 d1)
    (%point-array-row dd 2 d2)
    (vector-cross d1 d2 b)
    (setf num (vector-magnitude b))
    (setf den (vector-magnitude d1))
    (setf den (* (* den den) den))
    (when (float-op-is-bad num den +division+)
      (nurbs-error :num-err 'crv-planar-offset-get-deriv))
    (setf kap (/ num den))
    (setf den (vector-dot b nv))
    (setf sgn (if (> den 0d0) -1d0 1d0))
    (setf kpp (crv-get-curvature-deriv cur u ufl))
    ;; offset derivatives
    (setf f1 (+ 1d0 (* (* sgn kap) d)))
    (vector-scale d1 f1 o)
    (setf (aref od 1 0) (aref o 0) (aref od 1 1) (aref o 1) (aref od 1 2) (aref o 2))
    (setf f2 (* (* sgn kpp) d))
    (vector-combine f1 d2 f2 d1 o)
    (setf (aref od 2 0) (aref o 0) (aref od 2 1) (aref o 1) (aref od 2 2) (aref o 2))
    od))
