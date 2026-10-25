;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvAdv.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------- control point rows (L0-g7)

;;; NLib applies N_CopyCPt, N_Combine2CPts and N_VectorBlendCPt to elements of
;;; NL_CPOINT arrays; here those arrays are cpoint-arrays, so these do the same
;;; arithmetic on rows.  Each reads every input coordinate before writing.

(declaim (inline %copy-c-pt-row %combine-2-c-pt-rows %vector-blend-c-pt-row))

(defun %copy-c-pt-row (a ai b bi)
  "Row BI of B = row AI of A, as N_CopyCPt."
  (declare (type cpoint-array a b) (type fixnum ai bi))
  (let ((x (aref a ai 0)) (y (aref a ai 1)) (z (aref a ai 2)) (w (aref a ai 3)))
    (setf (aref b bi 0) x (aref b bi 1) y (aref b bi 2) z (aref b bi 3) w)
    b))

(defun %combine-2-c-pt-rows (alpha a ai beta b bi c ci)
  "Row CI of C = ALPHA * row AI of A + BETA * row BI of B, honoring +no-z+ and
+no-w+ as N_Combine2CPts does."
  (declare (type double-float alpha beta) (type cpoint-array a b c) (type fixnum ai bi ci))
  (let ((ax (aref a ai 0)) (ay (aref a ai 1)) (az (aref a ai 2)) (aw (aref a ai 3))
        (bx (aref b bi 0)) (by (aref b bi 1)) (bz (aref b bi 2)) (bw (aref b bi 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref c ci 0) (+ (* alpha ax) (* beta bx))
          (aref c ci 1) (+ (* alpha ay) (* beta by))
          (aref c ci 2) (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ (* alpha az) (* beta bz)))
                              ((/= az +no-z+) (* alpha az))
                              ((/= bz +no-z+) (* beta bz))
                              (t +no-z+))
          (aref c ci 3) (if (and (/= aw +no-w+) (/= bw +no-w+)) (+ (* alpha aw) (* beta bw)) +no-w+))
    c))

(defun %vector-blend-c-pt-row (alpha a ai b bi)
  "Row BI of B = row BI of B + ALPHA * row AI of A, honoring +no-z+ and +no-w+
of A as N_VectorBlendCPt does."
  (declare (type double-float alpha) (type cpoint-array a b) (type fixnum ai bi))
  (let ((ax (aref a ai 0)) (ay (aref a ai 1)) (az (aref a ai 2)) (aw (aref a ai 3))
        (bx (aref b bi 0)) (by (aref b bi 1)) (bz (aref b bi 2)) (bw (aref b bi 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref b bi 0) (+ bx (* alpha ax))
          (aref b bi 1) (+ by (* alpha ay))
          (aref b bi 2) (if (/= az +no-z+) (+ bz (* alpha az)) +no-z+)
          (aref b bi 3) (if (/= aw +no-w+) (+ bw (* alpha aw)) +no-w+))
    b))

;;; --------------------------------------------- Bezier decomposition (L0-g7)

(defun %next-bez-cpoints (nqw)
  "NQW, the next Bezier segment's control points, checked: where the C would
write through its still-NULL pointer (an unclamped single-span curve) and
crash, signal instead."
  (if (typep nqw 'cpoint-array)
      nqw
      (nurbs-error :cal-err '%crv-decompose-bez)))

(defun %crv-decompose-bez (cur)
  "A simple-vector of the Bezier segments of curve CUR, one degree-p curve
with knots (a ... a b ... b) per non-zero span, extracted without knot
refinement.  A port of N_CrvDecomposeBez (NL_CrvTool.cpp:811), which the
manifest does not list (its definition line carries a comment between the name
and the argument list); it belongs in crv-tool.lisp as crv-decompose-bez."
  (declare (type curve cur))
  (let* ((p (curve-degree cur))
         (pw (curve-cpoints cur))
         (up (curve-knots cur))
         (m (curve-m cur))
         (nsp (basis-get-span-count up p))
         (cura (make-array nsp)))
    (declare (type fixnum p m nsp) (type cpoint-array pw) (type real-vector up)
             (type simple-vector cura))
    (dotimes (i nsp)
      (setf (svref cura i) (%make-curve p (make-real-vector (+ p p 2)) (make-cpoint-array (1+ p)))))
    (fast
      (let ((alfs (make-array (1+ p) :element-type 'double-float :initial-element 0d0))
            (omas (make-array (1+ p) :element-type 'double-float :initial-element 0d0))
            (is p)
            (ie (1+ p))
            (iq -1)
            (nqw nil))
        (declare (dynamic-extent alfs omas) (type fixnum is ie iq))
        ;; skip over first knot over-multiplicity (bad data correction)
        (do () ((not (and (< ie m) (= (aref up (1- ie)) (aref up ie)))))
          (incf is)
          (incf ie))
        (let ((qw (curve-cpoints (svref cura 0))))
          (declare (type cpoint-array qw))
          (dotimes (i (1+ p))
            (%copy-c-pt-row pw (- (+ i is) p) qw i)))
        (do () ((not (< ie m)))
          (incf iq)
          (when (>= iq nsp)
            ;; NL_QUIT: the C fails without setting an error code
            (nurbs-error :cal-err '%crv-decompose-bez))
          (let ((qw (curve-cpoints (svref cura iq)))
                (uq (curve-knots (svref cura iq))))
            (declare (type cpoint-array qw) (type real-vector uq))
            (when (< iq (1- nsp))
              (setf nqw (curve-cpoints (svref cura (1+ iq)))))
            (let ((i ie))
              (declare (type fixnum i))
              (do () ((not (and (< ie m) (= (aref up ie) (aref up (1+ ie))))))
                (incf ie))
              (let* ((mlt (+ (- ie i) 1))
                     (r (- p mlt)))
                (declare (type fixnum mlt r))
                ;; insert the knot on the rhs only
                (when (< mlt p)
                  (let ((num (- (aref up ie) (aref up is))))
                    (declare (type double-float num))
                    (do ((i p (1- i))) ((<= i mlt))
                      (declare (type fixnum i))
                      (setf (aref alfs (- i mlt 1)) (/ num (- (aref up (+ is i)) (aref up is))))
                      (setf (aref omas (- i mlt 1)) (- 1d0 (aref alfs (- i mlt 1))))))
                  (do ((i 1 (1+ i))) ((> i r))
                    (declare (type fixnum i))
                    (let ((s (+ mlt i))
                          (save (- r i)))
                      (declare (type fixnum s save))
                      (do ((j p (1- j))) ((< j s))
                        (declare (type fixnum j))
                        (%combine-2-c-pt-rows (aref alfs (- j s)) qw j
                                              (aref omas (- j s)) qw (1- j)
                                              qw j))
                      (when (< ie m)
                        (%copy-c-pt-row qw p (%next-bez-cpoints nqw) save)))))
                ;; the Bezier segment's knot vector
                (dotimes (i (1+ p))
                  (setf (aref uq i) (aref up is)
                        (aref uq (+ i p 1)) (aref up ie)))
                ;; the next segment's remaining control points
                (when (< ie m)
                  (when (< r 0) (setf r 0))
                  (do ((i r (1+ i))) ((> i p))
                    (declare (type fixnum i))
                    (%copy-c-pt-row pw (+ (- ie p) i) (%next-bez-cpoints nqw) i)))
                (setf is ie
                      ie (1+ ie)))))))
      cura)))

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
         (bez (%crv-decompose-bez cur))
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
