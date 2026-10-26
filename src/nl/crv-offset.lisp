;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvOffset.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------- L6-g0

(defun crv-offset-pt-sampling (cur-p b d deg tan par tol)
  "An approximation of the offset by distance D of the planar, G1 curve CUR-P
along B x T (B the normal of its plane, not changed; T the unit tangent).  A
line (to *mtol*) is offset exactly, and so is a circle, made of the kind DEG
picks (4 quartic, 5 quintic, else quadratic) with radius r + D.  Otherwise the
offset is sampled at points spread over each Bezier segment by the bound on
its second derivative and tolerance TOL, interpolated by a curve of degree DEG
at CUR-P's parameters (PAR +inherited+) or ones computed from the points
(+uniform+, +chordlength+, +centripetal+), with the end tangent directions
of CUR-P kept when TAN is +yes+, and its knots removed within TOL.  A 2D
CUR-P gives a 2D offset.  Returns the fresh offset curve.  Signals :num-err
when 1/TOL is bad, and what crv-get-type, vector-normalize-ref,
crv-offset-get-max-2nd-deriv, crv-eval-tangent, fit-calc-crv-param-values,
crv-derivs, the interpolators and crv-remove-knots-deriv-constraints signal.
NLib: N_CrvOffsetPtSampling."
  (declare (type curve cur-p) (type point b) (type double-float d tol)
           (type fixnum deg tan par))
  (let ((q (make-point)) (v (make-point)) (x (make-point)) (y (make-point))
        (dv (make-point))
        (twod (if (crv-is-3d cur-p) +no+ +yes+)))
    (declare (dynamic-extent q v x y dv) (type fixnum twod))
    ;; the curve type
    (multiple-value-bind (r al typ) (crv-get-type cur-p *mtol* q v x y)
      (declare (type double-float r al) (type fixnum typ))
      ;; offset line
      (when (= typ +nline+)
        (vector-cross b v dv)
        (vector-normalize-ref dv)
        (vector-blend-pt d dv q)
        (let ((cur-q (crv-line-from-pt-and-vector q v)))
          (when (= twod +yes+)
            (crv-3d-to-2d cur-q))
          (return-from crv-offset-pt-sampling cur-q)))
      ;; offset circle
      (when (= typ +ncircle+)
        (let* ((ctp (cond ((= deg 4) +quartic+)
                          ((= deg 5) +quintic+)
                          (t +quadratic+)))
               (cur-q (create-circ-arc q x y (+ r d) 0.0d0 al ctp)))
          (declare (type fixnum ctp))
          (when (= twod +yes+)
            (crv-3d-to-2d cur-q))
          (return-from crv-offset-pt-sampling cur-q))))
    ;; offset free-form curve
    (let* ((up (curve-knots cur-p))
           (mp (curve-m cur-p))
           (p (curve-degree cur-p))
           (bn (make-point))
           (expo (if (= deg 1) 0.5d0 0.34d0))
           (f1 0d0) (f2 0d0))
      (declare (type real-vector up) (type fixnum mp p) (dynamic-extent bn)
               (type double-float expo f1 f2))
      (when (float-op-is-bad 1.0d0 tol +division+)
        (nurbs-error :num-err 'crv-offset-pt-sampling))
      (setf f1 (/ 1.0d0 tol))
      (setf f1 (c-pow f1 expo))
      (setf f2 (* 2.0d0 (float deg 1d0)))
      ;; B is passed by value in C
      (vector-copy b bn)
      (vector-normalize-ref bn)
      ;; the Bezier segments
      (multiple-value-bind (bez k) (crv-decompose-bez cur-p)
        (declare (type simple-vector bez) (type fixnum k))
        ;; the number of sampling points
        (let ((nu (make-array (1+ k) :element-type 'fixnum :initial-element 0))
              (mu 0)
              (unitspan (make-real-vector 2)))
          (declare (type (simple-array fixnum (*)) nu) (type fixnum mu))
          (create-interval unitspan 0d0 1d0)
          (dotimes (i (1+ k))
            (crv-reparam-to-interval (svref bez i) unitspan)
            (let* ((muu (crv-offset-get-max-2nd-deriv (svref bez i) d bn))
                   (num (* 0.125d0 muu))
                   (gro (nl-max f2 (c-sqrt num)))
                   (del (* f1 gro)))
              (declare (type double-float muu num gro del))
              (setf del (fceiling del))
              (setf (aref nu i) (truncate del))
              (incf mu (aref nu i))))
          ;; the parameters
          (let ((u (make-real-vector (1+ mu)))
                (ul (aref up p))
                (ur 0d0)
                (i (1+ p)) (j 0) (l 0))
            (declare (type real-vector u) (type double-float ul ur) (type fixnum i j l))
            (setf (aref u 0) ul)
            (do () ((not (< i mp)))
              (do () ((not (and (< i mp) (= (aref up i) (aref up (1+ i))))))
                (incf i))
              (setf ur (aref up i))
              (let ((uinc (/ (- ur ul) (float (aref nu j) 1d0))))
                (declare (type double-float uinc))
                (do ((kk 1 (1+ kk))) ((>= kk (aref nu j)))
                  (declare (type fixnum kk))
                  (setf (aref u (incf l)) (+ ul (* (float kk 1d0) uinc)))))
              (setf (aref u (incf l)) ur)
              (setf ul ur)
              (incf i)
              (incf j))
            ;; sample the offset
            (let ((pts (make-point-array (1+ mu)))
                  (ts (make-point)) (te (make-point)) (nv (make-point)) (pi* (make-point)))
              (declare (type point-array pts) (dynamic-extent ts te nv pi*))
              (dotimes (i (1+ mu))
                (crv-eval-tangent cur-p (aref u i) +left+ q ts)
                (vector-cross bn ts nv)
                (vector-pt-along-vector q d nv pi*)
                (%pt-row pi* pts i))
              ;; the knot vector
              (let* ((jj (if (= tan +yes+) (+ mu deg 3) (+ mu deg 1)))
                     (knu (make-real-vector (1+ jj)))
                     (cur-q nil))
                (declare (type fixnum jj) (type real-vector knu))
                (when (/= par +inherited+)
                  (fit-calc-crv-param-values pts mu +epoint+ par u))
                (if (= tan +yes+)
                    (fit-calc-knot-vector-end-derivs u mu deg knu)
                    (fit-crv-calc-knot-vector u mu deg knu))
                ;; interpolate the points
                (if (= tan +yes+)
                    (let ((cd (make-point-array 2)))
                      (declare (dynamic-extent cd))
                      ;; the end tangents
                      (crv-derivs cur-p (aref up 0) +left+ 1 cd)
                      (%row-pt cd 1 ts)
                      (crv-derivs cur-p (aref up mp) +left+ 1 cd)
                      (%row-pt cd 1 te)
                      (when (/= par +inherited+)
                        ;; adjust the magnitudes
                        (let ((del (dist-polygon pts mu)))
                          (declare (type double-float del))
                          (vector-normalize-ref ts)
                          (vector-normalize-ref te)
                          (vector-scale ts del ts)
                          (vector-scale te del te)))
                      (setf cur-q (fit-crv-knots-and-derivs pts mu +epoint+ u knu deg ts te)))
                    (setf cur-q (fit-crv-interp-given-params pts mu +epoint+ u knu deg)))
                (when (= twod +yes+)
                  (crv-3d-to-2d cur-q))
                ;; remove knots
                (setf cur-q (crv-remove-knots-deriv-constraints
                             cur-q tol +both+ (if (and (= tan +yes+) (= deg 2)) 1 0)))
                ;; compact the output curve
                (when (< (curve-n cur-q) mu)
                  (crv-compress cur-q))
                ;; rescale the knot vector if necessary
                (when (/= par +inherited+)
                  (let ((uq (curve-knots cur-q))
                        (mq (curve-m cur-q)))
                    (declare (type real-vector uq) (type fixnum mq))
                    (when (or (/= (aref up 0) (aref uq 0)) (/= (aref up mp) (aref uq mq)))
                      (let ((iv (make-real-vector 2)))
                        (create-interval iv (aref up 0) (aref up mp))
                        (crv-reparam-to-interval cur-q iv)))))
                cur-q))))))))
