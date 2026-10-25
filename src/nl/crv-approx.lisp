;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvApprox.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------- L6-g1

(defun approx-nurbs-with-non-rat-crv (curp tol deg par tn)
  "A fresh non-rational curve of degree DEG approximating curve CURP to
within TOL: each Bezier segment of CURP is sampled at a number of evenly
spaced parameters set by TOL and its second derivative bound, the samples are
interpolated (TN +yes+: with CURP's end tangents, scaled to its arc length
unless PAR is +inherited+), parametrized by PAR (+uniform+, +chordlength+,
+centripetal+, or +inherited+: CURP's own sample parameters), and then as
many knots are removed as TOL allows (keeping the end derivative knots for a
tangent-constrained quadratic).  Unless PAR is +inherited+, the result is
reparametrized onto CURP's domain.  A 2D CURP gives a 2D curve.  TOL must be
positive.  Signals :num-err when 1/TOL is out of range, and what the
decomposition, the sampling, the parametrization (:cal-err for another PAR),
the interpolation (:inp-err for TN +yes+ and DEG 1), the normalization of a
zero end tangent and the knot removal signal.
NLib: N_ApproxNurbsWithNonRatCrv."
  (declare (type curve curp) (type double-float tol) (type fixnum deg par tn))
  (let* ((mp (curve-m curp))
         (up (curve-knots curp))
         (p (curve-degree curp))
         (twod (not (crv-is-3d curp)))
         (ex (if (= deg 1) 0.5d0 0.34d0))
         (pct 0.01d0)
         (f1 0d0) (f2 0d0))
    (declare (type fixnum mp p) (type real-vector up) (type double-float ex pct f1 f2))
    (when (float-op-is-bad 1d0 tol +division+)
      (nurbs-error :num-err 'approx-nurbs-with-non-rat-crv))
    (setf f1 (/ 1d0 tol))
    (setf f1 (expt f1 ex))
    (setf f2 (float deg 1d0))
    ;; the Bezier segments
    (multiple-value-bind (bez k) (crv-decompose-bez curp)
      (declare (type simple-vector bez) (type fixnum k))
      ;; the number of sampling points per segment
      (let ((nu (make-index-vector (1+ k)))
            (mu 0)
            (unitspan (make-real-vector 2)))
        (declare (type index-vector nu) (type fixnum mu) (type real-vector unitspan))
        (create-interval unitspan 0d0 1d0)
        (dotimes (i (1+ k))
          (crv-reparam-to-interval (svref bez i) unitspan)
          (let* ((muu (crv-get-max-2nd-deriv (svref bez i)))
                 (num (* 0.125d0 muu))
                 (gro (nl-max f2 (sqrt (the (double-float 0d0) num))))
                 (del (* f1 gro)))
            (declare (type double-float muu num gro del))
            (setf del (fceiling del))
            (setf (aref nu i) (truncate del))
            (setf mu (+ mu (aref nu i)))))
        ;; the parameters
        (let ((u (make-real-vector (1+ mu)))
              (ul (aref up p))
              (ur 0d0)
              (i (1+ p))
              (j 0)
              (l 0))
          (declare (type real-vector u) (type double-float ul ur) (type fixnum i j l))
          (setf (aref u 0) ul)
          (do () ((not (< i mp)))
            (do () ((not (and (< i mp) (= (aref up i) (aref up (1+ i))))))
              (incf i))
            (setf ur (aref up i))
            (let ((uinc (/ (- ur ul) (float (aref nu j) 1d0))))
              (declare (type double-float uinc))
              (do ((kk 1 (1+ kk))) ((not (< kk (aref nu j))))
                (declare (type fixnum kk))
                (setf (aref u (incf l)) (+ ul (* (float kk 1d0) uinc))))
              (setf (aref u (incf l)) ur))
            (setf ul ur)
            (incf i)
            (incf j))
          ;; sample the curve
          (let ((pts (make-point-array (1+ mu)))
                (pt (make-point))
                (knu (make-real-vector (if (= tn +yes+) (+ mu deg 4) (+ mu deg 2))))
                (curq nil))
            (declare (type point-array pts) (type point pt) (type real-vector knu))
            (dotimes (ii (1+ mu))
              (crv-eval curp (aref u ii) +left+ pt)
              (%pt-row pt pts ii))
            ;; the knot vector
            (when (/= par +inherited+)
              (fit-calc-crv-param-values pts mu +epoint+ par u))
            (if (= tn +yes+)
                (fit-calc-knot-vector-end-derivs u mu deg knu)
                (fit-crv-calc-knot-vector u mu deg knu))
            ;; interpolate the points
            (cond
              ((= tn +yes+)
               ;; the end tangents
               (let ((cd (make-point-array 2))
                     (ts (make-point))
                     (te (make-point)))
                 (declare (type point-array cd) (type point ts te))
                 (crv-derivs curp (aref up 0) +left+ 1 cd)
                 (%row-pt cd 1 ts)
                 (crv-derivs curp (aref up mp) +left+ 1 cd)
                 (%row-pt cd 1 te)
                 (when (/= par +inherited+)
                   ;; adjust the magnitudes
                   (let ((del (crv-arc-length curp (aref up 0) (aref up mp) pct +relative+)))
                     (declare (type double-float del))
                     (vector-normalize-ref ts)
                     (vector-normalize-ref te)
                     (vector-scale ts del ts)
                     (vector-scale te del te)))
                 (setf curq (fit-crv-knots-and-derivs pts mu +epoint+ u knu deg ts te))))
              (t
               (setf curq (fit-crv-interp-given-params pts mu +epoint+ u knu deg))))
            (when twod
              (crv-3d-to-2d curq))
            ;; remove knots
            (let ((der (if (and (= tn +yes+) (= deg 2)) 1 0)))
              (declare (type fixnum der))
              (setf curq (crv-remove-knots-deriv-constraints curq tol +both+ der)))
            ;; compact the output curve
            (when (< (curve-n curq) mu)
              (crv-compress curq))
            ;; rescale the knot vector if necessary
            (when (/= par +inherited+)
              (let* ((uq (curve-knots curq))
                     (mq (1- (length uq))))
                (declare (type real-vector uq) (type fixnum mq))
                (when (or (/= (aref up 0) (aref uq 0)) (/= (aref up mp) (aref uq mq)))
                  (let ((iv (make-real-vector 2)))
                    (create-interval iv (aref up 0) (aref up mp))
                    (crv-reparam-to-interval curq iv)))))
            curq))))))
