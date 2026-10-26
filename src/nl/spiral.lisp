;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Spiral.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------- L4-g0

;;; NL_Spiral.cpp's file-private constant and variables.  The C keeps the
;;; spiral's data in file statics that N_CrvApproxSpiral sets before any
;;; static reads them; here CRV-APPROX-SPIRAL binds them, so nothing persists
;;; between calls (nothing in the C reads a value left by an earlier call).

(defconstant +spiral--twopi+ 6.2831853071795864d0
  "TWOPI, NL_Spiral.cpp.")

(defvar *spiral--rl* 1 "RL, NL_Spiral.cpp: 1 for a right hand spiral, -1 for a left.")
(defvar *spiral--n-spans* 8 "nSpans, NL_Spiral.cpp.")
(defvar *spiral--d-turns* 1d0 "dTurns, NL_Spiral.cpp: the number of turns.")
(defvar *spiral--h* 0.5d0 "h, NL_Spiral.cpp: the height.")
(defvar *spiral--r0* 1d0 "r0, NL_Spiral.cpp: the radius at z = 0.")
(defvar *spiral--r1* 1d0 "r1, NL_Spiral.cpp: the radius at z = h.")
(defvar *spiral--rdfun* 0d0 "Rdfun, NL_Spiral.cpp: the radius' derivative by angle.")

(declaim (type fixnum *spiral--rl* *spiral--n-spans*)
         (type double-float *spiral--d-turns* *spiral--h* *spiral--r0* *spiral--r1* *spiral--rdfun*))

(defun spiral--st-get-spiral-angle (theta)
  "The spiral's radius at angle THETA: r0 and r1 blended linearly over the
whole sweep 2 pi dTurns.
NLib: ST_GetSpiralAngle (NL_Spiral.cpp)."
  (declare (type double-float theta))
  (let ((r0 *spiral--r0*) (r1 *spiral--r1*) (dt *spiral--d-turns*))
    (/ (+ (* r0 (- (* +spiral--twopi+ dt) theta)) (* r1 theta))
       (* +spiral--twopi+ dt))))

(defun spiral--st-get-point-at-angle (theta xyz)
  "Write into point XYZ the spiral's point at angle THETA (clamped to the end
of the sweep) and return 1.
NLib: ST_GetPointAtAngle (NL_Spiral.cpp)."
  (declare (type double-float theta) (type point xyz))
  (let ((sweep (* +spiral--twopi+ *spiral--d-turns*)))
    (declare (type double-float sweep))
    (when (> theta sweep)
      (setf theta sweep))
    (let ((cos-theta (cos theta))
          (sin-theta (sin theta)))
      (declare (type double-float cos-theta sin-theta))
      (setf (aref xyz 0) (* (spiral--st-get-spiral-angle theta) cos-theta)
            (aref xyz 1) (* (* (float *spiral--rl* 1d0) (spiral--st-get-spiral-angle theta)) sin-theta)
            (aref xyz 2) (/ (* *spiral--h* theta) (* +spiral--twopi+ *spiral--d-turns*)))
      1)))

(defun spiral--st-get-point-at-angle-d (theta xyz)
  "Write into point XYZ the spiral's derivative by angle at THETA (clamped to
the end of the sweep) and return 1.
NLib: ST_GetPointAtAngleD (NL_Spiral.cpp)."
  (declare (type double-float theta) (type point xyz))
  (let ((sweep (* +spiral--twopi+ *spiral--d-turns*)))
    (declare (type double-float sweep))
    (when (> theta sweep)
      (setf theta sweep))
    (let ((cos-theta (cos theta))
          (sin-theta (sin theta))
          (rdfun *spiral--rdfun*))
      (declare (type double-float cos-theta sin-theta rdfun))
      (setf (aref xyz 0) (+ (* (- (spiral--st-get-spiral-angle theta)) sin-theta) (* rdfun cos-theta))
            (aref xyz 1) (* (float *spiral--rl* 1d0)
                            (+ (* (spiral--st-get-spiral-angle theta) cos-theta) (* rdfun sin-theta)))
            (aref xyz 2) (/ *spiral--h* (* +spiral--twopi+ *spiral--d-turns*)))
      1)))

(defun crv-approx-spiral (height rad0 rad1 num-turns ror-l num-points tol)
  "A cubic NURBS spiral around the z axis from (RAD0, 0, 0) to height HEIGHT,
where its radius is RAD1, over NUM-TURNS turns, right handed when ROR-L is 1
and left handed otherwise: the uniform cubic spline through equally spaced
points of the spiral, with its analytic end tangents.  With NUM-POINTS > 1
that many points are interpolated; otherwise the number of spans starts at
8 per turn and grows until the radial error at the span midpoints is below
TOL (raised to 10*mtol), up to *ccplim* - 6 spans.  Returns (values
num-points tol curve): the number of points the fit used (the C's count,
one span past the curve when the span limit stopped the growth), the error
reached and the fresh curve.  Signals :inp-err for a height or a number of
turns below mtol, or radii below -mtol or both below mtol.
NLib: N_CrvApproxSpiral."
  (declare (type double-float height rad0 rad1 num-turns tol) (type fixnum ror-l num-points))
  (let ((mtol *mtol*))
    (declare (type double-float mtol))
    (when (< height mtol)
      (nurbs-error :inp-err 'crv-approx-spiral))
    (when (or (< rad0 (- mtol)) (< rad1 (- mtol)))
      (nurbs-error :inp-err 'crv-approx-spiral))
    (when (< rad0 mtol)
      (setf rad0 0d0))
    (when (< rad1 mtol)
      (setf rad1 0d0))
    (when (< (+ rad0 rad1) mtol)
      (nurbs-error :inp-err 'crv-approx-spiral))
    (when (< num-turns mtol)
      (nurbs-error :inp-err 'crv-approx-spiral))
    (when (< tol (* 10d0 mtol))
      (setf tol (* 10d0 mtol))))
  (let* ((n-spans-max (- (the fixnum *ccplim*) 6))
         (*spiral--rl* (if (= ror-l 1) 1 -1))
         (*spiral--h* height)
         (*spiral--r0* rad0)
         (*spiral--r1* rad1)
         (*spiral--d-turns* num-turns)
         (*spiral--n-spans* (- num-points 1))
         (*spiral--rdfun* *spiral--rdfun*)
         (do-fit nil)
         (in-tol nil)
         (dmax 0d0)
         (crv nil)
         (deg 3)
         (dp0 (make-point))
         (dp1 (make-point))
         (c (make-point))
         (d (make-point))
         (pt (make-point)))
    (declare (type fixnum n-spans-max deg) (type double-float dmax) (type point dp0 dp1 c d pt))
    (when (<= *spiral--n-spans* 0)
      (setf do-fit t
            *spiral--n-spans* (truncate (+ (* 8.0d0 *spiral--d-turns*) 0.99d0))))
    (do () (in-tol)
      (let* ((n-spans *spiral--n-spans*)
             (d-turns *spiral--d-turns*)
             (pnts (make-point-array (1+ n-spans)))
             (d-theta (/ (* +spiral--twopi+ d-turns) (float n-spans 1d0)))
             (angle 0d0))
        (declare (type fixnum n-spans) (type double-float d-turns d-theta angle)
                 (type point-array pnts))
        (setf (aref pnts 0 0) *spiral--r0*
              (aref pnts 0 1) 0d0
              (aref pnts 0 2) 0d0)
        (do ((i 1 (1+ i))) ((> i n-spans))
          (declare (type fixnum i))
          (setf angle (+ angle d-theta))
          (spiral--st-get-point-at-angle angle pt)
          (%pt-row pt pnts i))
        (setf *spiral--rdfun* (/ (+ (- *spiral--r0*) *spiral--r1*) (* +spiral--twopi+ d-turns)))
        (spiral--st-get-point-at-angle-d 0d0 dp0)
        (spiral--st-get-point-at-angle-d (* +spiral--twopi+ d-turns) dp1)
        ;; the analytic tangents scaled uniformly to the fit's [0,1]
        (dotimes (k 3)
          (setf (aref dp0 k) (* (aref dp0 k) (* +spiral--twopi+ d-turns))))
        (dotimes (k 3)
          (setf (aref dp1 k) (* (aref dp1 k) (* +spiral--twopi+ d-turns))))
        ;; cubic interpolation at the knots
        (setf crv (fit-cubic-spline-interp pnts n-spans dp0 dp1 +uniform+))
        ;; goodness of fit: the radial error at each span's midpoint
        (let ((uu (curve-knots crv))
              (h *spiral--h*))
          (declare (type real-vector uu) (type double-float h))
          (setf dmax 0d0)
          (dotimes (i n-spans)
            (let ((u (* 0.5d0 (+ (aref uu (+ deg i)) (aref uu (+ deg i 1))))))
              (declare (type double-float u))
              (crv-eval crv u +left+ c)
              (let ((theta (/ (* (* (aref c 2) +spiral--twopi+) d-turns) h)))
                (declare (type double-float theta))
                (spiral--st-get-point-at-angle theta d)
                (let* ((dist-cxy (c-sqrt (+ (* (aref c 0) (aref c 0)) (* (aref c 1) (aref c 1)))))
                       (dist-dxy (c-sqrt (+ (* (aref d 0) (aref d 0)) (* (aref d 1) (aref d 1)))))
                       (dist (abs (- dist-cxy dist-dxy))))
                  (declare (type double-float dist-cxy dist-dxy dist))
                  (when (> dist dmax)
                    (setf dmax dist)))))))
        (setf in-tol (if do-fit (< dmax tol) t))
        ;; not in tolerance: more spans and repeat
        (unless in-tol
          (let ((n-spans0 n-spans)
                (factor (/ dmax tol))
                (delta-spans 0))
            (declare (type fixnum n-spans0 delta-spans) (type double-float factor))
            (when (> factor 1.0d0)
              (setf factor 1.0d0))
            (setf delta-spans (truncate (* (float n-spans 1d0) factor)))
            (when (< delta-spans 1)
              (setf delta-spans 1))
            (setf *spiral--n-spans* (+ n-spans delta-spans))
            (if (< n-spans0 n-spans-max)
                (when (> *spiral--n-spans* n-spans-max)
                  (setf *spiral--n-spans* n-spans-max))
                ;; the maximum number of spans is reached
                (setf in-tol t))))))
    (values (1+ *spiral--n-spans*) dmax crv)))
