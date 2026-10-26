;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-offset.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L9-g0

(defun %srf-offset-profile (plane &key line)
  "A random curve in the coordinate plane PLANE (:xz with x in [5, 45], off
the z axis, or :xy with z = 0): a line segment when LINE, else of degree 1 to
3, rational now and then."
  (let* ((p (if line 1 (1+ (rnd 3))))
         (n+1 (if line 2 (+ p 1 (rnd 4))))
         (rational (and (not line) (zerop (rnd 4))))
         (pw (make-cpoint-array n+1)))
    (multiple-value-bind (lo hi) (batten-random::random-domain)
      (dotimes (i n+1)
        (let ((w (if rational (+ 0.5d0 (rnd 2d0)) 1d0))
              (a (rnd-range 5d0 45d0)) (b (rnd-range -40d0 40d0)))
          (if (eq plane :xz)
              (setf (aref pw i 0) (* w a) (aref pw i 1) 0d0 (aref pw i 2) (* w b))
              (setf (aref pw i 0) (* w (- a 25d0)) (aref pw i 1) (* w b) (aref pw i 2) 0d0))
          (setf (aref pw i 3) (if rational w +no-w+))))
      (%make-curve p (random-knots p n+1 :lo lo :hi hi) pw))))

(defun %srf-offset-meridian ()
  "A profile in the xz-plane to revolve about the z axis: a general curve, a
line (a cone), a line parallel to the axis (a cylinder), or a circular arc
centred on the axis (a sphere) or off it (a torus)."
  (case (rnd 6)
    ((0 1) (%srf-offset-profile :xz))
    (2 (%srf-offset-profile :xz :line t))
    (3 (let ((c (%srf-offset-profile :xz :line t)))
         (setf (aref (curve-cpoints c) 1 0) (aref (curve-cpoints c) 0 0))
         c))
    (t (let* ((sphere (zerop (rnd 2)))
              (r (rnd-range 5d0 20d0))
              (cx (if sphere 0d0 (+ r (rnd-range 2d0 20d0))))
              (as (if sphere (rnd-range -85d0 -20d0) (rnd-range -170d0 0d0)))
              (ae (if sphere (rnd-range 20d0 85d0) (rnd-range 20d0 170d0))))
         (batten::with-c-float-environment
           (batten::create-circ-arc (make-point cx 0d0 (rnd-range -10d0 10d0))
                                    (make-point 1d0 0d0 0d0) (make-point 0d0 0d0 1d0)
                                    r as ae +quadratic+))))))

(defun %srf-offset-plane ()
  "A bilinear surface on a random plane, its corners a random quadrilateral in
it, now and then with two corners coincident."
  (let* ((o (random-point 50d0))
         (e1 (random-point 1d0)) (e2 (random-point 1d0))
         (a (rnd-range 5d0 40d0)) (b (rnd-range 5d0 40d0))
         (corner (lambda (s tt)
                   (make-point (+ (aref o 0) (* s (aref e1 0)) (* tt (aref e2 0)))
                               (+ (aref o 1) (* s (aref e1 1)) (* tt (aref e2 1)))
                               (+ (aref o 2) (* s (aref e1 2)) (* tt (aref e2 2))))))
         (p00 (funcall corner 0d0 0d0))
         (p10 (funcall corner a (rnd 5d0)))
         (p01 (funcall corner (rnd 5d0) b))
         (p11 (funcall corner (+ a (rnd 5d0)) (+ b (rnd 5d0)))))
    (case (rnd 8)
      (0 (setf p10 (copy-seq p00)))
      (1 (setf p01 (copy-seq p00))))
    (batten::with-c-float-environment (batten::create-srf-corner-pts p00 p10 p01 p11))))

(defun %srf-offset-scaled-surface ()
  "A random free-form surface, its control net scaled down now and then."
  (let* ((s (random-surface :max-degree 3 :max-extra 3 :two-d nil))
         (f (case (rnd 3) (0 1d0) (1 0.1d0) (t 0.03d0)))
         (net (surface-cpoints s)))
    (dotimes (i (array-dimension net 0) s)
      (dotimes (j (array-dimension net 1))
        (dotimes (c 3) (setf (aref net i j c) (* f (aref net i j c))))))))

(defun %srf-offset-input ()
  "A surface for N_SrfOffset: free-form, a plane, a surface of revolution (of
a general planar profile or a line: a cone or cylinder), or an extrusion of a
planar curve, perpendicular to its plane or oblique."
  (batten::with-c-float-environment
    (case (rnd 10)
      ((0 1 2 3) (%srf-offset-scaled-surface))
      (4 (%srf-offset-plane))
      ((5 6) (batten::create-revolved-srf (%srf-offset-meridian)
                                          (make-point 0d0 0d0 0d0) (make-point 0d0 0d0 1d0)
                                          (case (rnd 4) (0 360d0) (1 90d0) (t (rnd-range 20d0 350d0)))
                                          +quadratic+))
      (t (batten::create-srf-extrude-crv (%srf-offset-profile :xy)
                                         (if (zerop (rnd 4)) (make-point 0.3d0 0d0 1d0) (make-point 0d0 0d0 1d0))
                                         (rnd-range 5d0 40d0)
                                         (if (zerop (rnd 2)) +udir+ +vdir+))))))

(define-oracle-test srf-offset (:file "srf-offset" :trials 600)
  (let* ((s (%srf-offset-input))
         (d (* (if (zerop (rnd 2)) 1d0 -1d0) (rnd-range 0.1d0 4d0)))
         (pp (1+ (rnd 3))) (qq (1+ (rnd 3)))
         (par (case (rnd 6) ((0 1 2) +inherited+) (3 +uniform+) (4 +chordlength+) (t +centripetal+)))
         ;; a tiny TOL sends patches through the 1000-sample cap and the
         ;; 10000-point cap; a zero one is an error
         (tol (case (rnd 40) ((0 1 2 3) 1d-4) (4 0d0) (t (case (rnd 3) (0 0.05d0) (1 0.5d0) (t 2d0)))))
         ;; a large EPS makes nearly anything a plane
         (eps (case (rnd 5) (0 1d-6) (1 1d-4) (2 1d-3) (3 1d-2) (t (* 0.1d0 tol)))))
    (agree (list :sur-p s :d d :pp pp :qq qq :par par :eps eps :tol tol)
           (list (batten::srf-offset s d pp qq par eps tol))
           (with-foreign ((sp s :surface))
             (with-stacks (sg)
               (with-foreign-out ((sq :surface))
                 (let ((failed (c-flag-error-p (c-srf-offset sp d pp qq par eps tol sq sg))))
                   (values failed (unless failed (list (from-foreign sq :surface)))))))))))
