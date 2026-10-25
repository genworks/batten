;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-approx.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L6-g1

(defun %crv-approx-curve ()
  "A random curve for N_ApproxNurbsWithNonRatCrv: degree up to 5, a few
spans, rational or not, 2D or not, its control points within 10 of the
origin; now and then a single Bezier segment, a straight line (control
points on one), or a degenerate end (two equal end control points, a zero
end tangent)."
  (let* ((c (random-curve :max-degree 5 :max-extra 5))
         (pw (curve-cpoints c))
         (n+1 (array-dimension pw 0))
         (p (curve-degree c)))
    ;; scale down (homogeneous coordinates scale alike)
    (dotimes (i n+1)
      (dotimes (k 3)
        (unless (and (= k 2) (= (aref pw i 2) +no-z+))
          (setf (aref pw i k) (* 0.1d0 (aref pw i k))))))
    (case (rnd 12)
      (0 ;; a Bezier segment
       (let ((knots (curve-knots c)) (q (make-cpoint-array (1+ p))) (u (make-real-vector (* 2 (1+ p)))))
         (dotimes (i (1+ p))
           (setf (aref u i) (aref knots 0)
                 (aref u (+ p 1 i)) (aref knots (1- (length knots))))
           (dotimes (k 4) (setf (aref q i k) (aref pw i k))))
         (setf c (%make-curve p u q))))
      (1 ;; a straight line
       ;; from the origin to the first control point's Euclidean position
       (let* ((w0 (if (= (aref pw 0 3) +no-w+) 1d0 (aref pw 0 3)))
              (e (list (/ (aref pw 0 0) w0) (/ (aref pw 0 1) w0)
                       (if (= (aref pw 0 2) +no-z+) 0d0 (/ (aref pw 0 2) w0)))))
         (dotimes (i n+1)
           (let ((w (if (= (aref pw i 3) +no-w+) 1d0 (aref pw i 3)))
                 (s (/ (float i 1d0) (float (1- n+1) 1d0))))
             (dotimes (k 3)
               (unless (and (= k 2) (= (aref pw i 2) +no-z+))
                 (setf (aref pw i k) (* w (* s (elt e k))))))))))
      (2 ;; a zero end tangent; not for degree 1, whose first span would
         ;; then be a single point: its coincident samples give chord-length
         ;; parameters, so knots, of multiplicity above the degree, on which
         ;; N_CrvRemoveKnotMaxErr indexes outside its scratch array (Rw[p-s+2])
       (when (> p 1)
         (dotimes (k 4) (setf (aref pw 1 k) (aref pw 0 k))))))
    c))

(define-oracle-test approx-nurbs-with-non-rat-crv (:file "crv-approx" :trials 3000)
  ;; TOL positive (a negative one makes pow() a NaN the C casts to an index),
  ;; or now and then 0, which the division check rejects
  (let* ((c (%crv-approx-curve))
         (tol (case (rnd 50)
                (0 0d0)
                (t (case (rnd 6) (0 1d-3) (1 1d-2) (2 0.1d0) (3 1d0) (4 5d0) (t (rnd-range 1d-2 1d0))))))
         (deg (case (rnd 10) (0 1) (1 2) (2 2) (3 3) (4 3) (5 4) (6 5) (t (1+ (rnd 3)))))
         (par (case (rnd 25)
                (0 0)
                (1 7)
                (t (elt (list +uniform+ +chordlength+ +centripetal+ +inherited+) (rnd 4)))))
         (tn (if (zerop (rnd 2)) +yes+ +no+)))
    (agree (list :curve c :tol tol :deg deg :par par :tan tn)
           (list (batten::approx-nurbs-with-non-rat-crv c tol deg par tn))
           (with-foreign ((cp c :curve))
             (with-stacks (sg)
               (with-foreign-out ((q :curve))
                 (let ((failed (c-flag-error-p (c-approx-nurbs-with-non-rat-crv cp tol deg par tn q sg))))
                   (values failed (unless failed (list (from-foreign q :curve)))))))))))
