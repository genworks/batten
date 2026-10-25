;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-common.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test create-srf-corner-pts (:file "srf-common" :trials 1000)
  (let* ((pts (random-points 4))
         (corners (let (cs)
                    (dotimes (i 4 (nreverse cs))
                      (push (make-point (aref pts i 0) (aref pts i 1) (aref pts i 2)) cs)))))
    (destructuring-bind (p00 p10 p01 p11) corners
      (agree (list :corners corners)
             (list (batten::create-srf-corner-pts p00 p10 p01 p11))
             (with-foreign ((a p00 :point) (b p10 :point) (c p01 :point) (d p11 :point))
               (with-stacks (sg)
                 (with-foreign-out ((sur :surface))
                   (let ((failed (c-flag-error-p (c-create-srf-corner-pts a b c d sur sg))))
                     (values failed (unless failed (list (from-foreign sur :surface))))))))))))

;;; ------------------------------------------------------------- L2-g3

(define-oracle-test create-srf-extrude-crv (:file "srf-common" :trials 3000)
  (let* ((c (random-curve))
         (w (case (rnd 20)
              (0 (make-point 0d0 0d0 0d0))
              (1 (make-point 0d0 0d0 (rnd-range 0.1d0 3d0)))
              (t (random-point 5d0))))
         (d (case (rnd 10) (0 0d0) (1 (- (rnd 10d0))) (t (rnd 50d0))))
         (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+)))))
    (agree (list :curve c :w w :d d :dir dir)
           (let ((cc (%make-curve (curve-degree c) (copy-seq (curve-knots c))
                                  (let* ((a (curve-cpoints c))
                                         (b (make-array (array-dimensions a) :element-type 'double-float)))
                                    (dotimes (i (array-total-size a) b)
                                      (setf (row-major-aref b i) (row-major-aref a i))))))
                 (ww (copy-seq w)))
             (let ((s (batten::create-srf-extrude-crv cc ww d dir)))
               ;; W is the C's by-value argument: it must come back unchanged
               (list s cc ww)))
           (with-foreign ((cur c :curve) (wp w :point))
             (with-stacks (sg)
               (with-foreign-out ((sur :surface))
                 (let ((failed (c-flag-error-p (c-create-srf-extrude-crv cur wp d dir sur sg))))
                   (values failed (unless failed
                                    (list (from-foreign sur :surface) (from-foreign cur :curve)
                                          (from-foreign wp :point)))))))))))

;;; ------------------------------------------------------------- L4-g1

(defun %sc-copy-curve (c)
  "A fresh copy of curve C."
  (%make-curve (curve-degree c) (copy-seq (curve-knots c))
               (let* ((a (curve-cpoints c))
                      (b (make-array (array-dimensions a) :element-type 'double-float)))
                 (dotimes (i (array-total-size a) b)
                   (setf (row-major-aref b i) (row-major-aref a i))))))

(define-oracle-test create-ruled-srf-from-boundary-crvs (:file "srf-common" :trials 2000)
  (let* ((two-d (zerop (rnd 5)))
         (rat (rnd 4))                    ; 0 neither, 1 P, 2 Q, 3 both rational
         (cp (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (oddp rat)))
         (cq (case (rnd 4)
               ;; the same knots as P
               (0 (let ((c (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (>= rat 2))))
                    (if (and (= (curve-degree c) (curve-degree cp)))
                        c
                        (let ((d (%sc-copy-curve cp)))
                          (dotimes (i (1+ (curve-n d)) d)
                            (let ((w (if (>= rat 2) (+ 0.2d0 (rnd 4.8d0)) 1d0)))
                              (dotimes (k 3)
                                (unless (= (aref (curve-cpoints d) i k) +no-z+)
                                  (setf (aref (curve-cpoints d) i k) (* w (- (rnd 200d0) 100d0)))))
                              (setf (aref (curve-cpoints d) i 3) (if (>= rat 2) w +no-w+))))))))
               (t (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (>= rat 2)))))
         (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+)))))
    (agree (list :p cp :q cq :dir dir)
           (let ((lp (%sc-copy-curve cp)) (lq (%sc-copy-curve cq)))
             (list (batten::create-ruled-srf-from-boundary-crvs lp lq dir) lp lq))
           (with-foreign ((fp cp :curve) (fq cq :curve))
             (with-stacks (sg)
               (with-foreign-out ((sur :surface))
                 (let ((failed (c-flag-error-p (c-create-ruled-srf-from-boundary-crvs fp fq dir sur sg))))
                   (values failed (unless failed
                                    (list (from-foreign sur :surface) (from-foreign fp :curve)
                                          (from-foreign fq :curve)))))))))))
