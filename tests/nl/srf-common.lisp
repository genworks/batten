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
