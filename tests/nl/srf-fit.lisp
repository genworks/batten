;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-fit.lisp against NLib.

(in-package :batten-tests)

(defun %copy-real-matrix (a)
  (let ((b (make-real-matrix (array-dimension a 0) (array-dimension a 1))))
    (dotimes (i (array-total-size a) b)
      (setf (row-major-aref b i) (row-major-aref a i)))))

(defun %removable-knot (knots deg count)
  "A random knot index r of an interior knot with knots[r] /= knots[r+1], and
its multiplicity s, for a knot vector of degree DEG with control points
0..COUNT; NIL when there is no interior knot."
  (let ((rs nil))
    (do ((r (1+ deg) (1+ r))) ((> r count))
      (when (/= (aref knots r) (aref knots (1+ r))) (push r rs)))
    (when rs
      (let* ((r (nth (rnd (length rs)) rs))
             (s (do ((j r (1- j)) (c 0 (1+ c))) ((/= (aref knots j) (aref knots r)) c))))
        (values r s)))))

(define-oracle-test fit-srf-removal-boundary (:file "srf-fit" :trials 8000)
  (let ((dir (if (zerop (rnd 50)) 3 (if (zerop (rnd 2)) batten::+udir+ batten::+vdir+)))
        sur r s)
    ;; a surface with an interior knot in the direction of removal
    (do () (r)
      (setf sur (random-surface :two-d (zerop (rnd 5))))
      (multiple-value-setq (r s)
        (if (= dir batten::+vdir+)
            (%removable-knot (surface-knots-v sur) (surface-degree-v sur) (surface-m sur))
            (%removable-knot (surface-knots-u sur) (surface-degree-u sur) (surface-n sur)))))
    (let* ((n (surface-n sur)) (m (surface-m sur))
           (other (if (= dir batten::+vdir+) n m))
           (f (rnd (1+ other)))
           (l (if (zerop (rnd 10)) (max 0 (1- f)) (+ f (rnd (- (1+ other) f)))))
           (br (random-real-matrix (1+ n) (1+ m) 1d0))
           (mr (case (rnd 3) (0 batten::+bigd+) (1 -1d0) (t (rnd 100d0)))))
      (agree (list :surface sur :r r :s s :f f :l l :dir dir :mr mr)
             (let ((b (%copy-real-matrix br)))
               (list (batten::fit-srf-removal-boundary sur r s f l dir b mr) b))
             (with-foreign ((sp sur :surface) (bp br :real-matrix))
               (with-scalars ((mrp :double))
                 (setf (cffi:mem-ref mrp :double) mr)
                 (values (c-flag-error-p (c-fit-srf-removal-boundary sp r s f l dir bp mrp))
                         (list (cffi:mem-ref mrp :double)
                               (from-foreign bp :real-matrix (1+ n) (1+ m))))))))))
