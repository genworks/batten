;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Benchmarks of the Lisp routines, portable to any implementation.  The
;;;; curves and parameters come from batten/random's seeded generator, so a
;;;; row means the same work on every implementation.

(defpackage :batten-bench
  (:use :common-lisp :batten :batten-random)
  (:export #:run-bench #:bench-cases #:time-batten #:*cases*
           ;; defined by batten/bench-c when the oracle is loaded
           #:run-bench-c))

(in-package :batten-bench)

(defun seconds-since (t0)
  (/ (- (get-internal-real-time) t0) (float internal-time-units-per-second 1d0)))

(defun bytes-consed ()
  #+sbcl (sb-ext:get-bytes-consed)
  #+ccl (ccl::total-bytes-allocated)
  #-(or sbcl ccl) 0)

(defun curve-of (degree rational n+1)
  "The first random curve of the given shape from the seeded generator."
  (do ((c (random-curve) (random-curve)))
      ((and (= (curve-degree c) degree)
            (eq (and (curve-rational-p c) t) rational)
            (>= (1+ (curve-n c)) n+1)
            (/= (aref (curve-cpoints c) 0 2) +no-z+))
       c)))

(defparameter *cases* '((2 nil) (3 nil) (3 t) (5 nil) (5 t))
  "(degree rational) for each benchmark curve.")

(defun bench-cases (&key (seed 20261025))
  "A list of (degree rational curve params) for *cases*, from SEED."
  (reseed seed)
  (mapcar (lambda (spec)
            (destructuring-bind (p rat) spec
              (let ((c (curve-of p rat 10))
                    (params (make-array 1024 :element-type 'double-float)))
                (dotimes (i 1024) (setf (aref params i) (random-parameter (curve-knots c))))
                (list p rat c params))))
          *cases*))

(defparameter *repeats* 5
  "Each timing is the best of this many runs, to shed interference from other
work on the machine.")

(defun best-of (thunk)
  "The least of *repeats* timings of THUNK, in seconds."
  (let ((best most-positive-double-float))
    (dotimes (i *repeats* best)
      (let ((t0 (get-internal-real-time)))
        (funcall thunk)
        (setf best (min best (seconds-since t0)))))))

(defun time-batten (curve params evals der)
  "(values ns-per-eval bytes-per-eval) for curve-point! (DER nil) or
curve-derivatives! to order DER."
  (declare (type real-vector params) (type fixnum evals))
  (let* ((np (length params))
         (point (make-array 3 :element-type 'double-float))
         (ders (make-array (list (1+ (or der 0)) 3) :element-type 'double-float))
         (b0 (bytes-consed))
         (seconds (best-of
                   (if der
                       (lambda () (dotimes (i evals) (curve-derivatives! curve (aref params (mod i np)) der ders)))
                       (lambda () (dotimes (i evals) (curve-point! curve (aref params (mod i np)) point)))))))
    (values (* 1d9 (/ seconds evals))
            (/ (- (bytes-consed) b0) (float (* evals *repeats*) 1d0)))))

(defun run-bench (&key (evals 1000000) (seed 20261025))
  "Time every case; print a table and return its rows."
  (format t "~&~a ~a~%" (lisp-implementation-type) (lisp-implementation-version))
  (let (rows)
    (dolist (case (bench-cases :seed seed) (nreverse rows))
      (destructuring-bind (p rat c params) case
        (dolist (der '(nil 2))
          (multiple-value-bind (ns bytes) (time-batten c params evals der)
            (push (list :degree p :rational rat :der der :ns ns :bytes bytes) rows)
            (format t "~&degree ~d ~:[polynomial~;rational  ~]  ~12a ~7,1f ns  ~6,1f B/eval~%"
                    p rat (if der "derivs to 2" "point") ns bytes)))))))
