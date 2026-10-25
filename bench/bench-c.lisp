;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The benchmark cases beside NLib's C, through the oracle.  The C times
;;;; include one foreign call apiece.

(in-package :batten-bench)

(defun time-c (curve params evals der)
  (declare (type real-vector params) (type fixnum evals))
  (let ((np (length params)))
    (batten-oracle:with-foreign ((cur curve :curve))
      (cffi:with-foreign-object (out :double (* 3 (1+ (or der 0))))
        (* 1d9 (/ (best-of
                   (if der
                       (lambda () (dotimes (i evals) (batten-oracle::c-crv-derivs cur (aref params (mod i np)) 1 der out)))
                       (lambda () (dotimes (i evals) (batten-oracle::c-crv-eval cur (aref params (mod i np)) 1 out)))))
                  evals))))))

(defun run-bench-c (&key (evals 1000000) (seed 20261025))
  "Every case, Batten beside C; print a table and return its rows."
  (format t "~&~a ~a, against NLib C through CFFI~%" (lisp-implementation-type) (lisp-implementation-version))
  (let (rows)
    (dolist (case (bench-cases :seed seed) (nreverse rows))
      (destructuring-bind (p rat c params) case
        (dolist (der '(nil 2))
          (multiple-value-bind (ns bytes) (time-batten c params evals der)
            (let ((c-ns (time-c c params evals der)))
              (push (list :degree p :rational rat :der der :ns ns :c-ns c-ns :ratio (/ ns c-ns) :bytes bytes) rows)
              (format t "~&degree ~d ~:[polynomial~;rational  ~]  ~12a  Batten ~7,1f ns   C ~7,1f ns   ratio ~5,2f   ~5,1f B/eval~%"
                      p rat (if der "derivs to 2" "point") ns c-ns (/ ns c-ns) bytes))))))))
