;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; A digest of many evaluations, bit for bit: the same seed must give the
;;;; same digest on every implementation, since every routine performs the
;;;; same IEEE double operations in the same order.

(in-package :batten-random)

(defun mix-double (h x)
  "Fold the bits of double X into the 64-bit hash H (FNV-1a over the
significand, exponent and sign).  Zero is canonicalized first: the exponent
INTEGER-DECODE-FLOAT reports for it differs between implementations."
  (multiple-value-bind (sig exp sign)
      (if (zerop x)
          (values 0 0 (float-sign x))
          (integer-decode-float x))
    (dolist (v (list sig (ldb (byte 16 0) exp) (if (minusp sign) 1 0)) h)
      (setf h (ldb (byte 64 0) (* (logxor h (ldb (byte 64 0) v)) #x100000001B3))))))

(defun evaluation-digest (&key (trials 50000) (seed 20261025))
  "Evaluate random curves at random parameters -- points and derivatives to a
random order -- and return a 64-bit digest of every result bit."
  (reseed seed)
  (let ((h #xCBF29CE484222325) (errors 0))
    (dotimes (i trials)
      (let* ((c (random-curve))
             (u (random-parameter (curve-knots c)))
             (der (rnd (1+ +maxder+))))
        (handler-case
            (let ((p (curve-point c u))
                  (d (curve-derivatives c u der)))
              (dotimes (k 3) (setf h (mix-double h (aref p k))))
              (dotimes (k (array-total-size d)) (setf h (mix-double h (row-major-aref d k)))))
          (nurbs-error () (incf errors)))))
    (values (format nil "~16,'0x" h) errors)))
