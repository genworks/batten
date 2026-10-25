;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Random curves and parameters for the differential tests and benchmarks.

(defpackage :batten-random
  (:use :common-lisp :batten)
  (:export #:random-curve #:random-parameter #:*seed* #:reseed #:rnd
           #:evaluation-digest))

(in-package :batten-random)

;;; splitmix64, so that every implementation draws the same curves from the
;;; same seed: a failure reproduces anywhere, and benchmarks compare like
;;; with like.
(defvar *seed* 20261025)

(defun reseed (&optional (seed 20261025)) (setf *seed* seed))

(defun next-u64 ()
  (setf *seed* (ldb (byte 64 0) (+ *seed* #x9E3779B97F4A7C15)))
  (let ((z *seed*))
    (setf z (ldb (byte 64 0) (* (logxor z (ash z -30)) #xBF58476D1CE4E5B9))
          z (ldb (byte 64 0) (* (logxor z (ash z -27)) #x94D049BB133111EB)))
    (logxor z (ash z -31))))

(defun rnd (limit)
  "A random integer below an integer LIMIT, or a double in [0, LIMIT)."
  (if (integerp limit)
      (mod (next-u64) limit)
      (* limit (/ (float (ash (next-u64) -11) 1d0) 9007199254740992d0))))

(defun random-knots (p n+1 &key (lo 0d0) (hi 1d0))
  "A clamped knot vector for N+1 control points of degree P: P+1 equal knots at
each end, interior knots random, now and then repeated up to multiplicity P."
  (let* ((interior-count (- n+1 p 1))
         (interior (sort (let (ks)
                           (dotimes (i interior-count ks)
                             (push (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0)))) ks)))
                         #'<)))
    ;; repeat some interior knots, keeping the count
    (let ((v (coerce interior 'vector)))
      (dotimes (i (length v))
        (when (and (> i 0) (zerop (rnd 6)))
          (let ((mult 0))
            (do ((j (1- i) (1- j))) ((or (< j 0) (/= (aref v j) (aref v (1- i)))))
              (incf mult))
            (when (< mult p)
              (setf (aref v i) (aref v (1- i)))))))
      (setf interior (coerce v 'list)))
    (coerce (append (make-list (1+ p) :initial-element lo)
                    interior
                    (make-list (1+ p) :initial-element hi))
            'real-vector)))

(defun random-curve (&key (max-degree 7) (max-extra 12))
  "A random NURBS curve: degree 1..MAX-DEGREE, rational half the time, 2D a
fifth of the time, on a random parameter interval."
  (let* ((p (1+ (rnd max-degree)))
         (n+1 (+ p 1 (rnd max-extra)))
         (lo (if (zerop (rnd 3)) 0d0 (- (rnd 20d0) 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0)))
         (rational (zerop (rnd 2)))
         (two-d (zerop (rnd 5)))
         (cpoints (make-array (* 4 n+1) :element-type 'double-float)))
    (dotimes (i n+1)
      (let ((w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0))
            (b (* 4 i)))
        (setf (aref cpoints b) (* w (- (rnd 200d0) 100d0))
              (aref cpoints (+ b 1)) (* w (- (rnd 200d0) 100d0))
              (aref cpoints (+ b 2)) (if two-d +no-z+ (* w (- (rnd 200d0) 100d0)))
              (aref cpoints (+ b 3)) (if rational w +no-w+))))
    (%make-curve p (random-knots p n+1 :lo lo :hi hi) cpoints)))

(defun random-parameter (knots)
  "A parameter in the knot vector's domain: usually interior, now and then
exactly a knot or an end."
  (let ((lo (aref knots 0))
        (hi (aref knots (1- (length knots)))))
    (case (rnd 8)
      (0 lo)
      (1 hi)
      (2 (aref knots (rnd (length knots))))
      (t (+ lo (rnd (- hi lo)))))))
