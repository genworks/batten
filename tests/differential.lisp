;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests: each ported routine against NLib's C on random input.
;;;; The bar is bit-identical results, and the same inputs failing on both
;;;; sides.

(defpackage :batten-tests
  (:use :common-lisp :batten :batten-oracle :batten-random)
  (:export #:run-all))

(in-package :batten-tests)

(defstruct tally (name "") (trials 0) (identical 0) (both-failed 0) (mismatch 0)
  (max-rel 0d0) (first-mismatch nil))

(defun rel-diff (a b)
  (let ((d (abs (- a b))) (s (max (abs a) (abs b))))
    (if (zerop s) 0d0 (/ d s))))

(defun compare-arrays (a b)
  "(values identical-p max-relative-difference)"
  (let ((same t) (worst 0d0))
    (dotimes (i (array-total-size a))
      (let ((x (row-major-aref a i)) (y (row-major-aref b i)))
        (unless (eql x y)
          (setf same nil worst (max worst (rel-diff x y))))))
    (values same worst)))

(defun record (tally ok-cl ok-c same rel case)
  (incf (tally-trials tally))
  (cond ((and (not ok-cl) (not ok-c)) (incf (tally-both-failed tally)))
        ((and ok-cl ok-c same) (incf (tally-identical tally)))
        (t (incf (tally-mismatch tally))
           (setf (tally-max-rel tally) (max (tally-max-rel tally) (or rel 0d0)))
           (unless (tally-first-mismatch tally)
             (setf (tally-first-mismatch tally) case)))))

(defmacro trial (tally case cl-form c-form compare)
  "Run CL-FORM (signals nurbs-error on failure) and C-FORM (answers
(values result error-p)); record agreement in TALLY."
  `(let* ((cl-result nil) (ok-cl t))
     (handler-case (setf cl-result ,cl-form)
       (nurbs-error () (setf ok-cl nil)))
     (multiple-value-bind (c-result c-error) ,c-form
       (let ((ok-c (not c-error)))
         (if (and ok-cl ok-c)
             (multiple-value-bind (same rel) (funcall ,compare cl-result c-result)
               (record ,tally ok-cl ok-c same rel ,case))
             (record ,tally ok-cl ok-c nil nil ,case))))))

(defun test-find-span (trials)
  (let ((tally (make-tally :name "find-span (N_BasisFindSpan)")))
    (dotimes (i trials tally)
      (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
             (u (random-parameter k)) (side (if (zerop (rnd 2)) :left :right)))
        (trial tally (list :knots k :p p :u u :side side)
               (find-span k p u side)
               (c-find-span k p u side)
               (lambda (a b) (values (= a b) (if (= a b) 0d0 1d0))))))))

(defun test-basis-functions (trials)
  (let ((tally (make-tally :name "basis-functions! (N_BasisEval)")))
    (dotimes (i trials tally)
      (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
             (u (random-parameter k)) (side (if (zerop (rnd 2)) :left :right)))
        (trial tally (list :knots k :p p :u u :side side)
               (let ((n (make-array (1+ p) :element-type 'double-float)))
                 (cons (basis-functions! k p u side n) n))
               (multiple-value-bind (n spn err) (c-basis-functions k p u side)
                 (values (cons spn n) err))
               (lambda (a b)
                 (multiple-value-bind (same rel) (compare-arrays (cdr a) (cdr b))
                   (values (and same (= (car a) (car b))) rel))))))))

(defun test-basis-derivatives (trials)
  (let ((tally (make-tally :name "basis-derivatives! (N_BasisDerivs)")))
    (dotimes (i trials tally)
      (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
             (u (random-parameter k)) (side (if (zerop (rnd 2)) :left :right))
             (der (rnd (1+ +max-derivative+))))
        (trial tally (list :knots k :p p :u u :side side :der der)
               (let ((nd (make-array (list (1+ der) (1+ p)) :element-type 'double-float)))
                 (cons (basis-derivatives! k p u side der nd) nd))
               (multiple-value-bind (nd spn err) (c-basis-derivatives k p u side der)
                 (values (cons spn nd) err))
               (lambda (a b)
                 (multiple-value-bind (same rel) (compare-arrays (cdr a) (cdr b))
                   (values (and same (= (car a) (car b))) rel))))))))

(defun test-curve-point (trials)
  (let ((tally (make-tally :name "curve-point! (N_CrvEval)")))
    (dotimes (i trials tally)
      (let* ((c (random-curve)) (u (random-parameter (curve-knots c)))
             (side (if (zerop (rnd 2)) :left :right)))
        (trial tally (list :curve c :u u :side side)
               (curve-point c u :side side)
               (c-curve-point c u side)
               #'compare-arrays)))))

(defun test-curve-derivatives (trials)
  (let ((tally (make-tally :name "curve-derivatives! (N_CrvDerivs)")))
    (dotimes (i trials tally)
      (let* ((c (random-curve)) (u (random-parameter (curve-knots c)))
             (side (if (zerop (rnd 2)) :left :right))
             (der (rnd (1+ +max-derivative+))))
        (trial tally (list :curve c :u u :side side :der der)
               (curve-derivatives c u der :side side)
               (c-curve-derivatives c u der side)
               #'compare-arrays)))))

(defun report (tally)
  (format t "~&~38a ~7d trials  ~7d identical  ~5d both-failed  ~5d MISMATCH~@[  max rel ~,3e~]~%"
          (tally-name tally) (tally-trials tally) (tally-identical tally)
          (tally-both-failed tally) (tally-mismatch tally)
          (when (plusp (tally-mismatch tally)) (tally-max-rel tally)))
  (zerop (tally-mismatch tally)))

(defun run-all (&key (trials 20000) (seed 20261025))
  "Run every differential test; true when all agree with NLib to the bit."
  (reseed seed)
  (let ((tallies (list (test-find-span trials)
                       (test-basis-functions trials)
                       (test-basis-derivatives trials)
                       (test-curve-point trials)
                       (test-curve-derivatives trials))))
    (values (every #'identity (mapcar #'report tallies)) tallies)))
