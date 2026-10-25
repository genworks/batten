;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-basic.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test is-crv-rat (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (list (if (batten::is-crv-rat c) 1 0))
           (with-foreign ((cp c :curve))
             (values nil (list (c-is-crv-rat cp)))))))

(define-oracle-test crv-get-param-bounds (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (multiple-value-list (batten::crv-get-param-bounds c))
           (with-foreign ((cp c :curve))
             (with-scalars ((ul :double) (ur :double))
               (c-crv-get-param-bounds cp ul ur)
               (values nil (list (cffi:mem-ref ul :double) (cffi:mem-ref ur :double))))))))

(defun %random-cpoints-like (n+1 two-d rational)
  "A random cpoint-array of N+1 rows, 2D and rational as asked."
  (let ((a (make-cpoint-array n+1)))
    (dotimes (i n+1 a)
      (let ((w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
        (setf (aref a i 0) (* w (- (rnd 200d0) 100d0))
              (aref a i 1) (* w (- (rnd 200d0) 100d0))
              (aref a i 2) (if two-d +no-z+ (* w (- (rnd 200d0) 100d0)))
              (aref a i 3) (if rational w +no-w+))))))

(defun %degenerate-curve (&key (spread 0d0))
  "A random curve whose control points all sit near one point (within
SPREAD in each coordinate), rational or not, 2D or not."
  (let* ((c (random-curve))
         (pw (curve-cpoints c))
         (cx (- (rnd 200d0) 100d0)) (cy (- (rnd 200d0) 100d0)) (cz (- (rnd 200d0) 100d0)))
    (dotimes (i (array-dimension pw 0) c)
      (let ((w (if (= (aref pw i 3) +no-w+) 1d0 (aref pw i 3))))
        (setf (aref pw i 0) (* w (+ cx (- (rnd (* 2 spread)) spread)))
              (aref pw i 1) (* w (+ cy (- (rnd (* 2 spread)) spread))))
        (unless (= (aref pw i 2) +no-z+)
          (setf (aref pw i 2) (* w (+ cz (- (rnd (* 2 spread)) spread)))))))))

(defun %closed-curve ()
  "A random curve whose last control point is its first, give or take a
distance around *mtol*."
  (let* ((c (random-curve))
         (pw (curve-cpoints c))
         (n (1- (array-dimension pw 0)))
         (eps (case (rnd 4) (0 0d0) (1 5d-8) (2 1d-7) (t 2d-7))))
    (let* ((w0 (if (= (aref pw 0 3) +no-w+) 1d0 (aref pw 0 3)))
           (wn (if (= (aref pw n 3) +no-w+) 1d0 (aref pw n 3))))
      (setf (aref pw n 0) (* wn (+ (/ (aref pw 0 0) w0) eps))
            (aref pw n 1) (* wn (/ (aref pw 0 1) w0)))
      (unless (= (aref pw n 2) +no-z+)
        (setf (aref pw n 2) (* wn (/ (aref pw 0 2) w0)))))
    c))

(define-oracle-test crv-is-3d (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (list (if (batten::crv-is-3d c) 1 0))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-3d cp)))))))

(define-oracle-test crv-from-c-polygon-and-knot-vector (:file "crv-basic" :trials 500)
  (let ((c (random-curve)))
    (agree c
           (list (batten::crv-from-c-polygon-and-knot-vector
                  (%make-curve 0 (make-real-vector 0) (make-cpoint-array 0))
                  (curve-cpoints c) (curve-degree c) (curve-knots c)))
           (with-foreign ((pol (curve-cpoints c) :cpolygon) (knt (curve-knots c) :knots))
             (with-foreign-out ((cur :curve))
               (c-crv-from-c-polygon-and-knot-vector cur pol (curve-degree c) knt)
               (values nil (list (from-foreign cur :curve))))))))

(define-oracle-test crv-get-conic-data (:file "crv-basic" :trials 4000)
  (let* ((ctp (case (rnd 3) (0 batten::+quadratic+) (1 batten::+quartic+) (t batten::+quintic+)))
         (as (case (rnd 4)
               (0 (float (- (rnd 1440) 720) 1d0))
               (t (- (rnd 1440d0) 720d0))))
         (ae (case (rnd 6)
               (0 as)
               (1 (+ as 360d0))
               (2 (+ as (* 0.5d-6 (- (rnd 4d0) 2d0))))
               (3 (+ as 360d0 (* 1d-6 (- (rnd 4d0) 2d0))))
               (4 (+ as (float (* 30 (- (rnd 25) 12)) 1d0)))
               (t (- (rnd 2880d0) 1440d0)))))
    (agree (list :ctp ctp :as as :ae ae)
           (multiple-value-list (batten::crv-get-conic-data ctp as ae))
           (with-scalars ((aep :double) (n :long) (p :short) (m :long))
             (setf (cffi:mem-ref aep :double) ae)
             (c-crv-get-conic-data ctp as aep n p m)
             (values nil (list (cffi:mem-ref aep :double) (cffi:mem-ref n :long)
                               (cffi:mem-ref p :short) (cffi:mem-ref m :long)))))))

(define-oracle-test crv-set-c-pts-and-knots (:file "crv-basic" :trials 500)
  (let* ((c (random-curve))
         (p (curve-degree c))
         (n+1 (array-dimension (curve-cpoints c) 0))
         (pw (%random-cpoints-like n+1 (zerop (rnd 5)) (zerop (rnd 2))))
         (u (random-knots p n+1 :lo (- (rnd 4d0) 2d0) :hi (+ 3d0 (rnd 4d0)))))
    (agree (list :curve c :pw pw :u u)
           (list (batten::crv-set-c-pts-and-knots (batten::crv-copy c) pw u))
           (with-foreign ((cp c :curve) (fpw pw :cpoint-array) (fu u :real-vector))
             (c-crv-set-c-pts-and-knots cp fpw fu)
             (values nil (list (from-foreign cp :curve)))))))

(define-oracle-test crv-is-closed (:file "crv-basic" :trials 3000)
  (let ((c (if (zerop (rnd 2)) (random-curve) (%closed-curve))))
    (agree c
           (list (if (batten::crv-is-closed c) 1 0))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-closed cp)))))))

(define-oracle-test crv-is-degen (:file "crv-basic" :trials 3000)
  (let ((c (case (rnd 4)
             (0 (random-curve))
             (1 (%degenerate-curve))
             (2 (%degenerate-curve :spread 4d-8))
             (t (%degenerate-curve :spread 1d-7)))))
    (agree c
           (list (if (batten::crv-is-degen c) 1 0))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-degen cp)))))))

(define-oracle-test crv-copy (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (list (batten::crv-copy c))
           (with-stacks (s)
             (with-foreign ((cp c :curve))
               (with-foreign-out ((cq :curve))
                 (values (c-flag-error-p (c-crv-copy cp cq s))
                         (list (from-foreign cq :curve)))))))))

(define-oracle-test create-crv-from-num-and-denom (:file "crv-basic" :trials 2000)
  (let* ((num (random-curve :rational nil))
         (den (case (rnd 4)
                (0 nil)
                (1 (random-cfun))
                (t (%make-cfun (curve-degree num) (copy-seq (curve-knots num))
                               (random-real-vector (array-dimension (curve-cpoints num) 0) 5d0))))))
    (agree (list :num num :den den)
           (list (batten::create-crv-from-num-and-denom num den))
           (with-stacks (s)
             (with-foreign ((np num :curve))
               (with-foreign-out ((cur :curve))
                 (let ((rc (if den
                               (with-foreign ((dp den :cfun))
                                 (c-create-crv-from-num-and-denom np dp cur s))
                               (c-create-crv-from-num-and-denom np (cffi:null-pointer) cur s))))
                   (values (c-flag-error-p rc)
                           (unless (c-flag-error-p rc) (list (from-foreign cur :curve)))))))))))
