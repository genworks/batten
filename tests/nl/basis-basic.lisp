;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/basis-basic.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test basis-find-span (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)))
    (agree (list k p u flg)
           (list (batten::basis-find-span k p u flg))
           (with-foreign ((kp k :knots))
             (with-scalars ((spn :long))
               (values (c-flag-error-p (c-basis-find-span kp p u flg spn))
                       (list (cffi:mem-ref spn :long))))))))

(define-oracle-test basis-eval (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)))
    (agree (list k p u flg)
           (let ((n (make-real-vector (1+ p))))
             (list (batten::basis-eval k p u flg n) n))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((n :real-vector (1+ p)))
               (with-scalars ((spn :long))
                 (values (c-flag-error-p (c-basis-eval kp p u flg n spn))
                         (list (cffi:mem-ref spn :long) (from-foreign n :real-vector (1+ p))))))))))

(define-oracle-test basis-derivs (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)) (der (rnd (1+ +maxder+))))
    (agree (list k p u flg der)
           (let ((nd (make-real-matrix (1+ der) (1+ p))))
             (list (batten::basis-derivs k p u flg der nd) nd))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((nd :real-matrix (1+ der) (1+ p)))
               (with-scalars ((spn :long))
                 (values (c-flag-error-p (c-basis-derivs kp p u flg der nd spn))
                         (list (cffi:mem-ref spn :long) (from-foreign nd :real-matrix (1+ der) (1+ p))))))))))

(define-oracle-test basis-get-span-count (:file "basis-basic" :trials 2000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)))
    (agree (list k p)
           (list (batten::basis-get-span-count k p))
           (with-foreign ((kp k :knots))
             (with-scalars ((nsp :long))
               (c-basis-get-span-count kp p nsp)
               (values nil (list (cffi:mem-ref nsp :long))))))))

(define-oracle-test knot-vector-is-param-out-of-bounds (:file "basis-basic" :trials 2000)
  (let* ((k (curve-knots (random-curve)))
         (span (- (aref k (1- (length k))) (aref k 0)))
         (u (+ (aref k 0) (* span (rnd-range -0.2d0 1.2d0)))))
    (agree (list k u)
           (progn (batten::knot-vector-is-param-out-of-bounds k u 'test) (list t))
           (with-foreign ((kp k :knots))
             (values (c-flag-error-p (c-knot-vector-is-param-out-of-bounds kp u (cffi:null-pointer)))
                     (list t))))))

;;; ------------------------------------------------------------- L0-g2

(defun %spoil-knots (k p)
  "A copy of knot vector K, now and then made invalid for degree P: an end
multiplicity changed, a knot moved out of order, or an interior knot
repeated past P."
  (let* ((k (copy-seq k)) (m (1- (length k))))
    (case (rnd 10)
      (0 (setf (aref k p) (aref k (1+ p))))                       ; start multiplicity off
      (1 (setf (aref k (- m p)) (aref k (- m p 1))))              ; end multiplicity off
      (2 (when (> (- m p) (1+ p))                                 ; out of order
           (let ((i (+ p 1 (rnd (- m p p 1)))))
             (setf (aref k i) (+ (aref k i) (* 0.5d0 (- (aref k m) (aref k 0))))))))
      (3 (when (> (- m p) (1+ p))                                 ; interior multiplicity p+1
           (let ((i (+ p 1 (rnd (- m p p 1)))))
             (do ((j i (1+ j)) (c 0 (1+ c))) ((or (> c p) (>= j (- m p))))
               (setf (aref k j) (aref k i)))))))
    k))

(define-oracle-test knot-vector-is-valid (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (p0 (curve-degree c))
         (k (%spoil-knots (curve-knots c) p0))
         (p (case (rnd 16) (0 (1- p0)) (1 (1+ p0)) (t p0))))
    (agree (list k p)
           (progn (batten::knot-vector-is-valid k p 'test) (list t))
           (with-foreign ((kp k :knots))
             (values (c-flag-error-p (c-knot-vector-is-valid kp p (cffi:null-pointer)))
                     (list t))))))

(define-oracle-test knot-vector-from-real-array (:file "basis-basic" :trials 1000)
  (let* ((len (1+ (rnd 20))) (u (random-real-vector len)) (m (rnd len)))
    (agree (list u m)
           (list (batten::knot-vector-from-real-array u m))
           (with-foreign ((up u :real-vector))
             (with-foreign-out ((kp :knots))
               (c-knot-vector-from-real-array kp up m)
               (values nil (list (from-foreign kp :knots))))))))

(define-oracle-test knots-check (:file "basis-basic" :trials 2000)
  (let* ((fresh (zerop (rnd 4)))
         (k (if fresh nil (random-real-vector (1+ (rnd 20)))))
         (m (if fresh (rnd 20) (rnd (+ (length k) 2)))))
    ;; NLib's allocator leaves a fresh knot array uninitialized: compare only
    ;; its size there
    (agree (list k m)
           (let ((kv (batten::knots-check k m 'test)))
             (list (if fresh (length kv) kv)))
           (with-stacks (s)
             (cffi:with-foreign-object (kpp :pointer)
               (with-foreign ((kp (or k (make-real-vector 1)) :knots))
                 (setf (cffi:mem-ref kpp :pointer) (if fresh (cffi:null-pointer) kp))
                 (let ((err (c-flag-error-p (c-knots-check kpp m (cffi:null-pointer) s))))
                   (values err (unless err
                                 (let ((kv (from-foreign (cffi:mem-ref kpp :pointer) :knots)))
                                   (list (if fresh (length kv) kv))))))))))))

(defun %flag () (if (zerop (rnd 40)) 3 (random-side)))

(define-oracle-test basis-find-span-and-mult (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (%flag)))
    (agree (list k p u flg)
           (multiple-value-list (batten::basis-find-span-and-mult k p u flg))
           (with-foreign ((kp k :knots))
             (with-scalars ((spn :long) (mlt :long))
               (values (c-flag-error-p (c-basis-find-span-and-mult kp p u flg spn mlt))
                       (list (cffi:mem-ref spn :long) (cffi:mem-ref mlt :long))))))))

(define-oracle-test basis-find-index-node-array (:file "basis-basic" :trials 3000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (len (- (length k) 1 p)))
    (agree (list k p)
           (list (batten::basis-find-index-node-array k p (make-real-vector len)))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((tt :real-vector len))
               (c-basis-find-index-node-array kp p tt)
               (values nil (list (from-foreign tt :real-vector len))))))))

(defun %index-and-parameter (k p)
  "A basis function index (now and then out of range) and a parameter (now and
then outside the domain) for knot vector K of degree P."
  (let* ((n (- (length k) p 2))
         (i (case (rnd 30) (0 -1) (1 (1+ n)) (t (rnd (1+ n)))))
         (u (if (zerop (rnd 30))
                (+ (aref k (1- (length k))) 0.25d0)
                (random-parameter k))))
    (values i u)))

(define-oracle-test basis-i-eval (:file "basis-basic" :trials 10000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)) (flg (%flag)))
    (multiple-value-bind (i u) (%index-and-parameter k p)
      (agree (list k i p u flg)
             (list (batten::basis-i-eval k i p u flg))
             (with-foreign ((kp k :knots))
               (with-scalars ((n :double))
                 (values (c-flag-error-p (c-basis-i-eval kp i p u flg n))
                         (list (cffi:mem-ref n :double)))))))))

(define-oracle-test basis-i-derivs (:file "basis-basic" :trials 10000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)) (flg (%flag))
         (der (rnd (+ p 3))))
    (multiple-value-bind (i u) (%index-and-parameter k p)
      (agree (list k i p u flg der)
             (list (batten::basis-i-derivs k i p u flg der (make-real-vector (1+ der))))
             (with-foreign ((kp k :knots))
               (with-foreign-out ((nd :real-vector (1+ der)))
                 (values (c-flag-error-p (c-basis-i-derivs kp i p u flg der nd))
                         (list (from-foreign nd :real-vector (1+ der))))))))))
