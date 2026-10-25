;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/basis-adv.lisp against NLib.

(in-package :batten-tests)

(defun %some-knots ()
  "A random clamped knot vector and its degree, as (values knots p)."
  (let ((c (random-curve)))
    (values (curve-knots c) (curve-degree c))))

(defun %real-3d-to-foreign (dims)
  "An empty, zeroed NL_REAL*** of DIMS (a b c) inside the enclosing
WITH-FOREIGN, as (values pointer)."
  (destructuring-bind (a b c) dims
    (let* ((data (batten-oracle::falloc :double (* a b c)))
           (rows (batten-oracle::falloc :pointer (* a b)))
           (top (batten-oracle::falloc :pointer a)))
      (dotimes (k a)
        (dotimes (i b)
          (setf (cffi:mem-aref rows :pointer (+ (* k b) i))
                (cffi:inc-pointer data (* 8 (+ (* (+ (* k b) i) c))))))
        (setf (cffi:mem-aref top :pointer k) (cffi:inc-pointer rows (* 8 k b))))
      top)))

(defun %real-3d-from-foreign (ptr dims)
  "A fresh (simple-array double-float (a b c)) read from NL_REAL*** PTR."
  (destructuring-bind (a b c) dims
    (let ((out (make-array dims :element-type 'double-float)))
      (dotimes (k a out)
        (let ((plane (cffi:mem-aref ptr :pointer k)))
          (dotimes (i b)
            (let ((row (cffi:mem-aref plane :pointer i)))
              (dotimes (j c)
                (setf (aref out k i j) (cffi:mem-aref row :double j))))))))))

(define-oracle-test basis-get-longest-and-shortest-spans (:file "basis-adv" :trials 3000)
  (multiple-value-bind (knots p) (%some-knots)
    (agree (list :knots knots :p p)
           (multiple-value-list (batten::basis-get-longest-and-shortest-spans knots p))
           (with-foreign ((kp knots :knots))
             (with-scalars ((ds :double) (dl :double))
               (c-basis-get-longest-and-shortest-spans kp p ds dl)
               (values nil (list (cffi:mem-ref ds :double) (cffi:mem-ref dl :double))))))))

(define-oracle-test basis-increase-knot-mult (:file "basis-adv" :trials 3000)
  (multiple-value-bind (knots p) (%some-knots)
    (let ((tm (if (zerop (rnd 4)) (rnd (+ p 1)) (rnd (max 1 (- p 1))))))
      (agree (list :knots knots :p p :t tm)
             (list (batten::basis-increase-knot-mult knots p tm))
             (with-foreign ((kp knots :knots)
                            (kx (make-real-vector (* (length knots) (1+ tm))) :knots))
               (values (c-flag-error-p (c-basis-increase-knot-mult kp p tm kx))
                       (list (from-foreign kx :knots))))))))

(define-oracle-test basis-derivs-array (:file "basis-adv" :trials 3000)
  (multiple-value-bind (knots p) (%some-knots)
    (let* ((n (rnd 9))
           (flg (random-side))
           (der (rnd 7))
           (lo (aref knots 0))
           (hi (aref knots (1- (length knots))))
           (u (make-real-vector (1+ n))))
      ;; ascending parameters; for +left+ the inner ones stay below the end of
      ;; the domain, where the C's span walk would run off the knot vector
      ;; (n values below the end for +left+, one unrestricted: after sorting
      ;; only u[n] can be the end)
      (let ((us (sort (let ((l (list (random-parameter knots))))
                        (dotimes (j n l)
                          (let ((x (random-parameter knots)))
                            (push (if (and (= flg +left+) (>= x hi)) lo x) l))))
                      #'<)))
        (dotimes (j (1+ n)) (setf (aref u j) (pop us))))
      (let ((dims (list (1+ der) (1+ p) (1+ n))))
        (agree (list :knots knots :p p :u u :flg flg :der der)
               (let ((nd (make-array dims :element-type 'double-float :initial-element 0d0))
                     (spn (batten::make-index-vector (1+ n))))
                 (batten::basis-derivs-array knots p u n flg der nd spn)
                 (list nd spn))
               (with-foreign ((kp knots :knots) (up u :real-vector))
                 (let ((nd (%real-3d-to-foreign dims))
                       (spn (batten-oracle::falloc :long (1+ n))))
                   (values (c-flag-error-p (c-basis-derivs-array kp p up n flg der nd spn))
                           (list (%real-3d-from-foreign nd dims)
                                 (from-foreign spn :index-vector (1+ n)))))))))))
