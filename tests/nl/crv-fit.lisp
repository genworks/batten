;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-fit.lisp against NLib.

(in-package :batten-tests)

;;; NLib flag constants this file passes, unexported from BATTEN.
(defconstant +crv-fit-epoint+ batten::+epoint+)
(defconstant +crv-fit-hpoint+ batten::+hpoint+)
(defconstant +crv-fit-start+ batten::+start+)
(defconstant +crv-fit-end+ batten::+end+)

(defun %crv-fit-params (count)
  "COUNT ascending parameters: random in a random domain, now and then with
runs of repeated values."
  (let* ((lo (if (zerop (rnd 3)) 0d0 (- (rnd 20d0) 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0)))
         (v (sort (let (xs) (dotimes (i count xs) (push (rnd-range lo hi) xs))) #'<))
         (a (coerce v 'real-vector)))
    (when (zerop (rnd 4))
      (do ((i 1 (1+ i))) ((>= i count))
        (when (zerop (rnd 3)) (setf (aref a i) (aref a (1- i))))))
    (when (and (plusp count) (zerop (rnd 4)))
      (setf (aref a 0) 0d0 (aref a (1- count)) 1d0)
      (setf a (sort a #'<)))
    a))

(defun %crv-fit-knots-in (count)
  "A caller-allocated knot vector of COUNT knots holding random junk, so that
knots a routine does not write compare equal only if left alone on both sides."
  (random-real-vector count 1000d0))

(defmacro %crv-fit-knot-test (case lisp-call c-call u knt)
  "Run a knot-vector routine: LISP-CALL on fresh copies of U and KNT, C-CALL
with UU and KP bound to foreign copies; compare the knot vectors."
  `(let ((lisp-u (copy-seq ,u)) (lisp-knt (copy-seq ,knt)))
     (agree ,case
            (list (let ((u lisp-u) (knt lisp-knt)) (declare (ignorable u knt)) ,lisp-call))
            (with-foreign ((uu ,u :real-vector) (kp ,knt :knots))
              (let ((flag ,c-call))
                (values (c-flag-error-p flag) (list (from-foreign kp :knots))))))))

(define-oracle-test fit-calc-crv-param-values (:file "crv-fit" :trials 4000)
  (let* ((n (if (zerop (rnd 20)) 0 (1+ (rnd 15))))
         (ptp (case (rnd 12) (0 0) (1 3) (t (if (zerop (rnd 2)) +crv-fit-epoint+ +crv-fit-hpoint+))))
         (par (case (rnd 15) (0 0) (1 4) (t (1+ (rnd 3)))))
         (degenerate (zerop (rnd 8)))
         (scale (if (zerop (rnd 6)) 1d-7 100d0))
         (a (if (= ptp +crv-fit-hpoint+)
                (let ((pw (make-cpoint-array (1+ n)))
                      (two-d (zerop (rnd 4)))
                      (rational (case (rnd 3) (0 nil) (1 t) (t :random))))
                  (dotimes (i (1+ n) pw)
                    (let ((c (if (and degenerate (> i 0))
                                 (make-cpoint (aref pw 0 0) (aref pw 0 1) (aref pw 0 2) (aref pw 0 3))
                                 (random-cpoint :rational rational :two-d two-d :scale scale))))
                      (dotimes (j 4) (setf (aref pw i j) (aref c j))))))
                (let ((pts (random-points (1+ n) scale)))
                  (when degenerate
                    (do ((i 1 (1+ i))) ((> i n))
                      (dotimes (j 3) (setf (aref pts i j) (aref pts 0 j)))))
                  pts)))
         (u0 (random-real-vector (1+ n))))
    (agree (list :n n :ptp ptp :par par :a a)
           (list (batten::fit-calc-crv-param-values a n ptp par (copy-seq u0)))
           (with-foreign ((aa a (if (= ptp +crv-fit-hpoint+) :cpoint-array :point-array))
                          (uu u0 :real-vector))
             ;; with n = 0 and +uniform+ the C computes 1.0/0 (and never uses
             ;; it), which SBCL's floating point traps would turn into a signal
             (values (c-flag-error-p (sb-int:with-float-traps-masked (:divide-by-zero)
                                       (c-fit-calc-crv-param-values aa n ptp par uu)))
                     (list (from-foreign uu :real-vector (1+ n))))))))

(define-oracle-test fit-crv-calc-knot-vector (:file "crv-fit" :trials 4000)
  (let* ((n (1+ (rnd 20)))
         (p (1+ (rnd (min 9 (+ n 2)))))
         (u (%crv-fit-params (1+ n)))
         (knt (%crv-fit-knots-in (+ n p 2))))
    (%crv-fit-knot-test (list :n n :p p :u u)
                        (batten::fit-crv-calc-knot-vector u n p knt)
                        (c-fit-crv-calc-knot-vector uu n p kp)
                        u knt)))

(define-oracle-test fit-calc-knot-vector-end-derivs (:file "crv-fit" :trials 4000)
  (let* ((k (1+ (rnd 15)))
         (p (if (zerop (rnd 3)) 1 (1+ (rnd 8))))
         (u (%crv-fit-params (1+ k)))
         (knt (%crv-fit-knots-in (+ k p 4))))
    (%crv-fit-knot-test (list :k k :p p :u u)
                        (batten::fit-calc-knot-vector-end-derivs u k p knt)
                        (c-fit-calc-knot-vector-end-derivs uu k p kp)
                        u knt)))

(define-oracle-test fit-calc-knot-vector-deriv (:file "crv-fit" :trials 4000)
  (let* ((k (1+ (rnd 15)))
         (p (case (rnd 10) (0 (+ k 2 (rnd 3))) ((1 2 3) 1) (t (1+ (rnd (+ k 1))))))
         (whr (case (rnd 15) (0 0) (1 3) (t (if (zerop (rnd 2)) +crv-fit-start+ +crv-fit-end+))))
         (u (%crv-fit-params (1+ k)))
         (knt (%crv-fit-knots-in (+ k p 3))))
    (%crv-fit-knot-test (list :k k :p p :whr whr :u u)
                        (batten::fit-calc-knot-vector-deriv u k p whr knt)
                        (c-fit-calc-knot-vector-deriv uu k p whr kp)
                        u knt)))

(define-oracle-test fit-calc-knot-vector-crv-approx (:file "crv-fit" :trials 4000)
  (let* ((m (1+ (rnd 40)))
         (n (if (zerop (rnd 12)) (+ m 1 (rnd 3)) (rnd (1+ m))))
         (p (1+ (rnd (min 9 (+ n 2)))))
         (u (%crv-fit-params (1+ m)))
         (knt (%crv-fit-knots-in (+ n p 2))))
    (%crv-fit-knot-test (list :m m :n n :p p :u u)
                        (batten::fit-calc-knot-vector-crv-approx u m n p knt)
                        (c-fit-calc-knot-vector-crv-approx uu m n p kp)
                        u knt)))

(define-oracle-test fit-calc-knots-derivs (:file "crv-fit" :trials 4000)
  (let* ((k (1+ (rnd 14)))
         (p (case (rnd 6) (0 1) (1 (+ 4 (rnd 3))) ((2 3) 2) (t 3)))
         (u (%crv-fit-params (1+ k)))
         (knt (%crv-fit-knots-in (+ (* 2 k) p 3))))
    (%crv-fit-knot-test (list :k k :p p :u u)
                        (batten::fit-calc-knots-derivs u k p knt)
                        (c-fit-calc-knots-derivs uu k p kp)
                        u knt)))

(define-oracle-test fit-calc-knots-high-end-derivs (:file "crv-fit" :trials 4000)
  (let* ((n (1+ (rnd 15)))
         (p (1+ (rnd 8)))
         (k (rnd p))
         (l (rnd p))
         (u (%crv-fit-params (1+ n)))
         (knt (%crv-fit-knots-in (+ n k l p 2))))
    (%crv-fit-knot-test (list :n n :p p :k k :l l :u u)
                        (batten::fit-calc-knots-high-end-derivs u n p k l knt)
                        (c-fit-calc-knots-high-end-derivs uu n p k l kp)
                        u knt)))
