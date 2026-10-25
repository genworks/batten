;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_FuncsAdv.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun crv-func-refine (cfn-p knx)
  "Refine curve function CFN-P with the knots KNX (a real-vector, sorted,
strictly inside the domain).  Returns the refined curve function, fresh.
Signals :inp-err for an empty KNX and :par-err when its first or last knot is
on or outside an end of the domain.
NLib: N_CrvFuncRefine (the C's in-place case, cfnP = cfnQ, computes the same
function)."
  (declare (type cfun cfn-p) (type real-vector knx))
  (let* ((fp (cfun-vals cfn-p)) (n (1- (length fp)))
         (p (cfun-degree cfn-p))
         (up (cfun-knots cfn-p)) (m (1- (length up)))
         (x knx) (r (1- (length knx))))
    (declare (type real-vector fp up x) (type fixnum n p m r))
    ;; check input parameters
    (when (< r 0)
      (nurbs-error :inp-err 'crv-func-refine))
    (knot-vector-is-end-param up (aref x 0) 'crv-func-refine)
    (knot-vector-is-end-param up (aref x r) 'crv-func-refine)
    (let ((cfn-q (%make-cfun p (make-real-vector (+ m r 2)) (make-real-vector (+ n r 2)))))
      (multiple-value-bind (fq uq) (crv-func-cntrl-val-knots cfn-q)
        (declare (type real-vector fq uq))
        ;; knot spans
        (let* ((a (basis-find-span up p (aref x 0) +left+))
               (b (1+ (basis-find-span up p (aref x r) +left+)))
               (i 0) (k 0))
          (declare (type fixnum a b i k))
          (fast
            ;; output knot vector ends
            (do ((j 0 (1+ j))) ((> j a))
              (declare (type fixnum j))
              (setf (aref uq j) (aref up j)))
            (do ((j (+ b p) (1+ j))) ((> j m))
              (declare (type fixnum j))
              (setf (aref uq (+ j r 1)) (aref up j)))
            ;; unaltered control values
            (do ((j 0 (1+ j))) ((> j (- a p)))
              (declare (type fixnum j))
              (setf (aref fq j) (aref fp j)))
            (do ((j (1- b) (1+ j))) ((> j n))
              (declare (type fixnum j))
              (setf (aref fq (+ j r 1)) (aref fp j)))
            ;; refine the knot vector
            (setf i (+ b p -1)
                  k (+ b p r))
            (do ((j r (1- j))) ((< j 0))
              (declare (type fixnum j))
              (do () ((not (and (<= (aref x j) (aref up i)) (> i a))))
                (setf (aref fq (- k p 1)) (aref fp (- i p 1)))
                (setf (aref uq k) (aref up i))
                (decf k)
                (decf i))
              (setf (aref fq (- k p 1)) (aref fq (- k p)))
              (do ((l 1 (1+ l))) ((> l p))
                (declare (type fixnum l))
                (let ((tt (+ (- k p) l))
                      (alf (- (aref uq (+ k l)) (aref x j))))
                  (declare (type fixnum tt) (type double-float alf))
                  (if (<= (abs alf) 0d0)
                      (setf (aref fq (1- tt)) (aref fq tt))
                      (let* ((alf (/ alf (- (aref uq (+ k l)) (aref up (+ (- i p) l)))))
                             (oma (- 1d0 alf)))
                        (declare (type double-float alf oma))
                        (setf (aref fq (1- tt))
                              (+ (* oma (aref fq tt)) (* alf (aref fq (1- tt)))))))))
              (setf (aref uq k) (aref x j))
              (decf k)))))
      cfn-q)))

(defun crv-func-deriv (cfn d)
  "The D-th derivative function of curve function CFN, fresh: degree p-D on
the knots with D fewer at each end.  Signals :inp-err unless 0 <= D <= p,
:der-err for a zero-length span the differences divide by and :num-err for a
bad division.  NLib: N_CrvFuncDeriv."
  (declare (type cfun cfn) (type fixnum d))
  (let* ((fu (cfun-vals cfn)) (n (1- (length fu)))
         (p (cfun-degree cfn))
         (u (cfun-knots cfn)) (m (1- (length u))))
    (declare (type real-vector fu u) (type fixnum n p m))
    (when (or (< d 0) (> d p))
      (nurbs-error :inp-err 'crv-func-deriv))
    (let* ((nd (- n d))
           (pd (- p d))
           (md (+ nd pd 1))
           (der (%make-cfun pd (make-real-vector (1+ md)) (make-real-vector (1+ nd)))))
      (declare (type fixnum nd pd md))
      (multiple-value-bind (fd v) (crv-func-cntrl-val-knots der)
        (declare (type real-vector fd v))
        (let ((ft (make-real-vector (1+ n))))
          (declare (type real-vector ft))
          (fast
            (dotimes (i (1+ n))
              (setf (aref ft i) (aref fu i)))
            ;; control values of the derivative function
            (do ((k 1 (1+ k))) ((> k d))
              (declare (type fixnum k))
              (do ((i 0 (1+ i))) ((> i (- n k)))
                (declare (type fixnum i))
                (when (= (aref u (+ i p 1)) (aref u (+ i k)))
                  (nurbs-error :der-err 'crv-func-deriv))
                (let ((num (+ (- (float p 1d0) (float k 1d0)) 1d0))
                      (den (- (aref u (+ i p 1)) (aref u (+ i k)))))
                  (declare (type double-float num den))
                  (when (float-op-is-bad num den +division+)
                    (nurbs-error :num-err 'crv-func-deriv))
                  (let ((alf (/ num den)))
                    (declare (type double-float alf))
                    (setf (aref ft i) (* alf (- (aref ft (1+ i)) (aref ft i))))))))
            (dotimes (i (1+ nd))
              (setf (aref fd i) (aref ft i)))
            ;; the knot vector
            (let ((j -1))
              (declare (type fixnum j))
              (dotimes (i (1+ pd))
                (setf (aref v (incf j)) (aref u 0)))
              (do ((i (1+ p) (1+ i))) ((> i n))
                (declare (type fixnum i))
                (setf (aref v (incf j)) (aref u i)))
              (dotimes (i (1+ pd))
                (setf (aref v (incf j)) (aref u m)))))))
      der)))

(defun crv-func-multiply-constant (alpha cfn)
  "Multiply the control values of curve function CFN by ALPHA, in place.
Returns CFN.  NLib: N_CrvFuncMultiplyConstant."
  (declare (type double-float alpha) (type cfun cfn))
  (multiple-value-bind (n fu) (crv-func-cntrl-val cfn)
    (declare (type fixnum n) (type real-vector fu))
    (fast
      (dotimes (i (1+ n))
        (setf (aref fu i) (* alpha (aref fu i))))))
  cfn)
