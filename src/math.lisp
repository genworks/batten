;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Math.cpp (usd-brep source/SMLib/src).  Modified:
;;;; translated to Common Lisp.

(in-package :batten)

(defvar *bypass-float-checks* nil
  "When true, FLOAT-OP-BAD-P answers false without checking, as NLib does
when its NLIBCHKBYPASS is set.  NLib's default, and this one, is to check.")

(declaim (inline nl-sign))
(defun nl-sign (x)
  "1.0 for positive X, else -1.0 -- zero included.  NLib: NL_SIGN."
  (declare (type double-float x))
  (if (> x 0d0) 1d0 -1d0))

(declaim (inline float-op-bad-p))
(defun float-op-bad-p (x y op)
  "True when X OP Y would overflow or underflow a double.  OP is one of
:multiplication, :division, :addition, :subtraction.  NLib: N_FloatOpIsBad."
  (declare (type double-float x y))
  (cond
    (*bypass-float-checks* nil)
    ((eq op :division)
     (cond ((= y 0d0) t)
           ((= x 0d0) nil)
           (t (let ((x (abs x)) (y (abs y)))
                (cond ((and (<= y 1d0) (< x (* +biggest-double+ y))) nil)
                      ((and (>= y 1d0) (> x (* +smallest-double+ y))) nil)
                      (t t))))))
    ((eq op :multiplication)
     (cond ((= x 0d0) nil)
           ((= y 0d0) nil)
           (t (let ((x (abs x)) (y (abs y)))
                (cond ((and (>= x 1d0) (<= y 1d0)) nil)
                      ((and (<= x 1d0) (>= y 1d0)) nil)
                      ((>= x 1d0)
                       (not (or (<= y (/ +biggest-double+ x))
                                (<= x (/ +biggest-double+ y)))))
                      (t (not (> x (/ +smallest-double+ y)))))))))
    ((eq op :addition)
     (if (/= (nl-sign x) (nl-sign y))
         nil
         (not (< (abs y) (- +biggest-double+ (abs x))))))
    ((eq op :subtraction)
     (if (> (* (nl-sign x) (nl-sign y)) 0d0)
         nil
         (not (> (- (abs y)) (- (abs x) +biggest-double+)))))
    (t t)))

(declaim (inline homogeneous-to-euclid))
(defun homogeneous-to-euclid (wx wy wz w)
  "The Euclidean point of homogeneous point (WX WY WZ W), as three values,
honoring +no-z+ and +no-w+.  NLib: N_CPtToPtEuclid -- including its multiply
by the reciprocal of w, which rounds differently from dividing by w."
  (declare (type double-float wx wy wz w))
  (fast
    (if (/= w +no-w+)
        (let ((r (/ 1d0 w)))
          (values (* wx r) (* wy r) (if (/= wz +no-z+) (* wz r) 0d0)))
        (values wx wy (if (/= wz +no-z+) wz 0d0)))))

(declaim (inline cpoint-to-euclid))
(defun cpoint-to-euclid (cpoints i)
  "The Euclidean point of homogeneous control point I of CPOINTS, as three
values.  NLib: N_CPtToPtEuclid."
  (declare (type real-vector cpoints) (type index i))
  (fast
    (let ((base (* 4 i)))
      (homogeneous-to-euclid (aref cpoints base) (aref cpoints (+ base 1))
                             (aref cpoints (+ base 2)) (aref cpoints (+ base 3))))))
