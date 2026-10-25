;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvFit.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.
;;;;
;;;; Knot vector outputs (NL_KNOTVECTOR *knt, "MEMORY FOR knt IS ALLOCATED IN
;;;; THE CALLING ROUTINE") are real-vectors the caller allocates with exactly
;;;; the number of knots the C routine records in knt->m (its N_SetKnotIndex),
;;;; plus one; the routine writes the knots the C writes and leaves the others
;;;; as they were.

(in-package :batten)

;;; ------------------------------------------------- parameters for interpolation

(defun fit-calc-crv-param-values (a n ptp par u)
  "Write into U[0..N] parameter values for interpolating the points of A
(rows 0..N): a point-array when PTP is +epoint+, a cpoint-array when it is
+hpoint+.  PAR is +uniform+, +chordlength+ or +centripetal+.  U is allocated
by the caller; when the total chord length is not above *ptol*, U[1..N-1] are
left as they were.  Signals :cal-err for an unknown PTP or PAR (U[0] and U[N]
already written).  Returns U.
NLib: N_FitCalcCrvParamValues."
  (declare (type (simple-array double-float (* *)) a) (type fixnum n ptp par) (type real-vector u))
  (setf (aref u 0) 0d0
        (aref u n) 1d0)
  (when (= par +uniform+)
    ;; fact = 1.0/n; with n = 0 the C's infinite fact is never used
    (when (> n 0)
      (let ((fact (/ 1d0 (float n 1d0))))
        (declare (type double-float fact))
        (fast (do ((i 1 (1+ i))) ((>= i n))
                (declare (type fixnum i))
                (setf (aref u i) (* (float i 1d0) fact))))))
    (return-from fit-calc-crv-param-values u))
  (let ((dist (make-real-vector (1+ n)))
        (sum 0d0))
    (declare (type real-vector dist) (type double-float sum))
    (cond
      ((= ptp +epoint+)
       (let ((pts a) (p0 (make-point)) (p1 (make-point)))
         (declare (type point-array pts) (type point p0 p1) (dynamic-extent p0 p1))
         (flet ((dist-at (i)
                  (declare (type fixnum i))
                  (set-point p0 (aref pts (1- i) 0) (aref pts (1- i) 1) (aref pts (1- i) 2))
                  (set-point p1 (aref pts i 0) (aref pts i 1) (aref pts i 2))
                  (dist-pt-pt p0 p1)))
           (declare (inline dist-at))
           (cond ((= par +chordlength+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf sum (+ sum (aref dist i)))))
                 ((= par +centripetal+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf (aref dist i) (sqrt (the (double-float 0d0) (aref dist i))))
                    (setf sum (+ sum (aref dist i)))))
                 (t (nurbs-error :cal-err 'fit-calc-crv-param-values))))))
      ((= ptp +hpoint+)
       (let ((pw a) (q0 (make-cpoint)) (q1 (make-cpoint)))
         (declare (type cpoint-array pw) (type cpoint q0 q1) (dynamic-extent q0 q1))
         (flet ((dist-at (i)
                  (declare (type fixnum i))
                  (set-cpoint q0 (aref pw (1- i) 0) (aref pw (1- i) 1) (aref pw (1- i) 2) (aref pw (1- i) 3))
                  (set-cpoint q1 (aref pw i 0) (aref pw i 1) (aref pw i 2) (aref pw i 3))
                  (dist-cpt-cpt-homo q0 q1)))
           (declare (inline dist-at))
           (cond ((= par +chordlength+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf sum (+ sum (aref dist i)))))
                 ((= par +centripetal+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf (aref dist i) (sqrt (the (double-float 0d0) (aref dist i))))
                    (setf sum (+ sum (aref dist i)))))
                 (t (nurbs-error :cal-err 'fit-calc-crv-param-values))))))
      (t (nurbs-error :cal-err 'fit-calc-crv-param-values)))
    (when (> sum *ptol*)
      (fast (do ((i 1 (1+ i))) ((>= i n))
              (declare (type fixnum i))
              (setf (aref u i) (+ (aref u (1- i)) (/ (aref dist i) sum))))))
    u))

;;; ------------------------------------------------------------- knot vectors

(defun fit-crv-calc-knot-vector (u n p knt)
  "Write into KNT (N+P+2 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..N]: P+1 end knots
equal to U[0] and U[N], interior knots averages of P consecutive parameters.
Returns KNT.
NLib: N_FitCrvCalcKnotVector."
  (declare (type real-vector u knt) (type fixnum n p))
  (fast
    (dotimes (i (1+ p))
      (setf (aref knt i) (aref u 0)
            (aref knt (+ n i 1)) (aref u n)))
    (do ((i 1 (1+ i))) ((> i (- n p)))
      (declare (type fixnum i))
      (let ((sum 0d0))
        (declare (type double-float sum))
        (do ((j i (1+ j))) ((> j (+ i p -1)))
          (declare (type fixnum j))
          (setf sum (+ sum (aref u j))))
        (setf (aref knt (+ i p)) (/ sum (float p 1d0))))))
  knt)

(defun fit-calc-knot-vector-end-derivs (u k p knt)
  "Write into KNT (K+P+4 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..K] with both end
derivatives given.  Returns KNT.
NLib: N_FitCalcKnotVectorEndDerivs."
  (declare (type real-vector u knt) (type fixnum k p))
  (let ((n (+ k 2)))
    (declare (type fixnum n))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (if (= p 1)
          (progn
            (if (> k 1)
                (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1)))
                      (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k)))))
                (setf (aref knt 2) (+ (* 0.75d0 (aref u 0)) (* 0.25d0 (aref u 1)))
                      (aref knt 3) (+ (* 0.25d0 (aref u 0)) (* 0.75d0 (aref u 1)))))
            (do ((i 1 (1+ i))) ((>= i k))
              (declare (type fixnum i))
              (setf (aref knt (+ i 2)) (aref u i))))
          (do ((i 0 (1+ i))) ((> i (+ (- k p) 1)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p 1)) (/ sum (float p 1d0)))))))
    knt))

(defun fit-calc-knot-vector-deriv (u k p whr knt)
  "Write into KNT (K+P+3 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..K] with the end
derivative given at WHR, +start+ or +end+.  Signals :inp-err when P > K+1 or
WHR is neither.  Returns KNT.
NLib: N_FitCalcKnotVectorDeriv."
  (declare (type real-vector u knt) (type fixnum k p whr))
  (when (or (> p (1+ k)) (and (/= whr +start+) (/= whr +end+)))
    (nurbs-error :inp-err 'fit-calc-knot-vector-deriv))
  (let ((n (1+ k)))
    (declare (type fixnum n))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (when (= p 1)                     ; special case out
        (if (= k 1)
            (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1))))
            (progn
              (when (= whr +start+)
                (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1))))
                (do ((i 1 (1+ i))) ((>= i k))
                  (declare (type fixnum i))
                  (setf (aref knt (+ i 2)) (aref u i))))
              (when (= whr +end+)
                (setf (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k)))))
                (do ((i 1 (1+ i))) ((>= i k))
                  (declare (type fixnum i))
                  (setf (aref knt (+ i 1)) (aref u i)))))))
      (when (and (> p 1) (< p n))
        (when (= whr +start+)
          (do ((i 1 (1+ i))) ((> i (- k p)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p 1)) (/ sum (float p 1d0)))))
          (let ((sum 0d0))
            (declare (type double-float sum))
            (do ((j 0 (1+ j))) ((> j (1- p)))
              (declare (type fixnum j))
              (setf sum (+ sum (aref u j))))
            (setf (aref knt (1+ p)) (/ sum (float p 1d0)))))
        (when (= whr +end+)
          (do ((i 1 (1+ i))) ((> i (- k p)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p)) (/ sum (float p 1d0)))))
          (let ((sum 0d0))
            (declare (type double-float sum))
            (do ((j (1+ (- k p)) (1+ j))) ((> j k))
              (declare (type fixnum j))
              (setf sum (+ sum (aref u j))))
            (setf (aref knt n) (/ sum (float p 1d0)))))))
    knt))

(defun fit-calc-knot-vector-crv-approx (u m n p knt)
  "Write into KNT (N+P+2 knots, allocated by the caller) the knot vector for
global approximation of samples at ascending parameters U[0..M] by a curve of
degree P with control points 0..N: interior knots average P consecutive
representatives of the N+1 clusters of parameters.  Signals :inp-err when
N > M.  Returns KNT.
NLib: N_FitCalcKnotVectorCrvApprox."
  (declare (type real-vector u knt) (type fixnum m n p))
  (when (> n m)
    (nurbs-error :inp-err 'fit-calc-knot-vector-crv-approx))
  (fast
    (dotimes (i (1+ p))
      (setf (aref knt i) (aref u 0)
            (aref knt (+ n i 1)) (aref u m))))
  (let ((uk (make-real-vector (1+ n)))
        (d (/ (+ (float m 1d0) 1d0) (+ (float n 1d0) 1d0)))
        (il 0) (ih 0) (l -1d0))
    (declare (type real-vector uk) (type double-float d l) (type fixnum il ih))
    (fast
      (do ((i 0 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (setf l (+ l d))
        (setf ih (the fixnum (truncate (+ l 0.5d0)))) ; ROUND(l)
        (if (= il ih)
            (setf (aref uk i) (aref u il))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j il (1+ j))) ((> j ih))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref uk i) (/ sum (+ (- (float ih 1d0) (float il 1d0)) 1d0)))))
        (setf il (1+ ih)))
      (do ((i 1 (1+ i))) ((> i (- n p)))
        (declare (type fixnum i))
        (let ((sum 0d0))
          (declare (type double-float sum))
          (do ((j i (1+ j))) ((> j (+ i p -1)))
            (declare (type fixnum j))
            (setf sum (+ sum (aref uk j))))
          (setf (aref knt (+ i p)) (/ sum (float p 1d0))))))
    knt))

(defun fit-calc-knots-derivs (u k p knt)
  "Write into KNT (2K+P+3 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P (2 or 3; any other degree gets only its
end knots) at parameters U[0..K] with a derivative given at every point.
Returns KNT.
NLib: N_FitCalcKnotsDerivs."
  (declare (type real-vector u knt) (type fixnum k p))
  (let ((n (1+ (* 2 k)))
        (j (1+ p)))
    (declare (type fixnum n j))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (case p
        (2
         (setf (aref knt (1+ p)) (* 0.5d0 (+ (aref u 1) (aref u 0))))
         (do ((i 1 (1+ i))) ((>= i k))
           (declare (type fixnum i))
           (setf (aref knt (incf j)) (aref u i))
           (setf (aref knt (incf j)) (* 0.5d0 (+ (aref u i) (aref u (1+ i)))))))
        (3
         (when (> k 1)
           (setf (aref knt (1+ p)) (* 0.5d0 (+ (aref u 1) (aref u 0)))
                 (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k))))))
         (let ((alf (/ 1d0 3d0))
               (bet (/ 2d0 3d0)))
           (declare (type double-float alf bet))
           (do ((i 1 (1+ i))) ((>= i (1- k)))
             (declare (type fixnum i))
             (setf (aref knt (incf j)) (+ (* bet (aref u i)) (* alf (aref u (1+ i)))))
             (setf (aref knt (incf j)) (+ (* alf (aref u i)) (* bet (aref u (1+ i))))))))))
    knt))

(defun fit-calc-knots-high-end-derivs (u n p k l knt)
  "Write into KNT (N+K+L+P+2 knots, allocated by the caller) the knot vector
for global curve interpolation of degree P at parameters U[0..N] with the
first K derivatives given at the start and the first L at the end (K, L < P).
Returns KNT.
NLib: N_FitCalcKnotsHighEndDerivs."
  (declare (type real-vector u knt) (type fixnum n p k l))
  (let* ((m (+ n k l))
         (is (- 1 k))
         (ie (+ (- n p) l))
         (r p))
    (declare (type fixnum m is ie r))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ m i 1)) (aref u n)))
      (do ((i is (1+ i))) ((> i ie))
        (declare (type fixnum i))
        (let ((js (nl-max 0 i))
              (je (nl-min n (+ i p -1)))
              (sum 0d0))
          (declare (type fixnum js je) (type double-float sum))
          (do ((j js (1+ j))) ((> j je))
            (declare (type fixnum j))
            (setf sum (+ sum (aref u j))))
          (setf (aref knt (incf r)) (/ sum (+ (- (float je 1d0) (float js 1d0)) 1d0))))))
    knt))
