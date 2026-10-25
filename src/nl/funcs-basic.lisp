;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_FuncsBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; NL_CFUN and NL_SFUN hold their control values in an NL_CVALUE / NL_SVALUE
;;; object and their knots in NL_KNOTVECTOR objects; a cfun or sfun holds the
;;; values vector or matrix and the knot vectors directly.  So in these
;;; routines an NL_CVALUE is the real-vector of control values, an NL_SVALUE
;;; the real-matrix of them, and an NL_KNOTVECTOR the real-vector of knots.
;;; Where the C makes pointer assignments, these share the Lisp arrays.

(defun crv-is-func-sized (cfn nf mk rname)
  "Signal :sto-err in routine RNAME unless curve function CFN has room for
control values fu[0..NF] and knots U[0..MK].  NLib: N_CrvIsFuncSized."
  (declare (type cfun cfn) (type fixnum nf mk))
  (let ((n (1- (length (cfun-vals cfn))))
        (m (1- (length (cfun-knots cfn)))))
    (declare (type fixnum n m))
    (when (or (< n nf) (< m mk))
      (nurbs-error :sto-err rname))
    nil))

(defun c-value-from-array (fu n)
  "The control value object for function values FU[0..N]: FU itself when it
has exactly N+1 elements, else a fresh vector of its first N+1 (a Lisp vector
cannot share a prefix of another, as the C's pointer does).
NLib: N_CValueFromArray."
  (declare (type real-vector fu) (type fixnum n))
  (if (= (length fu) (1+ n))
      fu
      (subseq fu 0 (1+ n))))

(defun c-func-from-knot-vector (cfn cvl p knt)
  "Make curve function CFN, in place, from control values CVL, degree P and
knot vector KNT, sharing the arrays.  Return CFN.  NLib: N_CFuncFromKnotVector."
  (declare (type cfun cfn) (type real-vector cvl knt) (type fixnum p))
  (setf (cfun-vals cfn) cvl
        (cfun-degree cfn) p
        (cfun-knots cfn) knt)
  cfn)

(defun crv-func-cntrl-val (cfn)
  "The control values of curve function CFN: (values n fu), n the highest
index in fu.  NLib: N_CrvFuncCntrlVal."
  (declare (type cfun cfn))
  (let ((fu (cfun-vals cfn)))
    (values (1- (length fu)) fu)))

(defun crv-func-cntrl-val-knot-vector (cfn)
  "The control values, knot vector and knots of curve function CFN:
(values fu knt U), the knot vector and its knots being one array here.
NLib: N_CrvFuncCntrlValKnotVector."
  (declare (type cfun cfn))
  (values (cfun-vals cfn) (cfun-knots cfn) (cfun-knots cfn)))

(defun crv-func-cntrl-val-knots (cfn)
  "The control values and knots of curve function CFN: (values fu U).
NLib: N_CrvFuncCntrlValKnots."
  (declare (type cfun cfn))
  (values (cfun-vals cfn) (cfun-knots cfn)))

(defun crv-func-set-ptrs (cfn fu u)
  "Set the control values of curve function CFN to FU and its knots to U, in
place, sharing the arrays.  Return CFN.  NLib: N_CrvFuncSetPtrs."
  (declare (type cfun cfn) (type real-vector fu u))
  (setf (cfun-vals cfn) fu
        (cfun-knots cfn) u)
  cfn)

(defun srf-is-func-sized (sfn nf mf rk sk rname)
  "Signal :sto-err in routine RNAME unless surface function SFN has room for
control values fuv[0..NF][0..MF] and knots U[0..RK], V[0..SK].
NLib: N_SrfIsFuncSized."
  (declare (type sfun sfn) (type fixnum nf mf rk sk))
  (let ((n (1- (array-dimension (sfun-vals sfn) 0)))
        (m (1- (array-dimension (sfun-vals sfn) 1)))
        (r (1- (length (sfun-knots-u sfn))))
        (s (1- (length (sfun-knots-v sfn)))))
    (declare (type fixnum n m r s))
    (when (or (< n nf) (< m mf) (< r rk) (< s sk))
      (nurbs-error :sto-err rname))
    nil))

(defun s-value-from-array (fuv n m)
  "The control value object for function values FUV[0..N][0..M]: FUV itself
when it is exactly (N+1) x (M+1), else a fresh matrix of that corner of it.
NLib: N_SValueFromArray."
  (declare (type real-matrix fuv) (type fixnum n m))
  (if (and (= (array-dimension fuv 0) (1+ n)) (= (array-dimension fuv 1) (1+ m)))
      fuv
      (let ((out (make-real-matrix (1+ n) (1+ m))))
        (dotimes (i (1+ n) out)
          (dotimes (j (1+ m))
            (setf (aref out i j) (aref fuv i j)))))))

(defun s-func-from-knot-vectors (sfn svl p q knu knv)
  "Make surface function SFN, in place, from control values SVL, degrees P
and Q and knot vectors KNU and KNV, sharing the arrays.  Return SFN.
NLib: N_SFuncFromKnotVectors."
  (declare (type sfun sfn) (type real-matrix svl) (type fixnum p q) (type real-vector knu knv))
  (setf (sfun-vals sfn) svl
        (sfun-degree-u sfn) p
        (sfun-degree-v sfn) q
        (sfun-knots-u sfn) knu
        (sfun-knots-v sfn) knv)
  sfn)

(defun srf-func-cntrl-val (sfn)
  "The control values of surface function SFN: (values n m fuv), n and m the
highest indexes in fuv.  NLib: N_SrfFuncCntrlVal."
  (declare (type sfun sfn))
  (let ((fuv (sfun-vals sfn)))
    (values (1- (array-dimension fuv 0)) (1- (array-dimension fuv 1)) fuv)))

(defun s-func-get-knot-vectors (sfn)
  "The knot vectors of surface function SFN: (values knu knv).
NLib: N_SFuncGetKnotVectors."
  (declare (type sfun sfn))
  (values (sfun-knots-u sfn) (sfun-knots-v sfn)))

(defun s-func-get-array-and-knot-vectors (sfn)
  "The parts of surface function SFN: (values svl p q knu knv).
NLib: N_SFuncGetArrayAndKnotVectors."
  (declare (type sfun sfn))
  (values (sfun-vals sfn) (sfun-degree-u sfn) (sfun-degree-v sfn)
          (sfun-knots-u sfn) (sfun-knots-v sfn)))

;;; ------------------------------------------------------------ level 1, group 1

(defun c-func-from-knots (cfn fu n p u m)
  "Make curve function CFN, in place, from function values FU[0..N], degree
P and knots U[0..M], sharing the arrays where they are exactly that long, as
the C shares the pointers.  Return CFN.
NLib: N_CFuncFromKnots."
  (declare (type cfun cfn) (type real-vector fu u) (type fixnum n p m))
  (let ((cvl (c-value-from-array fu n))
        (knt (knot-vector-from-real-array u m)))
    (declare (type real-vector cvl knt))
    (c-func-from-knot-vector cfn cvl p knt)))

(defun crv-func-compact (cfn)
  "Give curve function CFN, in place, fresh control value and knot arrays
holding copies of its values and knots.  Return CFN.
NLib: N_CrvFuncCompact."
  (declare (type cfun cfn))
  (multiple-value-bind (fu u) (crv-func-cntrl-val-knots cfn)
    (declare (type real-vector fu u))
    (let* ((n (1- (length fu)))
           (m (1- (length u)))
           (fv (make-real-vector (1+ n)))
           (v (make-real-vector (1+ m))))
      (declare (type fixnum n m) (type real-vector fv v))
      (fast
        (dotimes (i (1+ n))
          (setf (aref fv i) (aref fu i)))
        (dotimes (i (1+ m))
          (setf (aref v i) (aref u i))))
      (crv-func-set-ptrs cfn fv v))))

(defun crv-func-copy (cfnp)
  "A fresh copy of curve function CFNP: its degree, control values and knots.
NLib: N_CrvFuncCopy."
  (declare (type cfun cfnp))
  (let* ((fp (cfun-vals cfnp))
         (up (cfun-knots cfnp))
         (n (1- (length fp)))
         (m (1- (length up)))
         (fq (make-real-vector (1+ n)))
         (uq (make-real-vector (1+ m))))
    (declare (type real-vector fp up fq uq) (type fixnum n m))
    (fast
      (dotimes (i (1+ n))
        (setf (aref fq i) (aref fp i)))
      (dotimes (i (1+ m))
        (setf (aref uq i) (aref up i))))
    (%make-cfun (cfun-degree cfnp) uq fq)))

(defun c-func-eval (cfn u flg)
  "The value of curve function CFN at parameter U, FLG (+left+ or +right+)
choosing the knot span at a knot.  Signals :par-err for U outside the
domain.
NLib: N_CFuncEval."
  (declare (type cfun cfn) (type double-float u) (type fixnum flg))
  (multiple-value-bind (fu knots) (crv-func-cntrl-val-knots cfn)
    (declare (type real-vector fu knots))
    (let ((p (cfun-degree cfn)))
      (declare (type fixnum p))
      (knot-vector-is-param-out-of-bounds knots u 'c-func-eval)
      (let ((n (make-array (1+ p) :element-type 'double-float :initial-element 0d0)))
        (declare (dynamic-extent n))
        (let ((spn (basis-eval knots p u flg n))
              (f 0d0))
          (declare (type fixnum spn) (type double-float f))
          (fast
            (dotimes (i (1+ p))
              (setf f (+ f (* (aref fu (+ (- spn p) i)) (aref n i))))))
          f)))))

(defun c-func-derivs (cfn u flg der fd)
  "Write the derivatives of curve function CFN at parameter U, orders 0 to
DER, into FD (a real-vector of at least DER+1, allocated by the caller):
FD[k] is the kth derivative.  FLG (+left+ or +right+) chooses the knot span
at a knot.  Signals :par-err for U outside the domain.  Return FD.
NLib: N_CFuncDerivs."
  (declare (type cfun cfn) (type double-float u) (type fixnum flg der) (type real-vector fd))
  (multiple-value-bind (fu knots) (crv-func-cntrl-val-knots cfn)
    (declare (type real-vector fu knots))
    (let ((p (cfun-degree cfn)))
      (declare (type fixnum p))
      (knot-vector-is-param-out-of-bounds knots u 'c-func-derivs)
      (let* ((nd (make-real-matrix (1+ der) (1+ p)))
             (spn (basis-derivs knots p u flg der nd)))
        (declare (type real-matrix nd) (type fixnum spn))
        (fast
          (do ((k 0 (1+ k))) ((> k der))
            (declare (type fixnum k))
            (setf (aref fd k) 0d0)
            (do ((j (- spn p) (1+ j))) ((> j spn))
              (declare (type fixnum j))
              (setf (aref fd k) (+ (aref fd k) (* (aref fu j) (aref nd k (+ (- j spn) p))))))))
        fd))))
