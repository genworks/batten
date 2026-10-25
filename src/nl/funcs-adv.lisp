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

;;; ------------------------------------------------------------------- L2-g3

(defun crv-func-eval-invert-pt (cfn f u0 tol itl)
  "The parameter u at which curve function CFN takes the value F, by Newton
iteration from the guess U0 with tolerance TOL and at most ITL iterations (a
degree-1 function is inverted directly, on the first span whose values
bracket F from below).  Signals :inp-err when F lies outside the values at
the ends of the domain and :con-err when the iteration fails or does not
converge.  Returns u.  For degree 1 with no bracketing span (possible only for
TOL <= 0) the C leaves u unset; this returns NIL there.
NLib: N_CrvFuncEvalInvertPt."
  (declare (type cfun cfn) (type double-float f u0 tol) (type fixnum itl))
  (let* ((n (1- (length (cfun-vals cfn))))
         (p (cfun-degree cfn))
         (uk (cfun-knots cfn))
         (m (1- (length uk))))
    (declare (type fixnum n p m) (type real-vector uk))
    ;; end values
    (let ((a (c-func-eval cfn (aref uk 0) +left+))
          (b (c-func-eval cfn (aref uk m) +left+)))
      (declare (type double-float a b))
      (when (< (abs (- f a)) tol)
        (return-from crv-func-eval-invert-pt (aref uk 0)))
      (when (< (abs (- b f)) tol)
        (return-from crv-func-eval-invert-pt (aref uk m)))
      (when (or (< f a) (> f b))
        (nurbs-error :inp-err 'crv-func-eval-invert-pt))
      ;; p = 1 as special case
      (when (= p 1)
        (do ((k p (1+ k))) ((> k n))
          (declare (type fixnum k))
          (setf b (c-func-eval cfn (aref uk (1+ k)) +left+))
          (when (and (< f b) (>= f a))
            (let ((alf (/ (- f a) (- b a))))
              (declare (type double-float alf))
              (return-from crv-func-eval-invert-pt
                (+ (* (- 1d0 alf) (aref uk k)) (* alf (aref uk (1+ k)))))))
          (setf a b))
        (return-from crv-func-eval-invert-pt nil)))
    ;; Newton iterations till convergence
    (let ((k 0)
          (unew u0)
          (fd (make-real-vector 2)))
      (declare (type fixnum k) (type double-float unew) (dynamic-extent fd))
      (do () ((not (<= k itl)))
        (c-func-derivs cfn unew +left+ 1 fd)
        (let ((num (- (aref fd 0) f)))
          (declare (type double-float num))
          (when (<= (abs num) tol)
            (return))
          (when (float-op-is-bad num (aref fd 1) +division+)
            (nurbs-error :con-err 'crv-func-eval-invert-pt))
          (let ((uold unew))
            (declare (type double-float uold))
            (setf unew (- uold (/ num (aref fd 1))))
            (when (< unew (aref uk 0))
              (setf unew (aref uk 0)))
            (when (> unew (aref uk m))
              (setf unew (aref uk m)))
            (when (<= (abs (* (- unew uold) (aref fd 1))) tol)
              (return))))
        (incf k))
      (when (> k itl)
        (nurbs-error :con-err 'crv-func-eval-invert-pt))
      unew)))

(defun crv-func-degree-elevate (cfn-p tt)
  "Curve function CFN-P with its degree raised by TT, fresh: each Bezier
segment is extracted by knot insertion, degree elevated, and the knots
between segments removed again.  Signals :deg-err unless 0 <= TT and the new
degree is at most *dmax*, :knt-err for an invalid knot vector, and the errors
of BEZ-FUNC-DEGREE-ELEVATE.
NLib: N_CrvFuncDegreeElevate (the C's in-place case, cfnP = cfnQ, computes the
same function)."
  (declare (type cfun cfn-p) (type fixnum tt))
  (let* ((fp (cfun-vals cfn-p)) (np (1- (length fp)))
         (p (cfun-degree cfn-p))
         (up (cfun-knots cfn-p)) (mp (1- (length up))))
    (declare (type real-vector fp up) (type fixnum np p mp))
    (when (or (< tt 0) (> (+ p tt) (the fixnum *dmax*)))
      (nurbs-error :deg-err 'crv-func-degree-elevate))
    (knot-vector-is-valid up p 'crv-func-degree-elevate)
    (let* ((s (basis-get-span-count up p))
           (nq (+ np (* tt s)))
           (q (+ p tt))
           (mq (+ mp (* tt (1+ s))))
           (fq (make-real-vector (1+ nq)))
           (uq (make-real-vector (1+ mq)))
           (fb (make-real-vector (1+ p)))
           (fn (make-real-vector (1+ p)))
           (fd (make-real-vector (1+ q)))
           (alfs (make-real-vector (1+ p)))
           (omas (make-real-vector (1+ p)))
           (a p) (b (1+ p)) (r -1) (find 1) (kind (1+ q)))
      (declare (type fixnum s nq q mq a b r find kind)
               (type real-vector fq uq fb fn fd alfs omas))
      ;; initialize
      (setf (aref fq 0) (aref fp 0))
      (do ((i 0 (1+ i))) ((> i q))
        (declare (type fixnum i))
        (setf (aref uq i) (aref up a)))
      (do ((i 0 (1+ i))) ((> i p))
        (declare (type fixnum i))
        (setf (aref fb i) (aref fp i)))
      ;; degree elevation matrix
      (let ((dm (make-rmatrix q p (bez-get-degree-elevation-matrix p tt) +mt-full+ p)))
        ;; extract each Bezier segment, elevate it, remove the knot between it
        ;; and the one before
        (do () ((not (< b mp)))
          (let* ((i0 b)
                 (mlt 0) (oldr r) (lbz 0) (rbz 0))
            (declare (type fixnum i0 mlt oldr lbz rbz))
            ;; multiplicity of the knot
            (do () ((not (and (< b mp) (= (aref up b) (aref up (1+ b))))))
              (incf b))
            (setf mlt (1+ (- b i0)))
            ;; insert the knot to get the Bezier segment
            (setf r (- p mlt))
            (setf lbz (if (> oldr 0) (truncate (+ oldr 2) 2) 1))
            (setf rbz (if (> r 0) (- q (truncate (1+ r) 2)) q))
            (when (> r 0)
              (let ((num (- (aref up b) (aref up a))))
                (declare (type double-float num))
                (fast
                  (do ((k p (1- k))) ((<= k mlt))
                    (declare (type fixnum k))
                    (setf (aref alfs (- k mlt 1)) (/ num (- (aref up (+ a k)) (aref up a))))
                    (setf (aref omas (- k mlt 1)) (- 1d0 (aref alfs (- k mlt 1)))))
                  (do ((j 1 (1+ j))) ((> j r))
                    (declare (type fixnum j))
                    (let ((save (- r j))
                          (ss (+ mlt j)))
                      (declare (type fixnum save ss))
                      (do ((k p (1- k))) ((< k ss))
                        (declare (type fixnum k))
                        (setf (aref fb k) (+ (* (aref alfs (- k ss)) (aref fb k))
                                             (* (aref omas (- k ss)) (aref fb (1- k))))))
                      (setf (aref fn save) (aref fb p)))))))
            ;; degree elevate the Bezier segment
            (bez-func-degree-elevate fb p tt dm lbz q fd)
            ;; remove the knot UP[a]
            (when (> oldr 1)
              (let* ((first (- kind 2))
                     (last kind)
                     (den (- (aref up b) (aref up a)))
                     (bet (/ (- (aref up b) (aref uq (1- kind))) den))
                     (omb (- 1d0 bet)))
                (declare (type fixnum first last) (type double-float den bet omb))
                (fast
                  (do ((k 1 (1+ k))) ((>= k oldr))
                    (declare (type fixnum k))
                    (let* ((i first)
                           (j last)
                           (l (1+ (- j kind))))
                      (declare (type fixnum i j l))
                      (do () ((not (> (- j i) k)))
                        (when (< i find)
                          (let* ((alf (/ (- (aref up b) (aref uq i)) (- (aref up a) (aref uq i))))
                                 (oma (- 1d0 alf)))
                            (declare (type double-float alf oma))
                            (setf (aref fq i) (+ (* alf (aref fq i)) (* oma (aref fq (1- i)))))))
                        (when (>= j lbz)
                          (if (<= (- j k) (+ (- kind q) oldr))
                              (let* ((gam (/ (- (aref up b) (aref uq (- j k))) den))
                                     (omg (- 1d0 gam)))
                                (declare (type double-float gam omg))
                                (setf (aref fd l) (+ (* gam (aref fd l)) (* omg (aref fd (1+ l))))))
                              (setf (aref fd l) (+ (* bet (aref fd l)) (* omb (aref fd (1+ l)))))))
                        (incf i)
                        (decf j)
                        (decf l)))
                    (decf first)
                    (incf last)))))
            ;; load the knot vector and control values
            (when (/= a p)
              (dotimes (i (- q oldr))
                (setf (aref uq kind) (aref up a))
                (incf kind)))
            (do ((i lbz (1+ i))) ((> i rbz))
              (declare (type fixnum i))
              (setf (aref fq find) (aref fd i))
              (incf find))
            ;; initialize for the next pass
            (if (< b mp)
                (progn
                  (dotimes (i r)
                    (setf (aref fb i) (aref fn i)))
                  (do ((i r (1+ i))) ((> i p))
                    (declare (type fixnum i))
                    (setf (aref fb i) (aref fp (+ (- b p) i))))
                  (setf a b)
                  (incf b))
                ;; the end knot
                (do ((i 0 (1+ i))) ((> i q))
                  (declare (type fixnum i))
                  (setf (aref uq (+ kind i)) (aref up b)))))))
      (%make-cfun q uq fq))))
