;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfNurb.cpp and NLI_Math.h (usd-brep source/SMLib).
;;;; Modified: translated to Common Lisp.  Scratch arrays keep the C's fixed
;;;; sizes ([NL_MAXDEG+1][NL_MAXDEG+1] and so on), flattened, on the stack.

(in-package :batten)

(deftype basis-derivatives () '(simple-array double-float (* * * *)))

(defconstant +deg+1+ (1+ +maxdeg+))
(defconstant +der+1+ (1+ +maxder+))

(declaim (inline srf-nurb--st-compute-basis-functions srf-nurb--st-compute-derivatives))

(defun srf-nurb--st-compute-basis-functions (p u knots ku ndu)
  "Fill NDU, a flat [NL_MAXDEG+1][NL_MAXDEG+1] array, with the triangular
table of basis functions and knot differences of degree P at U in span KU.
NLib: N_ComputeBasisFunctions (static in NL_SrfNurb.cpp)."
  (declare (type fixnum p ku) (type double-float u) (type real-vector knots ndu))
  (fast
    (let ((left (make-array +deg+1+ :element-type 'double-float))
          (right (make-array +deg+1+ :element-type 'double-float)))
      (declare (dynamic-extent left right))
      (macrolet ((ndu (r c) `(aref ndu (+ (the fixnum (* ,r +deg+1+)) ,c))))
        (setf (ndu 0 0) 1d0)
        (do ((j 1 (1+ j))) ((> j p))
          (declare (type fixnum j))
          (setf (aref left j) (- u (aref knots (- (+ ku 1) j)))
                (aref right j) (- (aref knots (+ ku j)) u))
          (let ((saved 0d0))
            (declare (type double-float saved))
            (dotimes (r j)
              (setf (ndu j r) (+ (aref right (1+ r)) (aref left (- j r))))
              (let ((temp (/ (ndu r (1- j)) (ndu j r))))
                (declare (type double-float temp))
                (setf (ndu r j) (+ saved (* (aref right (1+ r)) temp))
                      saved (* (aref left (- j r)) temp))))
            (setf (ndu j j) saved)))))))

(defun srf-nurb--st-compute-derivatives (p mder ndu du)
  "Fill DU[k][r], k = 1..MDER, a flat [NL_MAXDER+1][NL_MAXDEG+1] array, with
the unscaled basis function derivatives from the table NDU.
NLib: N_ComputeDerivatives / N_ComputeDerivativesImpl (static in NL_SrfNurb.cpp)."
  (declare (type fixnum p mder) (type real-vector ndu du))
  (fast
    (let ((a (make-array (* 2 +deg+1+) :element-type 'double-float)))
      (declare (dynamic-extent a))
      (macrolet ((ndu (r c) `(aref ndu (+ (the fixnum (* ,r +deg+1+)) ,c)))
                 (du (k r) `(aref du (+ (the fixnum (* ,k +deg+1+)) ,r)))
                 (a (s j) `(aref a (+ (the fixnum (* ,s +deg+1+)) ,j))))
        (dotimes (r (1+ p))
          (let ((s1 0) (s2 1))
            (declare (type fixnum s1 s2))
            (setf (a s1 0) 1d0)
            (do ((k 1 (1+ k))) ((> k mder))
              (declare (type fixnum k))
              (let* ((dd 0d0)
                     (rk (- r k))
                     (pk (- p k))
                     (j1 (if (>= rk -1) 1 (- rk)))
                     (j2 (if (<= (1- r) pk) (1- k) (- p r))))
                (declare (type double-float dd) (type fixnum rk pk j1 j2))
                (when (>= r k)
                  (setf (a s2 0) (/ (a s1 0) (ndu (1+ pk) rk))
                        dd (* (a s2 0) (ndu rk pk))))
                (do ((j j1 (1+ j))) ((> j j2))
                  (declare (type fixnum j))
                  (setf (a s2 j) (/ (- (a s1 j) (a s1 (1- j))) (ndu (1+ pk) (+ rk j))))
                  (incf dd (* (a s2 j) (ndu (+ rk j) pk))))
                (when (<= r pk)
                  (setf (a s2 k) (/ (- (a s1 (1- k))) (ndu (1+ pk) r)))
                  (incf dd (* (a s2 k) (ndu r pk))))
                (setf (du k r) dd)
                (rotatef s1 s2)))))))))

(defun %srf-basis-derivs (sur u v ufl vfl mfl udr vdr bd)
  "SRF-BASIS-DERIVS into a flat BD of dimensions (1+ UDR) (1+ VDR) (1+ p) (1+ q)
in row-major order.  Returns (values usp vsp)."
  (declare (type surface sur) (type double-float u v) (type fixnum ufl vfl mfl udr vdr)
           (type real-vector bd))
  (when (or (> udr +maxder+) (> vdr +maxder+))
    (nurbs-error :mxd-err 'srf-basis-derivs))
  (let* ((p (surface-degree-u sur)) (q (surface-degree-v sur))
         (knu (surface-knots-u sur)) (knv (surface-knots-v sur))
         (pw (surface-cpoints sur)))
    (declare (type fixnum p q) (type real-vector knu knv) (type cpoint-net pw))
    (fast
      (let* ((ku (basis-find-span knu p u ufl))
             (kv 0)
             (p+1 (1+ p)) (q+1 (1+ q)) (vdr+1 (1+ vdr))
             (ndu (make-array (* +deg+1+ +deg+1+) :element-type 'double-float))
             (du (make-array (* +der+1+ +deg+1+) :element-type 'double-float))
             (dv (make-array (* +der+1+ +deg+1+) :element-type 'double-float)))
        (declare (type fixnum ku kv p+1 q+1 vdr+1) (dynamic-extent ndu du dv))
        (macrolet ((du (k j) `(aref du (+ (the fixnum (* ,k +deg+1+)) ,j)))
                   (dv (k j) `(aref dv (+ (the fixnum (* ,k +deg+1+)) ,j)))
                   (ndu (r c) `(aref ndu (+ (the fixnum (* ,r +deg+1+)) ,c)))
                   (bd (k l i j) `(aref bd (+ (the fixnum (* (+ (the fixnum (* (+ (the fixnum (* ,k vdr+1)) ,l) p+1)) ,i) q+1)) ,j))))
          ;; u direction
          (let ((mder (nl-min p udr)))
            (declare (type fixnum mder))
            (do ((k (1+ p) (1+ k))) ((> k udr)) (declare (type fixnum k))
              (dotimes (j p+1) (setf (du k j) 0d0)))
            (srf-nurb--st-compute-basis-functions p u knu ku ndu)
            (dotimes (j p+1) (setf (du 0 j) (ndu j p)))
            (srf-nurb--st-compute-derivatives p mder ndu du)
            (let ((r p))
              (declare (type fixnum r))
              (do ((k 1 (1+ k))) ((> k mder)) (declare (type fixnum k))
                (dotimes (j p+1) (setf (du k j) (* (du k j) (float r 1d0))))
                (setf r (* r (- p k))))))
          ;; v direction
          (setf kv (basis-find-span knv q v vfl))
          (let ((mder (nl-min q vdr)))
            (declare (type fixnum mder))
            (do ((k (1+ q) (1+ k))) ((> k vdr)) (declare (type fixnum k))
              (dotimes (j q+1) (setf (dv k j) 0d0)))
            (srf-nurb--st-compute-basis-functions q v knv kv ndu)
            (dotimes (j q+1) (setf (dv 0 j) (ndu j q)))
            (srf-nurb--st-compute-derivatives q mder ndu dv)
            (let ((r q))
              (declare (type fixnum r))
              (do ((k 1 (1+ k))) ((> k mder)) (declare (type fixnum k))
                (dotimes (j q+1) (setf (dv k j) (* (dv k j) (float r 1d0))))
                (setf r (* r (- q k))))))
          (let ((upper-half (and (= mfl +true+) (= udr vdr))))
            (if (/= (aref pw 0 0 3) +no-w+)
                ;; rational
                (let ((tu (make-array +deg+1+ :element-type 'double-float))
                      (d (make-array (* +der+1+ +der+1+) :element-type 'double-float))
                      (s1 (- ku p)) (s2 (- kv q)))
                  (declare (dynamic-extent tu d) (type fixnum s1 s2))
                  (macrolet ((d (k l) `(aref d (+ (the fixnum (* ,k +der+1+)) ,l))))
                    (dotimes (l (1+ vdr))
                      (dotimes (i p+1)
                        (setf (aref tu i) 0d0)
                        (dotimes (j q+1)
                          (incf (aref tu i) (* (aref pw (+ s1 i) (+ s2 j) 3) (dv l j)))))
                      (dotimes (k (1+ udr))
                        (setf (d k l) 0d0)
                        (dotimes (i p+1)
                          (incf (d k l) (* (aref tu i) (du k i))))))
                    (when (float-op-is-bad 1d0 (d 0 0) +division+)
                      (nurbs-error :con-err 'srf-basis-derivs))
                    (let ((inv-w00 (/ 1d0 (d 0 0))))
                      (declare (type double-float inv-w00))
                      (dotimes (n p+1)
                        (dotimes (m q+1)
                          (let ((w (aref pw (+ s1 n) (+ s2 m) 3)))
                            (declare (type double-float w))
                            (dotimes (k (1+ udr))
                              (dotimes (l (1+ (if upper-half (- udr k) vdr)))
                                (let ((v1 (* (* w (du k n)) (dv l m))))
                                  (declare (type double-float v1))
                                  (do ((j 1 (1+ j))) ((> j l)) (declare (type fixnum j))
                                    (decf v1 (* (* (float (aref +pascal-tri+ l j) 1d0) (d 0 j))
                                                (bd k (- l j) n m))))
                                  (do ((i 1 (1+ i))) ((> i k)) (declare (type fixnum i))
                                    (decf v1 (* (* (float (aref +pascal-tri+ k i) 1d0) (d i 0))
                                                (bd (- k i) l n m)))
                                    (let ((v2 0d0))
                                      (declare (type double-float v2))
                                      (do ((j 1 (1+ j))) ((> j l)) (declare (type fixnum j))
                                        (incf v2 (* (* (float (aref +pascal-tri+ l j) 1d0) (d i j))
                                                    (bd (- k i) (- l j) n m))))
                                      (decf v1 (* (float (aref +pascal-tri+ k i) 1d0) v2))))
                                  (setf (bd k l n m) (* v1 inv-w00)))))))))))
                ;; polynomial
                (dotimes (k (1+ udr))
                  (dotimes (l (1+ (if upper-half (- udr k) vdr)))
                    (dotimes (i p+1)
                      (dotimes (j q+1)
                        (setf (bd k l i j) (* (du k i) (dv l j)))))))))
          (values ku kv))))))

(defun srf-basis-derivs (sur u v ufl vfl mfl udr vdr bd)
  "Fill BD, a basis-derivatives array of at least (1+ UDR) (1+ VDR) (1+ p) (1+ q),
with the non-vanishing (rational) basis functions of surface SUR at (U, V) and
their derivatives: BD[k][l][i][j] is the (k,l)th derivative of basis function
(usp-p+i, vsp-q+j).  With MFL +true+ and UDR = VDR only the upper half
k + l <= UDR is computed.  Returns (values usp vsp).  NLib: N_SrfBasisDerivs."
  (declare (type surface sur) (type basis-derivatives bd))
  (let* ((p+1 (1+ (surface-degree-u sur))) (q+1 (1+ (surface-degree-v sur)))
         (flat (make-real-vector (* (1+ udr) (1+ vdr) p+1 q+1))))
    (multiple-value-bind (usp vsp) (%srf-basis-derivs sur u v ufl vfl mfl udr vdr flat)
      (dotimes (k (1+ udr))
        (dotimes (l (1+ vdr))
          (dotimes (i p+1)
            (dotimes (j q+1)
              (setf (aref bd k l i j)
                    (aref flat (+ (* (+ (* (+ (* k (1+ vdr)) l) p+1) i) q+1) j)))))))
      (values usp vsp))))

(defun compute-derivative-matrix-spec (no-w no-z upper p q udr vdr pw l2 l3 bd sd)
  "Fill SD[k][l] (a point-net) with the surface derivatives from the flat basis
derivatives BD (dimensions (1+ UDR) (1+ VDR) (1+ P) (1+ Q)) and control points
PW from row L2, column L3.  NO-W and NO-Z say whether the points lack w and z;
with UPPER only k + l <= UDR.  NLib: N_ComputeDerivativeMatrixSpec."
  (declare (type fixnum p q udr vdr l2 l3) (type cpoint-net pw) (type real-vector bd)
           (type point-net sd))
  (fast
    (let ((pwx (make-array (* +deg+1+ +deg+1+) :element-type 'double-float))
          (pwy (make-array (* +deg+1+ +deg+1+) :element-type 'double-float))
          (pwz (make-array (* +deg+1+ +deg+1+) :element-type 'double-float))
          (p+1 (1+ p)) (q+1 (1+ q)) (vdr+1 (1+ vdr)))
      (declare (dynamic-extent pwx pwy pwz) (type fixnum p+1 q+1 vdr+1))
      (dotimes (i p+1)
        (dotimes (j q+1)
          (let ((wx (aref pw (+ i l2) (+ j l3) 0)) (wy (aref pw (+ i l2) (+ j l3) 1))
                (wz (aref pw (+ i l2) (+ j l3) 2)) (w (aref pw (+ i l2) (+ j l3) 3))
                (ix (+ (the fixnum (* i +deg+1+)) j)))
            (declare (type double-float wx wy wz w) (type fixnum ix))
            (if (not no-w)
                (let ((r (/ 1d0 w)))
                  (setf (aref pwx ix) (* wx r) (aref pwy ix) (* wy r)
                        (aref pwz ix) (if (not no-z) (* wz r) 0d0)))
                (setf (aref pwx ix) wx (aref pwy ix) wy
                      (aref pwz ix) (if (not no-z) wz 0d0))))))
      (dotimes (k (1+ udr))
        (dotimes (l (1+ (if upper (- udr k) vdr)))
          (let ((sx 0d0) (sy 0d0) (sz 0d0))
            (declare (type double-float sx sy sz))
            (dotimes (i p+1)
              (dotimes (j q+1)
                (let ((alpha (aref bd (+ (the fixnum (* (+ (the fixnum (* (+ (the fixnum (* k vdr+1)) l) p+1)) i) q+1)) j)))
                      (ix (+ (the fixnum (* i +deg+1+)) j)))
                  (declare (type double-float alpha) (type fixnum ix))
                  (setf sx (+ (* alpha (aref pwx ix)) sx)
                        sy (+ (* alpha (aref pwy ix)) sy))
                  (unless no-z
                    (setf sz (+ (* alpha (aref pwz ix)) sz))))))
            (setf (aref sd k l 0) sx (aref sd k l 1) sy (aref sd k l 2) sz))))
      sd)))

(defun compute-derivative-matrix (no-w no-z upper p q udr vdr pw l2 l3 bd sd)
  "NLib: N_ComputeDerivativeMatrix, which dispatches to the specialized one."
  (compute-derivative-matrix-spec (and no-w t) (and no-z t) (and upper t)
                                  p q udr vdr pw l2 l3 bd sd))

(defun srf-derivs (sur u v ufl vfl mfl udr vdr sd)
  "Write the point of surface SUR at (U, V) and its partial derivatives into
SD, a point-net of at least (1+ UDR) by (1+ VDR): SD[k][l] is the derivative
k times in u and l times in v.  With MFL +true+ and UDR = VDR only k + l <= UDR
is computed.  Returns SD.  NLib: N_SrfDerivs."
  (declare (type surface sur) (type double-float u v) (type fixnum ufl vfl mfl udr vdr)
           (type point-net sd))
  (let* ((p (surface-degree-u sur)) (q (surface-degree-v sur))
         (pw (surface-cpoints sur))
         (bd (make-array (* (1+ udr) (1+ vdr) (1+ p) (1+ q)) :element-type 'double-float)))
    (declare (type fixnum p q) (dynamic-extent bd))
    (multiple-value-bind (usp vsp) (%srf-basis-derivs sur u v ufl vfl mfl udr vdr bd)
      (declare (type fixnum usp vsp))
      (compute-derivative-matrix (= (aref pw 0 0 3) +no-w+) (= (aref pw 0 0 2) +no-z+)
                                 (and (= mfl +true+) (= udr vdr))
                                 p q udr vdr pw (- usp p) (- vsp q) bd sd))))

(defun srf-eval-pt (sur u v ufl vfl s)
  "Write the point of surface SUR at (U, V) into point S and return S.
NLib: N_SrfEvalPt (with N_SrfEvalPtImpl and NI_VectorBlendCPt_specialized)."
  (declare (type surface sur) (type double-float u v) (type fixnum ufl vfl) (type point s))
  (let* ((p (surface-degree-u sur)) (q (surface-degree-v sur))
         (knu (surface-knots-u sur)) (knv (surface-knots-v sur))
         (pw (surface-cpoints sur))
         (noz (= (aref pw 0 0 2) +no-z+))
         (now (= (aref pw 0 0 3) +no-w+)))
    (declare (type fixnum p q) (type real-vector knu knv) (type cpoint-net pw))
    (fast
      (let ((nu (make-array +deg+1+ :element-type 'double-float))
            (nv (make-array +deg+1+ :element-type 'double-float)))
        (declare (dynamic-extent nu nv))
        (let* ((usp (basis-eval knu p u ufl nu))
               (vsp (basis-eval knv q v vfl nv))
               (tmp1 (- usp p)) (tmp2 (- vsp q))
               (sx 0d0) (sy 0d0) (sz 0d0) (sw 0d0))
          (declare (type fixnum usp vsp tmp1 tmp2) (type double-float sx sy sz sw))
          (dotimes (i (1+ p))
            (let ((tx 0d0) (ty 0d0) (tz 0d0) (tw 0d0))
              (declare (type double-float tx ty tz tw))
              (dotimes (j (1+ q))
                (let ((alpha (aref nv j)) (r (+ tmp1 i)) (c (+ tmp2 j)))
                  (declare (type double-float alpha) (type fixnum r c))
                  (setf tx (+ tx (* alpha (aref pw r c 0)))
                        ty (+ ty (* alpha (aref pw r c 1)))
                        tz (if noz +no-z+ (+ tz (* alpha (aref pw r c 2))))
                        tw (if now +no-w+ (+ tw (* alpha (aref pw r c 3)))))))
              (let ((alpha (aref nu i)))
                (declare (type double-float alpha))
                (setf sx (+ sx (* alpha tx))
                      sy (+ sy (* alpha ty))
                      sz (if noz +no-z+ (+ sz (* alpha tz)))
                      sw (if now +no-w+ (+ sw (* alpha tw)))))))
          (if (not now)
              (let ((r (/ 1d0 sw)))
                (set-point s (* sx r) (* sy r) (if (not noz) (* sz r) 0d0)))
              (set-point s sx sy (if (not noz) sz 0d0))))))))

(defun srf-eval-pt-normal-at-pole (sur u v ufl vfl p su sv n sd)
  "Point, unit derivatives and unit normal of SUR at (U, V) where a first
derivative vanishes (a pole): the normal comes from the mixed derivative.
SD is a 2 by 2 point-net of scratch.  Signals :geo-err when both first
derivatives vanish or no normal exists.  NLib: N_SrfEvalPtNormalAtPole."
  (declare (type surface sur) (type double-float u v) (type point p su sv n) (type point-net sd))
  (let ((pole 0)
        (sd01 (make-point)) (sd10 (make-point)) (sd11 (make-point)))
    (declare (type fixnum pole))
    (srf-derivs sur u v ufl vfl +false+ 1 1 sd)
    (flet ((row (k l dst) (set-point dst (aref sd k l 0) (aref sd k l 1) (aref sd k l 2)))
           (normalize-ok (vec) (handler-case (progn (vector-normalize-ref vec) t)
                                 (nurbs-error () nil))))
      (row 0 0 p) (row 1 0 su) (row 0 1 sv)
      (row 0 1 sd01) (row 1 0 sd10) (row 1 1 sd11)
      (let* ((mag (vector-magnitude su)) (ok (normalize-ok su)))
        (when (or (<= mag (* 0.001d0 *mtol*)) (not ok))
          (setf pole 1)
          (set-point su 0d0 0d0 0d0)))
      (let* ((mag (vector-magnitude sv)) (ok (normalize-ok sv)))
        (when (or (<= mag (* 0.001d0 *mtol*)) (not ok))
          (set-point sv 0d0 0d0 0d0)
          (if (= pole 1)
              (nurbs-error :geo-err 'srf-eval-pt-normal-at-pole)
              (setf pole 2))))
      (case pole
        (0 (vector-cross su sv n))
        (1 (multiple-value-bind (us ue vs ve) (srf-get-parameter-bounds sur)
             (declare (ignore us ue vs))
             (when (= v ve) (vector-scale sd01 -1d0 sd01))
             (vector-cross sd11 sd01 n)))
        (2 (multiple-value-bind (us ue vs ve) (srf-get-parameter-bounds sur)
             (declare (ignore us vs ve))
             (when (= u ue) (vector-scale sd10 -1d0 sd10))
             (vector-cross sd10 sd11 n))))
      (when (< (vector-magnitude n) (* 0.001d0 *mtol*))
        (nurbs-error :geo-err 'srf-eval-pt-normal-at-pole))
      (vector-normalize-ref n)
      p)))

(defun srf-eval-pt-pt-deriv-normal-pole (sur u v ufl vfl p su sv n)
  "NLib: N_SrfEvalPtPtDerivNormalPole."
  (srf-eval-pt-normal-at-pole sur u v ufl vfl p su sv n (make-point-net 2 2)))

(defun srf-eval-pt-pt-deriv-normal-fast (sur u v ufl vfl p su sv n sd)
  "Write the point of SUR at (U, V) into P, the unit first derivatives into SU
and SV, and the unit normal into N, falling back to the pole routine where a
derivative or the normal vanishes.  SD is a 2 by 2 point-net of scratch.
NLib: N_SrfEvalPtPtDerivNormalFast."
  (declare (type surface sur) (type double-float u v) (type point p su sv n) (type point-net sd))
  (set-point p 0d0 0d0 0d0) (set-point su 0d0 0d0 0d0)
  (set-point sv 0d0 0d0 0d0) (set-point n 0d0 0d0 0d0)
  (srf-derivs sur u v ufl vfl +true+ 1 1 sd)
  (set-point p (aref sd 0 0 0) (aref sd 0 0 1) (aref sd 0 0 2))
  (set-point su (aref sd 1 0 0) (aref sd 1 0 1) (aref sd 1 0 2))
  (set-point sv (aref sd 0 1 0) (aref sd 0 1 1) (aref sd 0 1 2))
  (when (<= (vector-magnitude su) (* 0.001d0 *mtol*))
    (return-from srf-eval-pt-pt-deriv-normal-fast
      (srf-eval-pt-pt-deriv-normal-pole sur u v ufl vfl p su sv n)))
  (vector-normalize-ref su)
  (when (<= (vector-magnitude sv) (* 0.001d0 *mtol*))
    (return-from srf-eval-pt-pt-deriv-normal-fast
      (srf-eval-pt-pt-deriv-normal-pole sur u v ufl vfl p su sv n)))
  (vector-normalize-ref sv)
  (vector-cross su sv n)
  (when (<= (vector-magnitude n) *ltol*)
    (return-from srf-eval-pt-pt-deriv-normal-fast
      (srf-eval-pt-pt-deriv-normal-pole sur u v ufl vfl p su sv n)))
  (vector-normalize-ref n)
  p)

(defun srf-eval-pt-pt-deriv-normal (sur u v ufl vfl p su sv n)
  "Point, unit first derivatives and unit normal of SUR at (U, V) into P, SU,
SV and N.  NLib: N_SrfEvalPtPtDerivNormal."
  (srf-eval-pt-pt-deriv-normal-fast sur u v ufl vfl p su sv n (make-point-net 2 2)))

(defun srf-min-max-weight-pos-vectors (sur)
  "The smallest and largest control point weights of SUR (both +no-w+ when it
is not rational) and the smallest and largest magnitudes of its Euclidean
control points, as (values wmin wmax pmin pmax).
NLib: N_SrfMinMaxWeightPosVectors."
  (declare (type surface sur))
  (let* ((pw (surface-cpoints sur))
         (n (surface-n sur)) (m (surface-m sur))
         (wmin 0d0) (wmax 0d0) (pmin 0d0) (pmax 0d0)
         (cw (make-cpoint)) (p (make-point)))
    (declare (type cpoint-net pw) (type fixnum n m)
             (type double-float wmin wmax pmin pmax) (dynamic-extent cw p))
    (flet ((magnitude (i j)
             (setf (aref cw 0) (aref pw i j 0) (aref cw 1) (aref pw i j 1)
                   (aref cw 2) (aref pw i j 2) (aref cw 3) (aref pw i j 3))
             (c-pt-to-pt-euclid cw p)
             (vector-magnitude p)))
      (declare (inline magnitude))
      (if (is-srf-rat sur)
          (progn
            (setf wmin (aref pw 0 0 3) wmax (aref pw 0 0 3))
            (dotimes (i (1+ n))
              (dotimes (j (1+ m))
                (let ((w (aref pw i j 3)))
                  (declare (type double-float w))
                  (when (< w wmin) (setf wmin w))
                  (when (> w wmax) (setf wmax w))))))
          (setf wmin +no-w+ wmax +no-w+))
      (let ((mag (magnitude 0 0)))
        (declare (type double-float mag))
        (setf pmin mag pmax mag))
      (dotimes (i (1+ n))
        (dotimes (j (1+ m))
          (let ((mag (magnitude i j)))
            (declare (type double-float mag))
            (when (< mag pmin) (setf pmin mag))
            (when (> mag pmax) (setf pmax mag)))))
      (values wmin wmax pmin pmax))))

(defun srf-translate (sur tv)
  "Translate surface SUR in place by vector TV.  Returns SUR.
NLib: N_SrfTranslate."
  (declare (type surface sur) (type point tv))
  (let ((pw (surface-cpoints sur))
        (cw (make-cpoint)))
    (declare (type cpoint-net pw) (dynamic-extent cw))
    (dotimes (i (1+ (surface-n sur)))
      (dotimes (j (1+ (surface-m sur)))
        (setf (aref cw 0) (aref pw i j 0) (aref cw 1) (aref pw i j 1)
              (aref cw 2) (aref pw i j 2) (aref cw 3) (aref pw i j 3))
        (translate-c-pt cw tv cw)
        (setf (aref pw i j 0) (aref cw 0) (aref pw i j 1) (aref cw 1)
              (aref pw i j 2) (aref cw 2) (aref pw i j 3) (aref cw 3))))
    sur))

;;; ------------------------------------------------------------ level 2, group 2

(defun srf-get-b-box (sur box)
  "Write the min-max box of the Euclidean control net of surface SUR into BOX
(a real-vector of 6: xl xr yb yt zn zf) and return BOX.
NLib: N_SrfGetBBox."
  (declare (type surface sur) (type real-vector box))
  (let ((ntl (srf-get-e-net sur 0 (surface-n sur) 0 (surface-m sur))))
    (declare (type point-net ntl))
    (e-net-get-b-box ntl box)))

(defun srf-reparam (sur us ue vs ve)
  "Rescale the knot vectors of SUR in place to the rectangle [US,UE]x[VS,VE];
a direction whose end parameters already match is left alone.  Returns SUR.
NLib: N_SrfReparam."
  (declare (type surface sur) (type double-float us ue vs ve))
  (let ((r (make-real-vector 4)))
    (declare (type real-vector r))
    (create-rectangle r us ue vs ve)
    (multiple-value-bind (u1 u2 v1 v2) (srf-get-parameter-bounds sur)
      (declare (type double-float u1 u2 v1 v2))
      (let ((scaledirs (+ (if (or (/= us u1) (/= ue u2)) +udir+ 0)
                          (if (or (/= vs v1) (/= ve v2)) +vdir+ 0))))
        (declare (type fixnum scaledirs))
        (srf-reparam-to-interval sur r scaledirs)))))
