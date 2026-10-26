;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvNurb.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp; results written into caller arrays
;;;; as the C writes its NL_POINT output parameters.

(in-package :batten)

;;; Maybe-inline: the expansion is kept, but callers inline it only where
;;; they ask (the exported curve-point! and curve-derivatives!), so the API
;;; passes U unboxed without the port growing a copy at every call site.
(declaim (inline crv-eval crv-derivs))
(defun crv-eval (cur u flg c)
  "Write the point of curve CUR at parameter U into point C and return C.
FLG (+left+ or +right+) picks the span at a knot.  NLib: N_CrvEval."
  (declare (type curve cur) (type double-float u) (type fixnum flg) (type point c))
  (let* ((knots (curve-knots cur))
         (pw (curve-cpoints cur))
         (p (curve-degree cur)))
    (declare (type real-vector knots) (type cpoint-array pw) (type fixnum p))
    (knot-vector-is-param-out-of-bounds knots u 'crv-eval)
    (fast
      (let ((n (make-array +deg+1+ :element-type 'double-float)))
        (declare (dynamic-extent n))
        (let* ((spn (basis-eval knots p u flg n))
               (j (- spn p))
               (cx 0d0) (cy 0d0) (cz 0d0) (cw 0d0))
          (declare (type fixnum spn j) (type double-float cx cy cz cw))
          (do ((i j (1+ i))) ((> i spn))
            (declare (type fixnum i))
            (let ((alpha (aref n (- i j)))
                  (z (aref pw i 2))
                  (w (aref pw i 3)))
              (declare (type double-float alpha z w))
              (setf cx (+ cx (* alpha (aref pw i 0)))
                    cy (+ cy (* alpha (aref pw i 1)))
                    cz (if (/= z +no-z+) (+ cz (* alpha z)) +no-z+)
                    cw (if (/= w +no-w+) (+ cw (* alpha w)) +no-w+))))
          (with-euclid ((x y z) cx cy cz cw)
            (set-point c x y z)))))))

(defun crv-derivs (cur u flg der cd)
  "Write the point of curve CUR at U and its derivatives up to order DER into
CD, a point-array of at least (1+ DER) rows: CD[k] is the kth derivative.
Return CD.  DER is at most +maxder+.  NLib: N_CrvDerivs."
  (declare (type curve cur) (type double-float u) (type fixnum flg der) (type point-array cd))
  (let* ((p (curve-degree cur))
         (knots (curve-knots cur))
         (pw (curve-cpoints cur))
         (m (1- (length knots))))
    (declare (type fixnum p m) (type real-vector knots) (type cpoint-array pw))
    (when (or (< u (aref knots 0)) (> u (aref knots m)))
      (nurbs-error :par-err 'crv-derivs))
    (when (> der +maxder+)
      (nurbs-error :mxd-err 'crv-derivs))
    (fast
      (let* ((nsp (basis-find-span knots p u flg))
             (mder (min p der))
             (p+1 (1+ p))
             (der+1 (1+ der))
             ;; constant sizes (NL_MAXDEG bounds p, NL_MAXDER bounds der, both
             ;; checked above): stack-allocated on every implementation; the
             ;; macrolets below index with the true strides
             (bd (make-array +der+1-deg+1+ :element-type 'double-float))
             (ndu (make-array +deg+1-squared+ :element-type 'double-float))
             (left (make-array +deg+2+ :element-type 'double-float))
             (right (make-array +deg+2+ :element-type 'double-float))
             (a (make-array +2deg+2+ :element-type 'double-float))
             (dd (make-array +der+1+ :element-type 'double-float))
             (px (make-array +deg+1+ :element-type 'double-float))
             (py (make-array +deg+1+ :element-type 'double-float))
             (pz (make-array +deg+1+ :element-type 'double-float)))
        (declare (type fixnum nsp mder p+1 der+1)
                 (dynamic-extent bd ndu left right a dd px py pz))
        (macrolet ((bd (k j) `(aref bd (+ (the fixnum (* ,k p+1)) ,j)))
                   (ndu (r c) `(aref ndu (+ (the fixnum (* ,r p+1)) ,c)))
                   (a (r c) `(aref a (+ (the fixnum (* ,r p+1)) ,c))))
          (do ((k (1+ p) (1+ k))) ((> k der))
            (declare (type fixnum k))
            (dotimes (j p+1) (setf (bd k j) 0d0)))
          (setf (ndu 0 0) 1d0)
          (do ((j 1 (1+ j))) ((> j p))
            (declare (type fixnum j))
            (setf (aref left j) (- u (aref knots (- (+ nsp 1) j)))
                  (aref right j) (- (aref knots (+ nsp j)) u))
            (let ((saved 0d0))
              (declare (type double-float saved))
              (dotimes (r j)
                (setf (ndu j r) (+ (aref right (1+ r)) (aref left (- j r))))
                (when (float-op-is-bad (ndu r (1- j)) (ndu j r) +division+)
                  (nurbs-error :num-err 'crv-derivs))
                (let ((temp (/ (ndu r (1- j)) (ndu j r))))
                  (declare (type double-float temp))
                  (setf (ndu r j) (+ saved (* (aref right (1+ r)) temp))
                        saved (* (aref left (- j r)) temp))))
              (setf (ndu j j) saved)))
          (dotimes (j p+1)
            (setf (bd 0 j) (ndu j p)))
          (dotimes (r p+1)
            (let ((s1 0) (s2 1))
              (declare (type fixnum s1 s2))
              (setf (a 0 0) 1d0)
              (do ((k 1 (1+ k))) ((> k mder))
                (declare (type fixnum k))
                (let* ((d 0d0)
                       (rk (- r k))
                       (pk (- p k))
                       (j1 (if (>= rk -1) 1 (- rk)))
                       (j2 (if (<= (1- r) pk) (1- k) (- p r))))
                  (declare (type double-float d) (type fixnum rk pk j1 j2))
                  (when (>= r k)
                    (setf (a s2 0) (/ (a s1 0) (ndu (1+ pk) rk))
                          d (* (a s2 0) (ndu rk pk))))
                  (do ((j j1 (1+ j))) ((> j j2))
                    (declare (type fixnum j))
                    (setf (a s2 j) (/ (- (a s1 j) (a s1 (1- j))) (ndu (1+ pk) (+ rk j))))
                    (incf d (* (a s2 j) (ndu (+ rk j) pk))))
                  (when (<= r pk)
                    (setf (a s2 k) (/ (- (a s1 (1- k))) (ndu (1+ pk) r)))
                    (incf d (* (a s2 k) (ndu r pk))))
                  (setf (bd k r) d)
                  (rotatef s1 s2)))))
          (let ((r p))
            (declare (type fixnum r))
            (do ((k 1 (1+ k))) ((> k mder))
              (declare (type fixnum k))
              (dotimes (j p+1)
                (setf (bd k j) (* (bd k j) (float r 1d0))))
              (setf r (* r (- p k)))))
          (when (/= (aref pw 0 3) +no-w+)
            (dotimes (k der+1)
              (setf (aref dd k) 0d0)
              (do ((j (- nsp p) (1+ j))) ((> j nsp))
                (declare (type fixnum j))
                (incf (aref dd k) (* (aref pw j 3) (bd k (+ (- j nsp) p))))))
            (dotimes (i p+1)
              (let ((w (aref pw (+ (- nsp p) i) 3)))
                (declare (type double-float w))
                (dotimes (k der+1)
                  (let ((d (* w (bd k i))))
                    (declare (type double-float d))
                    (do ((j 1 (1+ j))) ((> j k))
                      (declare (type fixnum j))
                      (decf d (* (* (float (aref +pascal-tri+ k j) 1d0) (aref dd j))
                                 (bd (- k j) i))))
                    (setf (bd k i) (/ d (aref dd 0))))))))
          (let ((j (- nsp p)))
            (declare (type fixnum j))
            (do ((i j (1+ i))) ((> i nsp))
              (declare (type fixnum i))
              (with-euclid-row ((x y z) pw i)
                (setf (aref px (- i j)) x (aref py (- i j)) y (aref pz (- i j)) z))))
          (dotimes (k der+1)
            (let ((cx 0d0) (cy 0d0) (cz 0d0))
              (declare (type double-float cx cy cz))
              (dotimes (i p+1)
                (let ((d (bd k i)))
                  (declare (type double-float d))
                  (incf cx (* d (aref px i)))
                  (incf cy (* d (aref py i)))
                  (incf cz (* d (aref pz i)))))
              (setf (aref cd k 0) cx (aref cd k 1) cy (aref cd k 2) cz)))
          cd)))))
(declaim (notinline crv-eval crv-derivs))

(defun crv-eval-tangent (cur u flg p tv)
  "Write the point of curve CUR at parameter U into point P and the unit
tangent there into vector TV; return P.  FLG (+left+ or +right+) picks the
span at a knot.  Signals as crv-derivs does, and :geo-err for a zero
first derivative.
NLib: N_CrvEvalTangent."
  (declare (type curve cur) (type double-float u) (type fixnum flg) (type point p tv))
  (let ((cd (make-point-array 2))
        (d1 (make-point)))
    (declare (dynamic-extent cd d1))
    (crv-derivs cur u flg 1 cd)
    (%row-pt cd 0 p)
    (vector-normalize (%row-pt cd 1 d1) tv)
    p))

(defun crv-get-min-max-weights-and-pts (cur)
  "The minimum and maximum weights of curve CUR (+no-w+ for both when it is
not rational) and the minimum and maximum distances of its control points
from their Euclidean centroid, as (values wmin wmax pmin pmax).  PMAX is
always 0, as NLib sets it (see the C header).
NLib: N_CrvGetMinMaxWeightsAndPts."
  (declare (type curve cur))
  (let* ((pw (curve-cpoints cur))
         (n (curve-n cur))
         (wmin 0d0) (wmax 0d0) (pmin 0d0) (pmax 0d0)
         (p (make-point))
         (c (make-point)))
    (declare (type cpoint-array pw) (type fixnum n)
             (type double-float wmin wmax pmin pmax) (dynamic-extent p c))
    (cond ((is-crv-rat cur)
           (setf wmin (aref pw 0 3)
                 wmax (aref pw 0 3))
           (do ((i 1 (1+ i))) ((> i n))
             (declare (type fixnum i))
             (let ((w (aref pw i 3)))
               (declare (type double-float w))
               (when (< w wmin) (setf wmin w))
               (when (> w wmax) (setf wmax w)))))
          (t (setf wmin +no-w+
                   wmax +no-w+)))
    (fast
      (dotimes (i (1+ n))
        (with-euclid-row ((x y z) pw i)
          (set-point p x y z))
        (sum-2-pts p c c))
      (scale-pt (/ 1d0 (+ (float n 1d0) 1d0)) c c)
      (dotimes (i (1+ n))
        (with-euclid-row ((x y z) pw i)
          (set-point p x y z))
        (diff-2-pts p c p)
        (let ((mag (vector-magnitude p)))
          (declare (type double-float mag))
          (when (or (= i 0) (< mag pmin)) (setf pmin mag))
          (when (or (= i 0) (> mag pmax)) (setf pmax mag)))))
    (setf pmax 0d0)
    (values wmin wmax pmin pmax)))

(defun crv-line-from-pt-and-vector (p v)
  "A fresh degree-1 NURBS line from point P to P+V, with knots (0 0 1 1).
NLib: N_CrvLineFromPtAndVector."
  (declare (type point p v))
  (let ((pw (make-cpoint-array 2))
        (u (make-real-vector 4))
        (tt (make-point)))
    (declare (type cpoint-array pw) (type real-vector u) (dynamic-extent tt))
    (vector-sum p v tt)
    (setf (aref pw 0 0) (aref p 0) (aref pw 0 1) (aref p 1) (aref pw 0 2) (aref p 2)
          (aref pw 0 3) +no-w+
          (aref pw 1 0) (aref tt 0) (aref pw 1 1) (aref tt 1) (aref pw 1 2) (aref tt 2)
          (aref pw 1 3) +no-w+)
    (setf (aref u 0) 0d0 (aref u 1) 0d0 (aref u 2) 1d0 (aref u 3) 1d0)
    (%make-curve 1 u pw)))

(defun crv-get-curvature-deriv (cur u ufl)
  "The derivative of the curvature of curve CUR at parameter U, UFL (+left+
or +right+) picking the span at a knot.  Zero where |C' x C''| < *mtol*.
Signals as crv-derivs does, and :num-err when the final division would
overflow.
NLib: N_CrvGetCurvatureDeriv."
  (declare (type curve cur) (type double-float u) (type fixnum ufl))
  (let ((d (make-point-array 4))
        (d1 (make-point)) (d2 (make-point)) (d3 (make-point))
        (a (make-point)) (z (make-point)))
    (declare (dynamic-extent d d1 d2 d3 a z))
    (crv-derivs cur u ufl 3 d)
    (%row-pt d 1 d1)
    (%row-pt d 2 d2)
    (%row-pt d 3 d3)
    (vector-cross d1 d2 z)
    (let ((f1 (vector-magnitude z))
          (f2 0d0) (num 0d0) (den 0d0))
      (declare (type double-float f1 f2 num den))
      (when (< f1 (the double-float *mtol*))
        (return-from crv-get-curvature-deriv 0d0))
      (vector-normalize-ref z)
      (setf den (vector-magnitude d1))
      (vector-cross d1 d3 a)
      (setf f2 (vector-dot a z))
      (setf f1 (* den den))
      (setf num (* f1 f2))
      (setf f1 (vector-dot d1 d2))
      (vector-cross d1 d2 a)
      (setf f2 (vector-dot a z))
      (setf num (- num (* (* 3d0 f1) f2)))
      (setf den (* (* (* (* den den) den) den) den))
      (when (float-op-is-bad num den +division+)
        (nurbs-error :num-err 'crv-get-curvature-deriv))
      (/ num den))))

;;; ------------------------------------------------------------ level 1, group 4

(defun %reverse-knots (u m v)
  "Write the knots U[0..M] reversed into V[0..M] (V may be U):
V[i] = U[0] + U[M] - U[M-i], except that a knot equal to an end value takes
the other end value exactly.  Returns V.  NLib's inline NI_ReverseKnots
(NLI_Math.h)."
  (declare (type real-vector u v) (type fixnum m))
  (fast
    (let* ((u0 (aref u 0)) (um (aref u m)) (c (+ u0 um)))
      (declare (type double-float u0 um c))
      (do ((i 0 (1+ i))) ((> i (truncate m 2)) v)
        (declare (type fixnum i))
        (let ((a (aref u i)) (b (aref u (- m i))))
          (declare (type double-float a b))
          (setf (aref v i) (if (= b um) u0 (if (= b u0) um (- c b))))
          (setf (aref v (- m i)) (if (= a um) u0 (if (= a u0) um (- c a)))))))))

(defun crv-reverse (cur-p &optional cur-q)
  "Reverse the direction of curve CUR-P: the same curve traced out the other
way.  With CUR-Q the very object CUR-P the reversal is done in place;
otherwise (CUR-Q nil, NLib's NULL curve) a fresh reversed curve is made.
Returns the reversed curve.
NLib: N_CrvReverse."
  (declare (type curve cur-p))
  (let* ((n (curve-n cur-p)) (m (curve-m cur-p)) (p (curve-degree cur-p))
         (pw (curve-cpoints cur-p)) (up (curve-knots cur-p))
         (aw (make-cpoint)) (bw (make-cpoint)))
    (declare (type fixnum n m p) (type cpoint-array pw) (type real-vector up)
             (dynamic-extent aw bw))
    (cond
      ((eq cur-q cur-p)
       (let ((k (truncate n 2)))
         (declare (type fixnum k))
         (do ((i 0 (1+ i))) ((> i k))
           (declare (type fixnum i))
           (%row-c-pt pw i aw)
           (%row-c-pt pw (- n i) bw)
           (swap-c-pts aw bw)
           (%c-pt-row aw pw i)
           (%c-pt-row bw pw (- n i))))
       (%reverse-knots up m up)
       cur-p)
      (t
       (let ((qw (make-cpoint-array (1+ n))) (uq (make-real-vector (1+ m))))
         (declare (type cpoint-array qw) (type real-vector uq))
         (dotimes (i (1+ n))
           (%row-c-pt pw (- n i) aw)
           (copy-c-pt aw bw)
           (%c-pt-row bw qw i))
         (%reverse-knots up m uq)
         (%make-curve p uq qw))))))

;;; ------------------------------------------------------------ level 2, group 5

(defun crv-get-b-box (cur box)
  "Write the min-max box of the Euclidean control polygon of curve CUR into
BOX (a real-vector of 6: xl xr yb yt zn zf) and return BOX.
NLib: N_CrvGetBBox."
  (declare (type curve cur) (type real-vector box))
  (let ((ppl (crv-get-e-polygon cur 0 (curve-n cur))))
    (polygon-get-b-box ppl box)))

(defun crv-eval-frenet-frame (cur u flg p tv nv bv)
  "Write the point of curve CUR at parameter U into P and the unit tangent,
normal and binormal there into TV, NV and BV; return P.  FLG (+left+ or
+right+) picks the span at a knot.  Signals :deg-err for a degree below 2,
as crv-derivs does, and :geo-err for a zero first derivative or a zero
C' x C''.
NLib: N_CrvEvalFrenetFrame."
  (declare (type curve cur) (type double-float u) (type fixnum flg) (type point p tv nv bv))
  (when (<= (curve-degree cur) 1)
    (nurbs-error :deg-err 'crv-eval-frenet-frame))
  (let ((cd (make-point-array 3))
        (d1 (make-point)) (d2 (make-point)))
    (declare (dynamic-extent cd d1 d2))
    (crv-derivs cur u flg 2 cd)
    (%row-pt cd 0 p)
    (%row-pt cd 1 d1)
    (%row-pt cd 2 d2)
    (vector-normalize d1 tv)
    (vector-cross d1 d2 bv)
    (vector-normalize-ref bv)
    (vector-cross-ref bv tv nv)
    p))

(defun crv-transform (cur rma)
  "Transform curve CUR in place by the 4x4 matrix RMA (an rmatrix), as
transform-c-pt does each control point; a 2D curve is made 3D first.
Return CUR.
NLib: N_CrvTransform."
  (declare (type curve cur) (type rmatrix rma))
  (unless (crv-is-3d cur)
    (crv-2d-to-3d cur))
  (let ((pw (curve-cpoints cur))
        (cw (make-cpoint)))
    (declare (type cpoint-array pw) (dynamic-extent cw))
    (dotimes (i (1+ (curve-n cur)))
      (%row-c-pt pw i cw)
      (transform-c-pt cw rma cw)
      (%c-pt-row cw pw i)))
  cur)

(defun crv-translate (cur tv)
  "Translate curve CUR in place by vector TV, as translate-c-pt does each
control point; a 2D curve is made 3D first.  Return CUR.
NLib: N_CrvTranslate."
  (declare (type curve cur) (type point tv))
  (unless (crv-is-3d cur)
    (crv-2d-to-3d cur))
  (let ((pw (curve-cpoints cur))
        (cw (make-cpoint)))
    (declare (type cpoint-array pw) (dynamic-extent cw))
    (dotimes (i (1+ (curve-n cur)))
      (%row-c-pt pw i cw)
      (translate-c-pt cw tv cw)
      (%c-pt-row cw pw i)))
  cur)

(defun crv-reparam (cur ul ur)
  "Map the knots of curve CUR, in place, linearly onto [UL, UR], as
crv-reparam-to-interval does.  Return CUR.
NLib: N_CrvReparam."
  (declare (type curve cur) (type double-float ul ur))
  (let ((i (make-real-vector 2)))
    (declare (dynamic-extent i))
    (create-interval i ul ur)
    (crv-reparam-to-interval cur i)))

;;; ------------------------------------------------------------- L3-g6

(defun crv-eval-curvature (cur u ufl p tv cfl c nv)
  "The curvature of curve CUR at parameter U (the osculating circle's radius
is its reciprocal).  UFL +left+ means u is in [u_j, u_j+1) (right
derivatives), +right+ u in (u_j, u_j+1] (left derivatives).  With CFL +yes+
the point, unit tangent, center of the osculating circle and normal to its
plane (the binormal) are written into P, TV, C and NV, from the Frenet
frame; otherwise only P and TV are written, from the derivatives.  Signals as
crv-derivs and crv-eval-frenet-frame do, :num-err when |C'|^3 or the
curvature cannot be divided by, and :geo-err for a zero first derivative.
Returns the curvature.
NLib: N_CrvEvalCurvature."
  (declare (type curve cur) (type double-float u) (type fixnum ufl cfl)
           (type point p tv c nv))
  (let ((d (make-point-array 3))
        (d0 (make-point)) (d1 (make-point)) (d2 (make-point))
        (a (make-point)) (b (make-point))
        (num 0d0) (den 0d0) (kap 0d0) (rad 0d0))
    (declare (type point-array d) (type point d0 d1 d2 a b)
             (type double-float num den kap rad)
             (dynamic-extent d d0 d1 d2 a b))
    ;; derivatives
    (crv-derivs cur u ufl 2 d)
    (%row-pt d 0 d0)
    (%row-pt d 1 d1)
    (%row-pt d 2 d2)
    ;; curvature
    (vector-cross d1 d2 a)
    (setf num (vector-magnitude a))
    (setf den (vector-magnitude d1))
    (setf den (* (* den den) den))
    (when (float-op-is-bad num den +division+)
      (nurbs-error :num-err 'crv-eval-curvature))
    (setf kap (/ num den))
    ;; osculating circle
    (cond ((= cfl +yes+)
           (when (float-op-is-bad 1d0 kap +division+)
             (nurbs-error :num-err 'crv-eval-curvature))
           (setf rad (/ 1d0 kap))
           ;; the C hands its B as the Frenet normal and N as the binormal
           (crv-eval-frenet-frame cur u ufl p tv b nv)
           (vector-pt-along-vector d0 rad b c))
          (t
           (copy-pt d0 p)
           (vector-normalize d1 tv)))
    kap))

(defun crv-get-b-box-max-diag-dist (cur)
  "The length of the diagonal of the min-max box of curve CUR's control
polygon.
NLib: N_CrvGetBBoxMaxDiagDist."
  (declare (type curve cur))
  (let ((box (make-real-vector 6)))
    (declare (dynamic-extent box))
    (crv-get-b-box cur box)
    (b-box-get-diagonal box)))

;;; ------------------------------------------------------------- L4-g1

(defun crv-get-type (cur tol aa bb xx yy)
  "Classify curve CUR as a line, a circular arc or a free-form curve, to
within TOL (the largest distance from a line, or the largest variation of the
radius and distance from the circle's plane).  For a line, AA gets its start
point and BB its direction vector; for a circle AA gets the center, BB the
unit normal of its plane and XX, YY the unit axes in that plane; otherwise the
point outputs are left alone.  Returns the radius, the sweep angle in degrees
(both 0 unless a circle) and the type, +nline+, +ncircle+ or +nfreeform+.
Errors in evaluating the curve make it free-form; signals as
project-pt-plane, vector-normalize and vectors-angle do.
NLib: N_CrvGetType."
  (declare (type curve cur) (type double-float tol) (type point aa bb xx yy))
  (let* ((pw (curve-cpoints cur))
         (u (curve-knots cur))
         (p (curve-degree cur))
         (n (curve-n cur))
         (m (curve-m cur))
         (lin +yes+)
         (typ +nfreeform+)
         (rr 0d0) (al 0d0)
         (a (make-point)) (b (make-point)) (tv (make-point)) (q (make-point))
         (rpt (make-point)) (v (make-point)) (ca (make-point)) (na (make-point))
         (cw (make-cpoint))
         (dab 0d0))
    (declare (type cpoint-array pw) (type real-vector u) (type fixnum p n m lin typ)
             (type double-float rr al dab)
             (type point a b tv q rpt v ca na) (type cpoint cw)
             (dynamic-extent a b tv q rpt v ca na cw))
    (block out
      ;; check if line
      (c-pt-to-pt-euclid (%row-c-pt pw 0 cw) a)
      (c-pt-to-pt-euclid (%row-c-pt pw n cw) b)
      (vector-dir a b tv)
      (setf dab (vector-magnitude tv))
      (cond
        ((and (<= dab *mtol*) (<= dab *ptol*))
         (setf lin +no+))
        ((> n 1)
         (let ((dab2 (* dab dab))
               (kh (* 2 p))
               (tl (aref u p))
               (i (1+ p)))
           (declare (type double-float dab2 tl) (type fixnum kh i))
           (block spans
             (do () ((>= i m))
               (do () ((not (and (< i m) (= (aref u i) (aref u (1+ i))))))
                 (incf i))
               (let* ((tr (aref u i))
                      (ti (/ (- tr tl) (float kh 1d0))))
                 (declare (type double-float tr ti))
                 (when (>= i m)
                   (decf kh))
                 (do ((k 1 (1+ k))) ((> k kh))
                   (declare (type fixnum k))
                   (let ((tt (if (and (= k kh) (< i m)) tr (+ tl (* (float k 1d0) ti))))
                         (d 0d0))
                     (declare (type double-float tt d))
                     (handler-case (crv-eval cur tt +left+ q)
                       (nurbs-error () (return-from out)))
                     ;; project Q onto line [A,B] to get R
                     (vector-dir a q v)
                     (setf d (vector-dot v tv))
                     (setf d (/ d dab2))
                     (when (or (< d 0d0) (> d 1d0))
                       (setf lin +no+)
                       (return))
                     (vector-pt-along-vector a d tv rpt)
                     (setf d (dist-pt-pt q rpt))
                     (when (> d tol)
                       (setf lin +no+)
                       (return))
                     (setf d (dist-pt-pt q a))
                     (when (> d dab)
                       (setf lin +no+)
                       (return))))
                 (when (= lin +no+)
                   (return-from spans))
                 (setf tl tr)
                 (incf i)))))))
      (when (= lin +yes+)
        (copy-pt a aa)
        (copy-pt tv bb)
        (setf typ +nline+)
        (return-from out))
      ;; check if circle
      (let* ((d (dist-pt-pt a b))
             (closed (< d tol))
             (i (* (* 2 n) (1+ p)))
             (rd (make-real-vector (1+ i)))
             (c (make-point-array (1+ i)))
             (nn (make-point-array (1+ i)))
             (pp (make-point-array (1+ i)))
             (cj (make-point)) (nj (make-point)) (pj (make-point))
             (ra 0d0)
             (kh (* 2 p))
             (j -1)
             (tl (aref u p))
             (tol2 (* 2d0 tol))
             (maxrad 0d0)
             (minrad +bigd+)
             (kap 0d0)
             (fac 0d0)
             (pln (%make-plane (make-point) (make-point))))
        (declare (type double-float d ra tl tol2 maxrad minrad kap fac)
                 (type fixnum i kh j) (type real-vector rd)
                 (type point-array c nn pp) (type point cj nj pj)
                 (dynamic-extent cj nj pj))
        (setf i (1+ p))
        (set-point ca 0d0 0d0 0d0)
        (set-point na 0d0 0d0 0d0)
        ;; for all knot spans
        (do () ((>= i m))
          (let ((kl (if (= i (1+ p)) 0 1))
                (kr 0)
                (tr 0d0) (ti 0d0))
            (declare (type fixnum kl kr) (type double-float tr ti))
            (do () ((not (and (< i m) (= (aref u i) (aref u (1+ i))))))
              (incf i))
            (setf kr (if (and closed (= i m)) (1- kh) kh))
            (setf tr (aref u i)
                  ti (/ (- tr tl) (float kh 1d0)))
            (do ((k kl (1+ k))) ((> k kr))
              (declare (type fixnum k))
              (let ((tt (if (= k kh) tr (+ tl (* (float k 1d0) ti)))))
                (declare (type double-float tt))
                (incf j)
                (handler-case
                    (setf kap (crv-eval-curvature cur tt +left+ pj tv +yes+ cj nj))
                  (nurbs-error () (return-from out)))
                (%pt-row pj pp j)
                (%pt-row cj c j)
                (%pt-row nj nn j)
                (when (float-op-is-bad 1d0 kap +division+)
                  (return-from out))
                (setf (aref rd j) (/ 1d0 kap))
                (setf ra (+ ra (aref rd j)))
                (when (> minrad (aref rd j))
                  (setf minrad (aref rd j)))
                (when (< maxrad (aref rd j))
                  (setf maxrad (aref rd j)))
                (when (> (- maxrad minrad) tol2)
                  (return-from out))
                (sum-2-pts ca cj ca)
                (sum-2-pts na nj na)))
            (setf tl tr)
            (incf i)))
        ;; the sums to averages
        (setf fac (/ 1d0 (+ (float j 1d0) 1d0)))
        (scale-pt fac ca ca)
        (scale-pt fac na na)
        (create-plane-pt-normal pln ca na)
        ;; the radius from the sample positions
        (setf ra 0d0 maxrad 0d0 minrad +bigd+)
        (dotimes (i (1+ j))
          (%row-pt pp i pj)
          (project-pt-plane pln pj q)
          (setf d (dist-pt-pt pj ca))
          (setf (aref rd i) d)
          (setf ra (+ ra (aref rd i)))
          (when (> minrad (aref rd i))
            (setf minrad (aref rd i)))
          (when (< maxrad (aref rd i))
            (setf maxrad (aref rd i)))
          (setf d (dist-pt-pt pj q))
          (when (> d tol)
            (return-from out)))
        (setf fac (/ 1d0 (+ (float j 1d0) 1d0)))
        (setf ra (* fac ra))
        (when (> (abs (- maxrad ra)) tol)
          (return-from out))
        (when (> (abs (- ra minrad)) tol)
          (return-from out))
        ;; circle X axis = normalized (first control point - center)
        (copy-pt ca aa)
        (copy-pt na bb)
        (diff-2-pts a ca tv)
        (vector-normalize tv xx)
        ;; circle Y axis = normal x X axis
        (vector-cross na xx yy)
        (if closed
            (setf al 360d0)
            (progn
              (diff-2-pts b ca tv)
              (setf fac (vectors-angle xx tv))
              (setf d (vector-dot yy tv))
              (setf al (if (< d 0d0) (- 360d0 fac) fac))))
        (setf rr ra
              typ +ncircle+)))
    (values rr al typ)))
