;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Bezier.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.  NL_RMATRIX outputs are real-matrix
;;;; arrays: passed in when the caller has one, else allocated and returned,
;;;; as N_CheckMemRealMatrix does.

(in-package :batten)

;;; ------------------------------------------------- row helpers (not NLib)

;;; NLib applies N_CopyCPt, N_Combine2CPts, N_CopyPt and N_Combine2Pts to
;;; elements of NL_CPOINT and NL_POINT arrays; these do the same to rows of
;;; cpoint-array and point-array, with the same arithmetic, reading every
;;; input coordinate before writing, as C's by-value arguments do.


(defun %bezier-check-mem-real-matrix (mat n m rname)
  "MAT, or a fresh zeroed (1+ N) by (1+ M) real-matrix when MAT is nil;
signal :sto-err in RNAME when MAT is smaller, as N_CheckMemRealMatrix does."
  (declare (type fixnum n m))
  (cond ((null mat) (make-real-matrix (1+ n) (1+ m)))
        ((or (< (1- (array-dimension mat 0)) n) (< (1- (array-dimension mat 1)) m))
         (nurbs-error :sto-err rname))
        (t mat)))

;;; ------------------------------------------------------- matrices

(defun bez-get-degree-elevation-matrix (p tt &optional dm)
  "The matrix elevating a degree-P Bezier entity by TT degrees: DM[i][j],
0 <= i <= p+tt, 0 <= j <= p, written where non-zero.  DM, a real-matrix of
at least (1+ p+tt) by (1+ p), is filled when given (signaling :sto-err when
it is smaller); otherwise a fresh zeroed one is.  Returns DM.
NLib: N_BezGetDegreeElevationMatrix."
  (declare (type fixnum p tt))
  (let* ((q (+ p tt))
         (rm (%bezier-check-mem-real-matrix dm q p 'bez-get-degree-elevation-matrix))
         (r (truncate q 2)))
    (declare (type fixnum q r) (type real-matrix rm))
    (fast
      (setf (aref rm 0 0) 1d0
            (aref rm q p) 1d0)
      (do ((i 1 (1+ i))) ((> i r))
        (declare (type fixnum i))
        (let ((inv (/ 1d0 (float (aref +pascal-tri+ q i) 1d0)))
              (k (if (> (- i tt) 0) (- i tt) 0))
              (l (if (> i p) p i)))
          (declare (type double-float inv) (type fixnum k l))
          (do ((j k (1+ j))) ((> j l))
            (declare (type fixnum j))
            (setf (aref rm i j) (* (* inv (float (aref +pascal-tri+ p j) 1d0))
                                   (float (aref +pascal-tri+ tt (- i j)) 1d0))))))
      (do ((i (1+ r) (1+ i))) ((>= i q))
        (declare (type fixnum i))
        (let ((k (if (> (- i tt) 0) (- i tt) 0))
              (l (if (> i p) p i)))
          (declare (type fixnum k l))
          (do ((j k (1+ j))) ((> j l))
            (declare (type fixnum j))
            (setf (aref rm i j) (aref rm (- q i) (- p j)))))))
    rm))

(defun bez-func-multiply-bez-crv-matrix (p q &optional pm)
  "The matrix for the product of a degree-P and a degree-Q Bezier entity:
PM[i][j], 0 <= i <= p+q, 0 <= j <= p, written where non-zero.  PM, a
real-matrix of at least (1+ p+q) by (1+ p), is filled when given (signaling
:sto-err when it is smaller); otherwise a fresh zeroed one is.  Returns PM.
NLib: N_BezFuncMultiplyBezCrvMatrix."
  (declare (type fixnum p q))
  (let* ((n (+ p q))
         (m p)
         (a (%bezier-check-mem-real-matrix pm n m 'bez-func-multiply-bez-crv-matrix))
         (r (truncate n 2)))
    (declare (type fixnum n m r) (type real-matrix a))
    (fast
      (setf (aref a 0 0) 1d0
            (aref a n m) 1d0)
      (do ((i 1 (1+ i))) ((> i r))
        (declare (type fixnum i))
        (let ((inv (/ 1d0 (float (aref +pascal-tri+ n i) 1d0)))
              (jl (nl-max 0 (- i q)))
              (jh (nl-min p i)))
          (declare (type double-float inv) (type fixnum jl jh))
          (do ((j jl (1+ j))) ((> j jh))
            (declare (type fixnum j))
            (setf (aref a i j) (* (* inv (float (aref +pascal-tri+ p j) 1d0))
                                  (float (aref +pascal-tri+ q (- i j)) 1d0))))))
      (do ((i (1+ r) (1+ i))) ((>= i n))
        (declare (type fixnum i))
        (let ((jl (nl-max 0 (- i q)))
              (jh (nl-min p i)))
          (declare (type fixnum jl jh))
          (do ((j jl (1+ j))) ((> j jh))
            (declare (type fixnum j))
            (setf (aref a i j) (aref a (- n i) (- m j)))))))
    a))

;;; ------------------------------------------------- degree elevation, products

(defun %bezier-elevation-rm (dm p tt rname)
  "The storage of NL_RMATRIX DM, after filling it with the degree elevation
matrix for degree P raised by TT when DM is the NULL matrix, as the C's
N_BezGetDegreeElevationMatrix( p, t, dm, &S ) does."
  (declare (type rmatrix dm) (type fixnum p tt))
  (when (real-matrix-is-null dm)
    (%check-mem-real-matrix dm (+ p tt) p +mt-full+ p rname)
    (bez-get-degree-elevation-matrix p tt (rmatrix-rm dm)))
  (rmatrix-rm dm))

(defun %bezier-product-rm (pm p q rname)
  "The storage of NL_RMATRIX PM, after filling it with the product matrix for
degrees P and Q when PM is the NULL matrix, as the C's
N_BezFuncMultiplyBezCrvMatrix( p, q, pm, &S ) does."
  (declare (type rmatrix pm) (type fixnum p q))
  (when (real-matrix-is-null pm)
    (%check-mem-real-matrix pm (+ p q) p +mt-full+ p rname)
    (bez-func-multiply-bez-crv-matrix p q (rmatrix-rm pm)))
  (rmatrix-rm pm))

(defun bez-elevate-degree (pw p tt dm f l qw)
  "Elevate the degree-P Bezier curve with control points PW by TT degrees,
writing new control points F..L into QW (a cpoint-array the caller
allocates).  DM is the degree elevation matrix, an rmatrix: when it is the
NULL matrix it is computed (and DM given storage).  Returns QW.
NLib: N_BezElevateDegree."
  (declare (type cpoint-array pw qw) (type fixnum p tt f l) (type rmatrix dm))
  (let ((rm (%bezier-elevation-rm dm p tt 'bez-elevate-degree)))
    (declare (type real-matrix rm))
    (fast
      (do ((i f (1+ i))) ((> i l))
        (declare (type fixnum i))
        (let ((a (if (> (- i tt) 0) (- i tt) 0))
              (b (if (> i p) p i)))
          (declare (type fixnum a b))
          (setf (aref qw i 0) 0d0 (aref qw i 1) 0d0 (aref qw i 2) 0d0 (aref qw i 3) 0d0)
          (do ((j a (1+ j))) ((> j b))
            (declare (type fixnum j))
            (%vector-blend-c-pt-row (aref rm i j) pw j qw i))))))
  qw)

(defun bez-func-degree-elevate (f p tt dm s e g)
  "Elevate the degree-P Bezier function with control values F by TT degrees,
writing new control values S..E into G (allocated by the caller).  DM is the
degree elevation matrix, an rmatrix computed (and given storage) when it is
the NULL matrix.  Signals :inp-err for S < 0, E > P+TT, TT < 0 or S > E and
:deg-err when P+TT exceeds *dmax*.  Returns G.
NLib: N_BezFuncDegreeElevate."
  (declare (type real-vector f g) (type fixnum p tt s e) (type rmatrix dm))
  (when (or (< s 0) (> e (+ p tt)))
    (nurbs-error :inp-err 'bez-func-degree-elevate))
  (when (or (< tt 0) (> s e))
    (nurbs-error :inp-err 'bez-func-degree-elevate))
  (when (> (+ p tt) (the fixnum *dmax*))
    (nurbs-error :deg-err 'bez-func-degree-elevate))
  (let ((a (%bezier-elevation-rm dm p tt 'bez-func-degree-elevate)))
    (declare (type real-matrix a))
    (fast
      (do ((i s (1+ i))) ((> i e))
        (declare (type fixnum i))
        (setf (aref g i) 0d0)
        (let ((jl (nl-max 0 (- i tt)))
              (jh (nl-min p i)))
          (declare (type fixnum jl jh))
          (do ((j jl (1+ j))) ((> j jh))
            (declare (type fixnum j))
            (setf (aref g i) (+ (aref g i) (* (aref a i j) (aref f j)))))))))
  g)

(defun bez-func-multiply-bez-func (f p g q pm s e fg)
  "The control values S..E of the product of the degree-P Bezier function F
and the degree-Q Bezier function G, written into FG (allocated by the caller
to hold up to FG[P+Q]).  PM is the product matrix, an rmatrix computed (and
given storage) when it is the NULL matrix.  Signals :inp-err unless
0 <= S <= E <= P+Q.  Returns FG.
NLib: N_BezFuncMultiplyBezFunc."
  (declare (type real-vector f g fg) (type fixnum p q s e) (type rmatrix pm))
  (when (or (< s 0) (> e (+ p q)) (> s e))
    (nurbs-error :inp-err 'bez-func-multiply-bez-func))
  (let ((a (%bezier-product-rm pm p q 'bez-func-multiply-bez-func)))
    (declare (type real-matrix a))
    (fast
      (do ((i s (1+ i))) ((> i e))
        (declare (type fixnum i))
        (setf (aref fg i) 0d0)
        (let ((jl (nl-max 0 (- i q)))
              (jh (nl-min p i)))
          (declare (type fixnum jl jh))
          (do ((j jl (1+ j))) ((> j jh))
            (declare (type fixnum j))
            (setf (aref fg i) (+ (aref fg i) (* (* (aref a i j) (aref f j)) (aref g (- i j))))))))))
  fg)

(defun bez-func-multiply-bez-crv (f p pw q pm s e qw)
  "The control points S..E of the product of the degree-P Bezier function F
and the degree-Q Bezier curve with control points PW, written into QW (a
cpoint-array allocated by the caller to hold up to QW[P+Q]).  The function
scales the curve's x, y, z only.  PM is the product matrix, an rmatrix
computed (and given storage) when it is the NULL matrix.  Signals :inp-err
unless 0 <= S <= E <= P+Q.  Returns QW.
NLib: N_BezFuncMultiplyBezCrv."
  (declare (type real-vector f) (type cpoint-array pw qw) (type fixnum p q s e) (type rmatrix pm))
  (when (or (< s 0) (> e (+ p q)) (> s e))
    (nurbs-error :inp-err 'bez-func-multiply-bez-crv))
  (let ((a (%bezier-product-rm pm p q 'bez-func-multiply-bez-crv))
        (tw (make-array 4 :element-type 'double-float :initial-element 0d0))
        (pt (make-array 4 :element-type 'double-float :initial-element 0d0))
        (acc (make-array 4 :element-type 'double-float :initial-element 0d0))
        (czero (make-array 4 :element-type 'double-float :initial-element 0d0)))
    (declare (type real-matrix a) (type cpoint tw pt acc czero) (dynamic-extent tw pt acc czero))
    (fast
      (do ((i s (1+ i))) ((> i e))
        (declare (type fixnum i))
        (copy-c-pt czero acc)
        (let ((jl (nl-max 0 (- i q)))
              (jh (nl-min p i)))
          (declare (type fixnum jl jh))
          (do ((j jl (1+ j))) ((> j jh))
            (declare (type fixnum j))
            (%row-c-pt pw (- i j) pt)
            (scale-c-pt-xyz (aref f j) pt tw)
            (vector-blend-c-pt (aref a i j) tw acc)))
        (%c-pt-row acc qw i))))
  qw)

(defun bez-srf-elevate-degree (pw r tt dm dir f l roc qw)
  "Elevate the degree R of a Bezier surface with control net PW by TT
degrees in direction DIR (+udir+ or +vdir+), one row or column at a time:
new control points F..L of column (u) or row (v) ROC are written into QW, a
cpoint-net the caller allocates.  DM is the degree elevation matrix, an
rmatrix computed (and given storage) when it is the NULL matrix.  Signals
:cal-err for another DIR.  Returns QW.
NLib: N_BezSrfElevateDegree."
  (declare (type cpoint-net pw qw) (type fixnum r tt dir f l roc) (type rmatrix dm))
  (let ((rm (%bezier-elevation-rm dm r tt 'bez-srf-elevate-degree))
        (pt (make-array 4 :element-type 'double-float :initial-element 0d0))
        (acc (make-array 4 :element-type 'double-float :initial-element 0d0))
        (czero (make-array 4 :element-type 'double-float :initial-element 0d0)))
    (declare (type real-matrix rm) (type cpoint pt acc czero) (dynamic-extent pt acc czero))
    (cond
      ((= dir +udir+)
       (fast
         (do ((i f (1+ i))) ((> i l))
           (declare (type fixnum i))
           (let ((a (if (> (- i tt) 0) (- i tt) 0))
                 (b (if (> i r) r i)))
             (declare (type fixnum a b))
             (copy-c-pt czero acc)
             (do ((k a (1+ k))) ((> k b))
               (declare (type fixnum k))
               (%net-c-pt pw k roc pt)
               (vector-blend-c-pt (aref rm i k) pt acc))
             (%c-pt-net acc qw i roc)))))
      ((= dir +vdir+)
       (fast
         (do ((j f (1+ j))) ((> j l))
           (declare (type fixnum j))
           (let ((a (if (> (- j tt) 0) (- j tt) 0))
                 (b (if (> j r) r j)))
             (declare (type fixnum a b))
             (copy-c-pt czero acc)
             (do ((k a (1+ k))) ((> k b))
               (declare (type fixnum k))
               (%net-c-pt pw roc k pt)
               (vector-blend-c-pt (aref rm j k) pt acc))
             (%c-pt-net acc qw roc j)))))
      (t (nurbs-error :cal-err 'bez-srf-elevate-degree))))
  qw)

;;; ---------------------------------------------------------------- conics

(defun bez-set-conic-weight (p0 p1 p2)
  "The middle weight of the conic arc on control triangle P0 P1 P2 that best
approximates the circular arc on it (exact for an isosceles triangle); 1 when
two vertices coincide.  NLib: N_BezSetConicWeight."
  (declare (type point p0 p1 p2))
  (let ((d (dist-pt-pt p0 p2))
        (fl (dist-pt-pt p0 p1))
        (fr (dist-pt-pt p2 p1)))
    (declare (type double-float d fl fr))
    (if (= (* (* d fl) fr) 0d0)
        1d0
        (let* ((d (* d 0.5d0))
               (sl (/ d (+ d fl)))
               (sr (/ d (+ d fr)))
               (s (* 0.5d0 (+ sl sr))))
          (declare (type double-float d sl sr s))
          (/ s (- 1d0 s))))))

(defun bez-conic-split (p wp q r)
  "Split the conic arc with unweighted control points P (a point-array of 3)
and middle weight WP (end weights 1) at its shoulder point: the halves'
control points go into Q and R (point-arrays of 3), and their middle weights
are returned as (values wq wr), in normal form.  NLib: N_BezConicSplit."
  (declare (type point-array p q r) (type double-float wp))
  (let ((wq 0d0) (wr 0d0))
    (declare (type double-float wq wr))
    (%copy-pt-row p 0 q 0)
    (%copy-pt-row p 2 r 2)
    (when (/= wp 0d0)
      (let* ((alf (/ 1d0 (+ 1d0 wp)))
             (oma (- 1d0 alf)))
        (declare (type double-float alf oma))
        (%combine-2-pt-rows alf p 0 oma p 1 q 1)
        (%combine-2-pt-rows alf p 2 oma p 1 r 1)
        (%combine-2-pt-rows 0.5d0 q 1 0.5d0 r 1 q 2)
        (%copy-pt-row q 2 r 0)
        (setf wq (c-sqrt (* 0.5d0 (+ 1d0 wp)))
              wr wq)))
    (when (= wp 0d0)
      (%combine-2-pt-rows 1d0 p 0 1d0 p 1 q 1)
      (%combine-2-pt-rows 1d0 p 2 1d0 p 1 r 1)
      (%combine-2-pt-rows 0.5d0 q 1 0.5d0 r 1 q 2)
      (%copy-pt-row q 2 r 0)
      (setf wq (* 0.5d0 (sqrt 2d0))
            wr wq))
    (values wq wr)))

;;; ------------------------------------------------------------- splitting

(defun bez-split (pw p u qw rw)
  "Split the degree-P Bezier curve with control points PW at U in [0, 1]
(signaling :par-err outside it): the left segment's p+1 control points go
into QW and the right's into RW, cpoint-arrays the caller allocates.
Returns QW.  NLib: N_BezSplit."
  (declare (type cpoint-array pw qw rw) (type fixnum p) (type double-float u))
  (when (or (< u 0d0) (> u 1d0))
    (nurbs-error :par-err 'bez-split))
  (fast
    (cond
      ((= u 0d0)
       (dotimes (i (1+ p))
         (%copy-c-pt-row pw 0 qw i)
         (%copy-c-pt-row pw i rw i)))
      ((= u 1d0)
       (dotimes (i (1+ p))
         (%copy-c-pt-row pw i qw i)
         (%copy-c-pt-row pw p rw i)))
      (t
       (let ((aw (make-array (list (1+ p) 4) :element-type 'double-float))
             (omu (- 1d0 u)))
         (declare (dynamic-extent aw) (type double-float omu))
         (dotimes (i (1+ p))
           (%copy-c-pt-row pw i aw i))
         (%copy-c-pt-row pw 0 qw 0)
         (%copy-c-pt-row pw p rw p)
         (do ((k 1 (1+ k))) ((> k p))
           (declare (type fixnum k))
           (do ((i 0 (1+ i))) ((> i (- p k)))
             (declare (type fixnum i))
             (%combine-2-c-pt-rows omu aw i u aw (1+ i) aw i))
           (%copy-c-pt-row aw 0 qw k)
           (%copy-c-pt-row aw (- p k) rw (- p k)))))))
  qw)

;;; ------------------------------------------------------- reparametrization

(defun bez-reparam (pw p f q bin qw)
  "Reparametrize the degree-P Bezier curve with control points PW by the
degree-Q Bezier function with coefficients F (monotonic for the result to
make sense), writing the P*Q+1 control points of the result into QW (a
cpoint-array the caller allocates).  BIN is Pascal's triangle through row
P*Q, an index-matrix, or NIL to compute it.  Signals :deg-err when P*Q
exceeds *dmax*.  Returns QW.
NLib: N_BezReparam."
  (declare (type cpoint-array pw qw) (type fixnum p q) (type real-vector f)
           (type (or null index-matrix) bin))
  (let ((pq (* p q)))
    (declare (type fixnum pq))
    (when (> pq (the fixnum *dmax*))
      (nurbs-error :deg-err 'bez-reparam))
    (let ((bin (or bin (pascal-tri-row (make-array (list (1+ pq) (1+ pq)) :element-type 'fixnum
                                                                           :initial-element 0)
                                       pq)))
          (uw (make-cpoint-net (1+ p) (1+ pq)))
          (vw (make-cpoint-net (1+ p) (1+ pq)))
          (uarray nil)
          (a0 (make-array 4 :element-type 'double-float :initial-element 0d0))
          (a1 (make-array 4 :element-type 'double-float :initial-element 0d0))
          (aw (make-array 4 :element-type 'double-float :initial-element 0d0))
          (acc (make-array 4 :element-type 'double-float :initial-element 0d0))
          (czero (make-array 4 :element-type 'double-float :initial-element 0d0)))
      (declare (type index-matrix bin) (type cpoint-net uw vw) (type cpoint a0 a1 aw acc czero)
               (dynamic-extent a0 a1 aw acc czero))
      (dotimes (i (1+ p))
        (dotimes (k 4) (setf (aref uw i 0 k) (aref pw i k))))
      (fast
        (do ((s 1 (1+ s))) ((> s p))
          (declare (type fixnum s))
          (let ((qs (* q s))
                (qs1 (* q (1- s)))
                ;; the pass reads one net and writes the other
                (src (if uarray vw uw))
                (dst (if uarray uw vw)))
            (declare (type fixnum qs qs1) (type cpoint-net src dst))
            (do ((r 0 (1+ r))) ((> r qs))
              (declare (type fixnum r))
              (do ((i 0 (1+ i))) ((> i (- p s)))
                (declare (type fixnum i))
                (copy-c-pt czero acc)
                (let ((jmin (nl-max 0 (- r q)))
                      (jmax (nl-min r qs1)))
                  (declare (type fixnum jmin jmax))
                  (do ((j jmin (1+ j))) ((> j jmax))
                    (declare (type fixnum j))
                    (let* ((alf (aref f (- r j)))
                           (oma (- 1d0 alf))
                           (fac (* (float (aref bin qs1 j) 1d0) (float (aref bin q (- r j)) 1d0))))
                      (declare (type double-float alf oma fac))
                      (%net-c-pt src i j a0)
                      (%net-c-pt src (1+ i) j a1)
                      (combine-2-c-pts oma a0 alf a1 aw)
                      (vector-blend-c-pt fac aw acc))))
                (let ((fac (/ 1d0 (float (aref bin qs r) 1d0))))
                  (declare (type double-float fac))
                  (scale-c-pt fac acc acc))
                (%c-pt-net acc dst i r))))
          (setf uarray (not uarray))))
      (let ((src (if uarray vw uw)))
        (declare (type cpoint-net src))
        (dotimes (r (1+ pq))
          (dotimes (k 4) (setf (aref qw r k) (aref src 0 r k)))))))
  qw)

;;; ------------------------------------------------------------- L2-g7

(defun bez-conic-arc-from-pts-and-tangents (p0 t0 p2 t2 p p1)
  "The middle control point and weight of the conic arc from P0 (tangent T0)
to P2 (tangent T2) through point P.  Writes the middle control point into P1
and returns the weight.  When the tangents are parallel the middle point is at
infinity: P1 receives its direction (a multiple of T0) and the weight is 0.
Signals as isect-line-line does, :cal-err when a needed intersection is
missing or its parameter makes the division fail, and :num-err when the
weight's division would fail.  NLib: N_BezConicArcFromPtsAndTangents."
  (declare (type point p0 t0 p2 t2 p p1))
  ;; the C takes the points by value: copy them, so P1 may share an input
  (let ((cp0 (make-point)) (ct0 (make-point)) (cp2 (make-point)) (ct2 (make-point))
        (cp (make-point)) (m (make-point))
        (v0 (make-point)) (v1 (make-point)) (v2 (make-point))
        (l0 (%make-lineseg (make-point) (make-point) 0))
        (l2 (%make-lineseg (make-point) (make-point) 0))
        (w1 0d0))
    (declare (type point cp0 ct0 cp2 ct2 cp m v0 v1 v2) (type lineseg l0 l2)
             (type double-float w1) (dynamic-extent cp0 ct0 cp2 ct2 cp m v0 v1 v2))
    (vector-copy p0 cp0)
    (vector-copy t0 ct0)
    (vector-copy p2 cp2)
    (vector-copy t2 ct2)
    (vector-copy p cp)
    ;; middle control point
    (create-line-start-dir-vector l0 cp0 ct0 +unbounded+)
    (create-line-start-dir-vector l2 cp2 ct2 +unbounded+)
    (let ((ints (isect-line-line l0 l2 p1)))
      (declare (type fixnum ints))
      ;; middle weight
      (when (= ints +true+)
        (create-line-pt-pt l0 cp0 cp2 +bounded+)
        (create-line-pt-pt l2 p1 cp +unbounded+)
        (multiple-value-bind (flg tt0) (isect-line-line l0 l2 m)
          (declare (type fixnum flg) (type double-float tt0))
          (setf ints flg)
          (when (= ints +false+)
            (nurbs-error :cal-err 'bez-conic-arc-from-pts-and-tangents))
          (when (float-op-is-bad tt0 (- 1d0 tt0) +division+)
            (nurbs-error :cal-err 'bez-conic-arc-from-pts-and-tangents))
          (let* ((a (c-sqrt (/ tt0 (- 1d0 tt0))))
                 (u (/ a (+ 1d0 a))))
            (declare (type double-float a u))
            (vector-diff cp cp0 v0)
            (vector-diff p1 cp v1)
            (vector-diff cp cp2 v2)
            (let* ((alf (vector-dot v0 v1))
                   (bet (vector-dot v1 v2))
                   (gam (vector-dot v1 v1))
                   (a (* (- 1d0 u) (- 1d0 u)))
                   (b (* u u))
                   (c (* (* 2d0 u) (- 1d0 u)))
                   (num (+ (* a alf) (* b bet)))
                   (den (* c gam)))
              (declare (type double-float alf bet gam a b c num den))
              (when (float-op-is-bad num den +division+)
                (nurbs-error :num-err 'bez-conic-arc-from-pts-and-tangents))
              (setf w1 (/ num den))))))
      ;; direction vector
      (when (= ints +false+)
        (create-line-start-dir-vector l0 cp ct0 +unbounded+)
        (create-line-pt-pt l2 cp0 cp2 +bounded+)
        (multiple-value-bind (flg tt0 tt2) (isect-line-line l0 l2 m)
          (declare (type fixnum flg) (type double-float tt0 tt2))
          (when (= flg +false+)
            (nurbs-error :cal-err 'bez-conic-arc-from-pts-and-tangents))
          (when (float-op-is-bad tt2 (- 1d0 tt2) +division+)
            (nurbs-error :cal-err 'bez-conic-arc-from-pts-and-tangents))
          (let* ((a (c-sqrt (/ tt2 (- 1d0 tt2))))
                 (u (/ a (+ 1d0 a)))
                 (b (* (* 2d0 u) (- 1d0 u)))
                 (num (* (- tt0) (- 1d0 b))))
            (declare (type double-float a u b num))
            (when (float-op-is-bad num b +division+)
              (nurbs-error :num-err 'bez-conic-arc-from-pts-and-tangents))
            (setf b (/ num b)
                  w1 0d0)
            (vector-scale ct0 b p1)))))
    w1))

;;; ------------------------------------------------------------- L3-g5

(defun bezier--st-bez-calc-arc-length (pw p eps len depth)
  "The recursive body of bez-calc-arc-length: add to LEN the arc length of
the degree-P Bezier curve PW, to tolerance EPS, at recursion depth DEPTH.
A piece whose control polygon length UB and chord length LB differ by no more
than EPS contributes UB + LB (the caller halves the total); otherwise it is
split in half and each half recursed on with half the tolerance.  Returns
(values eps len): the tolerance actually achieved and the accumulated length.
Signals :quit beyond *recurselim* levels (the C returns 1 with no error set).
NLib: ST_BezCalcArcLength."
  (declare (type cpoint-array pw) (type fixnum p depth) (type double-float eps len))
  (when (> depth (the fixnum *recurselim*))
    (nurbs-error :quit 'bez-calc-arc-length))
  (let* ((ub (dist-c-polygon pw p))
         (lb (let ((a (make-cpoint)) (b (make-cpoint)))
               (declare (dynamic-extent a b))
               (%row-c-pt pw 0 a)
               (%row-c-pt pw p b)
               (dist-cpt-cpt a b)))
         (del (abs (- ub lb))))
    (declare (type double-float ub lb del))
    (if (<= del eps)
        (values del (+ len (+ ub lb)))
        (let ((qw (make-cpoint-array (1+ p)))
              (rw (make-cpoint-array (1+ p))))
          (declare (type cpoint-array qw rw))
          (bez-split pw p 0.5d0 qw rw)
          (let ((epl (/ eps 2d0))
                (epr 0d0))
            (declare (type double-float epl epr))
            (multiple-value-setq (epl len) (bezier--st-bez-calc-arc-length qw p epl len (1+ depth)))
            (setf epr (- eps epl))
            (multiple-value-setq (epr len) (bezier--st-bez-calc-arc-length rw p epr len (1+ depth)))
            (values (+ epl epr) len))))))

(defun bez-calc-arc-length (pw p eps len)
  "Accumulate into LEN the arc length estimate of the degree-P Bezier curve
with control points PW, to tolerance EPS, by recursive subdivision (see
bezier--st-bez-calc-arc-length): each piece adds its control polygon length
plus its chord length, so the caller halves the sum.  Returns (values eps len):
EPS the tolerance actually achieved and LEN the accumulated sum.  Signals
:quit when subdivision needs more than *recurselim* levels.
NLib: N_BezCalcArcLength."
  (declare (type cpoint-array pw) (type fixnum p) (type double-float eps len))
  (bezier--st-bez-calc-arc-length pw p eps len 0))

;;; ------------------------------------------------------------- L4-g0

(defun bez-calc-length (pw p tol)
  "The arc length of the degree-P Bezier curve with control points PW, to
tolerance TOL: BEZ-CALC-ARC-LENGTH with tolerance 2*TOL from a zero sum,
halved.  Signals :quit when subdivision needs more than *recurselim* levels.
NLib: N_BezCalcLength."
  (declare (type cpoint-array pw) (type fixnum p) (type double-float tol))
  (let ((eps (* 2.0d0 tol)))
    (declare (type double-float eps))
    (multiple-value-bind (e total) (bez-calc-arc-length pw p eps 0.0d0)
      (declare (ignore e) (type double-float total))
      (* 0.5d0 total))))

;;; ------------------------------------------------------------- L10-g0

(defun bez-eval-one-basis (p i u)
  "The value at U in [0, 1] of the I-th Bernstein polynomial of degree P, by
the recursive definition.  Signals :ind-err for I outside 0..P and :par-err
for U outside [0, 1].
NLib: N_BezEvalOneBasis."
  (declare (type fixnum p i) (type double-float u))
  (when (or (< i 0) (> i p))
    (nurbs-error :ind-err 'bez-eval-one-basis))
  (when (or (< u 0d0) (> u 1d0))
    (nurbs-error :par-err 'bez-eval-one-basis))
  (cond
    ((= u 0d0) (if (= i 0) 1d0 0d0))
    ((= u 1d0) (if (= i p) 1d0 0d0))
    (t
     (let ((tmp (make-array (1+ p) :element-type 'double-float :initial-element 0d0))
           (omu (- 1d0 u)))
       (declare (type real-vector tmp) (dynamic-extent tmp) (type double-float omu))
       (setf (aref tmp (- p i)) 1d0)
       (fast
         (do ((k 1 (1+ k))) ((> k p))
           (declare (type fixnum k))
           (do ((j p (1- j))) ((< j k))
             (declare (type fixnum j))
             (setf (aref tmp j) (+ (* omu (aref tmp j)) (* u (aref tmp (1- j))))))))
       (aref tmp p)))))

(defun bez-degree-reduce-coefs (p alf oma bet omb)
  "Write the coefficients for reducing a degree-P Bezier curve by one into
the caller's real-vectors: ALF[1..r] and OMA[1..r] (1 - ALF) from the left,
BET[r..p-2] and OMB[r..p-2] (1 - BET) from the right, r = (p-1)/2 truncated.
Returns ALF.
NLib: N_BezDegreeReduceCoefs."
  (declare (type fixnum p) (type real-vector alf oma bet omb))
  (let ((r (truncate (- p 1) 2))
        (fp (float p 1d0)))
    (declare (type fixnum r) (type double-float fp))
    (do ((i 1 (1+ i))) ((> i r))
      (declare (type fixnum i))
      (setf (aref alf i) (/ fp (- fp (float i 1d0)))
            (aref oma i) (- 1d0 (aref alf i))))
    (do ((i (- p 2) (1- i))) ((< i r))
      (declare (type fixnum i))
      (setf (aref bet i) (/ fp (+ (float i 1d0) 1d0))
            (aref omb i) (- 1d0 (aref bet i)))))
  alf)

;;; ------------------------------------------------------------- L11-g0

(defun bez-reduce-degree (pw p alf oma bet omb qw)
  "Reduce the degree-P Bezier curve with control points PW by one, writing
the P new control points into the caller's QW, with the reduction
coefficients ALF OMA BET OMB precomputed by bez-degree-reduce-coefs.  Returns
the maximum error of the reduction.  Signals as bez-eval-one-basis does.
NLib: N_BezReduceDegree."
  (declare (type cpoint-array pw qw) (type fixnum p) (type real-vector alf oma bet omb))
  (let ((pl (make-cpoint)) (pr (make-cpoint)) (aw (make-cpoint)) (bw (make-cpoint))
        (r (truncate (- p 1) 2))
        (fp (float p 1d0)))
    (declare (type cpoint pl pr aw bw) (dynamic-extent pl pr aw bw)
             (type fixnum r) (type double-float fp))
    ;; end points
    (%copy-c-pt-row pw 0 qw 0)
    (%copy-c-pt-row pw p qw (- p 1))
    (cond
      ((oddp p)
       ;; from the left
       (do ((i 1 (1+ i))) ((> i (- r 1)))
         (declare (type fixnum i))
         (%combine-2-c-pt-rows (aref alf i) pw i (aref oma i) qw (- i 1) qw i))
       ;; from the right
       (do ((i (- p 2) (1- i))) ((< i (+ r 1)))
         (declare (type fixnum i))
         (%combine-2-c-pt-rows (aref bet i) pw (+ i 1) (aref omb i) qw (+ i 1) qw i))
       ;; middle control point
       (combine-2-c-pts (aref alf r) (%row-c-pt pw r aw) (aref oma r) (%row-c-pt qw (- r 1) bw) pl)
       (combine-2-c-pts (aref bet r) (%row-c-pt pw (+ r 1) aw) (aref omb r) (%row-c-pt qw (+ r 1) bw) pr)
       (combine-2-c-pts 0.5d0 pl 0.5d0 pr aw)
       (%c-pt-row aw qw r)
       ;; the error
       (let* ((u (* 0.5d0 (- 1d0 (c-sqrt (/ 1d0 fp)))))
              (b (bez-eval-one-basis p r u))
              (b1 (bez-eval-one-basis p (+ r 1) u))
              (dw (dist-cpt-cpt-homo pl pr))
              (a (/ (* 0.5d0 (- fp (float r 1d0))) fp))
              (bb (abs (- b b1))))
         (declare (type double-float u b b1 dw a bb))
         (* (* a bb) dw)))
      (t
       ;; from the left
       (do ((i 1 (1+ i))) ((> i r))
         (declare (type fixnum i))
         (%combine-2-c-pt-rows (aref alf i) pw i (aref oma i) qw (- i 1) qw i))
       ;; from the right
       (do ((i (- p 2) (1- i))) ((< i (+ r 1)))
         (declare (type fixnum i))
         (%combine-2-c-pt-rows (aref bet i) pw (+ i 1) (aref omb i) qw (+ i 1) qw i))
       ;; the error
       (let* ((u (/ (+ (float r 1d0) 1d0) fp))
              (b1 (bez-eval-one-basis p (+ r 1) u)))
         (declare (type double-float u b1))
         (combine-2-c-pts 0.5d0 (%row-c-pt qw r aw) 0.5d0 (%row-c-pt qw (+ r 1) bw) pl)
         (* b1 (dist-cpt-cpt-homo (%row-c-pt pw (+ r 1) aw) pl)))))))

(defun bez-center-radius-from-3-c-pts (p0w p1w p2w c ps pe)
  "The center, radius and start and end angles of the 2D circular arc whose
Bezier control points are P0W P1W P2W (assumed to make a circle).  The
Euclidean end points go into PS and PE and the center into C; returns
(values r alf bet), the radius and the start and end angles in degrees.
When the control points are collinear (P1 within *mtol* of segment P0 P2) C
is (+infinite+ +infinite+ +infinite+), r +infinite+ and the angles
+undefined+.  Signals :inp-err when the radius' division fails, and as
dist-pt-line-seg, vector-normalize-ref and vector-directed-angle do.
NLib: N_BezCenterRadiusFrom3CPts."
  (declare (type cpoint p0w p1w p2w) (type point c ps pe))
  (let ((p0 (make-point)) (p1 (make-point)) (p2 (make-point)) (m (make-point))
        (n (make-point)) (v (make-point)) (rr (make-point))
        (r 0d0) (alf 0d0) (bet 0d0))
    (declare (type point p0 p1 p2 m n v rr) (dynamic-extent p0 p1 p2 m n v rr)
             (type double-float r alf bet))
    ;; Euclidean points and the end points
    (c-pt-to-pt-euclid p0w p0)
    (c-pt-to-pt-euclid p1w p1)
    (c-pt-to-pt-euclid p2w p2)
    (copy-pt p0 ps)
    (copy-pt p2 pe)
    ;; is the arc a line?
    (let ((d (dist-pt-line-seg p1 p0 p2)))
      (declare (type double-float d))
      (when (< d *mtol*)
        (pt-from-xyz +infinite+ +infinite+ +infinite+ c)
        (return-from bez-center-radius-from-3-c-pts
          (values +infinite+ +undefined+ +undefined+))))
    ;; radius, center and angles
    (combine-2-pts 0.5d0 p0 0.5d0 p2 m)
    (let ((a (dist-pt-pt p0 m))
          (b (dist-pt-pt m p1))
          (cc (dist-pt-pt p1 p0)))
      (declare (type double-float a b cc))
      (when (float-op-is-bad (* a cc) b +division+)
        (nurbs-error :inp-err 'bez-center-radius-from-3-c-pts))
      (setf r (/ (* a cc) b)))
    (vector-diff p0 p1 v)
    (vector-perpendicular v n)
    (vector-diff p2 p0 v)
    (when (< (vector-dot n v) 0d0)
      (vector-reverse-in-place n))
    (vector-normalize-ref n)
    (vector-pt-along-vector p0 r n c)
    (vector-create 1d0 0d0 0d0 rr)
    (vector-diff p0 c v)
    (setf alf (vector-directed-angle rr v))
    (vector-diff p2 c v)
    (setf bet (vector-directed-angle rr v))
    (values r alf bet)))

;;; ------------------------------------------------------------- L14-g0

(defun bez-eval-basis (p u b)
  "Write the values at U in [0, 1] of all P+1 Bernstein polynomials of degree
P into B[0..P] (a real-vector the caller allocates), by the triangular
scheme.  Signals :par-err for U outside [0, 1].  Returns B.
NLib: N_BezEvalBasis."
  (declare (type fixnum p) (type double-float u) (type real-vector b))
  (when (or (< u 0d0) (> u 1d0))
    (nurbs-error :par-err 'bez-eval-basis))
  (cond
    ((= u 0d0)
     (setf (aref b 0) 1d0)
     (do ((i 1 (1+ i))) ((> i p))
       (declare (type fixnum i))
       (setf (aref b i) 0d0)))
    ((= u 1d0)
     (do ((i 0 (1+ i))) ((>= i p))
       (declare (type fixnum i))
       (setf (aref b i) 0d0))
     (setf (aref b p) 1d0))
    (t
     (let ((omu (- 1d0 u)))
       (declare (type double-float omu))
       (setf (aref b 0) 1d0)
       (fast
         (do ((i 1 (1+ i))) ((> i p))
           (declare (type fixnum i))
           (let ((saved 0d0))
             (declare (type double-float saved))
             (do ((j 0 (1+ j))) ((>= j i))
               (declare (type fixnum j))
               (let ((temp (aref b j)))
                 (declare (type double-float temp))
                 (setf (aref b j) (+ saved (* omu temp))
                       saved (* u temp))))
             (setf (aref b i) saved)))))))
  b)

;;; ------------------------------------------------------------- L15-g0

(defun bez-interp-4-pts (q pw)
  "Write into PW[0..3] (a cpoint-array the caller allocates) the control
points of the non-rational cubic Bezier curve interpolating the points
Q[0..3] (a point-array) at chord-length parameters (1/3 and 2/3 when a
division by the total chord length fails).  Returns PW.  Signals as
bez-eval-basis, real-matrix-lu-decompose and real-matrix-for-back do.
NLib: N_BezInterp4Pts."
  (declare (type point-array q) (type cpoint-array pw))
  (let ((u (make-real-vector 4)) (dist (make-real-vector 4)) (b (make-real-vector 4))
        (sum 0d0) (pa (make-point)) (pb (make-point)))
    (declare (type real-vector u dist b) (type double-float sum)
             (dynamic-extent u dist b pa pb))
    ;; parameters
    (setf (aref u 0) 0d0
          (aref u 3) 1d0)
    (do ((i 1 (1+ i))) ((> i 3))
      (declare (type fixnum i))
      (%row-pt q (1- i) pa)
      (%row-pt q i pb)
      (setf (aref dist i) (dist-pt-pt pa pb))
      (setf sum (+ sum (aref dist i))))
    (if (or (float-op-is-bad (aref dist 1) sum +division+)
            (float-op-is-bad (aref dist 2) sum +division+))
        (setf (aref u 1) (/ 1d0 3d0)
              (aref u 2) (/ 2d0 3d0))
        (progn
          (setf (aref u 1) (+ (aref u 0) (/ (aref dist 1) sum)))
          (setf (aref u 2) (+ (aref u 1) (/ (aref dist 2) sum)))))
    ;; coefficient matrix
    (let* ((cm (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3))
           (a (rmatrix-rm cm)))
      (declare (type real-matrix a))
      (dotimes (j 4)
        (setf (aref a 3 j) 0d0
              (aref a 0 j) 0d0))
      (setf (aref a 0 0) 1d0
            (aref a 3 3) 1d0)
      (do ((i 1 (1+ i))) ((>= i 3))
        (declare (type fixnum i))
        (bez-eval-basis 3 (aref u i) b)
        (dotimes (j 4)
          (setf (aref a i j) (aref b j))))
      ;; control points
      (real-matrix-lu-decompose cm)
      (real-matrix-for-back cm q +epoint+ pw))
    pw))
