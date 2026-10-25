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

(declaim (inline %bezier-copy-c-row %bezier-combine-2-c-rows %bezier-copy-row
                 %bezier-combine-2-rows))

(defun %bezier-copy-c-row (a i b j)
  "Row J of cpoint-array B := row I of A, as copy-c-pt."
  (declare (type cpoint-array a b) (type index i j))
  (let ((x (aref a i 0)) (y (aref a i 1)) (z (aref a i 2)) (w (aref a i 3)))
    (setf (aref b j 0) x (aref b j 1) y (aref b j 2) z (aref b j 3) w)))

(defun %bezier-combine-2-c-rows (alpha a ia beta b ib c ic)
  "Row IC of C := ALPHA * row IA of A + BETA * row IB of B, as combine-2-c-pts."
  (declare (type double-float alpha beta) (type cpoint-array a b c) (type index ia ib ic))
  (let ((ax (aref a ia 0)) (ay (aref a ia 1)) (az (aref a ia 2)) (aw (aref a ia 3))
        (bx (aref b ib 0)) (by (aref b ib 1)) (bz (aref b ib 2)) (bw (aref b ib 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref c ic 0) (+ (* alpha ax) (* beta bx))
          (aref c ic 1) (+ (* alpha ay) (* beta by))
          (aref c ic 2) (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ (* alpha az) (* beta bz)))
                              ((/= az +no-z+) (* alpha az))
                              ((/= bz +no-z+) (* beta bz))
                              (t +no-z+))
          (aref c ic 3) (if (and (/= aw +no-w+) (/= bw +no-w+)) (+ (* alpha aw) (* beta bw)) +no-w+))))

(defun %bezier-copy-row (a i b j)
  "Row J of point-array B := row I of A, as copy-pt."
  (declare (type point-array a b) (type index i j))
  (let ((x (aref a i 0)) (y (aref a i 1)) (z (aref a i 2)))
    (setf (aref b j 0) x (aref b j 1) y (aref b j 2) z)))

(defun %bezier-combine-2-rows (alpha a ia beta b ib c ic)
  "Row IC of point-array C := ALPHA * row IA of A + BETA * row IB of B, as
combine-2-pts."
  (declare (type double-float alpha beta) (type point-array a b c) (type index ia ib ic))
  (let ((ax (aref a ia 0)) (ay (aref a ia 1)) (az (aref a ia 2))
        (bx (aref b ib 0)) (by (aref b ib 1)) (bz (aref b ib 2)))
    (declare (type double-float ax ay az bx by bz))
    (setf (aref c ic 0) (+ (* alpha ax) (* beta bx))
          (aref c ic 1) (+ (* alpha ay) (* beta by))
          (aref c ic 2) (+ (* alpha az) (* beta bz)))))

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
    (%bezier-copy-row p 0 q 0)
    (%bezier-copy-row p 2 r 2)
    (when (/= wp 0d0)
      (let* ((alf (/ 1d0 (+ 1d0 wp)))
             (oma (- 1d0 alf)))
        (declare (type double-float alf oma))
        (%bezier-combine-2-rows alf p 0 oma p 1 q 1)
        (%bezier-combine-2-rows alf p 2 oma p 1 r 1)
        (%bezier-combine-2-rows 0.5d0 q 1 0.5d0 r 1 q 2)
        (%bezier-copy-row q 2 r 0)
        (setf wq (sqrt (the (double-float 0d0) (* 0.5d0 (+ 1d0 wp))))
              wr wq)))
    (when (= wp 0d0)
      (%bezier-combine-2-rows 1d0 p 0 1d0 p 1 q 1)
      (%bezier-combine-2-rows 1d0 p 2 1d0 p 1 r 1)
      (%bezier-combine-2-rows 0.5d0 q 1 0.5d0 r 1 q 2)
      (%bezier-copy-row q 2 r 0)
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
         (%bezier-copy-c-row pw 0 qw i)
         (%bezier-copy-c-row pw i rw i)))
      ((= u 1d0)
       (dotimes (i (1+ p))
         (%bezier-copy-c-row pw i qw i)
         (%bezier-copy-c-row pw p rw i)))
      (t
       (let ((aw (make-array (list (1+ p) 4) :element-type 'double-float))
             (omu (- 1d0 u)))
         (declare (dynamic-extent aw) (type double-float omu))
         (dotimes (i (1+ p))
           (%bezier-copy-c-row pw i aw i))
         (%bezier-copy-c-row pw 0 qw 0)
         (%bezier-copy-c-row pw p rw p)
         (do ((k 1 (1+ k))) ((> k p))
           (declare (type fixnum k))
           (do ((i 0 (1+ i))) ((> i (- p k)))
             (declare (type fixnum i))
             (%bezier-combine-2-c-rows omu aw i u aw (1+ i) aw i))
           (%bezier-copy-c-row aw 0 qw k)
           (%bezier-copy-c-row aw (- p k) rw (- p k)))))))
  qw)
