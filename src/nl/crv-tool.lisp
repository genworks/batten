;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvTool.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.  Output curves the C sizes or
;;;; allocates on a stack are returned as fresh curves.

(in-package :batten)

;;; ------------------------------------------------- row helpers (not NLib)

;;; N_CopyCPt and N_Combine2CPts on elements of NL_CPOINT arrays, done on rows
;;; of cpoint-arrays with the same arithmetic, every input coordinate read
;;; before any output is written, as C's by-value arguments are.

(declaim (inline %crv-tool-copy-c-row %crv-tool-combine-2-c-rows %crv-tool-row-to-cpoint))

(defun %crv-tool-copy-c-row (a i b j)
  "Row J of cpoint-array B := row I of A, as copy-c-pt."
  (declare (type cpoint-array a b) (type index i j))
  (let ((x (aref a i 0)) (y (aref a i 1)) (z (aref a i 2)) (w (aref a i 3)))
    (setf (aref b j 0) x (aref b j 1) y (aref b j 2) z (aref b j 3) w)))

(defun %crv-tool-combine-2-c-rows (alpha a ia beta b ib c ic)
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

(defun %crv-tool-row-to-cpoint (a i cp)
  "Copy row I of cpoint-array A into cpoint CP and return CP."
  (declare (type cpoint-array a) (type index i) (type cpoint cp))
  (set-cpoint cp (aref a i 0) (aref a i 1) (aref a i 2) (aref a i 3)))

(defun %crv-tool-knot-vector-is-end-param (knots u rname)
  "Signal :par-err in RNAME when U is on or outside the ends of KNOTS.
N_KnotVectorIsEndParam (NL_BasisBasic.cpp), which has no Lisp port yet."
  (declare (type real-vector knots) (type double-float u))
  (when (or (<= u (aref knots 0)) (>= u (aref knots (1- (length knots)))))
    (nurbs-error :par-err rname)))

;;; ------------------------------------------------------- decomposition

(defun crv-decompose-continuity (curp)
  "Split curve CURP at its interior knots of multiplicity at least its degree,
into pieces each at least C1 inside.  Returns (values curves k): CURVES a
simple-vector of the k+1 fresh pieces (one copy of CURP when it has no such
knot, k = 0).  NLib: N_CrvDecomposeContinuity."
  (declare (type curve curp))
  (let* ((n (curve-n curp))
         (pw (curve-cpoints curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (r (1+ (truncate n p)))
         (ks (make-index-vector (1+ r)))
         (ke (make-index-vector (1+ r)))
         (ns (make-index-vector (1+ r)))
         (ne (make-index-vector (1+ r)))
         (ie (1+ p)))
    (declare (type fixnum n p m r ie) (type cpoint-array pw) (type real-vector up)
             (type index-vector ks ke ns ne))
    (setf r -1
          (aref ns 0) 0
          (aref ks 0) ie)
    (do () ((not (< ie m)))
      (let ((i ie))
        (declare (type fixnum i))
        (do () ((not (and (< ie m) (= (aref up ie) (aref up (1+ ie))))))
          (incf ie))
        (when (>= (+ (- ie i) 1) p)
          (incf r)
          (setf (aref ne r) (1- i)
                (aref ke r) i)
          (when (< ie m)
            (setf (aref ns (1+ r)) (aref ne r)
                  (aref ks (1+ r)) (1+ ie))))
        (incf ie)))
    (let ((cura (make-array (1+ r))))
      (dotimes (i (1+ r))
        (let* ((ni (- (aref ne i) (aref ns i)))
               (mi (+ (- (aref ke i) (aref ks i)) (* 2 p) 1))
               (qw (make-cpoint-array (1+ ni)))
               (uq (make-real-vector (1+ mi)))
               (l 0))
          (declare (type fixnum ni mi l))
          (dotimes (j (1+ ni))
            (%crv-tool-copy-c-row pw (+ (aref ns i) j) qw j))
          (dotimes (j (1+ p))
            (setf (aref uq j) (aref up (1- (aref ks i)))))
          (do ((j (1+ p) (1+ j))) ((> j ni))
            (declare (type fixnum j))
            (setf (aref uq j) (aref up (+ (aref ks i) l)))
            (incf l))
          (do ((j (1+ ni) (1+ j))) ((> j mi))
            (declare (type fixnum j))
            (setf (aref uq j) (aref up (aref ke i))))
          (setf (svref cura i) (%make-curve p uq qw))))
      (values cura r))))

;;; ---------------------------------------------------------- refinement

(defun crv-refine (curp knx)
  "Curve CURP with the knots of KNX (a non-decreasing real-vector strictly
inside CURP's domain) inserted, as a fresh curve.  Signals :inp-err when KNX
is empty and :par-err when its first or last knot is not inside the domain.
NLib: N_CrvRefine."
  (declare (type curve curp) (type real-vector knx))
  (let* ((n (curve-n curp))
         (pw (curve-cpoints curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (r (1- (length knx)))
         (x knx))
    (declare (type fixnum n p m r) (type cpoint-array pw) (type real-vector up x))
    (when (< r 0)
      (nurbs-error :inp-err 'crv-refine))
    (%crv-tool-knot-vector-is-end-param up (aref x 0) 'crv-refine)
    (%crv-tool-knot-vector-is-end-param up (aref x r) 'crv-refine)
    (let* ((qw (make-cpoint-array (+ n r 2)))
           (uq (make-real-vector (+ m r 2)))
           (a (basis-find-span up p (aref x 0) +left+))
           (b (1+ (basis-find-span up p (aref x r) +left+))))
      (declare (type cpoint-array qw) (type real-vector uq) (type fixnum a b))
      (fast
        (do ((j 0 (1+ j))) ((> j a))
          (declare (type fixnum j))
          (setf (aref uq j) (aref up j)))
        (do ((j (+ b p) (1+ j))) ((> j m))
          (declare (type fixnum j))
          (setf (aref uq (+ j r 1)) (aref up j)))
        (do ((j 0 (1+ j))) ((> j (- a p)))
          (declare (type fixnum j))
          (%crv-tool-copy-c-row pw j qw j))
        (do ((j (1- b) (1+ j))) ((> j n))
          (declare (type fixnum j))
          (%crv-tool-copy-c-row pw j qw (+ j r 1)))
        (let ((i (+ b p -1))
              (k (+ b p r)))
          (declare (type fixnum i k))
          (do ((j r (1- j))) ((< j 0))
            (declare (type fixnum j))
            (do () ((not (and (<= (aref x j) (aref up i)) (> i a))))
              (%crv-tool-copy-c-row pw (- i p 1) qw (- k p 1))
              (setf (aref uq k) (aref up i))
              (decf k)
              (decf i))
            (%crv-tool-copy-c-row qw (- k p) qw (- k p 1))
            (do ((l 1 (1+ l))) ((> l p))
              (declare (type fixnum l))
              (let ((tt (+ (- k p) l))
                    (alf (- (aref uq (+ k l)) (aref x j))))
                (declare (type fixnum tt) (type double-float alf))
                (if (<= (abs alf) 0d0)
                    (%crv-tool-copy-c-row qw tt qw (1- tt))
                    (let* ((alf (/ alf (- (aref uq (+ k l)) (aref up (+ (- i p) l)))))
                           (oma (- 1d0 alf)))
                      (declare (type double-float alf oma))
                      (%crv-tool-combine-2-c-rows alf qw (1- tt) oma qw tt qw (1- tt))))))
            (setf (aref uq k) (aref x j))
            (decf k))))
      (%make-curve p uq qw))))

;;; ------------------------------------------------------------ knot removal

(defun crv-remove-knot-max-err (cur r s)
  "The maximum error of removing knot U[R], of multiplicity S, once from curve
CUR.  U[R] must be an interior knot with U[R] /= U[R+1] and 0 < S <= p; this
is not checked.  NLib: N_CrvRemoveKnotMaxErr."
  (declare (type curve cur) (type fixnum r s))
  (let* ((pw (curve-cpoints cur))
         (p (curve-degree cur))
         (u (curve-knots cur))
         (rw (make-array (list (1+ (* 2 p)) 4) :element-type 'double-float :initial-element 0d0))
         (ca (make-cpoint)) (cb (make-cpoint)) (cc (make-cpoint))
         (i (- r p))
         (j (- r s))
         (ii 1)
         (jj (+ (- p s) 1)))
    (declare (type cpoint-array pw rw) (type fixnum p i j ii jj) (type real-vector u)
             (dynamic-extent rw ca cb cc))
    (fast
      (%crv-tool-copy-c-row pw (- r p 1) rw 0)
      (%crv-tool-copy-c-row pw (+ (- r s) 1) rw (+ (- p s) 2))
      (do () ((not (> (- j i) 0)))
        (let* ((alf (/ (- (aref u (+ i p 1)) (aref u i)) (- (aref u r) (aref u i))))
               (oma (- 1d0 alf))
               (bet (/ (- (aref u (+ j p 1)) (aref u j)) (- (aref u (+ j p 1)) (aref u r))))
               (omb (- 1d0 bet)))
          (declare (type double-float alf oma bet omb))
          (%crv-tool-combine-2-c-rows alf pw i oma rw (1- ii) rw ii)
          (%crv-tool-combine-2-c-rows bet pw j omb rw (1+ jj) rw jj)
          (incf i)
          (decf j)
          (incf ii)
          (decf jj))))
    (if (< (- j i) 0)
        (dist-cpt-cpt-homo (%crv-tool-row-to-cpoint rw (1- ii) ca)
                           (%crv-tool-row-to-cpoint rw (1+ jj) cb))
        (let* ((alf (/ (- (aref u r) (aref u i)) (- (aref u (+ i p 1)) (aref u i))))
               (oma (- 1d0 alf)))
          (declare (type double-float alf oma))
          (combine-2-c-pts alf (%crv-tool-row-to-cpoint rw (1+ ii) ca)
                           oma (%crv-tool-row-to-cpoint rw (1- ii) cb) cc)
          (dist-cpt-cpt-homo (%crv-tool-row-to-cpoint pw i ca) cc)))))
