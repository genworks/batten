;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Tessellate.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.  A triangulation's database, the C's
;;;; NL_INDEX **DT and NL_INDEX *hd, is a simple-vector of index-vectors (DT[i],
;;;; the neighbors of point i, or DT[i][0] = -1 for a point not in the
;;;; triangulation) and an index-vector (hd[i], the highest index in DT[i]).
;;;; Index arrays the C allocates are returned as fresh index-vectors.

(in-package :batten)

;;; ------------------------------------------------------------- L14-g0

(defun tess-get-triangle-edges (dt hd m)
  "The edges of the triangulation DT, HD over points 0..M: every pair i < j
with j in the list DT[i] of a point in the triangulation, in order of i and
then of the list.  Returns (values alf bet n): fresh index-vectors of the
edges' start and end points and N, their highest index (-1 for none).  As in
the C, the edges are gathered in arrays of 4M+1; more than that is outside
the C's domain (it overruns them) and signals an error here.
NLib: N_TessGetTriangleEdges."
  (declare (type simple-vector dt) (type index-vector hd) (type fixnum m))
  (let ((al (make-index-vector (1+ (* 4 m))))
        (bl (make-index-vector (1+ (* 4 m))))
        (tt -1))
    (declare (type index-vector al bl) (type fixnum tt))
    ;; find all edges
    (do ((i 0 (1+ i))) ((> i m))
      (declare (type fixnum i))
      (let ((dti (svref dt i)))
        (declare (type index-vector dti))
        (when (/= (aref dti 0) -1)
          (do ((r 0 (1+ r))) ((> r (aref hd i)))
            (declare (type fixnum r))
            (let ((j (aref dti r)))
              (declare (type fixnum j))
              (when (< i j)
                (incf tt)
                (setf (aref al tt) i
                      (aref bl tt) j)))))))
    ;; output edges
    (let ((ao (make-index-vector (1+ tt)))
          (bo (make-index-vector (1+ tt))))
      (declare (type index-vector ao bo))
      (fast
        (do ((i 0 (1+ i))) ((> i tt))
          (declare (type fixnum i))
          (setf (aref ao i) (aref al i)
                (aref bo i) (aref bl i))))
      (values ao bo tt))))

(defun tess-get-triangles (dt hd m)
  "The triangles of the triangulation DT, HD over points 0..M: every i < j,
i < k with j and k consecutive in the list DT[i] of a point in the
triangulation, counterclockwise as the lists run.  Returns (values alf bet
gam n): fresh index-vectors of the triangles' first, second and third points
and N, their highest index (-1 for none).  As in the C, the triangles are
gathered in one array of 3k+1, k = 2M+1, the second points from k on and the
third from 2k on: a (k+1)-th triangle overwrites the first one's second and
third points, as the C's does, and more is outside the C's domain (it
overruns the array) and signals an error here.
NLib: N_TessGetTriangles."
  (declare (type simple-vector dt) (type index-vector hd) (type fixnum m))
  (let* ((k (+ (* 2 m) 1))
         (al (make-index-vector (1+ (* 3 k))))
         (bl k)
         (gl (+ bl k))
         (tt -1))
    (declare (type fixnum k bl gl tt) (type index-vector al))
    ;; find all triangles
    (do ((i 0 (1+ i))) ((> i m))
      (declare (type fixnum i))
      (let ((dti (svref dt i)))
        (declare (type index-vector dti))
        (when (/= (aref dti 0) -1)
          (do ((r 0 (1+ r))) ((>= r (aref hd i)))
            (declare (type fixnum r))
            (let ((j (aref dti r))
                  (kk (aref dti (1+ r))))
              (declare (type fixnum j kk))
              (when (and (< i j) (< i kk))
                (incf tt)
                (setf (aref al tt) i
                      (aref al (+ bl tt)) j
                      (aref al (+ gl tt)) kk)))))))
    ;; output triangles
    (let ((ao (make-index-vector (1+ tt)))
          (bo (make-index-vector (1+ tt)))
          (gm (make-index-vector (1+ tt))))
      (declare (type index-vector ao bo gm))
      (fast
        (do ((i 0 (1+ i))) ((> i tt))
          (declare (type fixnum i))
          (setf (aref ao i) (aref al i)
                (aref bo i) (aref al (+ bl i))
                (aref gm i) (aref al (+ gl i)))))
      (values ao bo gm tt))))

(defun tess-srf-area (sur)
  "An estimate of surface SUR's area and of its average lengths in u and v,
for sizing a tessellation: the product of the longest control polygon of a
column (u) and of a row (v) of its Euclidean control points, and the mean
lengths of those polygons.  Returns (values a lu lv).
NLib: N_TessSrfArea."
  (declare (type surface sur))
  (let* ((n (surface-n sur)) (m (surface-m sur))
         (pw (surface-cpoints sur))
         (p (make-point-array (1+ (nl-max n m))))
         (cw (make-cpoint)) (a (make-point)) (b (make-point))
         (maxu 0d0) (sumu 0d0) (maxv 0d0) (sumv 0d0))
    (declare (type fixnum n m) (type cpoint-net pw) (type point-array p)
             (type cpoint cw) (type point a b) (dynamic-extent cw a b)
             (type double-float maxu sumu maxv sumv))
    ;; maximum polygon lengths in u and v
    (flet ((euclid (i j k)
             ;; P[k] = the Euclidean point of Pw[i][j]
             (c-pt-to-pt-euclid (%net-c-pt pw i j cw) a)
             (%pt-row a p k))
           (dist (k)
             ;; |P[k] - P[k-1]|
             (dist-pt-pt (%row-pt p k a) (%row-pt p (1- k) b))))
      (do ((j 0 (1+ j))) ((> j m))
        (declare (type fixnum j))
        (euclid 0 j 0)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((i 1 (1+ i))) ((> i n))
            (declare (type fixnum i))
            (euclid i j i)
            (setf len (+ len (dist i))))
          (setf sumu (+ sumu len))
          (when (> len maxu) (setf maxu len))))
      (do ((i 0 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (euclid i 0 0)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((j 1 (1+ j))) ((> j m))
            (declare (type fixnum j))
            (euclid i j j)
            (setf len (+ len (dist j))))
          (setf sumv (+ sumv len))
          (when (> len maxv) (setf maxv len)))))
    ;; area and lengths
    (values (* maxu maxv)
            (/ sumu (+ (float m 1d0) 1d0))
            (/ sumv (+ (float n 1d0) 1d0)))))
