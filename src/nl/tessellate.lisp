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

;;; ------------------------------------------------------------- L18-g0
;;;
;;; N_TessTrimmedSrf and its 31 file-private helpers.  Representations:
;;;   NL_REAL ***x, ***y   a simple-vector (outer loop i) of simple-vectors
;;;                        (polygon j: 0 the outer one, 1..ht[i] the inner
;;;                        ones) of real-vectors (vertices)
;;;   NL_INDEX *ht         an index-vector
;;;   NL_INDEX **hp, a polygon list p, a triangulation DT
;;;                        a simple-vector of index-vectors
;;;   NL_REAL **UL ... VT  a simple-vector of real-vectors
;;;   NL_FLAG *R, *T and other flag arrays   index-vectors of +yes+ / +no+ ...
;;;   NL_INDEX **TM, NL_FLAG **SM (rectangular)   index-matrices
;;;   NL_CURVE ***cuo      a simple-vector (outer loop i) of simple-vectors of
;;;                        curves (its segments)
;;;   NL_CURVE ****cui     a simple-vector (outer loop i) of simple-vectors
;;;                        (inner loop j) of simple-vectors of curves
;;; An array the C passes as T ** and reallocates is passed in and returned,
;;; possibly fresh, among the values in C argument order, together with the
;;; scalar in/out arguments.  NLib's reallocations (N_Realloc1dRealArray and
;;; the like) are fresh arrays holding the old elements 0..min(old, new); the
;;; rest, unset in the C, are 0 here and never read before written.

(defun %tess-realloc-real-vector (vec n m)
  "N_Realloc1dRealArray(&VEC, N, M): a fresh real-vector of highest index M
holding VEC[0..min(N, M)], the rest 0."
  (declare (type real-vector vec) (type fixnum n m))
  (let ((w (make-real-vector (1+ m))))
    (replace w vec :end2 (min (1+ (min n m)) (length vec)))
    w))

(defun %tess-realloc-index-vector (vec n m)
  "N_Realloc1dIntArray / N_Realloc1dFlagArray(&VEC, N, M): a fresh
index-vector of highest index M holding VEC[0..min(N, M)], the rest 0."
  (declare (type index-vector vec) (type fixnum n m))
  (let ((w (make-index-vector (1+ m))))
    (replace w vec :end2 (min (1+ (min n m)) (length vec)))
    w))

(defun %tess-c-index (x)
  "(NL_INDEX)X as the C compiles it on x86-64 (cvttsd2si): X truncated, or,
for a NaN, the 'integer indefinite' LONG_MIN -- here most-negative-fixnum,
which compares the same against every fixnum (as do the infinities and
anything beyond 2^62, which the C's arrays could never hold).  Used only
where the C's result feeds comparisons: ST_TessBBoxXTrim's polygon size
estimates, which a degenerate trimming polygon box makes 0 * inf."
  (declare (type double-float x))
  (if (and (= x x) (< -4.6d18 x 4.6d18))
      (values (truncate x))
      most-negative-fixnum))

(defun tess-trimmed-srf (sur cuo cui m ho hi hs epc eps tol)
  "Triangulate the trimmed NURBS surface SUR.  The trimmed domain is bounded
by the outer loops CUO[0..M] (CUO[i] a simple-vector of the curves of loop
i, highest index HO[i]; counterclockwise) and, inside outer loop i, the inner
loops CUI[i][0..HI[i]] (CUI[i][j] a simple-vector of curves, highest index
HS[i][j]; clockwise), all 2D curves in SUR's parameter space, disjoint and
closed.  EPC is the polygonization tolerance of the trimming curves, EPS the
tessellation tolerance, TOL the parameter space tolerance (lowered to
0.001 EPS when larger).  SUR's knot vectors are rescaled in place.  Returns
(values u v n dt hd): fresh real-vectors U, V of the vertex parameters
0..N (they may be longer than N+1), and the triangulation, DT a simple-vector
of index-vectors (DT[i] the counterclockwise neighbors of vertex i, DT[i][0]
= -1 for a vertex out of the domain) and HD an index-vector of their highest
indexes.  Signals :inp-err for a surface of no area, :num-err for a
degenerate domain, :tol-err and :geo-err where the polygonization or the
triangulation fails, and as the routines it calls do.
NLib: N_TessTrimmedSrf."
  (declare (type surface sur) (type simple-vector cuo cui hs) (type fixnum m)
           (type index-vector ho hi) (type double-float epc eps tol))
  (let ((side 0) (scor 0) (ecor 0)
        (ucl +no+) (vcl +no+) (uin +no+) (vin +no+)
        (np 0) (nt 0) (is 0) (kxy 0) (kk 0)
        (seg (vector))
        (ut (make-real-vector 0)) (vt (make-real-vector 0))
        (tmp 0d0) (lu 0d0) (lv 0d0) (usl 0d0) (vsl 0d0) (off 0d0)
        (xl 0d0) (xr 0d0) (yb 0d0) (yt 0d0)
        (d1 0d0) (d2 0d0) (d3 0d0) (d4 0d0)
        (ax 0d0) (ay 0d0) (bx 0d0) (by 0d0)
        (nxy 20)
        (pa (make-point)))
    (declare (type fixnum side scor ecor ucl vcl uin vin np nt is kxy kk nxy)
             (type simple-vector seg) (type real-vector ut vt)
             (type double-float tmp lu lv usl vsl off xl xr yb yt d1 d2 d3 d4 ax ay bx by)
             (type point pa))
    ;; check parameter space is small compared to 3D tol
    (when (> tol (* eps 0.001d0))
      (setf tol (* 0.001d0 eps)))
    ;; get some constants
    (multiple-value-setq (xl xr yb yt) (srf-get-parameter-bounds sur))
    (multiple-value-setq (tmp lu lv) (tess-srf-area sur))
    (when (<= tmp (* (the double-float *mtol*) (the double-float *mtol*)))
      (nurbs-error :inp-err 'tess-trimmed-srf))
    (when (float-op-is-bad lu (- xr xl) +division+)
      (nurbs-error :num-err 'tess-trimmed-srf))
    (when (float-op-is-bad lv (- yt yb) +division+)
      (nurbs-error :num-err 'tess-trimmed-srf))
    (setf usl (/ lu (- xr xl))
          vsl (/ lv (- yt yb))
          off (* 30.0d0 tol))
    ;; polygonize trimming curves: allocate the arrays
    (let ((ht (make-index-vector (1+ m)))
          (hp (make-array (1+ m)))
          (x (make-array (1+ m)))
          (y (make-array (1+ m))))
      (declare (type index-vector ht) (type simple-vector hp x y))
      (dotimes (i (1+ m))
        (setf (aref ht i) (+ (aref hi i) 1)))
      (dotimes (i (1+ m))
        (setf (svref hp i) (make-index-vector (1+ (aref ht i)))
              (svref x i) (make-array (1+ (aref ht i)))
              (svref y i) (make-array (1+ (aref ht i))))
        (dotimes (j (1+ (aref ht i)))
          ;; N_AllocReal1dArray(nxy + 1): highest index nxy + 1
          (setf (svref (svref x i) j) (make-real-vector (+ nxy 2))
                (svref (svref y i) j) (make-real-vector (+ nxy 2)))))
      ;; for each outer loop, polygonize inner and outer loops
      (dotimes (i (1+ m))
        (let ((cuoi (svref cuo i))
              (xi (svref x i))
              (yi (svref y i))
              (hpi (svref hp i)))
          (declare (type simple-vector cuoi xi yi) (type index-vector hpi))
          (symbol-macrolet ((x0 (the real-vector (svref xi 0)))
                            (y0 (the real-vector (svref yi 0))))
            (flet ((on-boundary (cur)
                     ;; ST_TessIsOnBoundary; side, scor, ecor keep their old
                     ;; values where the C leaves them unset
                     (multiple-value-bind (bnd sd sc ec)
                         (tessellate--st-tess-is-on-boundary cur tol xl xr yb yt)
                       (when sd (setf side sd))
                       (when sc (setf scor sc))
                       (when ec (setf ecor ec))
                       bnd))
                   (put (k xv yv)
                     ;; x[i][0][k] = xv, y[i][0][k] = yv
                     (setf (aref x0 k) xv
                           (aref y0 k) yv)))
              (declare (inline put))
              ;; initialize first point of polygon
              (with-euclid-row ((px py pz) (curve-cpoints (svref cuoi 0)) 0)
                (setf (aref x0 0) px
                      (aref y0 0) py))
              (multiple-value-setq (seg is) (crv-decompose-continuity (svref cuoi 0)))
              (when (> is 0)
                ;; eliminate degenerate segment at start
                (let* ((s0 (svref seg 0))
                       (len (dist-c-polygon (curve-cpoints s0) (curve-n s0))))
                  (declare (type curve s0) (type double-float len))
                  (when (<= len (* 0.1d0 tol))
                    (do ((j 1 (1+ j))) ((> j is))
                      (declare (type fixnum j))
                      (setf (svref seg (- j 1)) (svref seg j)))
                    (decf is))))
              (when (on-boundary (svref seg 0))
                (cond ((= side +left+)
                       (when (= scor +upperleft+)
                         (put 0 (- (aref x0 0) off) (+ (aref y0 0) off))))
                      ((= side +bottom+)
                       (when (= scor +lowerleft+)
                         (put 0 (- (aref x0 0) off) (- (aref y0 0) off))))
                      ((= side +right+)
                       (when (= scor +lowerright+)
                         (put 0 (+ (aref x0 0) off) (- (aref y0 0) off))))
                      ((= side +top+)
                       (when (= scor +upperright+)
                         (put 0 (+ (aref x0 0) off) (+ (aref y0 0) off))))))
              ;; polygonize each outer segment
              (setf np 0
                    kxy nxy)
              (do ((j 0 (1+ j))) ((> j (aref ho i)))
                (declare (type fixnum j))
                (when (> j 0)
                  (multiple-value-setq (seg is) (crv-decompose-continuity (svref cuoi j))))
                (do ((ks 0 (1+ ks))) ((> ks is))
                  (declare (type fixnum ks))
                  (let ((sk (svref seg ks)))
                    (declare (type curve sk))
                    (if (on-boundary sk)
                        ;; segment lies along the boundary
                        (let ((pw (curve-cpoints sk))
                              (nc (curve-n sk)))
                          (declare (type cpoint-array pw) (type fixnum nc))
                          (with-euclid-row ((px py pz) pw 0)
                            (setf ax px ay py))
                          (with-euclid-row ((px py pz) pw nc)
                            (setf bx px by py))
                          (when (> (+ np 3) kxy)
                            (setf kk (* 3 (+ (- is ks) 1)))
                            (setf (svref xi 0) (%tess-realloc-real-vector x0 np (+ kxy kk 1)))
                            (setf (svref yi 0) (%tess-realloc-real-vector y0 np (+ kxy kk 1)))
                            (setf kxy (+ kxy kk)))
                          (cond
                            ((= side +left+)
                             (cond ((and (= scor +upperleft+) (= ecor +lowerleft+))
                                    (when (>= (aref x0 np) (- xl tol))
                                      (put (+ np 1) (- ax off) ay)
                                      (incf np))
                                    (put (+ np 1) (- bx off) (- by off))
                                    (incf np))
                                   ((= scor +upperleft+)
                                    (when (>= (aref x0 np) (- xl tol))
                                      (put (+ np 1) (- ax off) ay)
                                      (incf np))
                                    (put (+ np 1) (- bx off) by)
                                    (put (+ np 2) bx by)
                                    (incf np 2))
                                   ((= ecor +lowerleft+)
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) (- ax off))
                                                 (= (aref y0 (- np 1)) ay))
                                        (decf np 2)))
                                    (put (+ np 1) (- ax off) ay)
                                    (put (+ np 2) (- bx off) (- by off))
                                    (incf np 2))
                                   (t
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) (- ax off))
                                                 (= (aref y0 (- np 1)) ay))
                                        (decf np 2)))
                                    (put (+ np 1) (- ax off) ay)
                                    (put (+ np 2) (- bx off) by)
                                    (put (+ np 3) bx by)
                                    (incf np 3))))
                            ((= side +bottom+)
                             (cond ((and (= scor +lowerleft+) (= ecor +lowerright+))
                                    (when (>= (aref y0 np) (- yb tol))
                                      (put (+ np 1) ax (- ay off))
                                      (incf np))
                                    (put (+ np 1) (+ bx off) (- by off))
                                    (incf np))
                                   ((= scor +lowerleft+)
                                    (when (>= (aref y0 np) (- yb tol))
                                      (put (+ np 1) ax (- ay off))
                                      (incf np))
                                    (put (+ np 1) bx (- by off))
                                    (put (+ np 2) bx by)
                                    (incf np 2))
                                   ((= ecor +lowerright+)
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) ax)
                                                 (= (aref y0 (- np 1)) (- ay off)))
                                        (decf np 2)))
                                    (put (+ np 1) ax (- ay off))
                                    (put (+ np 2) (+ bx off) (- by off))
                                    (incf np 2))
                                   (t
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) ax)
                                                 (= (aref y0 (- np 1)) (- ay off)))
                                        (decf np 2)))
                                    (put (+ np 1) ax (- ay off))
                                    (put (+ np 2) bx (- by off))
                                    (put (+ np 3) bx by)
                                    (incf np 3))))
                            ((= side +right+)
                             (cond ((and (= scor +lowerright+) (= ecor +upperright+))
                                    (when (<= (aref x0 np) (+ xr tol))
                                      (put (+ np 1) (+ ax off) ay)
                                      (incf np))
                                    (put (+ np 1) (+ bx off) (+ by off))
                                    (incf np))
                                   ((= scor +lowerright+)
                                    (when (<= (aref x0 np) (+ xr tol))
                                      (put (+ np 1) (+ ax off) ay)
                                      (incf np))
                                    (put (+ np 1) (+ bx off) by)
                                    (put (+ np 2) bx by)
                                    (incf np 2))
                                   ((= ecor +upperright+)
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) (+ ax off))
                                                 (= (aref y0 (- np 1)) ay))
                                        (decf np 2)))
                                    (put (+ np 1) (+ ax off) ay)
                                    (put (+ np 2) (+ bx off) (+ by off))
                                    (incf np 2))
                                   (t
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) (+ ax off))
                                                 (= (aref y0 (- np 1)) ay))
                                        (decf np 2)))
                                    (put (+ np 1) (+ ax off) ay)
                                    (put (+ np 2) (+ bx off) by)
                                    (put (+ np 3) bx by)
                                    (incf np 3))))
                            ((= side +top+)
                             (cond ((and (= scor +upperright+) (= ecor +upperleft+))
                                    (when (<= (aref y0 np) (+ yt tol))
                                      (put (+ np 1) ax (+ ay off))
                                      (incf np))
                                    (put (+ np 1) (- bx off) (+ by off))
                                    (incf np))
                                   ((= scor +upperright+)
                                    (when (<= (aref y0 np) (+ yt tol))
                                      (put (+ np 1) ax (+ ay off))
                                      (incf np))
                                    (put (+ np 1) bx (+ by off))
                                    (put (+ np 2) bx by)
                                    (incf np 2))
                                   ((= ecor +upperleft+)
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) ax)
                                                 (= (aref y0 (- np 1)) (+ ay off)))
                                        (decf np 2)))
                                    (put (+ np 1) ax (+ ay off))
                                    (put (+ np 2) (- bx off) (+ by off))
                                    (incf np 2))
                                   (t
                                    (when (> np 0)
                                      (when (and (= (aref x0 (- np 1)) ax)
                                                 (= (aref y0 (- np 1)) (+ ay off)))
                                        (decf np 2)))
                                    (put (+ np 1) ax (+ ay off))
                                    (put (+ np 2) bx (+ by off))
                                    (put (+ np 3) bx by)
                                    (incf np 3))))))
                        ;; not on the boundary: polygonize segment
                        (progn
                          (multiple-value-setq (ut vt nt)
                            (tessellate--st-tess-polygonize-trim-crv sur lu lv sk epc))
                          (when (< nt 0)
                            (nurbs-error :geo-err 'tess-trimmed-srf))
                          (when (= nt 0)
                            (when (= ks is)
                              (when (< np 3)
                                (nurbs-error :geo-err 'tess-trimmed-srf))
                              (multiple-value-setq (d1 d2) (crv-get-param-bounds sk))
                              (crv-eval sk d2 +left+ pa)
                              (setf (aref ut 1) (aref pa 0)
                                    (aref vt 1) (aref pa 1))
                              (decf np)
                              (setf nt 1)))
                          (when (> (+ np nt) kxy)
                            (setf kk (* nt (+ (- is ks) 1)))
                            (setf (svref xi 0) (%tess-realloc-real-vector x0 np (+ kxy kk 1)))
                            (setf (svref yi 0) (%tess-realloc-real-vector y0 np (+ kxy kk 1)))
                            (setf kxy (+ kxy kk)))
                          (do ((k 1 (1+ k))) ((> k nt))
                            (declare (type fixnum k))
                            (put (+ np k) (aref ut k) (aref vt k)))
                          (setf np (+ np nt)))))))
              (when (<= np 2)
                (nurbs-error :tol-err 'tess-trimmed-srf))
              ;; check for closure of outer loop
              (setf d1 (- (aref x0 0) (aref x0 np))
                    d2 (- (aref y0 0) (aref y0 np)))
              (when (> (c-sqrt (+ (* d1 d1) (* d2 d2))) tol)
                ;; close it up
                (incf np)
                (put np (aref x0 0) (aref y0 0)))
              (setf (aref hpi 0) np)))
          ;; now polygonize all inner loops residing inside the i-th outer loop
          (do ((j 0 (1+ j))) ((> j (aref hi i)))
            (declare (type fixnum j))
            (let ((cuij (svref (the simple-vector (svref cui i)) j))
                  (hsij (aref (the index-vector (svref hs i)) j)))
              (declare (type simple-vector cuij) (type fixnum hsij))
              (symbol-macrolet ((xj (the real-vector (svref xi (+ j 1))))
                                (yj (the real-vector (svref yi (+ j 1)))))
                (setf np 0
                      kxy nxy)
                (with-euclid-row ((px py pz) (curve-cpoints (svref cuij 0)) 0)
                  (setf (aref xj 0) px
                        (aref yj 0) py))
                (do ((k 0 (1+ k))) ((> k hsij))
                  (declare (type fixnum k))
                  (multiple-value-setq (ut vt nt)
                    (tessellate--st-tess-polygonize-trim-crv sur lu lv (svref cuij k) epc))
                  (when (<= nt 0)
                    (nurbs-error :geo-err 'tess-trimmed-srf))
                  (when (> (+ np nt) kxy)
                    (setf kk (* nt (+ (- hsij k) 1)))
                    (setf (svref xi (+ j 1)) (%tess-realloc-real-vector xj np (+ kxy kk)))
                    (setf (svref yi (+ j 1)) (%tess-realloc-real-vector yj np (+ kxy kk)))
                    (setf kxy (+ kxy kk)))
                  (do ((l 1 (1+ l))) ((> l nt))
                    (declare (type fixnum l))
                    (setf (aref xj (+ np l)) (aref ut l)
                          (aref yj (+ np l)) (aref vt l)))
                  (setf np (+ np nt)))
                (when (<= np 2)
                  (nurbs-error :tol-err 'tess-trimmed-srf))
                (setf (aref hpi (+ j 1)) np))))))
      ;; validate outer loop polygons
      (dotimes (i (1+ m))
        (let ((x0 (svref (the simple-vector (svref x i)) 0))
              (y0 (svref (the simple-vector (svref y i)) 0))
              (h (aref (the index-vector (svref hp i)) 0))
              (tol7 (* 7.0d0 tol)))
          (declare (type real-vector x0 y0) (type fixnum h) (type double-float tol7))
          (do ((k 0 (1+ k))) ((> k h))
            (declare (type fixnum k))
            (when (<= (abs (- (aref x0 k) xl)) tol7) (setf (aref x0 k) xl))
            (when (<= (abs (- (aref x0 k) xr)) tol7) (setf (aref x0 k) xr))
            (when (<= (abs (- (aref y0 k) yb)) tol7) (setf (aref y0 k) yb))
            (when (<= (abs (- (aref y0 k) yt)) tol7) (setf (aref y0 k) yt)))
          ;; a run of vertices on a domain side whose neighbors both lie
          ;; beyond it is moved to their mean
          (macrolet ((straighten (arr edge beyond)
                       `(do ((k 1 (1+ k))) ((>= k h))
                          (declare (type fixnum k))
                          (when (and (= (aref ,arr k) ,edge) (,beyond (aref ,arr (- k 1)) ,edge))
                            (let ((kk (+ k 1)))
                              (declare (type fixnum kk))
                              (do () ((not (and (< kk h) (= (aref ,arr kk) ,edge))))
                                (incf kk))
                              (when (,beyond (aref ,arr kk) ,edge)
                                (setf d1 (* 0.5d0 (+ (aref ,arr (- k 1)) (aref ,arr kk))))
                                (do ((kkk k (1+ kkk))) ((>= kkk kk))
                                  (declare (type fixnum kkk))
                                  (setf (aref ,arr kkk) d1))
                                (setf k kk)))))))
            (straighten x0 xl <)
            (straighten x0 xr >)
            (straighten y0 yb <)
            (straighten y0 yt >))))
      ;; rescale the polygons and subdivide the surface
      (let ((ulb (make-array (1+ m)))
            (urb (make-array (1+ m)))
            (vbb (make-array (1+ m)))
            (vtb (make-array (1+ m))))
        (declare (type simple-vector ulb urb vbb vtb))
        (dotimes (i (1+ m))
          (setf (svref ulb i) (make-real-vector (1+ (aref ht i)))
                (svref urb i) (make-real-vector (1+ (aref ht i)))
                (svref vbb i) (make-real-vector (1+ (aref ht i)))
                (svref vtb i) (make-real-vector (1+ (aref ht i)))))
        (setf tol (* tol (nl-max lu lv)))
        (setf d1 (+ lu (* 0.001d0 tol))
              d2 (+ lv (* 0.001d0 tol)))
        (dotimes (i (1+ m))
          (let ((uli (svref ulb i)) (uri (svref urb i)) (vbi (svref vbb i)) (vti (svref vtb i))
                (hpi (svref hp i)))
            (declare (type real-vector uli uri vbi vti) (type index-vector hpi))
            (dotimes (j (1+ (aref ht i)))
              (let ((xij (svref (the simple-vector (svref x i)) j))
                    (yij (svref (the simple-vector (svref y i)) j)))
                (declare (type real-vector xij yij))
                (setf (aref vbi j) +bigd+
                      (aref uli j) +bigd+)
                (setf (aref vti j) (- +bigd+)
                      (aref uri j) (- +bigd+))
                (do ((k 0 (1+ k))) ((> k (aref hpi j)))
                  (declare (type fixnum k))
                  (setf (aref xij k) (* usl (- (aref xij k) xl)))
                  (setf (aref yij k) (* vsl (- (aref yij k) yb)))
                  (when (< (aref xij k) 0d0) (setf (aref xij k) (* -3.1d0 tol)))
                  (when (> (aref xij k) d1) (setf (aref xij k) (+ lu (* 3.1d0 tol))))
                  (when (< (aref yij k) 0d0) (setf (aref yij k) (* -3.1d0 tol)))
                  (when (> (aref yij k) d2) (setf (aref yij k) (+ lv (* 3.1d0 tol))))
                  (when (< (aref xij k) (aref uli j)) (setf (aref uli j) (aref xij k)))
                  (when (> (aref xij k) (aref uri j)) (setf (aref uri j) (aref xij k)))
                  (when (< (aref yij k) (aref vbi j)) (setf (aref vbi j) (aref yij k)))
                  (when (> (aref yij k) (aref vti j)) (setf (aref vti j) (aref yij k))))))))
        (multiple-value-bind (u v rf tf n tm sm r s)
            (tessellate--st-tess-subdivide-srf sur eps x y ht hp m ulb urb vbb vtb tmp lu lv tol)
          (declare (type real-vector u v) (type index-vector rf tf) (type fixnum n r s)
                   (type index-matrix tm sm))
          ;; for closed surfaces, match subdivision points along common boundary
          (when (srf-is-closed sur +udir+) (setf ucl +yes+))
          (when (srf-is-closed sur +vdir+) (setf vcl +yes+))
          (when (or (= ucl +yes+) (= vcl +yes+))
            ;; see how many points to add
            (let ((ul u) (vl v) (nl n) (ni 0))
              (declare (type real-vector ul vl) (type fixnum nl ni))
              (when (= ucl +yes+)
                (dotimes (j (1+ s))
                  (let ((a (aref tm 0 j)) (b (aref tm r j)))
                    (declare (type fixnum a b))
                    (when (or (and (< a 0) (>= b 0)) (and (>= a 0) (< b 0)))
                      (incf ni)
                      (setf uin +yes+)))))
              (when (= vcl +yes+)
                (dotimes (i (1+ r))
                  (let ((a (aref tm i 0)) (b (aref tm i s)))
                    (declare (type fixnum a b))
                    (when (or (and (< a 0) (>= b 0)) (and (>= a 0) (< b 0)))
                      (incf ni)
                      (setf vin +yes+)))))
              (when (> ni 0)
                ;; reallocate memory
                (setf ul (%tess-realloc-real-vector ul nl (+ nl ni))
                      vl (%tess-realloc-real-vector vl nl (+ nl ni))
                      rf (%tess-realloc-index-vector rf nl (+ nl ni))
                      tf (%tess-realloc-index-vector tf nl (+ nl ni)))
                ;; insert new points
                (when (and (= ucl +yes+) (= uin +yes+))
                  (dotimes (j (1+ s))
                    (let ((a (aref tm 0 j)) (b (aref tm r j)))
                      (declare (type fixnum a b))
                      (when (and (< a 0) (>= b 0))
                        (setf (aref tm 0 j) (+ nl 1)
                              (aref ul (+ nl 1)) 0d0
                              (aref vl (+ nl 1)) (aref vl b)
                              (aref sm 0 j) +nostatus+
                              (aref rf (+ nl 1)) +no+
                              (aref tf (+ nl 1)) +no+)
                        (incf nl))
                      (when (and (>= a 0) (< b 0))
                        (setf (aref tm r j) (+ nl 1)
                              (aref ul (+ nl 1)) lu
                              (aref vl (+ nl 1)) (aref vl a)
                              (aref sm r j) +nostatus+
                              (aref rf (+ nl 1)) +no+
                              (aref tf (+ nl 1)) +no+)
                        (incf nl)))))
                (when (and (= vcl +yes+) (= vin +yes+))
                  (dotimes (i (1+ r))
                    (let ((a (aref tm i 0)) (b (aref tm i s)))
                      (declare (type fixnum a b))
                      (when (and (< a 0) (>= b 0))
                        (setf (aref tm i 0) (+ nl 1)
                              (aref ul (+ nl 1)) (aref ul b)
                              (aref vl (+ nl 1)) 0d0
                              (aref sm i 0) +nostatus+
                              (aref rf (+ nl 1)) +no+
                              (aref tf (+ nl 1)) +no+)
                        (incf nl))
                      (when (and (>= a 0) (< b 0))
                        (setf (aref tm i s) (+ nl 1)
                              (aref ul (+ nl 1)) (aref ul a)
                              (aref vl (+ nl 1)) lv
                              (aref sm i s) +nostatus+
                              (aref rf (+ nl 1)) +no+
                              (aref tf (+ nl 1)) +no+)
                        (incf nl))))))
              (setf u ul
                    v vl
                    n nl)))
          ;; now triangulate the subdivided domain
          (multiple-value-setq (xl xr yb yt) (srf-get-parameter-bounds sur))
          (multiple-value-bind (u v n dt hd)
              (tessellate--st-tess-region-of-rectangles u v n rf tf tm sm r s x y ht hp m
                                                        ulb urb vbb vtb xl xr yb yt tol)
            (declare (type real-vector u v) (type fixnum n) (type simple-vector dt)
                     (type index-vector hd))
            ;; check that tessellation points are not out of the valid region
            (when (= m 0)
              (setf tmp (* 1.2d0 tol))
              (setf d1 (- (aref (the real-vector (svref ulb 0)) 0) tmp)
                    d2 (+ (aref (the real-vector (svref urb 0)) 0) tmp)
                    d3 (- (aref (the real-vector (svref vbb 0)) 0) tmp)
                    d4 (+ (aref (the real-vector (svref vtb 0)) 0) tmp))
              (dotimes (i (1+ n))
                (when (/= (aref (the index-vector (svref dt i)) 0) -1)
                  (when (or (< (aref u i) d1) (> (aref u i) d2))
                    (nurbs-error :geo-err 'tess-trimmed-srf))
                  (when (or (< (aref v i) d3) (> (aref v i) d4))
                    (nurbs-error :geo-err 'tess-trimmed-srf)))))
            (values u v n dt hd)))))))

(defconstant +tessellate--prc+ 0.01d0
  "NL_PRIVATE NL_REAL prc of NL_Tessellate.cpp: the fraction of a trimming
curve's box diagonal its polygonization tolerance is held below.")

;;; ------------------------------------------------------------------ ST_TessSubdivideSrf

(defun tessellate--st-tess-subdivide-srf (sur eps x y ht hp m uls urs vbs vts area lu lv tol)
  "Subdivide surface SUR into rectangles each within EPS of the surface.  The
knots of SUR are first rescaled in place to [0, LU] x [0, LV].  X, Y, HT, HP,
M are the trimming polygons (see the file header), ULS, URS, VBS, VTS their
bounding boxes, AREA the surface's approximate area, TOL the parameter space
tolerance.  Rectangles classified +outp+ by tessellate--st-tess-is-inside-trim
are not subdivided.  Returns (values u v rf tf n tm sm r s): U, V fresh
real-vectors of the rectangle corners 0..N; RF, TF fresh index-vectors of
the right and top neighbor flags; TM the topology matrix (TM[i][j] >= 0 the
index of the corner (us[i], vs[j]), -1 none) and SM the status matrix, both
fresh index-matrices of (R+1) x (S+1).  When the subdivision is a single
rectangle, the polygon vertices X, Y are clamped to the domain and
near-duplicates dropped, in place, with HP updated in place.  Signals
:tol-err for EPS below *mtol*, :par-err for a degenerate rectangle, :num-err
when a corner is not found in the hash tables, and as srf-extract-patch and
update-real-hash do; srf-is-flat's errors are ignored, as in the C.
NLib: ST_TessSubdivideSrf (file-private, NL_Tessellate.cpp)."
  (declare (type surface sur) (type double-float eps area lu lv tol)
           (type simple-vector x y hp uls urs vbs vts) (type index-vector ht)
           (type fixnum m))
  (when (< eps (the double-float *mtol*))
    (nurbs-error :tol-err 'tessellate--st-tess-subdivide-srf))
  ;; rescale knot rectangle
  (let ((rct (make-real-vector 4)))
    (declare (dynamic-extent rct))
    (create-rectangle rct 0d0 lu 0d0 lv)
    (srf-reparam-to-interval sur rct +uvdir+))
  (let* ((ptol *ptol*)
         (fac (/ area eps))
         (nuv 0) (nh 0) (ns 0)
         (ctol 0d0))
    (declare (type double-float ptol fac ctol) (type fixnum nuv nh ns))
    ;; estimate various indexes
    (when (> fac 500000.0d0)
      (setf fac 500000.0d0))
    (setf nuv (nl-max 6 (the fixnum (truncate (* 0.25d0 fac)))))
    (when (< nuv 20000)
      (setf nuv (* nuv 2)))
    (setf nh (nl-max 20 (the fixnum (truncate (* 4.0d0 (c-sqrt fac))))))
    (setf ns (nl-max 4 (the fixnum (truncate (* 4.0d0 (the double-float (log fac)))))))
    (setf ns (* ns 2))
    (setf ctol (* 20.0d0 tol))
    ;; (1) allocate memory
    (let ((up (make-real-vector (1+ nuv)))
          (vp (make-real-vector (1+ nuv)))
          (rp (make-index-vector (1+ nuv)))
          (tp (make-index-vector (1+ nuv)))
          (rs (make-index-vector (1+ nuv)))
          (s-lb (make-index-vector (1+ ns)))
          (s-rb (make-index-vector (1+ ns)))
          (s-lt (make-index-vector (1+ ns)))
          (s-rt (make-index-vector (1+ ns)))
          (s-rs (make-index-vector (1+ ns)))
          (status 0) (bstat 0)
          ;; the C's flt and dir live across iterations; srf-is-flat leaves
          ;; them alone where N_SrfIsFlat does not set them
          (flt +no+) (dir 0)
          (top 0) (np 3)
          (minmax 0) (dm 0d0))
      (declare (type real-vector up vp)
               (type index-vector rp tp rs s-lb s-rb s-lt s-rt s-rs)
               (type fixnum status bstat flt dir top np minmax)
               (type double-float dm))
      (fast
        (do ((i 0 (1+ i))) ((> i nuv))
          (declare (type fixnum i))
          (setf (aref rs i) +nostatus+)))
      ;; (2) initialize stacks and output point set
      (setf status (tessellate--st-tess-is-inside-trim x y ht hp m uls urs vbs vts
                                                       0d0 lu 0d0 lv tol))
      (setf (aref s-lb 0) 0
            (aref s-rb 0) 1
            (aref s-lt 0) 2
            (aref s-rt 0) 3)
      (setf (aref s-rs 0) status
            (aref rs 0) status)
      (setf (aref up 0) 0d0 (aref vp 0) 0d0
            (aref up 1) lu  (aref vp 1) 0d0
            (aref up 2) 0d0 (aref vp 2) lv
            (aref up 3) lu  (aref vp 3) lv)
      (setf (aref rp 0) +yes+ (aref tp 0) +yes+
            (aref rp 1) +no+  (aref tp 1) +yes+
            (aref rp 2) +yes+ (aref tp 2) +no+
            (aref rp 3) +no+  (aref tp 3) +no+)
      (setf top 0
            np 3)
      ;; (3) subdivide while stack is not empty
      (do () ((not (>= top 0)))
        (block iteration
          ;; (3.1) pop rectangle stacks
          (let* ((lb (aref s-lb top))
                 (rb (aref s-rb top))
                 (lt (aref s-lt top))
                 (rt (aref s-rt top))
                 (ul (aref up lb))
                 (ur (aref up rt))
                 (vb (aref vp lb))
                 (vt (aref vp rt))
                 (um 0d0) (vm 0d0))
            (declare (type fixnum lb rb lt rt) (type double-float ul ur vb vt um vm))
            (when (or (< (- ur ul) ptol) (< (- vt vb) ptol))
              (nurbs-error :par-err 'tessellate--st-tess-subdivide-srf))
            (setf bstat (aref s-rs top))
            (decf top)
            ;; (3.2) extract surface
            (when (= bstat +outp+)
              (return-from iteration))
            (let ((sura (srf-extract-patch sur ul ur vb vt)))
              (declare (type surface sura))
              ;; (3.3) check flatness; the C ignores N_SrfIsFlat's error, which
              ;; leaves *flt alone when it comes from a closed surface's iso
              ;; curve, and +yes+ otherwise, and *dir alone either way
              (handler-case
                  (multiple-value-bind (f d) (srf-is-flat sura eps)
                    (setf flt f)
                    (when d (setf dir d)))
                (nurbs-error ()
                  (unless (or (srf-is-closed sura +udir+) (srf-is-closed sura +vdir+))
                    (setf flt +yes+)))))
            (when (= flt +yes+)
              (return-from iteration))
            ;; (3.4) surface must be subdivided
            ;; (3.4.1) collect subdivision points
            (when (> (+ np 2) nuv)
              (let ((nnew (+ nuv nuv)))
                (declare (type fixnum nnew))
                (let ((a (make-real-vector (1+ nnew)))) (replace a up :end2 (1+ nuv)) (setf up a))
                (let ((a (make-real-vector (1+ nnew)))) (replace a vp :end2 (1+ nuv)) (setf vp a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a rp :end2 (1+ nuv)) (setf rp a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a tp :end2 (1+ nuv)) (setf tp a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a rs :end2 (1+ nuv)) (setf rs a))
                (do ((i (1+ nuv) (1+ i))) ((> i nnew))
                  (declare (type fixnum i))
                  (setf (aref rs i) +nostatus+))
                (setf nuv nnew)))
            (setf um (* 0.5d0 (+ ul ur)))
            (setf vm (* 0.5d0 (+ vb vt)))
            (if (= dir +udir+)
                (progn
                  (when (> (- vt vm) ctol)
                    (let ((closest lv))
                      (declare (type double-float closest))
                      ;; don't subdivide too close to boundary of trim curve boxes
                      (do ((ii 0 (1+ ii))) ((> ii m))
                        (declare (type fixnum ii))
                        (let ((vbi (svref vbs ii)) (vti (svref vts ii)))
                          (declare (type real-vector vbi vti))
                          (do ((jj 0 (1+ jj))) ((> jj (aref ht ii)))
                            (declare (type fixnum jj))
                            (when (< (abs (- vm (aref vbi jj))) closest)
                              (setf dm (aref vbi jj))
                              (setf closest (abs (- vm dm)))
                              (setf minmax 1))
                            (when (< (abs (- vm (aref vti jj))) closest)
                              (setf dm (aref vti jj))
                              (setf closest (abs (- vm dm)))
                              (setf minmax 2)))))
                      (when (< closest ctol)
                        (if (= minmax 1)
                            (when (> (- dm ctol) vb)
                              (setf vm (- dm ctol)))
                            (when (< (+ dm ctol) vt)
                              (setf vm (+ dm ctol)))))))
                  (setf (aref up (+ np 1)) ul
                        (aref vp (+ np 1)) vm
                        (aref up (+ np 2)) ur
                        (aref vp (+ np 2)) vm)
                  (setf (aref rp (+ np 1)) +yes+
                        (aref tp (+ np 1)) +yes+
                        (aref rp (+ np 2)) +no+
                        (aref tp (+ np 2)) +yes+))
                (progn
                  (when (> (- ur um) ctol)
                    (let ((closest lu))
                      (declare (type double-float closest))
                      ;; don't subdivide too close to boundary of trim curve boxes
                      (do ((ii 0 (1+ ii))) ((> ii m))
                        (declare (type fixnum ii))
                        (let ((uli (svref uls ii)) (uri (svref urs ii)))
                          (declare (type real-vector uli uri))
                          (do ((jj 0 (1+ jj))) ((> jj (aref ht ii)))
                            (declare (type fixnum jj))
                            (when (< (abs (- um (aref uli jj))) closest)
                              (setf dm (aref uli jj))
                              (setf closest (abs (- um dm)))
                              (setf minmax 1))
                            (when (< (abs (- um (aref uri jj))) closest)
                              (setf dm (aref uri jj))
                              (setf closest (abs (- um dm)))
                              (setf minmax 2)))))
                      (when (< closest ctol)
                        (if (= minmax 1)
                            (when (> (- dm ctol) ul)
                              (setf um (- dm ctol)))
                            (when (< (+ dm ctol) ur)
                              (setf um (+ dm ctol)))))))
                  (setf (aref up (+ np 1)) um
                        (aref vp (+ np 1)) vb
                        (aref up (+ np 2)) um
                        (aref vp (+ np 2)) vt)
                  (setf (aref rp (+ np 1)) +yes+
                        (aref tp (+ np 1)) +yes+
                        (aref rp (+ np 2)) +yes+
                        (aref tp (+ np 2)) +no+)))
            ;; (3.4.2) put new parameter squares on the stack
            (when (> (+ top 2) ns)
              (let ((nnew (+ ns ns)))
                (declare (type fixnum nnew))
                (let ((a (make-index-vector (1+ nnew)))) (replace a s-lb :end2 (1+ ns)) (setf s-lb a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a s-rb :end2 (1+ ns)) (setf s-rb a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a s-lt :end2 (1+ ns)) (setf s-lt a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a s-rt :end2 (1+ ns)) (setf s-rt a))
                (let ((a (make-index-vector (1+ nnew)))) (replace a s-rs :end2 (1+ ns)) (setf s-rs a))
                (setf ns nnew)))
            (if (= dir +udir+)
                (setf (aref s-lb (+ top 1)) lb
                      (aref s-rb (+ top 1)) rb
                      (aref s-lt (+ top 1)) (+ np 1)
                      (aref s-rt (+ top 1)) (+ np 2)
                      (aref s-lb (+ top 2)) (+ np 1)
                      (aref s-rb (+ top 2)) (+ np 2)
                      (aref s-lt (+ top 2)) lt
                      (aref s-rt (+ top 2)) rt)
                (setf (aref s-lb (+ top 1)) lb
                      (aref s-rb (+ top 1)) (+ np 1)
                      (aref s-lt (+ top 1)) lt
                      (aref s-rt (+ top 1)) (+ np 2)
                      (aref s-lb (+ top 2)) (+ np 1)
                      (aref s-rb (+ top 2)) rb
                      (aref s-lt (+ top 2)) (+ np 2)
                      (aref s-rt (+ top 2)) rt))
            ;; (3.4.3) classify new rectangles
            (if (= bstat +inp+)
                (setf (aref s-rs (+ top 1)) +inp+
                      (aref s-rs (+ top 2)) +inp+
                      (aref rs lb) +inp+
                      (aref rs (+ np 1)) +inp+)
                (if (= dir +udir+)
                    (progn
                      (setf status (tessellate--st-tess-is-inside-trim
                                    x y ht hp m uls urs vbs vts ul ur vb vm tol))
                      (setf (aref s-rs (+ top 1)) status
                            (aref rs lb) status)
                      (setf status (tessellate--st-tess-is-inside-trim
                                    x y ht hp m uls urs vbs vts ul ur vm vt tol))
                      (setf (aref s-rs (+ top 2)) status
                            (aref rs (+ np 1)) status))
                    (progn
                      (setf status (tessellate--st-tess-is-inside-trim
                                    x y ht hp m uls urs vbs vts ul um vb vt tol))
                      (setf (aref s-rs (+ top 1)) status
                            (aref rs lb) status)
                      (setf status (tessellate--st-tess-is-inside-trim
                                    x y ht hp m uls urs vbs vts um ur vb vt tol))
                      (setf (aref s-rs (+ top 2)) status
                            (aref rs (+ np 1)) status))))
            (setf np (+ np 2))
            (setf top (+ top 2)))))
      ;; hash vertices of subdivision rectangles into hash tables
      (let ((hu (init-real-hash nh))
            (hv (init-real-hash nh))
            (mu -1) (mv -1))
        (declare (type simple-vector hu hv) (type fixnum mu mv))
        (do ((k 0 (1+ k))) ((> k np))
          (declare (type fixnum k))
          (setf mu (update-real-hash hu nh 0d0 (aref up k) lu mu))
          (setf mv (update-real-hash hv nh 0d0 (aref vp k) lv mv)))
        ;; retrieve subdivision parameters
        (let ((us (make-real-vector (1+ mu)))
              (vs (make-real-vector (1+ mv)))
              (tmp (make-array (list (1+ mu) (1+ mv)) :element-type 'fixnum :initial-element 0))
              (smp (make-array (list (1+ mu) (1+ mv)) :element-type 'fixnum :initial-element 0)))
          (declare (type real-vector us vs) (type index-matrix tmp smp))
          (retrieve-real-hash hu nh us)
          (retrieve-real-hash hv nh vs)
          ;; build topology matrix
          (let ((uo (make-real-vector (1+ np)))
                (vo (make-real-vector (1+ np)))
                (ro (make-index-vector (1+ np)))
                (to (make-index-vector (1+ np)))
                (a -1))
            (declare (type real-vector uo vo) (type index-vector ro to) (type fixnum a))
            (fast
              (do ((i 0 (1+ i))) ((> i mu))
                (declare (type fixnum i))
                (do ((j 0 (1+ j))) ((> j mv))
                  (declare (type fixnum j))
                  (setf (aref tmp i j) -1)
                  (setf (aref smp i j) +nostatus+))))
            (do ((k 0 (1+ k))) ((> k np))
              (declare (type fixnum k))
              (let ((i (search-insert-real us mu (aref up k)))
                    (j (search-insert-real vs mv (aref vp k))))
                (declare (type fixnum i j))
                (when (or (< i 0) (< j 0))
                  (nurbs-error :num-err 'tessellate--st-tess-subdivide-srf))
                (let ((l (aref tmp i j)))
                  (declare (type fixnum l))
                  (if (>= l 0)
                      (progn
                        (when (and (= (aref ro l) +no+) (= (aref rp k) +yes+))
                          (setf (aref ro l) +yes+))
                        (when (and (= (aref to l) +no+) (= (aref tp k) +yes+))
                          (setf (aref to l) +yes+)))
                      (progn
                        (incf a)
                        (setf (aref tmp i j) a)
                        (setf (aref uo a) (aref up k)
                              (aref vo a) (aref vp k)
                              (aref ro a) (aref rp k)
                              (aref to a) (aref tp k)))))
                (when (/= (aref rs k) +nostatus+)
                  (setf (aref smp i j) (aref rs k)))))
            ;; compact output
            (setf uo (subseq uo 0 (1+ a))
                  vo (subseq vo 0 (1+ a))
                  ro (subseq ro 0 (1+ a))
                  to (subseq to 0 (1+ a)))
            ;; handle 1 rectangle as special case (+over+)
            (when (= a 3)
              (setf (aref smp 0 0) +over+)
              (let ((tolsq (* tol tol)))
                (declare (type double-float tolsq))
                ;; ensure that all polygon vertices are within bounds
                (do ((i 0 (1+ i))) ((> i m))
                  (declare (type fixnum i))
                  (let ((xi (svref x i)) (yi (svref y i)) (hpi (svref hp i)))
                    (declare (type simple-vector xi yi) (type index-vector hpi))
                    (do ((j 0 (1+ j))) ((> j (aref ht i)))
                      (declare (type fixnum j))
                      (let ((xij (svref xi j)) (yij (svref yi j)))
                        (declare (type real-vector xij yij))
                        (do ((k 0 (1+ k))) ((> k (aref hpi j)))
                          (declare (type fixnum k))
                          (cond ((< (aref xij k) 0d0) (setf (aref xij k) 0d0))
                                ((> (aref xij k) lu) (setf (aref xij k) lu)))
                          (cond ((< (aref yij k) 0d0) (setf (aref yij k) 0d0))
                                ((> (aref yij k) lv) (setf (aref yij k) lv)))
                          ;; check for and eliminate duplicates
                          (when (and (> k 0) (> (aref hpi j) 3))
                            (let* ((d1 (- (aref xij k) (aref xij (1- k))))
                                   (d2 (- (aref yij k) (aref yij (1- k))))
                                   (dist (+ (* d1 d1) (* d2 d2))))
                              (declare (type double-float d1 d2 dist))
                              (when (<= dist tolsq)
                                (do ((l k (1+ l))) ((>= l (aref hpi j)))
                                  (declare (type fixnum l))
                                  (setf (aref xij l) (aref xij (1+ l))
                                        (aref yij l) (aref yij (1+ l))))
                                (decf (aref hpi j))
                                (decf k)))))))))))
            ;; create the output
            (values uo vo ro to a tmp smp mu mv)))))))

;;; ------------------------------------------------------------------ ST_TessIsInsideTrim

(defun tessellate--st-tess-is-inside-trim (u v ht hp n ul ur vb vt xl xr yb yt tol)
  "Classify the rectangle [XL, XR] x [YB, YT] against the trimmed domain of
loops 0..N: U, V the polygons (U[i][0] the outer one of loop i, U[i][1..HT[i]]
its inner ones; see the file header), HP their highest vertex indexes, UL,
UR, VB, VT their bounding boxes, TOL the parameter space tolerance.  Returns
+onandover+, +over+, +inp+, +onp+ or +outp+.
NLib: ST_TessIsInsideTrim (file-private, NL_Tessellate.cpp)."
  (declare (type simple-vector u v hp ul ur vb vt) (type index-vector ht)
           (type fixnum n) (type double-float xl xr yb yt tol))
  (let ((status 0) (on +no+) (over +no+) (in +no+))
    (declare (type fixnum status on over in))
    ;; check box with respect to all loops
    (do ((i 0 (1+ i))) ((> i n))
      (declare (type fixnum i))
      (let ((ui (svref u i)) (vi (svref v i)) (hpi (svref hp i))
            (uli (svref ul i)) (uri (svref ur i)) (vbi (svref vb i)) (vti (svref vt i)))
        (declare (type simple-vector ui vi) (type index-vector hpi)
                 (type real-vector uli uri vbi vti))
        (setf status (tessellate--st-tess-is-inside-polygon
                      (svref ui 0) (svref vi 0) (aref hpi 0)
                      (aref uli 0) (aref uri 0) (aref vbi 0) (aref vti 0)
                      xl xr yb yt tol))
        (cond ((= status +inp+) (setf in +yes+))
              ((= status +onp+) (setf on +yes+))
              ((= status +over+) (setf over +yes+)))
        (do ((j 1 (1+ j))) ((> j (aref ht i)))
          (declare (type fixnum j))
          (setf status (tessellate--st-tess-is-inside-polygon
                        (svref ui j) (svref vi j) (aref hpi j)
                        (aref uli j) (aref uri j) (aref vbi j) (aref vti j)
                        xl xr yb yt tol))
          (when (= status +inp+)
            (return-from tessellate--st-tess-is-inside-trim +outp+))
          (cond ((= status +onp+) (setf on +yes+))
                ((= status +over+) (setf over +yes+))))
        (when (= in +yes+)
          (return))))
    ;; determine box status
    (cond ((= on +yes+)
           (if (= over +yes+) +onandover+ +onp+))
          ((= over +yes+) +over+)
          ((= in +yes+) +inp+)
          (t +outp+))))

;;; ------------------------------------------------------------------ ST_TessIsInsideLoop

(defun tessellate--st-tess-is-inside-loop (u v h n ul ur vb vt xl xr yb yt tol)
  "Classify the rectangle [XL, XR] x [YB, YT] against one trimming loop: U,
V simple-vectors of the polygons' real-vectors (U[0] the outer polygon,
U[1..N] the inner ones), H an index-vector of their highest vertex indexes,
UL, UR, VB, VT real-vectors of their bounding boxes, TOL the parameter space
tolerance.  Returns +onandover+, +over+, +inp+, +onp+ or +outp+.
NLib: ST_TessIsInsideLoop (file-private, NL_Tessellate.cpp)."
  (declare (type simple-vector u v) (type index-vector h) (type fixnum n)
           (type real-vector ul ur vb vt) (type double-float xl xr yb yt tol))
  (let ((status 0) (on +no+) (over +no+))
    (declare (type fixnum status on over))
    ;; check box with respect to the loop
    (setf status (tessellate--st-tess-is-inside-polygon
                  (svref u 0) (svref v 0) (aref h 0)
                  (aref ul 0) (aref ur 0) (aref vb 0) (aref vt 0)
                  xl xr yb yt tol))
    (when (= status +outp+)
      (return-from tessellate--st-tess-is-inside-loop +outp+))
    (when (= status +onp+)
      (setf on +yes+))
    (when (= status +over+)
      (setf over +yes+))
    (do ((i 1 (1+ i))) ((> i n))
      (declare (type fixnum i))
      (setf status (tessellate--st-tess-is-inside-polygon
                    (svref u i) (svref v i) (aref h i)
                    (aref ul i) (aref ur i) (aref vb i) (aref vt i)
                    xl xr yb yt tol))
      (when (= status +inp+)
        (return-from tessellate--st-tess-is-inside-loop +outp+))
      (when (= status +onp+)
        (setf on +yes+))
      (when (= status +over+)
        (setf over +yes+)))
    ;; determine box status
    (cond ((= on +yes+)
           (if (= over +yes+) +onandover+ +onp+))
          ((= over +yes+) +over+)
          (t +inp+))))

;;; ------------------------------------------------------------------ ST_TessIsInsidePolygon

(defun tessellate--st-tess-is-inside-polygon (u v n ul ur vb vt xl xr yb yt tol)
  "Classify the rectangle [XL, XR] x [YB, YT] against the closed polygon
U[0..N], V[0..N] (real-vectors) with bounding box UL, UR, VB, VT, to
parameter space tolerance TOL.  Returns +over+ (the polygon is inside the
rectangle), +inp+ (the rectangle is inside the polygon), +onp+ (they
intersect) or +outp+.
NLib: ST_TessIsInsidePolygon (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum n)
           (type double-float ul ur vb vt xl xr yb yt tol))
  (fast
    (let* ((xlm (- xl tol))
           (xlp (+ xl tol))
           (xrm (- xr tol))
           (xrp (+ xr tol))
           (ybm (- yb tol))
           (ybp (+ yb tol))
           (ytm (- yt tol))
           (ytp (+ yt tol))
           (ovr +yes+)
           (its +no+)
           (klb 0) (krb 0) (klt 0) (krt 0)
           (lb 5) (rb 5) (lt 5) (rt 5))
      (declare (type double-float xlm xlp xrm xrp ybm ybp ytm ytp)
               (type fixnum ovr its klb krb klt krt lb rb lt rt))
      ;; quick box test
      (when (or (< xrp (- ul tol)) (> xlm (+ ur tol)) (< ytp (- vb tol)) (> ybm (+ vt tol)))
        (return-from tessellate--st-tess-is-inside-polygon +outp+))
      ;; check overlapping
      (setf ovr +yes+)
      (do ((i 0 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (when (or (< (aref u i) xlp) (> (aref u i) xrm) (< (aref v i) ybp) (> (aref v i) ytm))
          (setf ovr +no+)))
      (when (= ovr +yes+)
        (return-from tessellate--st-tess-is-inside-polygon +over+))
      ;; check each vertex with respect to the polygon
      (do ((i 1 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (let ((u0 (aref u (1- i))) (u1 (aref u i))
              (v0 (aref v (1- i))) (v1 (aref v i)))
          (declare (type double-float u0 u1 v0 v1))
          ;; check vertical and horizontal line segments
          (when (and (< (abs (- v0 yb)) tol) (< (abs (- v1 yb)) tol))
            (when (or (and (> xl (- u0 tol)) (< xl (+ u1 tol)))
                      (and (> xl (- u1 tol)) (< xl (+ u0 tol))))
              (setf lb +onp+))
            (when (or (and (> xr (- u0 tol)) (< xr (+ u1 tol)))
                      (and (> xr (- u1 tol)) (< xr (+ u0 tol))))
              (setf rb +onp+)))
          (when (and (< (abs (- u0 xr)) tol) (< (abs (- u1 xr)) tol))
            (when (or (and (> yb (- v0 tol)) (< yb (+ v1 tol)))
                      (and (> yb (- v1 tol)) (< yb (+ v0 tol))))
              (setf rb +onp+))
            (when (or (and (> yt (- v0 tol)) (< yt (+ v1 tol)))
                      (and (> yt (- v1 tol)) (< yt (+ v0 tol))))
              (setf rt +onp+)))
          (when (and (< (abs (- v0 yt)) tol) (< (abs (- v1 yt)) tol))
            (when (or (and (> xl (- u0 tol)) (< xl (+ u1 tol)))
                      (and (> xl (- u1 tol)) (< xl (+ u0 tol))))
              (setf lt +onp+))
            (when (or (and (> xr (- u0 tol)) (< xr (+ u1 tol)))
                      (and (> xr (- u1 tol)) (< xr (+ u0 tol))))
              (setf rt +onp+)))
          (when (and (< (abs (- u0 xl)) tol) (< (abs (- u1 xl)) tol))
            (when (or (and (> yb (- v0 tol)) (< yb (+ v1 tol)))
                      (and (> yb (- v1 tol)) (< yb (+ v0 tol))))
              (setf lb +onp+))
            (when (or (and (> yt (- v0 tol)) (< yt (+ v1 tol)))
                      (and (> yt (- v1 tol)) (< yt (+ v0 tol))))
              (setf lt +onp+)))
          ;; check in/out status by ray shooting
          (when (or (/= lb +onp+) (/= rb +onp+))
            (when (or (and (> v1 ybp)
                           (or (< v0 ybm) (< (abs (- v0 yb)) tol)))
                      (and (> v0 ybp)
                           (or (< v1 ybm) (< (abs (- v1 yb)) tol))))
              (let ((xb (+ u1 (/ (* (- yb v1) (- u1 u0)) (- v1 v0)))))
                (declare (type double-float xb))
                (when (/= lb +onp+)
                  (cond ((< (abs (- xb xl)) tol) (setf lb +onp+))
                        ((> xb xl) (incf klb))))
                (when (/= rb +onp+)
                  (cond ((< (abs (- xb xr)) tol) (setf rb +onp+))
                        ((> xb xr) (incf krb)))))))
          (when (or (/= lt +onp+) (/= rt +onp+))
            (when (or (and (> v1 ytp)
                           (or (< v0 ytm) (< (abs (- v0 yt)) tol)))
                      (and (> v0 ytp)
                           (or (< v1 ytm) (< (abs (- v1 yt)) tol))))
              (let ((xt (+ u1 (/ (* (- yt v1) (- u1 u0)) (- v1 v0)))))
                (declare (type double-float xt))
                (when (/= lt +onp+)
                  (cond ((< (abs (- xt xl)) tol) (setf lt +onp+))
                        ((> xt xl) (incf klt))))
                (when (/= rt +onp+)
                  (cond ((< (abs (- xt xr)) tol) (setf rt +onp+))
                        ((> xt xr) (incf krt)))))))))
      ;; define status of corners
      (when (/= lb +onp+)
        (setf lb (if (/= (rem klb 2) 0) +inp+ +outp+)))
      (when (/= rb +onp+)
        (setf rb (if (/= (rem krb 2) 0) +inp+ +outp+)))
      (when (/= lt +onp+)
        (setf lt (if (/= (rem klt 2) 0) +inp+ +outp+)))
      (when (/= rt +onp+)
        (setf rt (if (/= (rem krt 2) 0) +inp+ +outp+)))
      (when (or (= lb +onp+) (= rb +onp+) (= lt +onp+) (= rt +onp+))
        (return-from tessellate--st-tess-is-inside-polygon +onp+))
      ;; check intersection
      (setf its +no+)
      (do ((i 1 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (block segment
          (let ((u0 (aref u (1- i))) (u1 (aref u i))
                (v0 (aref v (1- i))) (v1 (aref v i)))
            (declare (type double-float u0 u1 v0 v1))
            ;; see if segment might intersect box
            (when (or (< (+ (nl-max u0 u1) tol) xlm) (> (- (nl-min u0 u1) tol) xrp))
              (return-from segment))
            (when (or (< (+ (nl-max v0 v1) tol) ybm) (> (- (nl-min v0 v1) tol) ytp))
              (return-from segment))
            ;; see if segment straddles the sides of the box
            (when (or (and (> u1 xlm) (< u0 xlp)) (and (> u0 xlm) (< u1 xlp)))
              (when (or (and (< (abs (- u0 xl)) tol) (and (< v0 ytp) (> v0 ybm)))
                        (and (< (abs (- u1 xl)) tol) (and (< v1 ytp) (> v1 ybm)))
                        (line-segs-2d-are-intersecting u0 v0 u1 v1 xl yb xl yt tol))
                (setf its +yes+)
                (return)))
            (when (or (and (> u1 xrm) (< u0 xrp)) (and (> u0 xrm) (< u1 xrp)))
              (when (or (and (< (abs (- u0 xr)) tol) (and (< v0 ytp) (> v0 ybm)))
                        (and (< (abs (- u1 xr)) tol) (and (< v1 ytp) (> v1 ybm)))
                        (line-segs-2d-are-intersecting u0 v0 u1 v1 xr yb xr yt tol))
                (setf its +yes+)
                (return)))
            (when (or (and (> v1 ybm) (< v0 ybp)) (and (> v0 ybm) (< v1 ybp)))
              (when (or (and (< (abs (- v0 yb)) tol) (and (< u0 xrp) (> u0 xlm)))
                        (and (< (abs (- v1 yb)) tol) (and (< u1 xrp) (> u1 xlm)))
                        (line-segs-2d-are-intersecting u0 v0 u1 v1 xl yb xr yb tol))
                (setf its +yes+)
                (return)))
            (when (or (and (> v1 ytm) (< v0 ytp)) (and (> v0 ytm) (< v1 ytp)))
              (when (or (and (< (abs (- v0 yt)) tol) (and (< u0 xrp) (> u0 xlm)))
                        (and (< (abs (- v1 yt)) tol) (and (< u1 xrp) (> u1 xlm)))
                        (line-segs-2d-are-intersecting u0 v0 u1 v1 xl yt xr yt tol))
                (setf its +yes+)
                (return))))))
      (cond ((= its +yes+) +onp+)
            ((and (= lb +inp+) (= rb +inp+) (= lt +inp+) (= rt +inp+)) +inp+)
            (t +outp+)))))

;;; ------------------------------------------------------------------ ST_TessPolygonizeTrimCrv

(defun tessellate--st-tess-polygonize-trim-crv (sur lu lv trm eps)
  "Polygonize the trimming curve TRM, in the parameter domain of surface SUR
(whose average lengths in u and v are LU and LV), so that no leg deviates
from the curve's image on SUR by more than EPS (EPS is lowered to about 1% of
TRM's box diagonal, scaled by the domain, for small curves).  A degree-1 TRM
on a bilinear SUR whose corner edges are perpendicular gives TRM's control
points.  Returns (values u v n): U, V fresh real-vectors of the polygon's
parameters, N their highest index.  As in the C, U and V may be longer than
N+1 (when returning from inside the subdivision); after the final compaction
they are N+1 long but never shorter than 2, the caller's nt = 0 case writing
index 1.  Signals :tol-err for EPS below *mtol*, :num-err when a segment's
chord cannot be checked and is too long, and as crv-eval,
srf-eval-pt-crv-on-srf, crv-decompose-continuity and
tessellate--st-tess-check-trim do.
NLib: ST_TessPolygonizeTrimCrv (file-private, NL_Tessellate.cpp)."
  (declare (type surface sur) (type double-float lu lv eps) (type curve trm))
  (when (< eps (the double-float *mtol*))
    (nurbs-error :tol-err 'tessellate--st-tess-polygonize-trim-crv))
  (let* ((kt (curve-knots trm))
         (m (curve-m trm))
         (p (curve-degree trm))
         (epl eps)
         (mtol *mtol*)
         (ptol *ptol*))
    (declare (type real-vector kt) (type fixnum m p) (type double-float epl mtol ptol))
    ;; treat degree 1 curve and rectangular plane as special case
    (when (= p 1)
      (let ((nn (surface-n sur)) (mm (surface-m sur))
            (pw (surface-cpoints sur)))
        (declare (type fixnum nn mm) (type cpoint-net pw))
        (when (and (= nn 1) (= mm 1))
          (let ((v1 (make-point)) (v2 (make-point))
                (cw1 (make-cpoint)) (cw2 (make-cpoint))
                (dot 0d0) (d1 0d0) (d2 0d0))
            (declare (dynamic-extent v1 v2 cw1 cw2) (type double-float dot d1 d2))
            (vector-diff-c-pts (%net-c-pt pw 0 0 cw1) (%net-c-pt pw 1 0 cw2) v1)
            (vector-diff-c-pts (%net-c-pt pw 0 0 cw1) (%net-c-pt pw 0 1 cw2) v2)
            (setf dot (vector-dot v1 v2))
            (setf d1 (vector-dot v1 v1))
            (setf d2 (vector-dot v2 v2))
            (setf d1 (c-sqrt (* d1 d2)))
            (unless (float-op-is-bad dot d1 +division+)
              (setf d2 (/ dot d1))     ; cosine
              (when (< (abs d2) 0.002d0) ; consider them perpendicular
                ;; just load the control points of the trim curve
                (let* ((cw (curve-cpoints trm))
                       (nc (curve-n trm))
                       (us (make-real-vector (1+ nc)))
                       (vs (make-real-vector (1+ nc))))
                  (declare (type cpoint-array cw) (type fixnum nc) (type real-vector us vs))
                  (do ((i 0 (1+ i))) ((> i nc))
                    (declare (type fixnum i))
                    (with-euclid-row ((cx cy cz) cw i)
                      (setf (aref us i) cx
                            (aref vs i) cy)))
                  (return-from tessellate--st-tess-polygonize-trim-crv
                    (values us vs nc)))))))))
    ;; check trimming tolerance
    (let ((len 0d0) (fac 0d0) (epa 0d0) (eplh 0d0)
          (cur nil) (k 0)
          (is 0) (ie 0) (tp 0d0) (dt 0d0) (dtg 0d0) (dist 0d0)
          (nh 0) (nhg 0)
          (ps (make-point)) (pe (make-point)) (cpt (make-point)))
      (declare (type double-float len fac epa eplh tp dt dtg dist)
               (type fixnum k is ie nh nhg) (type point ps pe cpt)
               (dynamic-extent ps pe cpt))
      (let ((box (make-real-vector 6)))
        (declare (dynamic-extent box))
        (crv-get-b-box trm box)
        (setf len (b-box-get-diagonal box)))
      (multiple-value-bind (xl xr yb yt) (srf-get-parameter-bounds sur)
        (declare (type double-float xl xr yb yt))
        (setf fac (/ (nl-max (- xr xl) (- yb yt)) (nl-max lu lv))))
      (setf epa (* +tessellate--prc+ len))
      (when (> (* epl fac) epa)
        (setf epl (/ epa fac)))
      (when (< epl mtol)
        (setf epl mtol))
      (setf eplh (/ epl 2.0d0))
      ;; decompose 2d trimming curve into continuous segments
      (multiple-value-setq (cur k) (crv-decompose-continuity trm))
      (let ((cur cur))
        (declare (type simple-vector cur))
        ;; estimate the number of vertices
        (setf is p
              ie (+ p 1)
              tp (aref kt 0)
              len 0d0)
        (srf-eval-pt-crv-on-srf sur trm tp ps)
        (do () ((not (< ie m)))
          (do () ((not (and (< ie m) (= (aref kt ie) (aref kt (1+ ie))))))
            (incf ie))
          (setf dt (/ (- (aref kt ie) (aref kt is)) (float p 1d0)))
          (do ((i 1 (1+ i))) ((> i p))
            (declare (type fixnum i))
            (setf tp (+ tp dt))
            (when (> tp (aref kt m))
              (setf tp (aref kt m)))
            (srf-eval-pt-crv-on-srf sur trm tp pe)
            (setf dist (dist-pt-pt ps pe))
            (setf len (+ len dist))
            (copy-pt pe ps))
          (setf is ie)
          (incf ie))
        (setf nh (nl-max 1 (the fixnum (truncate (/ len (c-sqrt epl))))))
        (setf dtg (/ (- (aref kt m) (aref kt 0)) (float nh 1d0)))
        (if (< nh 10)
            (setf nh (* nh 3))
            (setf nh (* nh 2)))
        (setf nhg (* 2 nh))
        (setf dt dtg)
        (let ((us (make-real-vector (1+ nh)))
              (vs (make-real-vector (1+ nh)))
              (next-seg nil)
              (ts (aref kt 0))
              (te 0d0)
              (mindt 0d0)
              (l 0))
          (declare (type real-vector us vs) (type double-float ts te mindt) (type fixnum l))
          ;; for each curve segment compute polygonal approximation
          (crv-eval trm ts +left+ cpt)
          (multiple-value-bind (cx cy) (pt-to-xyz cpt)
            (declare (type double-float cx cy))
            (multiple-value-bind (cu cv) (clamp-srf-at-params sur cx cy)
              (declare (type double-float cu cv))
              (setf (aref us 0) cu
                    (aref vs 0) cv)))
          (do ((i 0 (1+ i))) ((> i k))
            (declare (type fixnum i))
            (let* ((curi (svref cur i))
                   (ku (curve-knots curi))
                   (r (curve-m curi)))
              (declare (type curve curi) (type real-vector ku) (type fixnum r))
              (setf mindt dtg)
              (when (> i 0)
                (when (and (< (/ dist eplh) 0.001d0) (< (/ dt dtg) 0.001d0))
                  (setf dt dtg)))
              (setf te (+ ts dtg))
              (when (> te (aref ku r))
                (setf te (aref ku r)))
              (do () ((not (and (<= te (aref ku r)) (< ts te))))
                ;; step ahead if guess point is too close
                (block step-ahead
                  (do () (nil)
                    (handler-case
                        (setf dist (tessellate--st-tess-check-trim sur curi ts te))
                      (nurbs-error (c)
                        (cond ((< te (aref ku r))
                               (setf dist 0d0))
                              ((> ts (aref ku 0))
                               (error c))
                              (t
                               ;; the C ignores both evaluations' errors
                               (handler-case (srf-eval-pt-crv-on-srf sur curi ts ps)
                                 (nurbs-error () nil))
                               (handler-case (srf-eval-pt-crv-on-srf sur curi te pe)
                                 (nurbs-error () nil))
                               (setf dist (dist-pt-pt ps pe))
                               (when (> dist epl)
                                 (nurbs-error :num-err 'tessellate--st-tess-polygonize-trim-crv))
                               (when (= k 0)
                                 (return-from tessellate--st-tess-polygonize-trim-crv
                                   (values us vs l)))
                               (if (= i k)
                                   (setf dist (* 0.5d0 (+ epl eplh)))
                                   (progn
                                     (setf next-seg t)
                                     (return-from step-ahead)))))))
                    (when (or (>= dist eplh) (= te (aref ku r)))
                      (return-from step-ahead))
                    (setf te (+ te dt))
                    (when (> te (aref ku r))
                      (setf te (aref ku r)))))
                ;; find point so that epl/2 <= |Ps, Pe| < epl
                (setf dt (- te ts))
                (when (< dt mindt)
                  (setf mindt dt))
                (let ((count 0))
                  (declare (type fixnum count))
                  (block refine
                    (do () (next-seg)
                      (incf count)
                      (if (or (< dist epl) (> count 5000))
                          (if (or (>= dist eplh) (= te (aref ku r)) (> count 5000))
                              (progn
                                (incf l)
                                (when (> l nh)
                                  (let ((a (make-real-vector (+ nh nhg 1))))
                                    (replace a us :end2 (1+ nh))
                                    (setf us a))
                                  (let ((a (make-real-vector (+ nh nhg 1))))
                                    (replace a vs :end2 (1+ nh))
                                    (setf vs a))
                                  (setf nh (+ nh nhg)))
                                (when (<= (- (aref ku r) te) (* 0.01d0 mindt))
                                  (setf te (aref ku r)))
                                (crv-eval trm te +left+ cpt)
                                (multiple-value-bind (cx cy) (pt-to-xyz cpt)
                                  (declare (type double-float cx cy))
                                  (multiple-value-bind (cu cv) (clamp-srf-at-params sur cx cy)
                                    (declare (type double-float cu cv))
                                    (setf (aref us l) cu
                                          (aref vs l) cv)))
                                (when (< te (aref ku r))
                                  (when (and (< (abs (- (aref us l) (aref us (1- l)))) ptol)
                                             (< (abs (- (aref vs l) (aref vs (1- l)))) ptol))
                                    (decf l)))
                                (when (and (= ts (aref ku 0)) (= te (aref ku r)))
                                  (when (and (< (abs (- (aref us l) (aref us (1- l)))) ptol)
                                             (< (abs (- (aref vs l) (aref vs (1- l)))) ptol))
                                    (when (= k 0)
                                      (return-from tessellate--st-tess-polygonize-trim-crv
                                        (values us vs (1- l))))
                                    (decf l)
                                    (when (= i k)
                                      (setf (aref us l) (aref us (1+ l))
                                            (aref vs l) (aref vs (1+ l))))
                                    (return-from refine)))
                                (setf dt (- te ts))
                                (when (< dt mindt)
                                  (setf mindt dt))
                                (return-from refine))
                              (progn
                                (setf dt (/ dt 2.0d0))
                                (when (< dt mindt)
                                  (setf mindt dt))
                                (setf te (+ te dt))
                                (when (> te (aref ku r))
                                  (setf te (aref ku r)))))
                          (progn
                            (setf dt (/ dt 2.0d0))
                            (when (< dt mindt)
                              (setf mindt dt))
                            (setf te (- te dt))
                            (when (< te (aref ku 0))
                              (setf te (aref ku 0)))))
                      (setf dist (tessellate--st-tess-check-trim sur curi ts te)))))
                (setf ts te)
                (setf te (+ te dt))
                (when (> te (aref ku r))
                  (setf te (aref ku r)))
                (setf next-seg nil))))
          ;; compact output arrays
          (when (< l nh)
            (let ((size (1+ (nl-max l 1))))
              (declare (type fixnum size))
              (setf us (subseq us 0 size)
                    vs (subseq vs 0 size))))
          (values us vs l))))))

;;; ------------------------------------------------------------------ ST_TessCheckTrim

(defun tessellate--st-tess-check-trim (sur trm ts te)
  "The deviation of the segment [TS, TE] of trimming curve TRM, mapped onto
surface SUR, from its chord: the larger of the distances of the inner
control points of the cubic Bezier interpolating four equally spaced samples
from the chord of the end samples, and of the distance between SUR at the
midpoint of the domain chord and the cubic's midpoint.  Returns the
deviation.  Signals as crv-eval, srf-eval-pt and bez-interp-4-pts do;
dist-pt-line-seg's error is ignored, as in the C.
NLib: ST_TessCheckTrim (file-private, NL_Tessellate.cpp)."
  (declare (type surface sur) (type curve trm) (type double-float ts te))
  (let ((tv (make-real-vector 4))
        (q (make-point-array 4))
        (pw (make-cpoint-array 4))
        (pp (make-point)) (rp (make-point))
        (q0 (make-point)) (q1 (make-point)) (q2 (make-point)) (q3 (make-point))
        (qi (make-point)) (cw (make-cpoint))
        (dist 0d0) (um 0d0) (vm 0d0) (dev 0d0))
    (declare (type real-vector tv) (type point-array q) (type cpoint-array pw)
             (type point pp rp q0 q1 q2 q3 qi) (type cpoint cw)
             (dynamic-extent tv q pw pp rp q0 q1 q2 q3 qi cw)
             (type double-float dist um vm dev))
    ;; sample segment
    (setf (aref tv 0) ts)
    (setf (aref tv 1) (/ (+ (* 2.0d0 ts) te) 3.0d0))
    (setf (aref tv 2) (/ (+ ts (* 2.0d0 te)) 3.0d0))
    (setf (aref tv 3) te)
    (do ((i 0 (1+ i))) ((> i 3))
      (declare (type fixnum i))
      (crv-eval trm (aref tv i) +left+ pp)
      (multiple-value-bind (u v) (pt-to-xyz pp)
        (declare (type double-float u v))
        (multiple-value-setq (u v) (clamp-srf-at-params sur u v))
        (srf-eval-pt sur u v +left+ +left+ qi)
        (%pt-row qi q i)
        (cond ((= i 0)
               (setf um u
                     vm v))
              ((= i 3)
               (setf um (* 0.5d0 (+ u um))
                     vm (* 0.5d0 (+ v vm)))))))
    ;; get 4-point cubic and check deviation
    (bez-interp-4-pts q pw)
    (c-pt-to-pt-euclid (%row-c-pt pw 0 cw) q0)
    (c-pt-to-pt-euclid (%row-c-pt pw 3 cw) q3)
    (setf dev 0d0)
    (c-pt-to-pt-euclid (%row-c-pt pw 1 cw) q1)
    (handler-case (setf dist (dist-pt-line-seg q1 q0 q3))
      (nurbs-error () nil))
    (when (> dist dev)
      (setf dev dist))
    (c-pt-to-pt-euclid (%row-c-pt pw 2 cw) q2)
    (handler-case (setf dist (dist-pt-line-seg q2 q0 q3))
      (nurbs-error () nil))
    (when (> dist dev)
      (setf dev dist))
    ;; check 3D image of midpoint of uv-chord against midpoint of 3D chord
    (srf-eval-pt sur um vm +left+ +left+ pp)
    (let ((f18 (/ 1d0 8d0))
          (f38 (/ 3d0 8d0)))
      (declare (type double-float f18 f38))
      (combine-4-pts f18 q0 f38 q1 f38 q2 f18 q3 rp))
    (setf dist (dist-pt-pt pp rp))
    (when (> dist dev)
      (setf dev dist))
    dev))

;;; ------------------------------------------------------------------ ST_TessIsOnBoundary

(defun tessellate--st-tess-is-on-boundary (cur tol ul ur vb vt)
  "Whether curve CUR is a line, to within TOL, along a side of the parameter
rectangle [UL, UR] x [VB, VT].  Returns (values bnd side scor ecor): BND
true or false; SIDE +left+, +right+, +bottom+ or +top+ (the last side
matched), or NIL where the C leaves *side unset (BND false); SCOR and ECOR
the corner (+upperleft+, +lowerleft+, +lowerright+, +upperright+ or
+nocorner+) at which the line starts and ends, or NIL where the C leaves
them unset (CUR not a line).
NLib: ST_TessIsOnBoundary (file-private, NL_Tessellate.cpp)."
  (declare (type curve cur) (type double-float tol ul ur vb vt))
  (unless (crv-is-line cur tol)
    (return-from tessellate--st-tess-is-on-boundary (values nil nil nil nil)))
  (let* ((pw (curve-cpoints cur))
         (n (curve-n cur))
         (ax 0d0) (ay 0d0) (bx 0d0) (by 0d0)
         (bnd nil) (side nil)
         (scor +nocorner+) (ecor +nocorner+))
    (declare (type cpoint-array pw) (type fixnum n)
             (type double-float ax ay bx by) (type fixnum scor ecor))
    (with-euclid-row ((x y z) pw 0)
      (setf ax x ay y))
    (with-euclid-row ((x y z) pw n)
      (setf bx x by y))
    (when (and (< (abs (- ax ul)) tol) (< (abs (- bx ul)) tol))
      (setf bnd t
            side +left+)
      (when (< (abs (- ay vt)) tol)
        (setf scor +upperleft+))
      (when (< (abs (- by vb)) tol)
        (setf ecor +lowerleft+)))
    (when (and (< (abs (- ax ur)) tol) (< (abs (- bx ur)) tol))
      (setf bnd t
            side +right+)
      (when (< (abs (- ay vb)) tol)
        (setf scor +lowerright+))
      (when (< (abs (- by vt)) tol)
        (setf ecor +upperright+)))
    (when (and (< (abs (- ay vb)) tol) (< (abs (- by vb)) tol))
      (setf bnd t
            side +bottom+)
      (when (< (abs (- ax ul)) tol)
        (setf scor +lowerleft+))
      (when (< (abs (- bx ur)) tol)
        (setf ecor +lowerright+)))
    (when (and (< (abs (- ay vt)) tol) (< (abs (- by vt)) tol))
      (setf bnd t
            side +top+)
      (when (< (abs (- ax ur)) tol)
        (setf scor +upperright+))
      (when (< (abs (- bx ul)) tol)
        (setf ecor +upperleft+)))
    (values bnd side scor ecor)))

;;; ------------------------------------------------------------------ ST_TessRegionOfRectangles
;;;
;;; The triangulation of the subdivided, trimmed rectangles and its small
;;; helpers.  Its vertex lists DL (NL_INDEX ***, 4 x (nh+1) lists) are a
;;; (simple-array t (* *)) of index-vectors or NIL.

(defconstant %tess-eia 2147483647
  "NL_EIA (NL_DataStruct.h), NL_BIGI = INT_MAX: the end-of-array marker that
closes each of ST_TessRegionOfRectangles' vertex lists.")

(defun tessellate--st-tess-region-of-rectangles
    (u v n rf tf tm sm r s x y ht hp m ul-box ur-box vb-box vt-box
     srf-u-min srf-u-max srf-v-min srf-v-max tol)
  "Triangulate the trimmed region subdivided into rectangles by
ST_TessSubdivideSrf: every rectangle (lower left (i,j) of the topology matrix
TM, 0 <= i < R, 0 <= j < S, status SM[i][j] not +outp+) of every outer loop LP
in 0..M is cut against the trimming polygons X[LP], Y[LP] (HP[LP] their
highest vertex indexes, HT[LP] their highest polygon index, UL-BOX ...
VT-BOX[LP] their bounds) by ST_TessBBoxXTrim, its polygons are triangulated
by an advancing front, and the per-rectangle vertex lists are merged into
the triangulation DT, HD.  U, V (real-vectors, highest index N) are the
points; RF, TF the right and top flags (C's R, T); SRF-U-MIN ... SRF-V-MAX
the surface's parameter bounds; TOL the parameter-space tolerance.
U, V and N are the C's in-out arguments, extended by ST_TessBBoxXTrim.
Returns (values u v n dt hd): the final point arrays and highest index, and
the fresh triangulation -- DT a simple-vector of index-vectors (DT[i][0] = -1
for a point not in it), HD an index-vector of their highest indexes.
Signals :geo-err or :mem-err where the C does, and lets the callees' errors
through.
NLib: ST_TessRegionOfRectangles (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum n r s m)
           (type index-vector rf tf ht) (type index-matrix tm sm)
           (type simple-vector x y hp ul-box ur-box vb-box vt-box)
           (type double-float srf-u-min srf-u-max srf-v-min srf-v-max tol))
  (macrolet ((cfor ((var init test step) &body body)
               ;; C's for (var = init; test; step) on an existing variable,
               ;; which keeps its final value as in C.  (return) is C's
               ;; break, (skip) its continue.
               `(progn (setf ,var ,init)
                       (do () ((not ,test))
                         (block %skip ,@body)
                         ,step)))
             (cwhile (test &body body)
               `(do () ((not ,test))
                  (block %skip ,@body)))
             (skip () '(return-from %skip))
             (err (code)
               `(nurbs-error ,code 'tessellate--st-tess-region-of-rectangles))
             (dlv (lis pt)
               ;; DL[lis][pt]
               `(the index-vector (aref dl ,lis ,pt))))
    (let* ((ns 30) (np 10) (sinc 20) (pin 5) (top -1) (nh n)
           (bi -1) (bs 75) (bn 30)
           (ul u) (vl v)
           (a 0) (b 0) (c 0) (i 0) (j 0) (k -1) (l -1) (lp 0)
           (ind 0) (nl 0) (lo 0) (hid 0) (ni 0) (dpl 0)
           (alf 0) (bet 0) (gam 0) (outer 0) (tt -1)
           (lia 0) (lib 0) (lig 0) (lis 0)
           (way 0) (left 0) (start 0) (join 0)
           (prev 0) (curr 0) (next 0) (cnt -1) (touch 0)
           (inner +no+) (dupl +no+) (lsd nil) (rsd nil)
           (cmin 0d0) (csa 0d0) (xs 0d0) (ys 0d0) (xe 0d0) (ye 0d0)
           (xf 0d0) (yf 0d0) (xl 0d0) (xr 0d0) (yb 0d0) (yt 0d0)
           (xes 0d0) (yes 0d0) (xef 0d0) (yef 0d0) (xsf 0d0) (ysf 0d0)
           (p (vector)) (g (make-index-vector 0)) (ins (make-index-vector 0))
           (idup (make-index-vector 0)) (lid-off 0) (ld-off 0)
           (buf (make-index-vector 0))
           (dl (make-array '(0 0) :initial-element nil))
           (hi (make-array '(0 0) :element-type 'fixnum :initial-element 0))
           (li (make-index-vector 0))
           (used (make-index-vector 0))
           (s-a (make-index-vector 0))
           (s-b (make-index-vector 0)))
      (declare (type fixnum ns np sinc pin top nh bi bs bn a b c i j k l lp
                     ind nl lo hid ni dpl alf bet gam outer tt lia lib lig lis
                     way left start join prev curr next cnt touch inner dupl
                     lid-off ld-off)
               (type real-vector ul vl)
               (type double-float cmin csa xs ys xe ye xf yf xl xr yb yt
                     xes yes xef yef xsf ysf)
               ;; nil after a box ST_TessBBoxXTrim finds empty (tt = -1)
               (type (or null simple-vector) p)
               (type (or null index-vector) g ins)
               (type index-vector idup buf li used s-a s-b)
               (type (simple-array t (* *)) dl)
               (type index-matrix hi))
      (labels ((realloc-iv (vec old new)
                 ;; N_Realloc1dIntArray(&vec, old, new): a fresh array of
                 ;; highest index NEW holding vec[0 .. min(old, new)]
                 (declare (type index-vector vec) (type fixnum old new))
                 (let ((w (make-index-vector (1+ new))))
                   (declare (type index-vector w))
                   (fast
                     (do ((q 0 (1+ q))) ((> q (nl-min old new)))
                       (declare (type fixnum q))
                       (setf (aref w q) (aref vec q))))
                   w))
               (grow (lis pt)
                 ;; HI[lis][pt]++; reallocate DL[lis][pt] at its end marker
                 (declare (type fixnum lis pt))
                 (incf (aref hi lis pt))
                 (setf hid (aref hi lis pt))
                 (when (= (aref (dlv lis pt) hid) %tess-eia)
                   (let ((w (realloc-iv (dlv lis pt) hid (+ hid pin))))
                     (declare (type index-vector w))
                     (setf (aref dl lis pt) w)
                     (tessellate--st-tess-init-array w hid (+ hid pin) 0)
                     (setf (aref w (+ hid pin)) %tess-eia))))
               (add (lis pt target new after)
                 ;; grow DL[lis][pt], find TARGET in it among 0..hid-1, and
                 ;; insert NEW just after it (AFTER) or just before it
                 (declare (type fixnum lis pt target new))
                 (grow lis pt)
                 (setf b -1)
                 (let ((lst (dlv lis pt)))
                   (declare (type index-vector lst))
                   (cfor (a 0 (< a hid) (incf a))
                     (when (= (aref lst a) target)
                       (setf b a)
                       (return)))
                   (when (= b -1)
                     (err :geo-err))
                   (if after
                       (progn
                         (cfor (a hid (>= a (+ b 2)) (decf a))
                           (setf (aref lst a) (aref lst (- a 1))))
                         (setf (aref lst (+ b 1)) new))
                       (progn
                         (cfor (a hid (>= a (+ b 1)) (decf a))
                           (setf (aref lst a) (aref lst (- a 1))))
                         (setf (aref lst b) new)))))
               (grow-stack (by)
                 ;; S_A, S_B reallocated from ns to ns + by
                 (declare (type fixnum by))
                 (setf s-a (realloc-iv s-a ns (+ ns by)))
                 (setf s-b (realloc-iv s-b ns (+ ns by)))
                 (setf ns (+ ns by)))
               (push-edge (ea eb)
                 (declare (type fixnum ea eb))
                 (setf (aref s-a (+ top 1)) ea
                       (aref s-b (+ top 1)) eb)
                 (incf top))
               (new-list (pt)
                 ;; the first free list c of point PT gets fresh memory
                 (declare (type fixnum pt))
                 (cfor (c 0 (<= c 3) (incf c))
                   (when (null (aref dl c pt))
                     (let ((w (make-index-vector (1+ np))))
                       (declare (type index-vector w))
                       (setf (aref dl c pt) w)
                       (tessellate--st-tess-init-array w 0 np 0)
                       (setf (aref w np) %tess-eia)
                       (setf (aref li pt) c))
                     (return)))))
        (declare (inline push-edge))
        ;; Allocate memory and initialize
        (setf c (+ 1 (truncate (nl-max r s) 2)))
        (cfor (a 0 (<= a m) (incf a))
          (cfor (b 0 (<= b (aref ht a)) (incf b))
            (setf nh (+ nh (* (aref (the index-vector (svref hp a)) b) c)))))
        (setf buf (make-index-vector (1+ bs))
              dl (make-array (list 4 (1+ nh)) :initial-element nil)
              hi (make-array (list 4 (1+ nh)) :element-type 'fixnum
                                              :initial-element 0)
              li (make-index-vector (1+ nh))
              used (make-index-vector (1+ nh))
              s-a (make-index-vector (1+ ns))
              s-b (make-index-vector (1+ ns)))
        (cfor (b 0 (<= b nh) (incf b))
          (setf (aref li b) -1
                (aref used b) +no+)
          (cfor (a 0 (<= a 3) (incf a))
            (setf (aref dl a b) nil
                  (aref hi a b) -1)))
        ;; For all rectangles in all loops do
        (cfor (lp 0 (<= lp m) (incf lp))
          (cfor (i 0 (< i r) (incf i))
            (cfor (j 0 (< j s) (incf j))
              ;; Find opposite corner
              (multiple-value-setq (k l)
                (tessellate--st-tess-find-b-box-indices tm rf tf i j))
              (when (and (>= k 0) (>= l 0) (/= (aref sm i j) +outp+))
                ;; See if box is inside the current outer loop
                (when (> m 0)
                  (setf xl (aref ul (aref tm i j))
                        xr (aref ul (aref tm k l))
                        yb (aref vl (aref tm i j))
                        yt (aref vl (aref tm k l)))
                  (when (= (the fixnum
                                (tessellate--st-tess-is-inside-loop
                                 (svref x lp) (svref y lp) (svref hp lp) (aref ht lp)
                                 (svref ul-box lp) (svref ur-box lp)
                                 (svref vb-box lp) (svref vt-box lp)
                                 xl xr yb yt tol))
                           +outp+)
                    (skip)))
                ;; Get trimming polygons
                (multiple-value-setq (u v n buf bi bs p g ins tt)
                  (tessellate--st-tess-b-box-x-trim
                   u v n
                   tm r s i j k l
                   (aref sm i j) (svref x lp) (svref y lp) (svref hp lp) (aref ht lp)
                   (svref ul-box lp) (svref ur-box lp) (svref vb-box lp) (svref vt-box lp)
                   srf-u-min srf-u-max srf-v-min srf-v-max
                   tol
                   buf bi bs bn))
                (setf ul u
                      vl v
                      nl n)
                (when (= tt -1)
                  (skip))
                (when (> nl nh)
                  (setf a (* 2 (- nl nh))
                        b (truncate nh 2)
                        ni (nl-max a b))
                  ;; N_Realloc2dIntPtrArray, N_Realloc2dIntArray: 4 x (nh+ni+1)
                  (let ((ndl (make-array (list 4 (+ nh ni 1)) :initial-element nil))
                        (nhi (make-array (list 4 (+ nh ni 1)) :element-type 'fixnum
                                                              :initial-element 0)))
                    (declare (type (simple-array t (* *)) ndl) (type index-matrix nhi))
                    (dotimes (q 4)
                      (dotimes (w (1+ nh))
                        (setf (aref ndl q w) (aref dl q w)
                              (aref nhi q w) (aref hi q w))))
                    (setf dl ndl
                          hi nhi))
                  (setf li (realloc-iv li nh (+ nh ni)))
                  (setf used (realloc-iv used nh (+ nh ni)))
                  (cfor (b (+ nh 1) (<= b (+ nh ni)) (incf b))
                    (setf (aref used b) +no+
                          (aref li b) -1)
                    (cfor (a 0 (<= a 3) (incf a))
                      (setf (aref dl a b) nil
                            (aref hi a b) -1)))
                  (setf nh (+ nh ni)))
                ;; Postprocess polygons
                (cfor (a 0 (<= a tt) (incf a))
                  (let ((pa (svref p a)))
                    (declare (type index-vector pa))
                    (when (= (aref pa 0) (aref pa (aref g a)))
                      (decf (aref g a)))
                    ;; remove repeated points to avoid an infinite loop later on
                    (cfor (b 0 (< b (aref g a)) (incf b))
                      (when (= (aref pa b) (aref pa (+ b 1)))
                        (cfor (c (+ b 1) (< c (aref g a)) (incf c))
                          (setf (aref pa c) (aref pa (+ c 1))))
                        (decf (aref g a))))))
                (setf inner +no+
                      dpl 0)
                (cfor (a 0 (<= a tt) (incf a))
                  (when (>= (aref ins a) 0)
                    (setf inner +yes+))
                  (when (> (aref g a) dpl)
                    (setf dpl (aref g a))))
                ;; idup, LID = &idup[dpl+1], ld = &LID[dpl+1]: one array, as in C
                (setf idup (make-index-vector (1+ (* 3 (+ dpl 1))))
                      lid-off (+ dpl 1)
                      ld-off (* 2 (+ dpl 1)))
                ;; For all disjoint outer polygons do
                (setf outer 0)
                (cwhile (and (<= outer tt) (= (aref ins outer) -1))
                  (when (<= (aref g outer) 1)
                    (err :geo-err))
                  ;; Allocate memory for outer polygon vertex lists and initialize
                  (setf b (aref g outer))
                  (let ((po (svref p outer)))
                    (declare (type index-vector po))
                    (cfor (a 0 (<= a b) (incf a))
                      (setf bet (aref po a))
                      ;; the C writes used[bet] before this check; an
                      ;; out-of-range bet errors either way
                      (when (or (< bet 0) (> bet nl))
                        (err :geo-err))
                      (setf (aref used bet) +no+)
                      (new-list bet)
                      (if (= a 0)
                          (setf alf (aref po b))
                          (setf alf (aref po (- a 1))))
                      (if (= a b)
                          (setf gam (aref po 0))
                          (setf gam (aref po (+ a 1))))
                      (when (or (< alf 0) (> alf nl))
                        (err :geo-err))
                      (when (or (< gam 0) (> gam nl))
                        (err :geo-err))
                      (when (> c 3)
                        (err :mem-err))
                      (setf (aref (dlv c bet) 0) gam
                            (aref (dlv c bet) 1) alf)
                      (incf (aref hi c bet) 2)))
                  ;; Allocate memory for inner pol vertex lists and initialize
                  (when (= inner +yes+)
                    (cfor (a 0 (<= a tt) (incf a))
                      (when (= (aref ins a) outer)
                        (let ((pa (svref p a)))
                          (declare (type index-vector pa))
                          (cfor (b 0 (<= b (aref g a)) (incf b))
                            (setf bet (aref pa b))
                            (when (or (< bet 0) (> bet nl))
                              (err :geo-err))
                            (setf (aref used bet) +no+)
                            (new-list bet)
                            (if (= b 0)
                                (setf alf (aref pa (aref g a)))
                                (setf alf (aref pa (- b 1))))
                            (if (= b (aref g a))
                                (setf gam (aref pa 0))
                                (setf gam (aref pa (+ b 1))))
                            (when (or (< alf 0) (> alf nl))
                              (err :geo-err))
                            (when (or (< gam 0) (> gam nl))
                              (err :geo-err))
                            (when (> c 3)
                              (err :mem-err))
                            (setf (aref (dlv c bet) 0) gam
                                  (aref (dlv c bet) 1) alf)
                            (incf (aref hi c bet) 2))))))
                  ;; Buffer indexes of multiple points
                  (setf dpl -1)
                  (let ((po (svref p outer)))
                    (declare (type index-vector po))
                    (cfor (a 0 (< a (aref g outer)) (incf a))
                      (cfor (b (+ a 1) (<= b (aref g outer)) (incf b))
                        (when (= (aref po a) (aref po b))
                          (incf dpl)
                          (setf (aref idup dpl) (aref po a))
                          (setf (aref idup (+ ld-off dpl)) a)
                          (setf (aref idup (+ lid-off dpl)) (aref li (aref po a)))))))
                  ;; Put outer polygon vertices on the stack
                  (when (> (+ top (aref g outer) 2) ns)
                    (grow-stack (+ (aref g outer) 2)))
                  (let ((po (svref p outer)))
                    (declare (type index-vector po))
                    (cfor (b (aref g outer) (>= b 0) (decf b))
                      (setf alf (aref po b))
                      (if (= b (aref g outer))
                          (setf bet (aref po 0))
                          (setf bet (aref po (+ b 1))))
                      (push-edge alf bet)))
                  ;; Put all inner polygon vertices on the stack
                  (when (= inner +yes+)
                    (cfor (a 0 (<= a tt) (incf a))
                      (when (= (aref ins a) outer)
                        (when (> (+ top (aref g a) 2) ns)
                          (grow-stack (+ (aref g a) 2)))
                        (let ((pa (svref p a)))
                          (declare (type index-vector pa))
                          (cfor (b 0 (<= b (aref g a)) (incf b))
                            (setf alf (aref pa b))
                            (if (= b (aref g a))
                                (setf bet (aref pa 0))
                                (setf bet (aref pa (+ b 1))))
                            (push-edge alf bet))))))
                  ;; While stack is not empty do
                  (setf (aref used (aref s-a top)) +yes+)
                  (setf (aref used (aref s-b top)) +yes+)
                  (setf cnt -1)
                  (cwhile (>= top 0)
                    ;; Get current edge -> pop stack
                    (setf alf (aref s-a top)
                          bet (aref s-b top)
                          xs (aref ul alf)
                          ys (aref vl alf)
                          xe (aref ul bet)
                          ye (aref vl bet)
                          xes (- xe xs)
                          yes (- ye ys))
                    ;; Search for the third point
                    (setf cmin +bigd+
                          gam -1)
                    (cfor (a 0 (<= a tt) (incf a))
                      (let ((pa (svref p a)))
                        (declare (type index-vector pa))
                        (cfor (b 0 (<= b (aref g a)) (incf b))
                          (setf ind (aref pa b))
                          (when (or (= ind alf) (= ind bet))
                            (skip))
                          (setf xf (aref ul ind)
                                yf (aref vl ind)
                                xsf (- xs xf)
                                ysf (- ys yf))
                          (when (= (the fixnum (tessellate--st-tess-side-test xsf ysf xes yes))
                                   +right+)
                            (skip))
                          (setf xef (- xe xf)
                                yef (- ye yf))
                          (setf csa (tessellate--st-tess-cosine xsf ysf xef yef))
                          (when (>= csa cmin)
                            (skip))
                          (when (tessellate--st-tess-line-x-polygon
                                 xs ys xf yf (- xsf) (- ysf) ul vl p g tt tol)
                            (skip))
                          (when (tessellate--st-tess-line-x-polygon
                                 xf yf xe ye xef yef ul vl p g tt tol)
                            (skip))
                          (setf cmin csa
                                gam ind))))
                    (when (or (< gam 0) (> gam nl))
                      (err :geo-err))
                    ;; Update point lists
                    ;; Determine touch case
                    (multiple-value-setq (touch top)
                      (tessellate--st-tess-edge-touch alf bet gam s-a s-b top used))
                    (setf lia (aref li alf)
                          lib (aref li bet)
                          lig (aref li gam))
                    (when (or (< lia 0) (> lia 3))
                      (err :mem-err))
                    (when (or (< lib 0) (> lib 3))
                      (err :mem-err))
                    (when (or (< lig 0) (> lig 3))
                      (err :mem-err))
                    ;; Handle multiple points
                    (when (>= dpl 0)
                      (setf dupl +no+)
                      (cfor (a 0 (<= a dpl) (incf a))
                        (setf ind (aref idup a))
                        (when (or (= alf ind) (= bet ind) (= gam ind))
                          (setf dupl +yes+
                                c (aref idup (+ lid-off a))
                                lo (aref idup (+ ld-off a)))
                          (return)))
                      (when (= dupl +yes+)
                        (setf b (aref g outer))
                        (let ((po (svref p outer)))
                          (declare (type index-vector po))
                          (if (= lo 0)
                              (setf prev (aref po b))
                              (setf prev (aref po (- lo 1))))
                          (setf curr (aref po lo))
                          (if (= lo b)
                              (setf next (aref po 0))
                              (setf next (aref po (+ lo 1)))))
                        (cond
                          ((= alf ind)
                           (setf lsd (tessellate--st-tess-pt-in-span ul vl prev curr next bet tol)
                                 rsd (tessellate--st-tess-pt-in-span ul vl prev curr next gam tol))
                           (if (and lsd rsd)
                               (setf lia (- c 1))
                               (setf lia c))
                           (when (or (< lia 0) (> lia 3))
                             (err :mem-err)))
                          ((= bet ind)
                           (setf lsd (tessellate--st-tess-pt-in-span ul vl prev curr next alf tol)
                                 rsd (tessellate--st-tess-pt-in-span ul vl prev curr next gam tol))
                           (if (and lsd rsd)
                               (setf lib (- c 1))
                               (setf lib c))
                           (when (or (< lib 0) (> lib 3))
                             (err :mem-err)))
                          ((= gam ind)
                           (setf lsd (tessellate--st-tess-pt-in-span ul vl prev curr next alf tol)
                                 rsd (tessellate--st-tess-pt-in-span ul vl prev curr next bet tol))
                           (if (and lsd rsd)
                               (setf lig (- c 1))
                               (setf lig c))
                           (when (or (< lig 0) (> lig 3))
                             (err :mem-err))))))
                    (setf (aref used gam) +yes+)
                    (cond
                      ((= touch +twotouch+)
                       ;; Skip to the next edge
                       (setf cnt -1)
                       (decf top)
                       (skip))
                      ((= touch +notouch+)
                       ;; Put edge to the bottom of the stack
                       (incf cnt)
                       (when (> cnt (+ top 3))
                         (err :geo-err))
                       (cfor (a top (>= a 1) (decf a))
                         (setf (aref s-a a) (aref s-a (- a 1))
                               (aref s-b a) (aref s-b (- a 1))))
                       (setf (aref s-a 0) alf
                             (aref s-b 0) bet)
                       (skip))
                      ((= touch +regular+)
                       (setf cnt -1)
                       ;; Add gam to alf's list after bet
                       (add lia alf bet gam t)
                       ;; Add gam to bet's list before alf
                       (add lib bet alf gam nil)
                       ;; Append alf and bet to gam's list
                       (incf (aref hi lig gam) 2)
                       (let ((lst (dlv lig gam)))
                         (declare (type index-vector lst))
                         (setf (aref lst 3) (aref lst 1))
                         (setf (aref lst 1) alf)
                         (setf (aref lst 2) bet)))
                      ((= touch +righttouch+)
                       (setf cnt -1)
                       ;; Add gam to bet's list before alf
                       (add lib bet alf gam nil)
                       ;; Add bet to gam's list after alf
                       (add lig gam alf bet t))
                      ((= touch +lefttouch+)
                       (setf cnt -1)
                       ;; Add gam to alf's list after bet
                       (add lia alf bet gam t)
                       ;; Add alf to gam's list before bet
                       (add lig gam bet alf nil)))
                    ;; Put non-boundary edges on the stack
                    (decf top)
                    (when (> (+ top 2) ns)
                      (grow-stack sinc))
                    (when (and (not (tessellate--st-tess-is-edge-on-boundary alf gam p g tt))
                               (/= touch +righttouch+))
                      (push-edge alf gam))
                    (when (and (not (tessellate--st-tess-is-edge-on-boundary gam bet p g tt))
                               (/= touch +lefttouch+))
                      (push-edge gam bet)))
                  (incf outer))))))
        ;; Triangulation of each rectangle is finished -> merge triangulations
        (let ((dout (make-array (1+ nl) :initial-element nil))
              (ho (make-index-vector (1+ nl))))
          (declare (type simple-vector dout) (type index-vector ho))
          (cfor (ind 0 (<= ind nl) (incf ind))
            ;; Get highest index
            (setf (aref ho ind) 0
                  way 0)
            (cfor (lis 0 (<= lis 3) (incf lis))
              (when (>= (aref hi lis ind) 0)
                (setf (aref ho ind) (+ (aref ho ind) (aref hi lis ind)))
                (incf way)))
            ;; Merge lists
            (let ((d (make-index-vector (1+ (aref ho ind)))))
              (declare (type index-vector d))
              (setf (svref dout ind) d)
              (when (= way 0)
                (setf (aref d 0) -1)
                (skip))
              ;; See if lists match up
              (setf start -1)
              (cfor (a 0 (< a way) (incf a))
                (setf left (- way 1)
                      hid (aref hi a ind)
                      join (aref (dlv a ind) hid))
                (cwhile (>= left 1)
                  (cfor (b 0 (< b way) (incf b))
                    (when (= a b)
                      (skip))
                    (when (= (aref (dlv b ind) 0) join)
                      (return)))
                  (when (= b way)
                    (return))
                  (setf hid (aref hi b ind)
                        join (aref (dlv b ind) hid))
                  (decf left))
                (when (< b way)
                  (setf start a)
                  (return)))
              ;; If no match, copy list
              (if (= start -1)
                  (progn
                    (when (> way 1)
                      (err :geo-err))
                    (setf hid (aref hi 0 ind))
                    (let ((lst (dlv 0 ind)))
                      (declare (type index-vector lst))
                      (cfor (b 0 (<= b hid) (incf b))
                        (setf (aref d b) (aref lst b)))))
                  ;; Match found -> hook up lists
                  (progn
                    (setf left (- way 1)
                          hid (aref hi start ind)
                          join (aref (dlv start ind) hid)
                          tt -1)
                    (let ((lst (dlv start ind)))
                      (declare (type index-vector lst))
                      (cfor (b 0 (<= b hid) (incf b))
                        (setf (aref d (incf tt)) (aref lst b))))
                    (cwhile (>= left 1)
                      (cfor (a 0 (< a way) (incf a))
                        (when (= (aref (dlv a ind) 0) join)
                          (return)))
                      (setf hid (aref hi a ind)
                            join (aref (dlv a ind) hid))
                      (let ((lst (dlv a ind)))
                        (declare (type index-vector lst))
                        (cfor (b 1 (<= b hid) (incf b))
                          (setf (aref d (incf tt)) (aref lst b))))
                      (decf left))))))
          (values u v n dout ho))))))

(defun tessellate--st-tess-find-b-box-indices (tm rf tf i j)
  "Given the topology matrix TM, the right and top flags RF and TF (C's R and
T) and the lower left corner (I, J) of a subdivision rectangle, find the
indexes (K, L) of its upper right corner: the first row K > I where TM[K][J]
is a point with a top neighbor and the first column L > J where TM[I][L] is
a point with a right neighbor.  Returns (values k l), both -1 when no box is
anchored at (I, J).
NLib: ST_TessFindBBoxIndices (file-private, NL_Tessellate.cpp)."
  (declare (type index-matrix tm) (type index-vector rf tf) (type fixnum i j))
  (let ((c (aref tm i j))
        (a 0) (b 0) (k -1))
    (declare (type fixnum c a b k))
    ;; Check if box exists at (i, j)
    (when (< c 0)
      (return-from tessellate--st-tess-find-b-box-indices (values -1 -1)))
    (when (or (= (aref rf c) +no+) (= (aref tf c) +no+))
      (return-from tessellate--st-tess-find-b-box-indices (values -1 -1)))
    ;; Locate indexes
    (setf a (+ i 1)
          b j)
    (fast
      (do () ((not (or (< (aref tm a b) 0) (= (aref tf (aref tm a b)) +no+))))
        (incf a)))
    (setf k a)
    (setf a i
          b (+ j 1))
    (fast
      (do () ((not (or (< (aref tm a b) 0) (= (aref rf (aref tm a b)) +no+))))
        (incf b)))
    (values k b)))

(defun tessellate--st-tess-init-array (dl s e v)
  "Set DL[S..E] to V.  Returns DL.
NLib: ST_TessInitArray (file-private, NL_Tessellate.cpp)."
  (declare (type index-vector dl) (type fixnum s e v))
  (fast
    (do ((i s (1+ i))) ((> i e))
      (declare (type fixnum i))
      (setf (aref dl i) v)))
  dl)

(defun tessellate--st-tess-side-test (xsf ysf xes yes)
  "+right+ when the cross product XES*YSF - YES*XSF is >= 0, else +left+.
NLib: ST_TessSideTest (file-private, NL_Tessellate.cpp)."
  (declare (type double-float xsf ysf xes yes))
  (let ((crs (- (* xes ysf) (* yes xsf))))
    (declare (type double-float crs))
    (if (>= crs 0d0)
        +right+
        +left+)))

(defun tessellate--st-tess-pt-in-span (u v p c n k tol)
  "Whether point K of U, V lies in the span at vertex C between the previous
vertex P and the next vertex N: P is to the right of C->K (cross product <=
TOL) and N to its left (cross product >= -TOL).
NLib: ST_TessPtInSpan (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum p c n k) (type double-float tol))
  (let* ((wx (- (aref u k) (aref u c)))
         (wy (- (aref v k) (aref v c)))
         (vx (- (aref u p) (aref u c)))
         (vy (- (aref v p) (aref v c)))
         (crs (- (* vx wy) (* vy wx)))
         (sp (if (<= crs tol) +right+ +left+))
         (sn 0))
    (declare (type double-float wx wy vx vy crs) (type fixnum sp sn))
    (setf vx (- (aref u n) (aref u c))
          vy (- (aref v n) (aref v c))
          crs (- (* vx wy) (* vy wx)))
    (if (>= crs (- tol))
        (setf sn +left+)
        (setf sn +right+))
    (and (= sp +right+) (= sn +left+))))

(defun tessellate--st-tess-cosine (u12 v12 u32 v32)
  "The signed squared cosine of the angle between (U12, V12) and (U32, V32):
ab*|ab| / (|a|^2 |b|^2), ab their dot product.  Signals :num-err when the
division would overflow or underflow.  Returns c.
NLib: ST_TessCosine (file-private, NL_Tessellate.cpp)."
  (declare (type double-float u12 v12 u32 v32))
  (let ((ab (+ (* u32 u12) (* v32 v12)))
        (a (+ (* u32 u32) (* v32 v32)))
        (b (+ (* u12 u12) (* v12 v12)))
        (den 0d0))
    (declare (type double-float ab a b den))
    (if (< ab 0d0)
        (setf ab (- (* ab ab)))
        (setf ab (* ab ab)))
    (setf den (* a b))
    (when (float-op-is-bad ab den +division+)
      (nurbs-error :num-err 'tessellate--st-tess-cosine))
    (/ ab den)))

(defun tessellate--st-tess-is-edge-on-boundary (r s p g tt)
  "Whether the edge (R, S) is an edge, in either direction, of one of the
closed polygons P[0..TT] (index-vectors, highest indexes G).
NLib: ST_TessIsEdgeOnBoundary (file-private, NL_Tessellate.cpp)."
  (declare (type fixnum r s tt) (type simple-vector p) (type index-vector g))
  (fast
    (do ((i 0 (1+ i))) ((> i tt))
      (declare (type fixnum i))
      (let ((pi* (svref p i))
            (gi (aref g i)))
        (declare (type index-vector pi*) (type fixnum gi))
        (do ((j 0 (1+ j))) ((> j gi))
          (declare (type fixnum j))
          (let ((k (if (= j gi) 0 (+ j 1)))
                (l (if (= j 0) gi (- j 1))))
            (declare (type fixnum k l))
            (when (= r (aref pi* j))
              (when (or (= s (aref pi* k)) (= s (aref pi* l)))
                (return-from tessellate--st-tess-is-edge-on-boundary t))))))))
  nil)

(defun tessellate--st-tess-edge-touch (alf bet gam s-a s-b top used)
  "Classify the triangle (ALF, BET, GAM) against the edge stack S-A, S-B
(highest index TOP): an edge ALF-GAM on the stack is removed (a left edge
found: +righttouch+), as is an edge GAM-BET (+lefttouch+); both found is
+twotouch+; neither, +notouch+ when GAM is already USED, else +regular+.
S-A and S-B are modified in place.  Returns (values touch top), TOP the new
highest stack index.
NLib: ST_TessEdgeTouch (file-private, NL_Tessellate.cpp)."
  (declare (type fixnum alf bet gam top) (type index-vector s-a s-b used))
  (let ((touch 0) (ag -1) (gb -1))
    (declare (type fixnum touch ag gb))
    ;; Check for left touch
    (fast
      (do ((i 0 (1+ i))) ((> i top))
        (declare (type fixnum i))
        (when (or (and (= alf (aref s-a i)) (= gam (aref s-b i)))
                  (and (= alf (aref s-b i)) (= gam (aref s-a i))))
          (setf ag i)
          (return)))
      (when (>= ag 0)
        (do ((i (+ ag 1) (1+ i))) ((> i top))
          (declare (type fixnum i))
          (setf (aref s-a (- i 1)) (aref s-a i)
                (aref s-b (- i 1)) (aref s-b i)))
        (decf top))
      ;; Check for right touch
      (do ((i 0 (1+ i))) ((> i top))
        (declare (type fixnum i))
        (when (or (and (= gam (aref s-a i)) (= bet (aref s-b i)))
                  (and (= gam (aref s-b i)) (= bet (aref s-a i))))
          (setf gb i)
          (return)))
      (when (>= gb 0)
        (do ((i (+ gb 1) (1+ i))) ((> i top))
          (declare (type fixnum i))
          (setf (aref s-a (- i 1)) (aref s-a i)
                (aref s-b (- i 1)) (aref s-b i)))
        (decf top)))
    ;; Classify touch
    (cond ((or (>= ag 0) (>= gb 0))
           (cond ((and (>= ag 0) (>= gb 0)) (setf touch +twotouch+))
                 ((>= ag 0) (setf touch +righttouch+))
                 ((>= gb 0) (setf touch +lefttouch+))))
          ((= (aref used gam) +yes+)
           (setf touch +notouch+))
          (t
           (setf touch +regular+)))
    (values touch top)))

(defun tessellate--st-tess-line-x-polygon (x1 y1 x2 y2 x12 y12 u v p g tt tol)
  "Whether the segment (X1, Y1)-(X2, Y2), of direction (X12, Y12), crosses
an edge of one of the closed polygons P[0..TT] (index-vectors into the
points U, V, highest indexes G) properly, by more than TOL; a segment
shorter than 0.01 TOL crosses nothing.
NLib: ST_TessLineXPolygon (file-private, NL_Tessellate.cpp)."
  (declare (type double-float x1 y1 x2 y2 x12 y12 tol)
           (type real-vector u v) (type simple-vector p) (type index-vector g)
           (type fixnum tt))
  (let ((mag12 (c-sqrt (+ (* x12 x12) (* y12 y12)))))
    (declare (type double-float mag12))
    (when (<= mag12 (* 0.01d0 tol))
      (return-from tessellate--st-tess-line-x-polygon nil))
    (let* ((tol2 (* (* 2.0d0 tol) mag12))
           (mnxp (+ (nl-min x1 x2) tol))
           (mnyp (+ (nl-min y1 y2) tol))
           (mxxm (- (nl-max x1 x2) tol))
           (mxym (- (nl-max y1 y2) tol))
           (ltol *ltol*))
      (declare (type double-float tol2 mnxp mnyp mxxm mxym ltol))
      (fast
        (do ((i 0 (1+ i))) ((> i tt))
          (declare (type fixnum i))
          (let ((pi* (svref p i))
                (gi (aref g i)))
            (declare (type index-vector pi*) (type fixnum gi))
            (do ((j 0 (1+ j))) ((> j gi))
              (declare (type fixnum j))
              (block %next
                (let* ((a (aref pi* j))
                       (x3 (aref u a))
                       (y3 (aref v a))
                       (k (if (= j gi) 0 (+ j 1)))
                       (a4 (aref pi* k))
                       (x4 (aref u a4))
                       (y4 (aref v a4)))
                  (declare (type fixnum a k a4) (type double-float x3 y3 x4 y4))
                  (when (< (- (nl-max x3 x4) tol) mnxp)
                    (return-from %next))
                  (when (< mxxm (+ (nl-min x3 x4) tol))
                    (return-from %next))
                  (when (< (- (nl-max y3 y4) tol) mnyp)
                    (return-from %next))
                  (when (< mxym (+ (nl-min y3 y4) tol))
                    (return-from %next))
                  (let* ((x31 (- x1 x3))
                         (y31 (- y1 y3))
                         (x41 (- x1 x4))
                         (y41 (- y1 y4))
                         (cr3 (- (* x12 y31) (* y12 x31)))
                         (cr4 (- (* x12 y41) (* y12 x41))))
                    (declare (type double-float x31 y31 x41 y41 cr3 cr4))
                    (when (> (* cr3 cr4) 0d0)
                      (when (and (> (abs cr3) tol2) (> (abs cr4) tol2))
                        (return-from %next)))
                    (let* ((x43 (- x3 x4))
                           (y43 (- y3 y4))
                           (den (- (* y12 x43) (* x12 y43)))
                           (numb cr3))
                      (declare (type double-float x43 y43 den numb))
                      (when (< (abs den) ltol)
                        (return-from %next))
                      (if (> den 0d0)
                          (let ((denp (- den tol)))
                            (declare (type double-float denp))
                            (when (or (< numb tol) (> numb denp))
                              (return-from %next))
                            (let ((numa (- (* y43 x31) (* x43 y31))))
                              (declare (type double-float numa))
                              (when (or (< numa tol) (> numa denp))
                                (return-from %next))))
                          (let ((denm (+ den tol)))
                            (declare (type double-float denm))
                            (when (or (> numb (- tol)) (< numb denm))
                              (return-from %next))
                            (let ((numa (- (* y43 x31) (* x43 y31))))
                              (declare (type double-float numa))
                              (when (or (> numa (- tol)) (< numa denm))
                                (return-from %next)))))
                      (return-from tessellate--st-tess-line-x-polygon t)))))))))
      nil)))

;;; ------------------------------------------------------------------ ST_IsectLineWith2dPolygon
;;;
;;; Intersecting a side of a box with a trimming polygon, and finding the
;;; crossing points of the polygons with the box sides.

(defun tessellate--st-isect-line-with-2d-polygon (xs ys xe ye flg x y n xl xr yb yt
                                                  srfumin srfumax srfvmin srfvmax tol)
  "Intersect the directed side (XS,YS)-(XE,YE) of a box -- FLG +left+ (the
vertical side on the left, run downwards), +bottom+, +right+ or +top+ -- with
the closed polygon (X[i], Y[i]), i = 0..N (X[N] = X[0]), whose bounding box is
XL, XR, YB, YT, to parameter tolerance TOL; SRFUMIN .. SRFVMAX are the
surface's parameter bounds, at whose corners intersections beyond the side's
ends are dropped.  Vertex intersections and legs overlapping the side are
classified by side tests.  The intersections are sorted along the side (by
decreasing y on the left, increasing x on the bottom, increasing y on the
right, decreasing x on the top) and the side's fixed coordinate written into
them.  Returns (values u v l itp m): fresh real-vectors U, V of the
intersection points, a fresh index-vector L of the legs they lie on, a fresh
index-vector ITP of their types (+inp+ the polygon enters the half space,
+outp+ leaves it, +turnp+ turns on the side, +touchp+ touches it from the
right or overlaps it), and M, their highest index (-1 for none, the arrays
then empty; the C leaves its outputs unwritten).  Signals :cal-err for
another FLG (after the search, as the C does).  As in the C, a leg whose
intersection parameter is a NaN makes the search loop for ever, and a NaN
size estimate (a degenerate box and side) is outside the C's domain (it
casts the NaN to an integer) and signals here.
NLib: ST_IsectLineWith2dPolygon (file-private, NL_Tessellate.cpp)."
  (declare (type double-float xs ys xe ye xl xr yb yt srfumin srfumax srfvmin srfvmax tol)
           (type fixnum flg n) (type real-vector x y))
  (let* ((fac (/ (nl-max (abs (- xe xs)) (abs (- ye ys)))
                 (nl-max (- xr xl) (- yt yb))))
         ;; (NL_INDEX)(NL_MIN(n, fac*n)) + 5: the conditional is a double
         (ne (+ (the fixnum (values (truncate (the double-float
                                                   (nl-min (float n 1d0)
                                                           (* fac (float n 1d0)))))))
                5))
         (ul (make-real-vector (1+ ne)))
         (vl (make-real-vector (1+ ne)))
         (ll (make-index-vector (1+ ne)))
         (itl (make-index-vector (1+ ne)))
         (i 0)
         (mx (- n 1))
         (np -1)
         (xp (- (+ xe ys) ye))
         (yp (- (+ ye xe) xs))
         (xn (- (+ xs ye) ys))
         (yn (- (+ ys xs) xe))
         (lengthb (c-sqrt (+ (* (- xe xs) (- xe xs)) (* (- ye ys) (- ye ys)))))
         (sa 0) (sb 0) (sc 0) (sd 0)
         (output +no+)
         (xi 0d0) (yi 0d0) (a 0d0) (b 0d0)
         (boo nil))
    (declare (type double-float fac xp yp xn yn lengthb xi yi a b)
             (type fixnum ne i mx np sa sb sc sd output)
             (type real-vector ul vl) (type index-vector ll itl))
    (labels ((side (x1 y1 x2 y2 px py)
               (tessellate--st-tess-side-test-2d x1 y1 x2 y2 px py))
             (overlap (px py pos)
               (setf output (tessellate--st-tess-overlap px py xs ys xe ye tol flg pos output)))
             (xline (x1 y1 x2 y2)
               ;; ST_TessLineXLine2d with discarded outputs (xa, ya, c, d)
               (values (tessellate--st-tess-line-x-line-2d x1 y1 x2 y2 xs ys xe ye tol)))
             (add (pu pv li)
               ;; np++, grow the local arrays, store the point and leg
               (incf np)
               (multiple-value-setq (ul vl ll itl ne)
                 (tessellate--st-realloc-arrays ul vl ll itl ne np))
               (setf (aref ul np) pu
                     (aref vl np) pv
                     (aref ll np) li)))
      (declare (inline side overlap xline))
      ;; get all intersections
      (do () ((not (<= i mx)))
        ;; see if current leg intersects line
        (multiple-value-bind (ok x0 y0 a0 b0)
            (tessellate--st-tess-line-x-line-2d (aref x i) (aref y i) (aref x (+ i 1)) (aref y (+ i 1))
                                                xs ys xe ye tol)
          (declare (type double-float x0 y0 a0 b0))
          (setf boo ok)
          (when ok
            (setf xi x0 yi y0 a a0 b b0)))
        (when boo
          (setf b (* b lengthb))      ; modify b from normalized to true length
          ;; don't allow intersections at corners of domain that are
          ;; supposed to be outside
          (when (or (< b (* -0.01d0 tol)) (> (- b lengthb) (* 0.01d0 tol)))
            (cond ((= flg +bottom+)
                   (when (and (= xs srfumin) (< b tol)) (setf boo nil))
                   (when (and (= xe srfumax) (> b lengthb)) (setf boo nil)))
                  ((= flg +top+)
                   (when (and (= xs srfumax) (< b tol)) (setf boo nil))
                   (when (and (= xe srfumin) (> b lengthb)) (setf boo nil)))
                  ((= flg +left+)
                   (when (and (= ys srfvmax) (< b tol)) (setf boo nil))
                   (when (and (= ye srfvmin) (> b lengthb)) (setf boo nil)))
                  (t
                   (when (and (= ys srfvmin) (< b tol)) (setf boo nil))
                   (when (and (= ye srfvmax) (> b lengthb)) (setf boo nil))))))
        (if boo
            ;; classify intersection types
            (let ((lengtha (c-sqrt (+ (* (- (aref x i) (aref x (+ i 1)))
                                         (- (aref x i) (aref x (+ i 1))))
                                      (* (- (aref y i) (aref y (+ i 1)))
                                         (- (aref y i) (aref y (+ i 1))))))))
              (declare (type double-float lengtha))
              (setf a (* a lengtha))  ; modify a from normalized to true length
              (cond
                ((and (>= a tol) (<= a (- lengtha tol)))
                 ;; through intersection
                 (let ((x1 (aref x (+ i 1))) (y1 (aref y (+ i 1)))
                       (x0 (aref x i)) (y0 (aref y i)))
                   (declare (type double-float x1 y1 x0 y0))
                   (setf sa (side xs ys xe ye x1 y1))
                   (setf output +yes+)
                   (cond ((< b tol)
                          (when (= sa +left+)
                            (setf sc (side xs ys xn yn x1 y1))
                            (when (= sc +right+) (setf output +no+)))
                          (when (= sa +right+)
                            (setf sc (side xs ys xn yn x0 y0))
                            (when (= sc +right+) (setf output +no+)))
                          (overlap x1 y1 +start+))
                         ((> b (- lengthb tol))
                          (when (= sa +left+)
                            (setf sc (side xe ye xp yp x1 y1))
                            (when (= sc +right+) (setf output +no+)))
                          (when (= sa +right+)
                            (setf sc (side xe ye xp yp x0 y0))
                            (when (= sc +right+) (setf output +no+)))
                          (overlap x1 y1 +end+)))
                   (when (= output +yes+)
                     (add xi yi (rem i n))
                     (setf (aref itl np) (if (= sa +left+) +inp+ +outp+)))
                   (incf i)))
                ((and (< a tol) (= i 0))
                 ;; corner or overlap intersection at the beginning
                 (let ((x1 (aref x 1)) (y1 (aref y 1))
                       (xm (aref x (- n 1))) (ym (aref y (- n 1))))
                   (declare (type double-float x1 y1 xm ym))
                   (if (xline xm ym (aref x 0) (aref y 0))
                       (progn
                         ;; intersects previous segment -> corner intersection
                         (setf sa (side xs ys xe ye x1 y1))
                         (setf sb (side xs ys xe ye xm ym))
                         (setf output +yes+)
                         (cond ((< b tol)
                                (setf sc (side xs ys xn yn x1 y1))
                                (setf sd (side xs ys xn yn xm ym))
                                (when (and (= sc +right+) (= sd +right+)) (setf output +no+))
                                (overlap x1 y1 +start+))
                               ((> b (- lengthb tol))
                                (setf sc (side xe ye xp yp x1 y1))
                                (setf sd (side xe ye xp yp xm ym))
                                (when (and (= sc +right+) (= sd +right+)) (setf output +no+))
                                (overlap x1 y1 +end+)))
                         (when (= output +yes+)
                           (add xi yi 0)
                           (cond ((and (= sa +right+) (= sb +right+))
                                  (setf (aref itl np) +touchp+))
                                 ((and (>= b tol) (<= b (- lengthb tol)))
                                  (cond ((and (= sa +left+) (= sb +left+))
                                         (setf (aref itl np) +turnp+))
                                        ((= sa +left+)
                                         (setf (aref itl np) +inp+))
                                        ((= sb +left+)
                                         (setf (aref itl np) +outp+))))
                                 (t
                                  (if (= sa sb)
                                      (cond ((and (= sc +left+) (= sd +right+))
                                             (setf (aref itl np) +inp+))
                                            ((and (= sc +right+) (= sd +left+))
                                             (setf (aref itl np) +outp+))
                                            ((and (= sc +left+) (= sd +left+))
                                             (setf (aref itl np) +turnp+)))
                                      (if (and (= sa +left+) (= sc +left+))
                                          (setf (aref itl np) +inp+)
                                          (setf (aref itl np) +outp+)))))))
                       (progn
                         ;; does not intersect previous segment -> overlap case
                         (setf sa (side xs ys xe ye x1 y1))
                         (setf output +yes+)
                         (cond ((< b tol)
                                (setf sc (side xs ys xn yn x1 y1))
                                (when (= sc +right+) (setf output +no+))
                                (overlap x1 y1 +start+))
                               ((> b (- lengthb tol))
                                (setf sc (side xe ye xp yp x1 y1))
                                (when (= sc +right+) (setf output +no+))
                                (overlap x1 y1 +end+)))
                         (when (= output +yes+)
                           (add xi yi 0)
                           (cond ((= sa +left+) (setf (aref itl np) +inp+))
                                 ((= sa +right+) (setf (aref itl np) +touchp+))))))
                   (incf i)
                   (setf mx (- n 2))))
                ((and (< a tol) (/= i 0))
                 ;; overlap case not at the beginning
                 (let ((x1 (aref x (+ i 1))) (y1 (aref y (+ i 1))))
                   (declare (type double-float x1 y1))
                   (setf sa (side xs ys xe ye x1 y1))
                   (setf output +yes+)
                   (cond ((< b tol)
                          (setf sc (side xs ys xn yn x1 y1))
                          (when (= sc +right+) (setf output +no+))
                          (overlap x1 y1 +start+))
                         ((> b (- lengthb tol))
                          (setf sc (side xe ye xp yp x1 y1))
                          (when (= sc +right+) (setf output +no+))
                          (overlap x1 y1 +end+)))
                   (when (= output +yes+)
                     (add xi yi (rem i n))
                     (cond ((= sa +left+) (setf (aref itl np) +inp+))
                           ((= sa +right+) (setf (aref itl np) +touchp+))))
                   (incf i)))
                ((> a (- lengtha tol))
                 ;; corner or overlap not at the beginning
                 (let ((x0 (aref x i)) (y0 (aref y i)))
                   (declare (type double-float x0 y0))
                   (if (and (< i (- n 1))
                            (xline (aref x (+ i 1)) (aref y (+ i 1)) (aref x (+ i 2)) (aref y (+ i 2))))
                       (let ((x2 (aref x (+ i 2))) (y2 (aref y (+ i 2))))
                         (declare (type double-float x2 y2))
                         ;; intersects next segment -> corner intersection
                         (setf sa (side xs ys xe ye x0 y0))
                         (setf sb (side xs ys xe ye x2 y2))
                         (setf output +yes+)
                         (cond ((< b tol)
                                (setf sc (side xs ys xn yn x0 y0))
                                (setf sd (side xs ys xn yn x2 y2))
                                (when (and (= sc +right+) (= sd +right+)) (setf output +no+))
                                (overlap x2 y2 +start+))
                               ((> b (- lengthb tol))
                                (setf sc (side xe ye xp yp x0 y0))
                                (setf sd (side xe ye xp yp x2 y2))
                                (when (and (= sc +right+) (= sd +right+)) (setf output +no+))
                                (overlap x2 y2 +end+)))
                         (when (= output +yes+)
                           (add xi yi (rem (+ i 1) n))
                           (cond ((and (= sa +right+) (= sb +right+))
                                  (setf (aref itl np) +touchp+))
                                 ((and (>= b tol) (<= b (- lengthb tol)))
                                  (cond ((and (= sa +left+) (= sb +left+))
                                         (setf (aref itl np) +turnp+))
                                        ((= sa +left+)
                                         (setf (aref itl np) +outp+))
                                        ((= sb +left+)
                                         (setf (aref itl np) +inp+))))
                                 (t
                                  (if (= sa sb)
                                      (cond ((and (= sc +left+) (= sd +right+))
                                             (setf (aref itl np) +outp+))
                                            ((and (= sc +right+) (= sd +left+))
                                             (setf (aref itl np) +inp+))
                                            ((and (= sc +left+) (= sd +left+))
                                             (setf (aref itl np) +turnp+)))
                                      (if (and (= sb +left+) (= sd +left+))
                                          (setf (aref itl np) +inp+)
                                          (setf (aref itl np) +outp+)))))))
                       (progn
                         ;; does not intersect next segment -> overlap case
                         (setf sa (side xs ys xe ye x0 y0))
                         (setf output +yes+)
                         (cond ((< b tol)
                                (setf sc (side xs ys xn yn x0 y0))
                                (when (= sc +right+) (setf output +no+))
                                (overlap x0 y0 +start+))
                               ((> b (- lengthb tol))
                                (setf sc (side xe ye xp yp x0 y0))
                                (when (= sc +right+) (setf output +no+))
                                (overlap x0 y0 +end+)))
                         (when (= output +yes+)
                           (add xi yi (rem (+ i 1) n))
                           (when (= sa +left+) (setf (aref itl np) +outp+))
                           (when (= sa +right+) (setf (aref itl np) +touchp+)))))
                   (incf i 2)))))
            ;; no intersection -> check overlap case
            (let ((x0 (aref x i)) (y0 (aref y i)))
              (declare (type double-float x0 y0))
              (setf output +no+)
              (cond ((= flg +left+)
                     (when (< (abs (- x0 xs)) tol)
                       (when (and (< ye y0) (< y0 ys)) (setf output +yes+))))
                    ((= flg +bottom+)
                     (when (< (abs (- y0 ys)) tol)
                       (when (and (< xs x0) (< x0 xe)) (setf output +yes+))))
                    ((= flg +right+)
                     (when (< (abs (- x0 xs)) tol)
                       (when (and (< ys y0) (< y0 ye)) (setf output +yes+))))
                    ((= flg +top+)
                     (when (< (abs (- y0 ys)) tol)
                       (when (and (< xe x0) (< x0 xs)) (setf output +yes+)))))
              (when (= output +yes+)
                (add x0 y0 i)
                (if (xline (aref x (- n 1)) (aref y (- n 1)) (aref x 0) (aref y 0))
                    (progn
                      (setf sa (side xs ys xe ye (aref x (- n 1)) (aref y (- n 1))))
                      (if (= sa +left+)
                          (setf (aref itl np) +outp+)
                          (setf (aref itl np) +touchp+))
                      (setf mx (- n 2)))
                    (setf (aref itl np) +touchp+)))
              (incf i)))))
    ;; sort intersection points
    (macrolet ((sort-by (key test fixed fixed-value)
                 ;; insertion sort of KEY (with ll, itl) by TEST, then
                 ;; FIXED[0..np] = FIXED-VALUE
                 `(fast
                    (do ((i 1 (1+ i))) ((> i np))
                      (declare (type fixnum i))
                      (let ((kv (aref ,key i))
                            (k (aref ll i))
                            (ita (aref itl i))
                            (j i))
                        (declare (type double-float kv) (type fixnum k ita j))
                        (do () ((not (and (>= j 1) (,test (aref ,key (- j 1)) kv))))
                          (setf (aref ,key j) (aref ,key (- j 1))
                                (aref ll j) (aref ll (- j 1))
                                (aref itl j) (aref itl (- j 1)))
                          (decf j))
                        (setf (aref ,key j) kv
                              (aref ll j) k
                              (aref itl j) ita)))
                    (do ((i 0 (1+ i))) ((> i np))
                      (declare (type fixnum i))
                      (setf (aref ,fixed i) ,fixed-value)))))
      (cond ((= flg +left+) (sort-by vl < ul xs))      ; by decreasing y
            ((= flg +bottom+) (sort-by ul > vl ys))    ; by increasing x
            ((= flg +right+) (sort-by vl > ul xs))     ; by increasing y
            ((= flg +top+) (sort-by ul < vl ys))       ; by decreasing x
            (t (nurbs-error :cal-err 'tessellate--st-isect-line-with-2d-polygon))))
    ;; create output
    (values (subseq ul 0 (1+ np))
            (subseq vl 0 (1+ np))
            (subseq ll 0 (1+ np))
            (subseq itl 0 (1+ np))
            np)))

(defun tessellate--st-tess-side-test-2d (xs ys xe ye x y)
  "The side of the directed line (XS,YS)-(XE,YE) on which the point (X,Y)
lies: +left+ when the cross product of the line's direction and the vector
to the point is >= 0 (the point on the line included), else +right+ (a NaN
included).
NLib: ST_TessSideTest2d (file-private, NL_Tessellate.cpp)."
  (declare (type double-float xs ys xe ye x y))
  (let* ((vx (- xe xs))
         (vy (- ye ys))
         (wx (- x xs))
         (wy (- y ys))
         (crs (- (* vx wy) (* vy wx))))
    (declare (type double-float vx vy wx wy crs))
    (if (>= crs 0d0) +left+ +right+)))

(defun tessellate--st-tess-overlap (x y xs ys xe ye tol flg pos output)
  "The overlap test of ST_IsectLineWith2dPolygon: OUTPUT (the C's in-out
NL_FLAG), or +no+ when the point (X,Y) is within TOL, along the side FLG
(+left+, +bottom+, +right+ or +top+), of the side's start (POS +start+, the
point (XS,YS)) or end (POS +end+, (XE,YE)).  Returns the new OUTPUT.
NLib: ST_TessOverlap (file-private, NL_Tessellate.cpp)."
  (declare (type double-float x y xs ys xe ye tol) (type fixnum flg pos output))
  (when (= pos +start+)
    (cond ((= flg +left+)
           (when (< (abs (- y ys)) tol) (setf output +no+)))
          ((= flg +bottom+)
           (when (< (abs (- x xs)) tol) (setf output +no+)))
          ((= flg +right+)
           (when (< (abs (- y ys)) tol) (setf output +no+)))
          ((= flg +top+)
           (when (< (abs (- x xs)) tol) (setf output +no+)))))
  (when (= pos +end+)
    (cond ((= flg +left+)
           (when (< (abs (- y ye)) tol) (setf output +no+)))
          ((= flg +bottom+)
           (when (< (abs (- x xe)) tol) (setf output +no+)))
          ((= flg +right+)
           (when (< (abs (- y ye)) tol) (setf output +no+)))
          ((= flg +top+)
           (when (< (abs (- x xe)) tol) (setf output +no+)))))
  output)

(defun tessellate--st-realloc-arrays (u v l f ne np)
  "Room for index NP in the local arrays of ST_IsectLineWith2dPolygon: when
NP exceeds their highest index NE, real-vectors U, V and index-vectors L, F
are replaced by fresh arrays of highest index 2NE holding their elements
0..NE (the rest zero; uninitialized in the C) and NE doubles.  Returns
(values u v l f ne), the C's in-out arguments.
NLib: ST_ReallocArrays (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type index-vector l f) (type fixnum ne np))
  (let ((nl ne))
    (declare (type fixnum nl))
    (when (> np nl)
      (let ((ul (make-real-vector (1+ (+ nl nl))))
            (vl (make-real-vector (1+ (+ nl nl))))
            (ll (make-index-vector (1+ (+ nl nl))))
            (fl (make-index-vector (1+ (+ nl nl)))))
        (declare (type real-vector ul vl) (type index-vector ll fl))
        (fast
          (do ((k 0 (1+ k))) ((> k nl))
            (declare (type fixnum k))
            (setf (aref ul k) (aref u k)
                  (aref vl k) (aref v k)
                  (aref ll k) (aref l k)
                  (aref fl k) (aref f k))))
        (setf u ul v vl l ll f fl)
        (setf nl (+ nl nl))))
    (values u v l f nl)))

(defun tessellate--st-tess-line-x-line-2d (x1 y1 x2 y2 x3 y3 x4 y4 tol)
  "Intersect the 2D segments (X1,Y1)-(X2,Y2) and (X3,Y3)-(X4,Y4), to
parameter tolerance TOL (a parameter outside [0,1] by at most TOL in length
counts).  Segments within *ltol* of parallel (coincident ones included) do
not intersect.  Returns (values intersect-p x y alf bet): true when they
intersect, the intersection point (the mean of the two segments' points) and
the parameters, (x, y) = P1 + alf (P2 - P1) = P3 + bet (P4 - P3).  When
there is no intersection the outputs the C leaves unwritten are 0d0 (ALF is
written by the C, and returned, once it is computed).
NLib: ST_TessLineXLine2d (file-private, NL_Tessellate.cpp)."
  (declare (type double-float x1 y1 x2 y2 x3 y3 x4 y4 tol))
  (let ((alf 0d0) (bet 0d0))
    (declare (type double-float alf bet))
    (flet ((no () (return-from tessellate--st-tess-line-x-line-2d
                    (values nil 0d0 0d0 alf bet))))
      ;; do bounding box check
      (when (<= (+ (nl-max x3 x4) tol) (- (nl-min x1 x2) tol)) (no))
      (when (<= (+ (nl-max x1 x2) tol) (- (nl-min x3 x4) tol)) (no))
      (when (<= (+ (nl-max y3 y4) tol) (- (nl-min y1 y2) tol)) (no))
      (when (<= (+ (nl-max y1 y2) tol) (- (nl-min y3 y4) tol)) (no))
      ;; check intersection
      (let* ((ax (- x2 x1))
             (ay (- y2 y1))
             (bx (- x3 x4))
             (by (- y3 y4))
             (cx (- x1 x3))
             (cy (- y1 y3))
             (den (- (* ay bx) (* ax by))))
        (declare (type double-float ax ay bx by cx cy den))
        (when (< (abs den) (the double-float *ltol*)) (no))
        (let ((numa (- (* by cx) (* bx cy)))
              (numb (- (* ax cy) (* ay cx)))
              (d 0d0))
          (declare (type double-float numa numb d))
          (when (float-op-is-bad numa den +division+) (no))
          (setf alf (/ numa den))
          (setf d 0d0)
          (let ((lengtha (c-sqrt (+ (* ax ax) (* ay ay)))))
            (declare (type double-float lengtha))
            (cond ((< alf 0d0) (setf d (* (- alf) lengtha)))
                  ((> alf 1d0) (setf d (* (- alf 1d0) lengtha)))))
          (when (> d tol) (no))
          (when (float-op-is-bad numb den +division+) (no))
          (setf bet (/ numb den))
          (setf d 0d0)
          (let ((lengthb (c-sqrt (+ (* bx bx) (* by by)))))
            (declare (type double-float lengthb))
            (cond ((< bet 0d0) (setf d (* (- bet) lengthb)))
                  ((> bet 1d0) (setf d (* (- bet 1d0) lengthb)))))
          (when (> d tol) (no))
          (let ((xa (+ x1 (* alf ax)))
                (ya (+ y1 (* alf ay)))
                (xb (- x3 (* bet bx)))
                (yb (- y3 (* bet by))))
            (declare (type double-float xa ya xb yb))
            (values t (* 0.5d0 (+ xa xb)) (* 0.5d0 (+ ya yb)) alf bet)))))))

(defun tessellate--st-tess-find-crossing-pt (li lp ll ls kl bi bp bl bs kb ri rp rl rs kr
                                             ti tp tl ts kt p h l s)
  "The box-side intersection where trimming polygon P, run from its leg L,
crosses out of the box, other than the one with point index S.  The four
sides' lists (left, bottom, right, top: point indexes LI, polygons LP, legs
LL, types LS, highest index KL; and likewise B.., R.., T..) are searched in
that order first for an intersection on leg L; then all intersections on the
next leg, (L+1) mod H[P], are gathered: the first is taken when there is one
or they all have the same point index, else the first +turnp+ among them.
Returns (values side sin pin): the side (+left+ .. +top+), the index in that
side's lists and the point index.  Signals :geo-err when none is found (the
C returns NL_YES, on which its one caller sets NL_GEO_ERR).  More than four
intersections on the next leg overrun the C's arrays (outside its domain);
here they signal an error.
NLib: ST_TessFindCrossingPt (file-private, NL_Tessellate.cpp)."
  (declare (type index-vector li lp ll ls bi bp bl bs ri rp rl rs ti tp tl ts h)
           (type fixnum kl kb kr kt p l s))
  (let ((found nil) (side 0) (sin 0) (pin 0))
    (declare (type fixnum side sin pin))
    ;; first try to find on current leg
    (macrolet ((scan (ii pp lg k sd)
                 `(unless found
                    (fast
                      (do ((r 0 (1+ r))) ((> r ,k))
                        (declare (type fixnum r))
                        (when (and (= (aref ,pp r) p) (= (aref ,lg r) l))
                          (when (/= (aref ,ii r) s)
                            (setf sin r
                                  pin (aref ,ii r)
                                  side ,sd
                                  found t)
                            (return))))))))
      (scan li lp ll kl +left+)
      (scan bi bp bl kb +bottom+)
      (scan ri rp rl kr +right+)
      (scan ti tp tl kt +top+))
    ;; now try the next leg to account for turning points
    (unless found
      (let ((n -1)
            (psa (make-index-vector 4))
            (sda (make-index-vector 4))
            (pia (make-index-vector 4))
            (sia (make-index-vector 4)))
        (declare (type fixnum n) (type index-vector psa sda pia sia))
        (macrolet ((gather (ii pp lg st k sd)
                     `(do ((r 0 (1+ r))) ((> r ,k))
                        (declare (type fixnum r))
                        (when (and (= (aref ,pp r) p)
                                   (= (aref ,lg r) (rem (+ l 1) (aref h p))))
                          (when (/= (aref ,ii r) s)
                            (incf n)
                            (setf (aref sia n) r
                                  (aref pia n) (aref ,ii r)
                                  (aref sda n) ,sd
                                  (aref psa n) (aref ,st r)))))))
          (gather li lp ll ls kl +left+)
          (gather bi bp bl bs kb +bottom+)
          (gather ri rp rl rs kr +right+)
          (gather ti tp tl ts kt +top+))
        (setf found (>= n 0))
        (when found
          (let ((same t))
            (do ((r 0 (1+ r))) ((>= r n))
              (declare (type fixnum r))
              (when (/= (aref pia r) (aref pia (+ r 1)))
                (setf same nil)
                (return)))
            (if (or (= n 0) same)
                (setf sin (aref sia 0)
                      pin (aref pia 0)
                      side (aref sda 0))
                (let ((m -1))
                  (declare (type fixnum m))
                  (do ((r 0 (1+ r))) ((> r n))
                    (declare (type fixnum r))
                    (when (= (aref psa r) +turnp+)
                      (setf m r)
                      (return)))
                  (if (= m -1)
                      (setf found nil)
                      (setf sin (aref sia m)
                            pin (aref pia m)
                            side (aref sda m)))))))))
    (unless found
      (nurbs-error :geo-err 'tessellate--st-tess-find-crossing-pt))
    (values side sin pin)))

(defun tessellate--st-tess-find-unused (li ls lu lp ll kl bi bs bu bp bl kb ri rs ru rp rl kr
                                        ti ts tu tp tl kt)
  "The first box-side intersection of type +inp+ not yet used (state
+notused+), searching the left, bottom, right and top sides' lists in that
order (point indexes LI, types LS, use states LU, polygons LP, legs LL,
highest index KL; and likewise B.., R.., T..), or failing that the first
unused +turnp+.  Returns (values p l side sin pin): its polygon and leg, the
side (+left+ .. +top+), its index in that side's lists and its point index.
Signals :geo-err when there is none (the C returns NL_YES, on which its one
caller sets NL_GEO_ERR).
NLib: ST_TessFindUnused (file-private, NL_Tessellate.cpp)."
  (declare (type index-vector li ls lu lp ll bi bs bu bp bl ri rs ru rp rl ti ts tu tp tl)
           (type fixnum kl kb kr kt))
  (let ((found nil) (p 0) (l 0) (side 0) (sin 0) (pin 0))
    (declare (type fixnum p l side sin pin))
    (macrolet ((scan (state ii st us pp lg k sd)
                 `(unless found
                    (fast
                      (do ((r 0 (1+ r))) ((> r ,k))
                        (declare (type fixnum r))
                        (when (and (= (aref ,st r) ,state) (= (aref ,us r) +notused+))
                          (setf p (aref ,pp r)
                                l (aref ,lg r)
                                sin r
                                pin (aref ,ii r)
                                side ,sd
                                found t)
                          (return)))))))
      ;; search for an unused NL_INP first
      (scan +inp+ li ls lu lp ll kl +left+)
      (scan +inp+ bi bs bu bp bl kb +bottom+)
      (scan +inp+ ri rs ru rp rl kr +right+)
      (scan +inp+ ti ts tu tp tl kt +top+)
      ;; now search for an unused NL_TURNP
      (scan +turnp+ li ls lu lp ll kl +left+)
      (scan +turnp+ bi bs bu bp bl kb +bottom+)
      (scan +turnp+ ri rs ru rp rl kr +right+)
      (scan +turnp+ ti ts tu tp tl kt +top+))
    (unless found
      (nurbs-error :geo-err 'tessellate--st-tess-find-unused))
    (values p l side sin pin)))

;;; ------------------------------------------------------------------ ST_TessBBoxXTrim
;;;
;;; ST_TessBBoxXTrim and the small statics it alone calls.

(defun tessellate--st-tess-b-box-x-trim (u v vtx-count tm rtm stm i j k l sta
                                         x y h m pul pur pvb pvt
                                         srf-u-min srf-u-max srf-v-min srf-v-max
                                         tol buf bi bs bn)
  "Intersect the box of the subdivision with lower left corner TM[I][J] and
upper right corner TM[K][L] (RTM, STM the highest indexes of TM) with the
trimming polygons X[a], Y[a] (real-vectors, a = 0..M, X[0] the outer loop,
H[a] the highest index of polygon a, PUL PUR PVB PVT their bounding boxes).
STA is the box's status: +inp+, +over+, +onp+ or +onandover+.  U, V
(real-vectors) are the parameter points so far, VTX-COUNT their highest
index; points of polygon vertices and intersections are appended to them.
BUF (an index-vector), BI and BS are the intersection-point buffer, its
highest used index and its highest allocated index, BN its growth
increment.  Returns (values u v vtx-count buf bi bs p g ins t): the point
arrays (fresh when grown) and their highest index, the buffer and its two
indexes, then the output polygons P (a simple-vector of index-vectors of
point indexes, valid material on the left), G their highest indexes, INS
(-1 for an outer polygon, else the index of the polygon containing it) and
T, the highest polygon index.  When T is -1, P, G and INS are nil (the C
leaves them unset).  Signals :geo-err where the C does.
NLib: ST_TessBBoxXTrim (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum vtx-count rtm stm i j k l sta m bn bi bs)
           (type index-matrix tm) (type simple-vector x y) (type index-vector h buf)
           (type real-vector pul pur pvb pvt)
           (type double-float srf-u-min srf-u-max srf-v-min srf-v-max tol))
  (let* ((nh vtx-count)
         (nl nh)
         (ul u)
         (vl v)
         (i00 (aref tm i j))
         (i11 (aref tm k l))
         (bl (aref ul i00))
         (br (aref ul i11))
         (bb (aref vl i00))
         (bt (aref vl i11))
         (blt (+ bl tol))
         (brt (- br tol))
         (bbt (+ bb tol))
         (btt (- bt tol))
         (ni -1)
         (s 0)
         (pin 0))
    (declare (type fixnum nh nl i00 i11 ni s pin) (type real-vector ul vl)
             (type double-float bl br bb bt blt brt bbt btt))
    (do ((a 0 (1+ a))) ((> a m))
      (declare (type fixnum a))
      (setf ni (+ ni (+ (aref h a) 1))))
    (setf ni (* ni 2))
    (macrolet ((done (p g ins tt)
                 `(return-from tessellate--st-tess-b-box-x-trim
                    (values ul vl nl buf bi bs ,p ,g ,ins ,tt)))
               (append-pt (xf yf)
                 ;; ST_TessAppendPt(&ul, &vl, &nh, ni, xf, yf, &nl, SG)
                 `(multiple-value-setq (ul vl nh nl)
                    (tessellate--st-tess-append-pt ul vl nh ni ,xf ,yf nl))))
      (flet ((box-polygon ()
               ;; the box's own points, counterclockwise from the upper left
               ;; corner: (values index-vector highest-index)
               (let ((r -1))
                 (declare (type fixnum r))
                 (do ((b l (1- b))) ((<= b j))
                   (declare (type fixnum b))
                   (when (>= (aref tm i b) 0) (incf r)))
                 (do ((a i (1+ a))) ((>= a k))
                   (declare (type fixnum a))
                   (when (>= (aref tm a j) 0) (incf r)))
                 (do ((b j (1+ b))) ((>= b l))
                   (declare (type fixnum b))
                   (when (>= (aref tm k b) 0) (incf r)))
                 (do ((a k (1- a))) ((<= a i))
                   (declare (type fixnum a))
                   (when (>= (aref tm a l) 0) (incf r)))
                 (let ((q (make-index-vector (1+ r))))
                   (declare (type index-vector q))
                   (setf r -1)
                   (do ((b l (1- b))) ((<= b j))
                     (declare (type fixnum b))
                     (when (>= (aref tm i b) 0)
                       (setf (aref q (incf r)) (aref tm i b))))
                   (do ((a i (1+ a))) ((>= a k))
                     (declare (type fixnum a))
                     (when (>= (aref tm a j) 0)
                       (setf (aref q (incf r)) (aref tm a j))))
                   (do ((b j (1+ b))) ((>= b l))
                     (declare (type fixnum b))
                     (when (>= (aref tm k b) 0)
                       (setf (aref q (incf r)) (aref tm k b))))
                   (do ((a k (1- a))) ((<= a i))
                     (declare (type fixnum a))
                     (when (>= (aref tm a l) 0)
                       (setf (aref q (incf r)) (aref tm a l))))
                   (values q r))))
             (trim-polygon (a)
               ;; trimming polygon A, its vertices 0..H[a]-1 appended to the
               ;; points: an index-vector of their indexes
               (declare (type fixnum a))
               (let ((q (make-index-vector (aref h a)))
                     (xa (svref x a))
                     (ya (svref y a)))
                 (declare (type index-vector q) (type real-vector xa ya))
                 (do ((r 0 (1+ r))) ((>= r (aref h a)))
                   (declare (type fixnum r))
                   (append-pt (aref xa r) (aref ya r))
                   (setf (aref q r) nl))
                 q)))

        ;; Box is entirely inside the domain
        (when (= sta +inp+)
          (let ((po (make-array 1 :initial-element nil))
                (go (make-index-vector 1))
                (ino (make-index-vector 1)))
            (multiple-value-bind (q r) (box-polygon)
              (setf (svref po 0) q
                    (aref go 0) r
                    (aref ino 0) -1))
            (done po go ino 0)))

        ;; Box contains trimming polygons
        (when (= sta +over+)
          (let ((pfl (make-index-vector (1+ m))))
            (declare (type index-vector pfl))
            (setf s -1)
            (do ((a 0 (1+ a))) ((> a m))
              (declare (type fixnum a))
              (cond ((and (= rtm 1) (= stm 1))
                     ;; special case only 1 rectangle: every polygon in
                     (incf s)
                     (setf (aref pfl a) +yes+))
                    ((and (>= (aref pul a) blt) (<= (aref pur a) brt)
                          (>= (aref pvb a) bbt) (<= (aref pvt a) btt))
                     (incf s)
                     (setf (aref pfl a) +yes+))
                    (t (setf (aref pfl a) +no+))))
            (when (/= (aref pfl 0) +yes+)
              (incf s))
            (let ((po (make-array (1+ s) :initial-element nil))
                  (go (make-index-vector (1+ s)))
                  (ino (make-index-vector (1+ s))))
              (declare (type simple-vector po) (type index-vector go ino))
              ;; output the boundary first
              (setf s -1)
              (when (/= (aref pfl 0) +yes+)
                (incf s)
                (multiple-value-bind (q r) (box-polygon)
                  (setf (svref po s) q
                        (aref go s) r
                        (aref ino s) -1)))
              ;; output the inner loops next
              (do ((a 0 (1+ a))) ((> a m))
                (declare (type fixnum a))
                (when (= (aref pfl a) +yes+)
                  (incf s)
                  (setf (svref po s) (trim-polygon a))
                  (setf (aref go s) (- (aref h a) 1))
                  (if (= a 0)
                      (setf (aref ino s) -1)
                      (setf (aref ino s) 0))))
              (done po go ino s))))

        ;; Box intersects polygons and may contain some others
        (unless (or (= sta +onp+) (= sta +onandover+))
          ;; the C goes on to read its unset gl (NULL) here
          (error "ST_TessBBoxXTrim: box status ~d is outside the C's domain." sta))
        (let* ((aa (make-real-vector (1+ m)))
               (pfl (make-index-vector (1+ m)))
               (im (+ (+ (* 2 (- k i)) (* 2 (- l j))) 1))
               (ba (* (- br bl) (- bt bb)))
               (ri 0)
               (xl (make-array (1+ m) :initial-element nil))
               (yl (make-array (1+ m) :initial-element nil))
               (xr (make-array (1+ m) :initial-element nil))
               (yr (make-array (1+ m) :initial-element nil))
               (xb (make-array (1+ m) :initial-element nil))
               (yb (make-array (1+ m) :initial-element nil))
               (xt (make-array (1+ m) :initial-element nil))
               (yt (make-array (1+ m) :initial-element nil))
               (ml (make-index-vector (1+ m)))
               (mr (make-index-vector (1+ m)))
               (mb (make-index-vector (1+ m)))
               (mt (make-index-vector (1+ m)))
               (ll (make-array (1+ m) :initial-element nil))
               (lr (make-array (1+ m) :initial-element nil))
               (lb (make-array (1+ m) :initial-element nil))
               (lt (make-array (1+ m) :initial-element nil))
               (itl (make-array (1+ m) :initial-element nil))
               (itr (make-array (1+ m) :initial-element nil))
               (itb (make-array (1+ m) :initial-element nil))
               (itt (make-array (1+ m) :initial-element nil))
               (kl -1) (kb -1) (kr -1) (kt -1))
          (declare (type real-vector aa) (type index-vector pfl ml mr mb mt)
                   (type fixnum im ri kl kb kr kt) (type double-float ba)
                   (type simple-vector xl yl xr yr xb yb xt yt ll lr lb lt itl itr itb itt))
          ;; estimate the highest index in each polygon
          (setf s 0)
          (do ((a 0 (1+ a))) ((> a m))
            (declare (type fixnum a))
            (let ((fac (/ 1d0 (* (- (aref pur a) (aref pul a)) (- (aref pvt a) (aref pvb a))))))
              (declare (type double-float fac))
              (setf (aref aa a) (b-box-2d-calc-overlap bl br bb bt (aref pul a) (aref pur a)
                                                       (aref pvb a) (aref pvt a) tol))
              (when (>= (aref aa a) 0d0)
                (let ((b (%tess-c-index (+ (* (/ (- ba (aref aa a)) ba) (float im 1d0)) 1d0)))
                      (r (%tess-c-index (+ (* (* (aref aa a) fac) (float (aref h a) 1d0)) 1d0))))
                  (declare (type fixnum b r))
                  (when (> (nl-max r b) s)
                    (setf s (nl-max r b)))))))
          (setf ri (+ s 1))

          ;; get intersection with each side
          (do ((a 0 (1+ a))) ((> a m))
            (declare (type fixnum a))
            (setf (aref ml a) -1 (aref mb a) -1 (aref mr a) -1 (aref mt a) -1)
            (when (> (aref aa a) (- tol))
              (let ((xa (svref x a)) (ya (svref y a)))
                (multiple-value-bind (pu pv pl pf pm)
                    (tessellate--st-isect-line-with-2d-polygon
                     bl bt bl bb +left+ xa ya (aref h a) (aref pul a) (aref pur a) (aref pvb a) (aref pvt a)
                     srf-u-min srf-u-max srf-v-min srf-v-max tol)
                  (setf (svref xl a) pu (svref yl a) pv (svref ll a) pl (svref itl a) pf
                        (aref ml a) pm))
                (multiple-value-bind (pu pv pl pf pm)
                    (tessellate--st-isect-line-with-2d-polygon
                     bl bb br bb +bottom+ xa ya (aref h a) (aref pul a) (aref pur a) (aref pvb a) (aref pvt a)
                     srf-u-min srf-u-max srf-v-min srf-v-max tol)
                  (setf (svref xb a) pu (svref yb a) pv (svref lb a) pl (svref itb a) pf
                        (aref mb a) pm))
                (multiple-value-bind (pu pv pl pf pm)
                    (tessellate--st-isect-line-with-2d-polygon
                     br bb br bt +right+ xa ya (aref h a) (aref pul a) (aref pur a) (aref pvb a) (aref pvt a)
                     srf-u-min srf-u-max srf-v-min srf-v-max tol)
                  (setf (svref xr a) pu (svref yr a) pv (svref lr a) pl (svref itr a) pf
                        (aref mr a) pm))
                (multiple-value-bind (pu pv pl pf pm)
                    (tessellate--st-isect-line-with-2d-polygon
                     br bt bl bt +top+ xa ya (aref h a) (aref pul a) (aref pur a) (aref pvb a) (aref pvt a)
                     srf-u-min srf-u-max srf-v-min srf-v-max tol)
                  (setf (svref xt a) pu (svref yt a) pv (svref lt a) pl (svref itt a) pf
                        (aref mt a) pm)))
              (setf kl (+ kl (+ (aref ml a) 1))
                    kb (+ kb (+ (aref mb a) 1))
                    kr (+ kr (+ (aref mr a) 1))
                    kt (+ kt (+ (aref mt a) 1)))))

          ;; merge lists of intx and box points
          (setf kl (+ kl (+ (- l j) 1))
                kb (+ kb (+ (- k i) 1))
                kr (+ kr (+ (- l j) 1))
                kt (+ kt (+ (- k i) 1)))

          (let ((l-i (make-index-vector (1+ kl))) (l-p (make-index-vector (1+ kl)))
                (l-l (make-index-vector (1+ kl))) (l-s (make-index-vector (1+ kl)))
                (l-u (make-index-vector (1+ kl)))
                (b-i (make-index-vector (1+ kb))) (b-p (make-index-vector (1+ kb)))
                (b-l (make-index-vector (1+ kb))) (b-s (make-index-vector (1+ kb)))
                (b-u (make-index-vector (1+ kb)))
                (r-i (make-index-vector (1+ kr))) (r-p (make-index-vector (1+ kr)))
                (r-l (make-index-vector (1+ kr))) (r-s (make-index-vector (1+ kr)))
                (r-u (make-index-vector (1+ kr)))
                (t-i (make-index-vector (1+ kt))) (t-p (make-index-vector (1+ kt)))
                (t-l (make-index-vector (1+ kt))) (t-s (make-index-vector (1+ kt)))
                (t-u (make-index-vector (1+ kt)))
                (ip (make-index-vector (1+ m)))
                (pix -1)
                (ext 0d0))
            (declare (type index-vector l-i l-p l-l l-s l-u b-i b-p b-l b-s b-u
                           r-i r-p r-l r-s r-u t-i t-p t-l t-s t-u ip)
                     (type fixnum pix) (type double-float ext))
            (macrolet ((merge-side (&key kk ii pp lll ss uu mm cv xv yv lv itv cmp extreme
                                         box-coord sp-init sp-cont sp-step tm-at buffer-test)
                         ;; one side's list: the box points and the intersection
                         ;; points sorted along the side's direction
                         `(progn
                            (do ((a 0 (1+ a))) ((> a m))
                              (declare (type fixnum a))
                              (setf (aref ip a) 0))
                            (setf pix -1
                                  ext ,extreme)
                            (do ((a 0 (1+ a))) ((> a m))
                              (declare (type fixnum a))
                              (when (>= (aref ,mm a) 0)
                                (when (,cmp (aref (the real-vector (svref ,cv a)) (aref ip a)) ext)
                                  (setf ext (aref (the real-vector (svref ,cv a)) (aref ip a))
                                        pix a))))
                            (if (= pix -1)
                                ;; no intersection: output box points
                                (progn
                                  (setf ,kk -1)
                                  (do ((sp ,sp-init (,sp-step sp))) ((not ,sp-cont))
                                    (declare (type fixnum sp))
                                    (when (>= ,tm-at 0)
                                      (incf ,kk)
                                      (setf (aref ,ii ,kk) ,tm-at
                                            (aref ,pp ,kk) -1
                                            (aref ,lll ,kk) -1
                                            (aref ,ss ,kk) +nostatus+
                                            (aref ,uu ,kk) +notused+))))
                                ;; merge intersection points with box points
                                (let ((sp ,sp-init)
                                      (its t)
                                      (moved nil))
                                  (declare (type fixnum sp))
                                  (setf ,kk -1)
                                  (do () ((not ,sp-cont))
                                    (if (< ,tm-at 0)
                                        (setf sp (,sp-step sp))
                                        (progn
                                          (incf ,kk)
                                          (when (and its (< pix 0))
                                            (nurbs-error :geo-err 'tessellate--st-tess-b-box-x-trim))
                                          (cond ((and its
                                                      (< (abs (- (aref ,box-coord ,tm-at)
                                                                 (aref (the real-vector (svref ,cv pix)) (aref ip pix))))
                                                         tol))
                                                 ;; use intx point
                                                 (setf (aref ,ii ,kk) ,tm-at
                                                       (aref ,pp ,kk) pix
                                                       (aref ,lll ,kk) (aref (the index-vector (svref ,lv pix)) (aref ip pix))
                                                       (aref ,ss ,kk) (aref (the index-vector (svref ,itv pix)) (aref ip pix))
                                                       (aref ,uu ,kk) +notused+)
                                                 (setf sp (,sp-step sp))
                                                 (incf (aref ip pix))
                                                 (setf moved t))
                                                ((or (not its)
                                                     (,cmp (aref ,box-coord ,tm-at)
                                                           (aref (the real-vector (svref ,cv pix)) (aref ip pix))))
                                                 ;; box point
                                                 (setf (aref ,ii ,kk) ,tm-at
                                                       (aref ,pp ,kk) -1
                                                       (aref ,lll ,kk) -1
                                                       (aref ,ss ,kk) +nostatus+
                                                       (aref ,uu ,kk) +notused+)
                                                 (setf sp (,sp-step sp))
                                                 (setf moved nil))
                                                (t
                                                 (let ((fnd +no+)
                                                       (px (aref (the real-vector (svref ,xv pix)) (aref ip pix)))
                                                       (py (aref (the real-vector (svref ,yv pix)) (aref ip pix))))
                                                   (declare (type fixnum fnd) (type double-float px py))
                                                   (when ,buffer-test
                                                     ;; look for intx point in buffer
                                                     (multiple-value-bind (nbuf nbi nbs nfnd nk)
                                                         (tessellate--st-tess-buffer ul vl nl px py tol buf bi bs bn)
                                                       (setf buf nbuf bi nbi bs nbs fnd nfnd)
                                                       (when (= nfnd +yes+)
                                                         (setf pin nk))))
                                                   (when (= fnd +no+)
                                                     ;; append this intx point to list
                                                     (append-pt px py)
                                                     (setf pin nl)))
                                                 (setf (aref ,ii ,kk) pin
                                                       (aref ,pp ,kk) pix
                                                       (aref ,lll ,kk) (aref (the index-vector (svref ,lv pix)) (aref ip pix))
                                                       (aref ,ss ,kk) (aref (the index-vector (svref ,itv pix)) (aref ip pix))
                                                       (aref ,uu ,kk) +notused+)
                                                 (incf (aref ip pix))
                                                 (setf moved t)))
                                          (when moved
                                            (setf pix -1
                                                  ext ,extreme
                                                  its nil)
                                            (do ((a 0 (1+ a))) ((> a m))
                                              (declare (type fixnum a))
                                              (when (and (>= (aref ,mm a) 0) (<= (aref ip a) (aref ,mm a)))
                                                (when (,cmp (aref (the real-vector (svref ,cv a)) (aref ip a)) ext)
                                                  (setf ext (aref (the real-vector (svref ,cv a)) (aref ip a))
                                                        pix a))
                                                (setf its t))))))))))))
              ;; NL_LEFT: from the highest intersection down
              (merge-side :kk kl :ii l-i :pp l-p :lll l-l :ss l-s :uu l-u :mm ml
                          :cv yl :xv xl :yv yl :lv ll :itv itl :cmp > :extreme (- bb tol)
                          :box-coord vl :sp-init l :sp-cont (>= sp j) :sp-step 1-
                          :tm-at (aref tm i sp) :buffer-test (/= i 0))
              ;; NL_BOTTOM: from the leftmost intersection right
              (merge-side :kk kb :ii b-i :pp b-p :lll b-l :ss b-s :uu b-u :mm mb
                          :cv xb :xv xb :yv yb :lv lb :itv itb :cmp < :extreme (+ br tol)
                          :box-coord ul :sp-init i :sp-cont (<= sp k) :sp-step 1+
                          :tm-at (aref tm sp j) :buffer-test (/= j 0))
              ;; NL_RIGHT: from the lowest intersection up
              (merge-side :kk kr :ii r-i :pp r-p :lll r-l :ss r-s :uu r-u :mm mr
                          :cv yr :xv xr :yv yr :lv lr :itv itr :cmp < :extreme (+ bt tol)
                          :box-coord vl :sp-init j :sp-cont (<= sp l) :sp-step 1+
                          :tm-at (aref tm k sp) :buffer-test (/= k rtm))
              ;; NL_TOP: from the rightmost intersection left
              (merge-side :kk kt :ii t-i :pp t-p :lll t-l :ss t-s :uu t-u :mm mt
                          :cv xt :xv xt :yv yt :lv lt :itv itt :cmp > :extreme (- bl tol)
                          :box-coord ul :sp-init k :sp-cont (>= sp i) :sp-step 1-
                          :tm-at (aref tm sp l) :buffer-test (/= l stm)))

            ;; All lists are created.  Now trace lists and form polygons.
            (let ((sh -1) (ipu 0) (oh 0) (c 0))
              (declare (type fixnum sh ipu oh c))
              ;; find number of output polygons
              (flet ((count-in (ss kk)
                       (declare (type index-vector ss) (type fixnum kk))
                       (do ((b 0 (1+ b))) ((> b kk))
                         (declare (type fixnum b))
                         (when (or (= (aref ss b) +inp+) (= (aref ss b) +turnp+))
                           (incf sh)))))
                (count-in l-s kl)
                (count-in b-s kb)
                (count-in r-s kr)
                (count-in t-s kt))

              ;; avoid touch cases
              (when (= sh -1)
                ;; check if polygons are inside the box
                (setf pfl (make-index-vector (1+ m)))
                (setf s -1)
                (do ((a 0 (1+ a))) ((> a m))
                  (declare (type fixnum a))
                  (if (and (>= (aref pul a) blt) (<= (aref pur a) brt)
                           (>= (aref pvb a) bbt) (<= (aref pvt a) btt))
                      (progn (incf s)
                             (setf (aref pfl a) +yes+))
                      (setf (aref pfl a) +no+)))
                (if (and (> s -1) (< s m))
                    (setf c (+ s 1))
                    (setf c s))
                (when (>= s 0)
                  (let ((po (make-array (1+ c) :initial-element nil))
                        (go (make-index-vector (1+ c)))
                        (ino (make-index-vector (1+ c))))
                    (declare (type simple-vector po) (type index-vector go ino))
                    (setf c -1)
                    (when (/= (aref pfl 0) +yes+)
                      (incf c)
                      (multiple-value-bind (q r) (box-polygon)
                        (setf (svref po c) q
                              (aref go c) r
                              (aref ino c) -1)))
                    (do ((a 0 (1+ a))) ((> a m))
                      (declare (type fixnum a))
                      (when (= (aref pfl a) +yes+)
                        (incf c)
                        (setf (svref po c) (trim-polygon a))
                        (setf (aref go c) (- (aref h a) 1))
                        (if (= a 0)
                            (setf (aref ino c) -1)
                            (setf (aref ino c) 0))))
                    (done po go ino c)))

                ;; check if box is inside the polygon
                (let ((mx (* 0.5d0 (+ bl br)))
                      (my (* 0.5d0 (+ bb bt)))
                      (pinpol +inp+))
                  (declare (type double-float mx my) (type fixnum pinpol))
                  (do ((a 0 (1+ a))) ((> a m))
                    (declare (type fixnum a))
                    (let ((pst (tessellate--st-is-pt-in-polygon
                                mx my (svref x a) (svref y a) (aref h a)
                                (aref pul a) (aref pur a) (aref pvb a) (aref pvt a) tol)))
                      (declare (type fixnum pst))
                      (when (and (= a 0) (= pst +outp+))
                        (setf pinpol +outp+)
                        (return))
                      (when (and (> a 0) (= pst +inp+))
                        (setf pinpol +outp+)
                        (return))))
                  (when (= pinpol +outp+)
                    (done nil nil nil -1)))

                ;; switch to NL_INP case
                (let* ((r (- (+ kl kb kr kt) 1))
                       (po (make-array 1 :initial-element nil))
                       (q (make-index-vector (1+ r)))
                       (go (make-index-vector 1))
                       (ino (make-index-vector 1)))
                  (declare (type fixnum r) (type simple-vector po) (type index-vector q go ino))
                  (setf (svref po 0) q)
                  (setf r -1)
                  (do ((b 0 (1+ b))) ((>= b kl))
                    (declare (type fixnum b))
                    (setf (aref q (incf r)) (aref l-i b)))
                  (do ((a 0 (1+ a))) ((>= a kb))
                    (declare (type fixnum a))
                    (setf (aref q (incf r)) (aref b-i a)))
                  (do ((b 0 (1+ b))) ((>= b kr))
                    (declare (type fixnum b))
                    (setf (aref q (incf r)) (aref r-i b)))
                  (do ((a 0 (1+ a))) ((>= a kt))
                    (declare (type fixnum a))
                    (setf (aref q (incf r)) (aref t-i a)))
                  (setf (aref go 0) r
                        (aref ino 0) -1)
                  (done po go ino 0)))

              ;; adjust number of NL_INP or NL_TURNP
              (setf ipu (+ sh 1))
              (flet ((in-or-turn (f) (or (= f +inp+) (= f +turnp+))))
                (declare (inline in-or-turn))
                ;; (a = INP or TURNP) and (b = INP or TURNP): the C's four cases
                (when (and (in-or-turn (aref t-s kt)) (in-or-turn (aref l-s 0)))
                  (decf ipu))
                (when (and (in-or-turn (aref l-s kl)) (in-or-turn (aref b-s 0)))
                  (decf ipu))
                (when (and (in-or-turn (aref b-s kb)) (in-or-turn (aref r-s 0)))
                  (decf ipu))
                (when (and (in-or-turn (aref r-s kr)) (in-or-turn (aref t-s 0)))
                  (decf ipu)))

              ;; consider the NL_ONANDOVER case
              (setf oh 0)
              (when (= sta +onandover+)
                (do ((a 0 (1+ a))) ((> a m))
                  (declare (type fixnum a))
                  (if (and (>= (aref pul a) blt) (<= (aref pur a) brt)
                           (>= (aref pvb a) bbt) (<= (aref pvt a) btt))
                      (progn (incf oh)
                             (setf (aref pfl a) +yes+))
                      (setf (aref pfl a) +no+))))

              ;; allocate memory for output polygons
              (let ((pl (make-array (+ sh oh 2) :initial-element nil))
                    (gl (make-index-vector (+ sh oh 2)))
                    (inl (make-index-vector (+ sh oh 2)))
                    (rh (make-index-vector (+ sh oh 2)))
                    (r -1))
                (declare (type simple-vector pl) (type index-vector gl inl rh) (type fixnum r))
                (do ((a 0 (1+ a))) ((> a sh))
                  (declare (type fixnum a))
                  (setf (svref pl a) (make-index-vector (1+ ri)))
                  (setf (aref rh a) ri))

                ;; while not all input points are used, create output polygons
                (macrolet ((append-vertex (idx)
                             ;; ST_TessAppendPolygonVertex(&pl[s], &rh[s], ri, idx, &r, &SL)
                             `(multiple-value-bind (nq nqh nr)
                                  (tessellate--st-tess-append-polygon-vertex
                                   (svref pl s) (aref rh s) ri ,idx r)
                                (setf (svref pl s) nq
                                      (aref rh s) nqh
                                      r nr))))
                  (setf s -1)          ; high index of output polygons
                  (do () ((<= ipu 0))
                    (let ((li 0) (side 0) (sn 0)
                          (pst +nostatus+) (cross t) (cri 0) (sti 0)
                          (ploop t) (sturn nil) (a 0) (crx 0d0) (cry 0d0) (cnt -1))
                      (declare (type fixnum li side sn pst cri sti a cnt)
                               (type double-float crx cry))
                      ;; find an unused NL_INP or NL_TURNP
                      (multiple-value-setq (pix li side sn pin)
                        (handler-case
                            (tessellate--st-tess-find-unused
                             l-i l-s l-u l-p l-l kl b-i b-s b-u b-p b-l kb
                             r-i r-s r-u r-p r-l kr t-i t-s t-u t-p t-l kt)
                          (nurbs-error ()
                            (nurbs-error :geo-err 'tessellate--st-tess-b-box-x-trim))))
                      (cond ((= side +left+)
                             (setf (aref l-u sn) +used+)
                             (when (= sn kl) (setf (aref b-u 0) +used+)))
                            ((= side +bottom+)
                             (setf (aref b-u sn) +used+)
                             (when (= sn kb) (setf (aref r-u 0) +used+)))
                            ((= side +right+)
                             (setf (aref r-u sn) +used+)
                             (when (= sn kr) (setf (aref t-u 0) +used+)))
                            ((= side +top+)
                             (setf (aref t-u sn) +used+)
                             (when (= sn kt) (setf (aref l-u 0) +used+))))
                      (decf ipu)

                      ;; get start point
                      (setf pst +nostatus+
                            cross t)
                      (incf s)
                      (setf r -1     ; high index in pl[s] (r + 1 distinct pts in output polygon)
                            cri pin  ; current index
                            sti pin  ; start index (IN- or TURN- point)
                            ploop t
                            sturn nil)
                      (append-vertex pin)
                      (setf a (rem (+ li 1) (aref h pix)))
                      ;; coordinates of current trim polygon point: may be
                      ;; corner or intx point once in inner loop
                      (setf crx (aref (the real-vector (svref x pix)) a)
                            cry (aref (the real-vector (svref y pix)) a)
                            cnt -1)

                      (do () ((not (or (/= cri sti) ploop)))
                        (incf cnt)
                        (when (> cnt (* 4 (aref h pix)))
                          (nurbs-error :geo-err 'tessellate--st-tess-b-box-x-trim))
                        (cond
                          ;; Case 1: point is inside the box -> polygon point
                          ((tessellate--st-tess-is-pt-in-b-box crx cry bl br bb bt tol)
                           (append-pt crx cry)
                           (append-vertex nl)
                           (setf li (rem (+ li 1) (aref h pix)))
                           (setf a (rem (+ li 1) (aref h pix)))
                           (setf cri nl)
                           (setf crx (aref (the real-vector (svref x pix)) a)
                                 cry (aref (the real-vector (svref y pix)) a)))

                          ;; point is not in box -> cross side
                          (cross
                           (multiple-value-setq (side sn pin)
                             (handler-case
                                 (tessellate--st-tess-find-crossing-pt
                                  l-i l-p l-l l-s kl b-i b-p b-l b-s kb
                                  r-i r-p r-l r-s kr t-i t-p t-l t-s kt
                                  pix h li cri)
                               (nurbs-error ()
                                 (nurbs-error :geo-err 'tessellate--st-tess-b-box-x-trim))))
                           (cond ((= side +left+) (setf pst (aref l-s sn)))
                                 ((= side +bottom+) (setf pst (aref b-s sn)))
                                 ((= side +right+) (setf pst (aref r-s sn)))
                                 ((= side +top+) (setf pst (aref t-s sn))))
                           (setf crx (aref ul pin)
                                 cry (aref vl pin))
                           (when (and (= pst +turnp+) (= pix 0))
                             (cond ((= side +left+) (setf (aref l-u sn) +used+))
                                   ((= side +bottom+) (setf (aref b-u sn) +used+))
                                   ((= side +right+) (setf (aref r-u sn) +used+))
                                   ((= side +top+) (setf (aref t-u sn) +used+)))
                             (decf ipu))
                           (setf cross nil
                                 cri pin)
                           (append-vertex pin)
                           (when ploop
                             (when (or (/= cri sti) (= pst +turnp+))
                               (setf ploop nil))))

                          ;; point is on the side
                          ((or (not (tessellate--st-tess-is-pt-b-box-corner side sn kl kb kr kt))
                               (and (= pst +turnp+) (= pix 0))
                               (and (= pst +turnp+) sturn)
                               (= pst +inp+))
                           ;; point is NL_TURNP point or NL_INP -> go inside or go along the side
                           (let ((goin nil))
                             (cond ((= pst +inp+) (setf goin t))
                                   ((and (= pst +turnp+) (= pix 0)) (setf goin t))
                                   ((and (= pst +turnp+) sturn) (setf goin t)))
                             (if goin
                                 (progn
                                   (cond ((= side +left+)
                                          (setf pix (aref l-p sn)
                                                li (aref l-l sn))
                                          (when (/= cri sti) (setf (aref l-u sn) +used+)))
                                         ((= side +bottom+)
                                          (setf pix (aref b-p sn)
                                                li (aref b-l sn))
                                          (when (/= cri sti) (setf (aref b-u sn) +used+)))
                                         ((= side +right+)
                                          (setf pix (aref r-p sn)
                                                li (aref r-l sn))
                                          (when (/= cri sti) (setf (aref r-u sn) +used+)))
                                         ((= side +top+)
                                          (setf pix (aref t-p sn)
                                                li (aref t-l sn))
                                          (when (/= cri sti) (setf (aref t-u sn) +used+))))
                                   (when (/= cri sti)
                                     (when (or (/= pst +turnp+) (/= pix 0))
                                       (decf ipu)))
                                   (setf a (rem (+ li 1) (aref h pix)))
                                   (setf crx (aref (the real-vector (svref x pix)) a)
                                         cry (aref (the real-vector (svref y pix)) a))
                                   (setf cross t
                                         pst +nostatus+
                                         sturn nil))
                                 ;; point is side point -> march along side
                                 (progn
                                   (cond ((= side +left+)
                                          (setf cri (aref l-i (+ sn 1))
                                                pst (aref l-s (+ sn 1))))
                                         ((= side +bottom+)
                                          (setf cri (aref b-i (+ sn 1))
                                                pst (aref b-s (+ sn 1))))
                                         ((= side +right+)
                                          (setf cri (aref r-i (+ sn 1))
                                                pst (aref r-s (+ sn 1))))
                                         ((= side +top+)
                                          (setf cri (aref t-i (+ sn 1))
                                                pst (aref t-s (+ sn 1)))))
                                   (setf crx (aref ul cri)
                                         cry (aref vl cri))
                                   (incf sn)
                                   (when (= pst +turnp+)
                                     (setf sturn t))
                                   (when (/= cri sti)
                                     (append-vertex cri)
                                     (setf ploop nil))))))

                          ;; point is corner point -> change direction
                          ((tessellate--st-tess-is-pt-b-box-corner side sn kl kb kr kt)
                           (cond ((= side +left+)
                                  (if (= sn 0)
                                      (setf side +left+ cri (aref l-i 1) pst (aref l-s 1))
                                      (setf side +bottom+ cri (aref b-i 1) pst (aref b-s 1))))
                                 ((= side +bottom+)
                                  (if (= sn 0)
                                      (setf side +bottom+ cri (aref b-i 1) pst (aref b-s 1))
                                      (setf side +right+ cri (aref r-i 1) pst (aref r-s 1))))
                                 ((= side +right+)
                                  (if (= sn 0)
                                      (setf side +right+ cri (aref r-i 1) pst (aref r-s 1))
                                      (setf side +top+ cri (aref t-i 1) pst (aref t-s 1))))
                                 ((= side +top+)
                                  (if (= sn 0)
                                      (setf side +top+ cri (aref t-i 1) pst (aref t-s 1))
                                      (setf side +left+ cri (aref l-i 1) pst (aref l-s 1)))))
                           (setf sn 1)
                           (setf crx (aref ul cri)
                                 cry (aref vl cri))
                           (when (= pst +turnp+)
                             (setf sturn t))
                           (when (/= cri sti)
                             (append-vertex cri)
                             (setf ploop nil))))) ; end of inner while

                      (setf (aref gl s) r
                            (aref inl s) -1)))) ; end of outer while

                ;; add polygons that are inside the box
                (when (= sta +onandover+)
                  ;; append trimming polygon vertices to the polygon list
                  (setf c s)
                  (do ((a 0 (1+ a))) ((> a m))
                    (declare (type fixnum a))
                    (when (= (aref pfl a) +yes+)
                      (incf s)
                      (setf (svref pl s) (trim-polygon a))
                      (setf (aref gl s) (- (aref h a) 1))
                      (setf pin (aref (the index-vector (svref pl s)) 0))
                      (let ((mx (aref ul pin))
                            (my (aref vl pin)))
                        (declare (type double-float mx my))
                        (setf (aref inl s) -1)
                        (do ((b 0 (1+ b))) ((> b c))
                          (declare (type fixnum b))
                          (when (tessellate--st-is-pt-in-poly ul vl (svref pl b) (aref gl b) mx my tol)
                            (setf (aref inl s) b)))))))

                ;; eliminate degenerate polygons and create output
                (setf r -1)
                (do ((a 0 (1+ a))) ((> a s))
                  (declare (type fixnum a))
                  (when (> (aref gl a) 1)
                    (incf r)))
                (when (= r -1)
                  (done nil nil nil -1))
                (let ((po (make-array (1+ r) :initial-element nil))
                      (go (make-index-vector (1+ r)))
                      (ino (make-index-vector (1+ r))))
                  (declare (type simple-vector po) (type index-vector go ino))
                  (setf r -1)
                  (do ((a 0 (1+ a))) ((> a s))
                    (declare (type fixnum a))
                    (when (> (aref gl a) 1)
                      (incf r)
                      (setf (aref go r) (aref gl a)
                            (aref ino r) (aref inl a))
                      (let ((q (make-index-vector (1+ (aref go r))))
                            (pa (svref pl a)))
                        (declare (type index-vector q pa))
                        (setf (svref po r) q)
                        (do ((b 0 (1+ b))) ((> b (aref go r)))
                          (declare (type fixnum b))
                          (setf (aref q b) (aref pa b))))))
                  (done po go ino r))))))))))

(defun tessellate--st-tess-append-pt (u v nh ni x y n)
  "Append the point (X, Y) to the parameter lists U, V (real-vectors) whose
highest used index is N and highest allocated index NH: when they are full,
both are first grown by NI, into fresh arrays holding the old elements
0..NH.  Returns (values u v nh n): the lists (fresh when grown), their
highest allocated index and N+1, the new point's index.
NLib: ST_TessAppendPt (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum nh ni n) (type double-float x y))
  (let ((mh nh)
        (m n)
        (ul u)
        (vl v))
    (declare (type fixnum mh m) (type real-vector ul vl))
    (setf m (+ m 1))
    (when (> m mh)
      ;; N_Realloc1dRealArray(&ul, mh, mh + ni) and the same for vl
      (let* ((kl (nl-min mh (+ mh ni)))
             (nu (make-real-vector (1+ (+ mh ni))))
             (nv (make-real-vector (1+ (+ mh ni)))))
        (declare (type fixnum kl) (type real-vector nu nv))
        (replace nu ul :end2 (min (1+ kl) (length ul)))
        (replace nv vl :end2 (min (1+ kl) (length vl)))
        (setf ul nu
              vl nv))
      (setf mh (+ mh ni)))
    (setf (aref ul m) x
          (aref vl m) y)
    (values ul vl mh m)))

(defun tessellate--st-tess-buffer (u v n x y tol buf bi bs bn)
  "Look for the point (X, Y) among the intersection points in buffer BUF (an
index-vector of indexes into U, V, highest used index BI, highest allocated
index BS): a buffered point within TOL of it (squared distance below TOL^2)
is removed from the buffer and reported found; otherwise N+1, the index the
point gets when appended to U, V, is buffered (the buffer grown by BN, into
a fresh array, when full).  Returns (values buf bi bs fnd k): the buffer
(fresh when grown), its two highest indexes, FND +yes+ or +no+, and K, the
index of the point found -- -1 when not found (the C leaves *k unset).
NLib: ST_TessBuffer (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type fixnum n bi bs bn) (type double-float x y tol)
           (type index-vector buf))
  (let ((bfl buf)
        (bli bi)
        (bls bs)
        (t2 (* tol tol))
        (f -1)
        (k -1))
    (declare (type index-vector bfl) (type fixnum bli bls f k) (type double-float t2))
    (fast
      (do ((i 0 (1+ i))) ((> i bli))
        (declare (type fixnum i))
        (let* ((j (aref bfl i))
               (d2 (+ (* (- x (aref u j)) (- x (aref u j)))
                      (* (- y (aref v j)) (- y (aref v j))))))
          (declare (type fixnum j) (type double-float d2))
          (when (< d2 t2)
            (setf k j
                  f i)
            (return)))))
    (when (>= f 0)
      ;; delete point from buffer, then exit
      (fast
        (do ((i (+ f 1) (1+ i))) ((> i bli))
          (declare (type fixnum i))
          (setf (aref bfl (- i 1)) (aref bfl i))))
      (setf bli (- bli 1))
      (return-from tessellate--st-tess-buffer (values bfl bli bls +yes+ k)))
    ;; add point to buffer
    (setf bli (+ bli 1))
    (when (> bli bls)
      ;; N_Realloc1dIntArray(&bfl, bls, bls + bn)
      (let* ((kl (nl-min bls (+ bls bn)))
             (nb (make-index-vector (1+ (+ bls bn)))))
        (declare (type fixnum kl) (type index-vector nb))
        (replace nb bfl :end2 (min (1+ kl) (length bfl)))
        (setf bfl nb))
      (setf bls (+ bls bn)))
    (setf (aref bfl bli) (+ n 1))
    (values bfl bli bls +no+ k)))

(defun tessellate--st-tess-append-polygon-vertex (p rh ri i r)
  "Append the point index I to the polygon P (an index-vector) whose highest
used index is R and highest allocated index RH: when it is full, it is first
grown by RI, into a fresh array holding the old elements 0..RH.  Returns
(values p rh r): the polygon (fresh when grown), its highest allocated index
and R+1.
NLib: ST_TessAppendPolygonVertex (file-private, NL_Tessellate.cpp)."
  (declare (type index-vector p) (type fixnum rh ri i r))
  (let ((sh rh)
        (s r)
        (q p))
    (declare (type fixnum sh s) (type index-vector q))
    (setf s (+ s 1))
    (when (> s sh)
      ;; N_Realloc1dIntArray(&q, sh, sh + ri)
      (let* ((kl (nl-min sh (+ sh ri)))
             (nq (make-index-vector (1+ (+ sh ri)))))
        (declare (type fixnum kl) (type index-vector nq))
        (replace nq q :end2 (min (1+ kl) (length q)))
        (setf q nq))
      (setf sh (+ sh ri)))
    (setf (aref q s) i)
    (values q sh s)))

(defun tessellate--st-tess-is-pt-in-b-box (x y xl xr yb yt tol)
  "True when the point (X, Y) lies inside the box [XL, XR] x [YB, YT] by at
least TOL.
NLib: ST_TessIsPtInBBox (file-private, NL_Tessellate.cpp)."
  (declare (type double-float x y xl xr yb yt tol))
  (cond ((or (< x (+ xl tol)) (> x (- xr tol))) nil)
        ((or (< y (+ yb tol)) (> y (- yt tol))) nil)
        (t t)))

(defun tessellate--st-tess-is-pt-b-box-corner (side i kl kb kr kt)
  "True when entry I of the box side list SIDE (+left+, +bottom+, +right+ or
+top+; KL, KB, KR, KT the lists' highest indexes) is a corner of the box: its
first or last entry.
NLib: ST_TessIsPtBBoxCorner (file-private, NL_Tessellate.cpp)."
  (declare (type fixnum side i kl kb kr kt))
  (cond ((= side +left+) (or (= i 0) (= i kl)))
        ((= side +bottom+) (or (= i 0) (= i kb)))
        ((= side +right+) (or (= i 0) (= i kr)))
        ((= side +top+) (or (= i 0) (= i kt)))
        (t nil)))

(defun tessellate--st-is-pt-in-polygon (xp yp x y n bxl bxr byb byt tol)
  "Whether the point (XP, YP) lies inside the polygon X[0..N], Y[0..N]
(real-vectors) with bounding box BXL BXR BYB BYT: +outp+ when it is outside
the box widened by TOL, else by the parity of the legs crossing the ray to
its right, +inp+ or +outp+.
NLib: ST_IsPtInPolygon (file-private, NL_Tessellate.cpp)."
  (declare (type double-float xp yp bxl bxr byb byt tol) (type real-vector x y)
           (type fixnum n))
  (when (or (< xp (- bxl tol)) (> xp (+ bxr tol)) (< yp (- byb tol)) (> yp (+ byt tol)))
    (return-from tessellate--st-is-pt-in-polygon +outp+))
  (let ((k 0)
        (ytl (+ yp tol)))
    (declare (type fixnum k) (type double-float ytl))
    (fast
      (do ((i 1 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (when (or (and (> (aref y i) yp) (< (aref y (- i 1)) ytl))
                  (and (> (aref y (- i 1)) yp) (< (aref y i) ytl)))
          (let ((xi (+ (aref x i)
                       (/ (* (- yp (aref y i)) (- (aref x i) (aref x (- i 1))))
                          (- (aref y i) (aref y (- i 1)))))))
            (declare (type double-float xi))
            (when (> xi xp)
              (incf k))))))
    (if (/= (rem k 2) 0)
        +inp+
        +outp+)))

(defun tessellate--st-is-pt-in-poly (u v p n x y tol)
  "True when the point (X, Y) lies inside the closed polygon whose vertices
are the points U[P[0..N]], V[P[0..N]] (P an index-vector, the leg from P[N]
back to P[0] included), by the parity of the legs crossing the ray to its
right.
NLib: ST_IsPtInPoly (file-private, NL_Tessellate.cpp)."
  (declare (type real-vector u v) (type index-vector p) (type fixnum n)
           (type double-float x y tol))
  (let ((k 0)
        (ytl (+ y tol)))
    (declare (type fixnum k) (type double-float ytl))
    (fast
      (do ((i 1 (1+ i))) ((> i (+ n 1)))
        (declare (type fixnum i))
        (let* ((j (if (= i (+ n 1))
                      (aref p 0)
                      (aref p i)))
               (j1 (aref p (- i 1))))
          (declare (type fixnum j j1))
          (when (or (and (> (aref v j) y) (< (aref v j1) ytl))
                    (and (> (aref v j1) y) (< (aref v j) ytl)))
            (let ((xi (+ (aref u j)
                         (/ (* (- y (aref v j)) (- (aref u j) (aref u j1)))
                            (- (aref v j) (aref v j1))))))
              (declare (type double-float xi))
              (when (> xi x)
                (incf k)))))))
    (/= (rem k 2) 0)))
