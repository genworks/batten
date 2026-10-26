;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(declaim (inline is-srf-rat))
(defun is-srf-rat (sur)
  "True when surface SUR is rational.  NLib: N_IsSrfRat."
  (declare (type surface sur))
  (/= (aref (surface-cpoints sur) 0 0 3) +no-w+))

(defun srf-get-parameter-bounds (sur)
  "The parameter domain of SUR, (values ul ur vb vt).  NLib: N_SrfGetParameterBounds."
  (declare (type surface sur))
  (let ((ku (surface-knots-u sur)) (kv (surface-knots-v sur)))
    (values (aref ku 0) (aref ku (1- (length ku)))
            (aref kv 0) (aref kv (1- (length kv))))))

(defun srf-weights-are-valid (sur rname)
  "Signal :wei-err in routine RNAME when a weight of rational surface SUR lies
outside [*wmin*, *wmax*].  Returns SUR otherwise.  NLib: N_SrfWeightsAreValid."
  (declare (type surface sur))
  (when (is-srf-rat sur)
    (let ((pw (surface-cpoints sur))
          (wmin *wmin*) (wmax *wmax*))
      (declare (type cpoint-net pw) (type double-float wmin wmax))
      (dotimes (i (1+ (surface-n sur)))
        (dotimes (j (1+ (surface-m sur)))
          (let ((w (aref pw i j 3)))
            (declare (type double-float w))
            (when (or (< w wmin) (> w wmax))
              (nurbs-error :wei-err rname)))))))
  sur)

(defun srf-from-c-net-and-knot-vectors (sur net p q knu knv)
  "Make SUR the surface of control net NET (a cpoint-net), degrees P and Q
and knot vectors KNU and KNV, sharing them as the C shares its pointers.
Returns SUR.  NLib: N_SrfFromCNetAndKnotVectors."
  (declare (type surface sur) (type cpoint-net net) (type fixnum p q)
           (type real-vector knu knv))
  (setf (surface-cpoints sur) net
        (surface-degree-u sur) p
        (surface-degree-v sur) q
        (surface-knots-u sur) knu
        (surface-knots-v sur) knv)
  sur)

(defun srf-get-net-and-knot-vectors (sur)
  "The control net, degrees and knot vectors of SUR, the surface's own
objects: (values net p q knu knv).  NLib: N_SrfGetNetAndKnotVectors."
  (declare (type surface sur))
  (values (surface-cpoints sur) (surface-degree-u sur) (surface-degree-v sur)
          (surface-knots-u sur) (surface-knots-v sur)))

(defun srf-get-c-pts-knot-vector-and-knots (sur)
  "The control points, knot vectors and knot arrays of SUR, the surface's own
objects: (values pw knu knv u v).  A Batten knot vector is its knot array, so
KNU is U and KNV is V.  NLib: N_SrfGetCPtsKnotVectorAndKnots."
  (declare (type surface sur))
  (values (surface-cpoints sur) (surface-knots-u sur) (surface-knots-v sur)
          (surface-knots-u sur) (surface-knots-v sur)))

(defun srf-set-c-pts-and-knots (sur pw u v)
  "Give SUR the control points PW and the knot arrays U and V, sharing them as
the C shares its pointers.  Returns SUR.  NLib: N_SrfSetCPtsAndKnots."
  (declare (type surface sur) (type cpoint-net pw) (type real-vector u v))
  (setf (surface-cpoints sur) pw
        (surface-knots-u sur) u
        (surface-knots-v sur) v)
  sur)

(defun srf-is-closed (sur dir)
  "True when SUR is closed in direction DIR (+udir+ or +vdir+): the Euclidean
control points of its first and last rows (+udir+) or columns (+vdir+) lie
within *mtol* of each other.  False for any other DIR.  NLib: N_SrfIsClosed."
  (declare (type surface sur) (type fixnum dir))
  (let* ((pw (surface-cpoints sur))
         (n (surface-n sur)) (m (surface-m sur))
         (mtol *mtol*)
         (cw (make-cpoint))
         (ps (make-point)) (pe (make-point)))
    (declare (type cpoint-net pw) (type fixnum n m) (type double-float mtol)
             (dynamic-extent cw ps pe))
    (flet ((euclid (i j p)
             (setf (aref cw 0) (aref pw i j 0) (aref cw 1) (aref pw i j 1)
                   (aref cw 2) (aref pw i j 2) (aref cw 3) (aref pw i j 3))
             (c-pt-to-pt-euclid cw p)))
      (declare (inline euclid))
      (cond
        ((= dir +udir+)
         (dotimes (j (1+ m) t)
           (euclid 0 j ps)
           (euclid n j pe)
           (when (> (dist-pt-pt ps pe) mtol)
             (return nil))))
        ((= dir +vdir+)
         (dotimes (i (1+ n) t)
           (euclid i 0 ps)
           (euclid i m pe)
           (when (> (dist-pt-pt ps pe) mtol)
             (return nil))))
        (t nil)))))

(defun 1d-pt-set-is-closed (a m ptp tol)
  "True when the first and last (index M) points of A lie within TOL of each
other: A is a point-array for PTP +epoint+ (Euclidean distance) or a
cpoint-array for +hpoint+ (homogeneous distance).  False for any other PTP.
NLib: N_1dPtSetIsClosed."
  (declare (type fixnum m ptp) (type double-float tol))
  (cond
    ((= ptp +epoint+)
     (let ((a a) (p0 (make-point)) (pm (make-point)))
       (declare (type point-array a) (dynamic-extent p0 pm))
       (set-point p0 (aref a 0 0) (aref a 0 1) (aref a 0 2))
       (set-point pm (aref a m 0) (aref a m 1) (aref a m 2))
       (not (> (dist-pt-pt p0 pm) tol))))
    ((= ptp +hpoint+)
     (let ((a a) (p0 (make-cpoint)) (pm (make-cpoint)))
       (declare (type cpoint-array a) (dynamic-extent p0 pm))
       (set-cpoint p0 (aref a 0 0) (aref a 0 1) (aref a 0 2) (aref a 0 3))
       (set-cpoint pm (aref a m 0) (aref a m 1) (aref a m 2) (aref a m 3))
       (not (> (dist-cpt-cpt-homo p0 pm) tol))))
    (t nil)))

(defun srf-copy (sur-p)
  "A fresh copy of surface SUR-P.  NLib: N_SrfCopy."
  (declare (type surface sur-p))
  (let* ((n (surface-n sur-p)) (m (surface-m sur-p))
         (pw (surface-cpoints sur-p))
         (up (surface-knots-u sur-p)) (vp (surface-knots-v sur-p))
         (qw (make-cpoint-net (1+ n) (1+ m)))
         (uq (make-real-vector (length up)))
         (vq (make-real-vector (length vp))))
    (declare (type fixnum n m) (type cpoint-net pw qw) (type real-vector up vp uq vq))
    (dotimes (i (1+ n))
      (dotimes (j (1+ m))
        (dotimes (c 4) (setf (aref qw i j c) (aref pw i j c)))))
    (dotimes (i (length up)) (setf (aref uq i) (aref up i)))
    (dotimes (j (length vp)) (setf (aref vq j) (aref vp j)))
    (%make-surface (surface-degree-u sur-p) (surface-degree-v sur-p) uq vq qw)))

;;; --------------------------------------------------------------- construction

(defun srf-from-c-pts-and-knots (sur pw n m p q u v r s)
  "Make SUR the surface of control points PW[0..N][0..M], degrees P and Q and
knots U[0..R] and V[0..S], sharing the arrays as the C shares its pointers
(a corner or prefix is copied when an array is larger).  Returns SUR.
NLib: N_SrfFromCPtsAndKnots."
  (declare (type surface sur) (type cpoint-net pw) (type fixnum n m p q r s)
           (type real-vector u v))
  (let ((net (c-net-from-c-pts pw n m))
        (knu (knot-vector-from-real-array u r))
        (knv (knot-vector-from-real-array v s)))
    (srf-from-c-net-and-knot-vectors sur net p q knu knv)))

(defun srf-from-c-pt-coords-and-knots (sur wx wy wz w n m p q u v r s)
  "Make SUR the surface whose control points have coordinates WX, WY, WZ and
W (real-matrices, [0..N][0..M]; W may hold +no-w+), with degrees P and Q and
knots U[0..R] and V[0..S], sharing the knot arrays as the C shares its
pointers.  Returns SUR.  NLib: N_SrfFromCPtCoordsAndKnots."
  (declare (type surface sur) (type real-matrix wx wy wz w) (type fixnum n m p q r s)
           (type real-vector u v))
  (let* ((pw (xyz-to-2d-c-pt-array wx wy wz w n m))
         (net (c-net-from-c-pts pw n m))
         (knu (knot-vector-from-real-array u r))
         (knv (knot-vector-from-real-array v s)))
    (srf-from-c-net-and-knot-vectors sur net p q knu knv)))

(defun srf-compress (sur)
  "Give SUR fresh control point and knot arrays of exactly its sizes, copies
of its own, in place.  Returns SUR.  NLib: N_SrfCompress."
  (declare (type surface sur))
  (let* ((n (surface-n sur)) (m (surface-m sur))
         (pw (surface-cpoints sur))
         (u (surface-knots-u sur)) (v (surface-knots-v sur))
         (r (1- (length u))) (s (1- (length v)))
         (qw (make-cpoint-net (1+ n) (1+ m)))
         (ss (make-real-vector (1+ r)))
         (tt (make-real-vector (1+ s))))
    (declare (type fixnum n m r s) (type cpoint-net pw qw) (type real-vector u v ss tt))
    (fast
      (dotimes (i (1+ n))
        (dotimes (j (1+ m))
          (%copy-c-pt-net pw i j qw i j)))
      (dotimes (i (1+ r)) (setf (aref ss i) (aref u i)))
      (dotimes (j (1+ s)) (setf (aref tt j) (aref v j))))
    (srf-set-c-pts-and-knots sur qw ss tt)))

(defun swap-uv (sur)
  "Swap the roles of u and v in SUR, in place: S_new(v,u) = S_old(u,v), the
control net transposed into a fresh array and the degrees and knot vectors
exchanged.  Returns SUR.  NLib: N_SwapUV."
  (declare (type surface sur))
  (multiple-value-bind (net p q knu knv) (srf-get-net-and-knot-vectors sur)
    (declare (ignore net) (type fixnum p q))
    (let* ((n (surface-n sur)) (m (surface-m sur))
           (pw (surface-cpoints sur))
           (qw (make-cpoint-net (1+ m) (1+ n))))
      (declare (type fixnum n m) (type cpoint-net pw qw))
      (fast
        (dotimes (i (1+ n))
          (dotimes (j (1+ m))
            (%copy-c-pt-net pw i j qw j i))))
      (srf-from-c-net-and-knot-vectors sur (c-net-from-c-pts qw m n) q p knv knu))))

;;; ------------------------------------------------------------------ utilities

(defun srf-get-e-net (sur ku lu kv lv)
  "A fresh point net of the Euclidean images of control points
Pw[KU..LU][KV..LV] of SUR, entry [i-KU][j-KV].  Signals :ind-err unless
0 <= KU <= LU <= n and 0 <= KV <= LV <= m.  NLib: N_SrfGetENet."
  (declare (type surface sur) (type fixnum ku lu kv lv))
  (let ((n (surface-n sur)) (m (surface-m sur)) (pw (surface-cpoints sur)))
    (declare (type fixnum n m) (type cpoint-net pw))
    (when (or (> ku lu) (< ku 0) (> lu n))
      (nurbs-error :ind-err 'srf-get-e-net))
    (when (or (> kv lv) (< kv 0) (> lv m))
      (nurbs-error :ind-err 'srf-get-e-net))
    (let ((p (make-point-net (1+ (- lu ku)) (1+ (- lv kv)))))
      (declare (type point-net p))
      (fast
        (do ((i ku (1+ i))) ((> i lu))
          (declare (type fixnum i))
          (do ((j kv (1+ j))) ((> j lv))
            (declare (type fixnum j))
            (with-euclid ((x y z) (aref pw i j 0) (aref pw i j 1) (aref pw i j 2) (aref pw i j 3))
              (declare (type double-float x y z))
              (setf (aref p (- i ku) (- j kv) 0) x
                    (aref p (- i ku) (- j kv) 1) y
                    (aref p (- i ku) (- j kv) 2) z)))))
      (e-net-from-pts (- lu ku) (- lv kv) p))))

(defun srf-rat-to-non-rat (sur)
  "Make SUR non-rational in place: every weight set to +no-w+.  Returns SUR.
NLib: N_SrfRatToNonRat."
  (declare (type surface sur))
  (let ((pw (surface-cpoints sur)))
    (declare (type cpoint-net pw))
    (fast
      (dotimes (i (array-dimension pw 0))
        (dotimes (j (array-dimension pw 1))
          (setf (aref pw i j 3) +no-w+)))))
  sur)

(defun srf-non-rat-to-rat (sur)
  "Make SUR rational in place: every weight set to 1.0.  Returns SUR.
NLib: N_SrfNonRatToRat."
  (declare (type surface sur))
  (let ((pw (surface-cpoints sur)))
    (declare (type cpoint-net pw))
    (fast
      (dotimes (i (array-dimension pw 0))
        (dotimes (j (array-dimension pw 1))
          (setf (aref pw i j 3) 1d0)))))
  sur)

(defun srf-reparam-to-interval (sur r dir)
  "Rescale the knot vectors of SUR in place to the rectangle R (ul ur vb vt):
the u knots for DIR +udir+, the v knots for +vdir+, both for +uvdir+.  A knot
vector already on its interval is left alone.  Returns SUR.
NLib: N_SrfReparamToInterval."
  (declare (type surface sur) (type real-vector r) (type fixnum dir))
  (let* ((p (surface-degree-u sur)) (q (surface-degree-v sur))
         (u (surface-knots-u sur)) (v (surface-knots-v sur))
         (rr (1- (length u))) (s (1- (length v))))
    (declare (type fixnum p q rr s) (type real-vector u v))
    (multiple-value-bind (a b c d) (rectangle-get-data r)
      (declare (type double-float a b c d))
      (fast
        (when (or (= dir +udir+) (= dir +uvdir+))
          (when (or (/= a (aref u 0)) (/= b (aref u rr)))
            (let* ((u0 (aref u 0))
                   (fac (/ (- b a) (- (aref u rr) (aref u 0)))))
              (declare (type double-float u0 fac))
              (do ((i 0 (1+ i))) ((> i p))
                (declare (type fixnum i))
                (setf (aref u i) a))
              (do ((i (1+ p) (1+ i))) ((> i (- rr p 1)))
                (declare (type fixnum i))
                (setf (aref u i) (+ (* fac (- (aref u i) u0)) a)))
              (do ((i (- rr p) (1+ i))) ((> i rr))
                (declare (type fixnum i))
                (setf (aref u i) b)))))
        (when (or (= dir +vdir+) (= dir +uvdir+))
          (when (or (/= c (aref v 0)) (/= d (aref v s)))
            (let* ((v0 (aref v 0))
                   (fac (/ (- d c) (- (aref v s) (aref v 0)))))
              (declare (type double-float v0 fac))
              (do ((j 0 (1+ j))) ((> j q))
                (declare (type fixnum j))
                (setf (aref v j) c))
              (do ((j (1+ q) (1+ j))) ((> j (- s q 1)))
                (declare (type fixnum j))
                (setf (aref v j) (+ (* fac (- (aref v j) v0)) c)))
              (do ((j (- s q) (1+ j))) ((> j s))
                (declare (type fixnum j))
                (setf (aref v j) d))))))))
  sur)

;;; ------------------------------------------------------------ level 2, group 5

(defun make-srfs-rat-compatible (sur k)
  "Make the surfaces SUR[0..K] (a simple-vector of surfaces) all rational,
in place, when any one of them is rational; otherwise leave them alone.
Return SUR.
NLib: N_MakeSrfsRatCompatible."
  (declare (type simple-vector sur) (type fixnum k))
  (let ((allnra +true+))
    (declare (type fixnum allnra))
    (do ((i 0 (1+ i))) ((> i k))
      (declare (type fixnum i))
      (when (is-srf-rat (svref sur i))
        (setf allnra +false+)
        (return)))
    (when (= allnra +false+)
      (do ((i 0 (1+ i))) ((> i k))
        (declare (type fixnum i))
        (unless (is-srf-rat (svref sur i))
          (srf-non-rat-to-rat (svref sur i))))))
  sur)

(defun srf-is-flat-cheap (sur eps)
  "Whether surface SUR is flat to within EPS, cheaply: every Euclidean
control point within EPS of the plane through the corners' centroid whose
normal is the cross product of the vectors from that centroid to the
midpoints of the corner edges P00-P10 and P10-P11.  Returns (values flt dir): FLT +yes+ or +no+, DIR +udir+
(subdivide at a v value) or +vdir+ (at a u value).  A surface closed in u
or v, or one whose corner plane is degenerate, is not flat.  Signals as
plane-pt-normal-to-implicit does.
NLib: N_SrfIsFlatCheap."
  (declare (type surface sur) (type double-float eps))
  (when (srf-is-closed sur +udir+)
    (return-from srf-is-flat-cheap (values +no+ +vdir+)))
  (when (srf-is-closed sur +vdir+)
    (return-from srf-is-flat-cheap (values +no+ +udir+)))
  (let ((n (surface-n sur)) (m (surface-m sur))
        (pw (surface-cpoints sur))
        (cw (make-cpoint))
        (p00 (make-point)) (p10 (make-point)) (p01 (make-point)) (p11 (make-point))
        (p (make-point)) (a (make-point)) (b (make-point))
        (nv (make-point)) (cv (make-point)) (dv (make-point))
        (flt +yes+) (dir +udir+))
    (declare (type fixnum n m flt dir) (type cpoint-net pw)
             (dynamic-extent cw p00 p10 p01 p11 p a b nv cv dv))
    (c-pt-to-pt-euclid (%net-c-pt pw 0 0 cw) p00)
    (c-pt-to-pt-euclid (%net-c-pt pw n 0 cw) p10)
    (c-pt-to-pt-euclid (%net-c-pt pw 0 m cw) p01)
    (c-pt-to-pt-euclid (%net-c-pt pw n m cw) p11)
    (combine-4-pts 0.25d0 p00 0.25d0 p10 0.25d0 p01 0.25d0 p11 p)
    (vector-combine 0.5d0 p00 0.5d0 p10 a)
    (vector-combine 0.5d0 p10 0.5d0 p11 b)
    (vector-diff a p cv)
    (vector-diff b p dv)
    (let ((lv (vector-magnitude cv))
          (lu (vector-magnitude dv))
          (ln 0d0))
      (declare (type double-float lv lu ln))
      (vector-cross cv dv nv)
      (setf ln (vector-magnitude nv))
      (when (< ln (the double-float *mtol*))
        (return-from srf-is-flat-cheap
          (values +no+ (if (> lu lv) +vdir+ +udir+))))
      (multiple-value-bind (pa pb pc pd) (plane-pt-normal-to-implicit p nv)
        (declare (type double-float pa pb pc pd))
        (block scan
          (dotimes (i (1+ n))
            (dotimes (j (1+ m))
              (c-pt-to-pt-euclid (%net-c-pt pw i j cw) p)
              (let ((dis (dist-signed-pt-plane pa pb pc pd p)))
                (declare (type double-float dis))
                (when (> (abs dis) eps)
                  (setf flt +no+)
                  (return-from scan))))))
        (setf dir (if (> lu lv) +vdir+ +udir+))
        (values flt dir)))))

;;; ------------------------------------------------------------ level 14, group 0

(defun clamp-srf-at-params (sur u v)
  "U and V clamped to the knot rectangle of surface SUR: (values u v).
NLib: N_ClampSrfAtParams."
  (declare (type surface sur) (type double-float u v))
  (let* ((ku (surface-knots-u sur)) (kv (surface-knots-v sur))
         (r (1- (length ku))) (s (1- (length kv))))
    (declare (type real-vector ku kv) (type fixnum r s))
    (when (< u (aref ku 0)) (setf u (aref ku 0)))
    (when (> u (aref ku r)) (setf u (aref ku r)))
    (when (< v (aref kv 0)) (setf v (aref kv 0)))
    (when (> v (aref kv s)) (setf v (aref kv s)))
    (values u v)))

(defun srf-is-degen (sur)
  "True when surface SUR is degenerate to a point: every Euclidean control
point within *mtol* of their centroid.
NLib: N_SrfIsDegen."
  (declare (type surface sur))
  (let* ((n (surface-n sur)) (m (surface-m sur))
         (pw (surface-cpoints sur))
         (fac (/ 1d0 (* (+ (float n 1d0) 1d0) (+ (float m 1d0) 1d0))))
         (mtol *mtol*)
         (cw (make-cpoint)) (q (make-point)) (mm (make-point)))
    (declare (type fixnum n m) (type cpoint-net pw) (type double-float fac mtol)
             (type cpoint cw) (type point q mm) (dynamic-extent cw q mm))
    (copy-pt (make-point) mm)
    (dotimes (i (1+ n))
      (dotimes (j (1+ m))
        (c-pt-to-pt-euclid (%net-c-pt pw i j cw) q)
        (sum-2-pts mm q mm)))
    (scale-pt fac mm mm)
    (dotimes (i (1+ n) t)
      (dotimes (j (1+ m))
        (c-pt-to-pt-euclid (%net-c-pt pw i j cw) q)
        (when (> (dist-pt-pt q mm) mtol)
          (return-from srf-is-degen nil))))))

(defun srf-is-singular (sur side)
  "True when surface SUR has a pole on boundary SIDE: +left+ (u = umin),
+right+ (u = umax), +bottom+ (v = vmin) or +top+ (v = vmax).  The boundary's
first control point must lie within *mtol* of the one past its middle, and
every control point of the boundary within *mtol* of their centroid.  False
for any other SIDE.
NLib: N_SrfIsSingular."
  (declare (type surface sur) (type fixnum side))
  (let* ((n (surface-n sur)) (m (surface-m sur))
         (pw (surface-cpoints sur))
         (mtol *mtol*)
         (cw (make-cpoint)) (q (make-point)) (mm (make-point)))
    (declare (type fixnum n m) (type cpoint-net pw) (type double-float mtol)
             (type cpoint cw) (type point q mm) (dynamic-extent cw q mm))
    (flet ((pole (fixed k mid)
             ;; FIXED: the fixed index (a row for u, a column for v); K the
             ;; highest index along the boundary, MID the probe index
             (flet ((entry (s)
                      (if (or (= side +left+) (= side +right+))
                          (%net-c-pt pw fixed s cw)
                          (%net-c-pt pw s fixed cw))))
               (c-pt-to-pt-euclid (entry 0) mm)
               (c-pt-to-pt-euclid (entry mid) q)
               (when (> (dist-pt-pt mm q) mtol)
                 (return-from srf-is-singular nil))
               (let ((fac (/ 1d0 (+ (float k 1d0) 1d0))))
                 (declare (type double-float fac))
                 (do ((s 1 (1+ s))) ((> s k))
                   (declare (type fixnum s))
                   (c-pt-to-pt-euclid (entry s) q)
                   (sum-2-pts mm q mm))
                 (scale-pt fac mm mm))
               (do ((s 0 (1+ s))) ((> s k) t)
                 (declare (type fixnum s))
                 (c-pt-to-pt-euclid (entry s) q)
                 (when (> (dist-pt-pt q mm) mtol)
                   (return-from srf-is-singular nil))))))
      (cond ((= side +left+) (pole 0 m (1+ (truncate m 2))))
            ((= side +right+) (pole n m (1+ (truncate m 2))))
            ((= side +bottom+) (pole 0 n (1+ (truncate n 2))))
            ((= side +top+) (pole m n (1+ (truncate n 2))))
            (t nil)))))

;;; ------------------------------------------------------------ level 15, group 0

(defun srf-is-flat (sur eps)
  "Whether surface SUR can be approximated to within EPS by a planar
quadrilateral, and if not which way to subdivide it.  Returns (values flt
dir): FLT +yes+ or +no+; DIR +udir+ (subdivide at a v value), +vdir+ (at a u
value), or NIL where the C leaves *dir unset (always when FLT is +yes+, and
for a surface with two opposite poles whose corners are not pairwise
coincident).  A surface closed in u or v is flat only when its control net
collapses across the closure and its boundary iso-curve is a line.  Otherwise
every Euclidean control point must lie within EPS of the plane through the
corners' centroid normal to the averaged corner normals (from evaluated
normals near the corners when those cancel), and the interior control points
of the boundaries within EPS of the chords between the corners.  Signals as
srf-extract-iso-crv and dist-pt-line-seg do, and :geo-err when no plane
normal is found.
NLib: N_SrfIsFlat."
  (declare (type surface sur) (type double-float eps))
  (let* ((n (surface-n sur)) (m (surface-m sur))
         (pw (surface-cpoints sur))
         (mtol *mtol*)
         (cw (make-cpoint)) (cw2 (make-cpoint)))
    (declare (type fixnum n m) (type cpoint-net pw) (type double-float mtol)
             (dynamic-extent cw cw2))
    ;; closed in u
    (when (srf-is-closed sur +udir+)
      (let ((degen t))
        (block scan
          (dotimes (j (1+ m))
            (do ((i 1 (1+ i))) ((> i n))
              (declare (type fixnum i))
              (when (> (dist-cpt-cpt (%net-c-pt pw (1- i) j cw) (%net-c-pt pw i j cw2)) mtol)
                (setf degen nil)
                (return-from scan)))))
        (when degen
          (let ((us (srf-get-parameter-bounds sur)))
            (declare (type double-float us))
            (return-from srf-is-flat
              (if (crv-is-line (srf-extract-iso-crv sur us +vdir+) eps)
                  (values +yes+ nil)
                  (values +no+ +udir+)))))
        (return-from srf-is-flat (values +no+ +vdir+))))
    ;; closed in v
    (when (srf-is-closed sur +vdir+)
      (let ((degen t))
        (block scan
          (dotimes (i (1+ n))
            (do ((j 1 (1+ j))) ((> j m))
              (declare (type fixnum j))
              (when (> (dist-cpt-cpt (%net-c-pt pw i (1- j) cw) (%net-c-pt pw i j cw2)) mtol)
                (setf degen nil)
                (return-from scan)))))
        (when degen
          (multiple-value-bind (us ue vs) (srf-get-parameter-bounds sur)
            (declare (ignore us ue) (type double-float vs))
            (return-from srf-is-flat
              (if (crv-is-line (srf-extract-iso-crv sur vs +udir+) eps)
                  (values +yes+ nil)
                  (values +no+ +vdir+)))))
        (return-from srf-is-flat (values +no+ +udir+))))
    (let ((dist (make-real-matrix (1+ n) (1+ m)))
          (sd (make-point-net 2 2))
          (p00 (make-point)) (p10 (make-point)) (p01 (make-point)) (p11 (make-point))
          (p (make-point)) (q (make-point))
          (n00 (make-point)) (n10 (make-point)) (n01 (make-point)) (n11 (make-point))
          (nn (make-point))
          (va (make-point)) (vb (make-point)) (vc (make-point)) (vd (make-point))
          (flat t) (straight nil)
          (maxu 0d0) (maxv 0d0) (du 0d0) (dv 0d0)
          (dir nil))
      (declare (type real-matrix dist) (type point-net sd)
               (dynamic-extent p00 p10 p01 p11 p q n00 n10 n01 n11 nn va vb vc vd)
               (type double-float maxu maxv du dv))
      ;; the four corner points
      (c-pt-to-pt-euclid (%net-c-pt pw 0 0 cw) p00)
      (c-pt-to-pt-euclid (%net-c-pt pw n 0 cw) p10)
      (c-pt-to-pt-euclid (%net-c-pt pw 0 m cw) p01)
      (c-pt-to-pt-euclid (%net-c-pt pw n m cw) p11)
      ;; degenerate boundaries
      (when (or (and (srf-is-singular sur +left+) (srf-is-singular sur +right+))
                (and (srf-is-singular sur +bottom+) (srf-is-singular sur +top+)))
        (let ((a (dist-pt-pt p00 p10)) (b (dist-pt-pt p01 p11))
              (c (dist-pt-pt p00 p01)) (d (dist-pt-pt p10 p11)))
          (declare (type double-float a b c d))
          (return-from srf-is-flat
            (cond ((and (< a mtol) (< b mtol) (< c mtol) (< d mtol))
                   (if (srf-is-degen sur)
                       (values +yes+ nil)
                       (values +no+ +udir+)))
                  ((and (< a mtol) (< b mtol))
                   (values +no+ +udir+))
                  ((and (< c mtol) (< d mtol))
                   (values +no+ +vdir+))
                  (t (values +no+ nil))))))
      ;; the plane of the four corner points
      (combine-4-pts 0.25d0 p00 0.25d0 p10 0.25d0 p01 0.25d0 p11 p)
      (vector-diff p10 p00 va)
      (vector-diff p01 p00 vb)
      (vector-cross va vb n00)
      (vector-diff p11 p10 vc)
      (vector-cross va vc n10)
      (vector-diff p11 p01 vd)
      (vector-cross vd vb n01)
      (vector-cross vd vc n11)
      (dotimes (i 2)
        (let* ((a (vector-magnitude n00)) (b (vector-magnitude n10))
               (c (vector-magnitude n01)) (d (vector-magnitude n11))
               (mx a))
          (declare (type double-float a b c d mx))
          ;; A is a copy of the longest corner normal
          (vector-copy n00 va)
          (when (> b mx) (setf mx b) (vector-copy n10 va))
          (when (> c mx) (setf mx c) (vector-copy n01 va))
          (when (> d mx) (setf mx d) (vector-copy n11 va)))
        (when (< (vector-dot va n00) 0d0) (vector-scale n00 -1d0 n00))
        (when (< (vector-dot va n10) 0d0) (vector-scale n10 -1d0 n10))
        (when (< (vector-dot va n01) 0d0) (vector-scale n01 -1d0 n01))
        (when (< (vector-dot va n11) 0d0) (vector-scale n11 -1d0 n11))
        (combine-4-pts 0.25d0 n00 0.25d0 n10 0.25d0 n01 0.25d0 n11 nn)
        (if (= i 1)
            (vector-normalize-ref nn)   ; its error is the C's NL_OUT
            (when (handler-case (progn (vector-normalize-ref nn) t)
                    (nurbs-error () nil))
              (return)))
        ;; corner normals evaluated just inside the corners; the C ignores
        ;; their errors
        (multiple-value-bind (us ue vs ve) (srf-get-parameter-bounds sur)
          (declare (type double-float us ue vs ve))
          (let ((du (* 0.01d0 (- ue us)))
                (dv (* 0.01d0 (- ve vs))))
            (declare (type double-float du dv))
            (setf us (+ us du) ue (- ue du)
                  vs (+ vs dv) ve (- ve dv))
            (flet ((normal (u v nv)
                     (handler-case
                         (srf-eval-pt-pt-deriv-normal-fast sur u v +left+ +left+ q vb vc nv sd)
                       (nurbs-error () nil))))
              (normal us vs n00)
              (normal us ve n01)
              (normal ue vs n10)
              (normal ue ve n11)))))
      (let ((pd (- (vector-dot p nn))))
        (declare (type double-float pd))
        (multiple-value-bind (a b c) (pt-to-xyz nn)
          (declare (type double-float a b c))
          ;; flatness of the control points
          (dotimes (i (1+ n))
            (dotimes (j (1+ m))
              (c-pt-to-pt-euclid (%net-c-pt pw i j cw) p)
              (setf (aref dist i j) (dist-signed-pt-plane a b c pd p))
              (setf (aref dist i j) (abs (aref dist i j)))
              (when (> (aref dist i j) eps)
                (setf flat nil))))))
      ;; straightness of the boundary curves
      (when flat
        (setf maxu 0d0 maxv 0d0
              straight t)
        (do ((i 1 (1+ i))) ((>= i n))
          (declare (type fixnum i))
          (c-pt-to-pt-euclid (%net-c-pt pw i 0 cw) p)
          (setf du (dist-pt-line-seg p p00 p10))
          (when (> du maxu) (setf maxu du))
          (c-pt-to-pt-euclid (%net-c-pt pw i m cw) p)
          (setf du (dist-pt-line-seg p p01 p11))
          (when (> du maxu) (setf maxu du)))
        (do ((j 1 (1+ j))) ((>= j m))
          (declare (type fixnum j))
          (c-pt-to-pt-euclid (%net-c-pt pw 0 j cw) p)
          (setf dv (dist-pt-line-seg p p00 p01))
          (when (> dv maxv) (setf maxv dv))
          (c-pt-to-pt-euclid (%net-c-pt pw n j cw) p)
          (setf dv (dist-pt-line-seg p p10 p11))
          (when (> dv maxv) (setf maxv dv)))
        (when (or (> maxu eps) (> maxv eps))
          (setf straight nil)))
      ;; which direction to subdivide
      (if (or (not flat) (not straight))
          (progn
            (if (not flat)
                (progn
                  (setf du 0d0)
                  (dotimes (j (1+ m))
                    (let ((mn (aref dist 0 j)) (mx (aref dist 0 j)))
                      (declare (type double-float mn mx))
                      (dotimes (i (1+ n))
                        (when (< (aref dist i j) mn) (setf mn (aref dist i j)))
                        (when (> (aref dist i j) mx) (setf mx (aref dist i j))))
                      (setf du (+ du (- mx mn)))))
                  (setf du (/ du (+ (float m 1d0) 1d0)))
                  (setf dv 0d0)
                  (dotimes (i (1+ n))
                    (let ((mn (aref dist i 0)) (mx (aref dist i 0)))
                      (declare (type double-float mn mx))
                      (dotimes (j (1+ m))
                        (when (< (aref dist i j) mn) (setf mn (aref dist i j)))
                        (when (> (aref dist i j) mx) (setf mx (aref dist i j))))
                      (setf dv (+ dv (- mx mn)))))
                  (setf dv (/ dv (+ (float n 1d0) 1d0)))
                  (setf dir (if (> du dv) +vdir+ +udir+)))
                (setf dir (if (> maxu maxv) +vdir+ +udir+)))
            (values +no+ dir))
          (values +yes+ dir)))))
