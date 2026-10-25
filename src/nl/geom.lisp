;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Geom.cpp (usd-brep source/SMLib/src).  Modified:
;;;; translated to Common Lisp.  Vector arguments, passed by value in C, are
;;;; arrays here; every routine reads its inputs before writing its output.

(in-package :batten)

(declaim (inline vector-dot vector-magnitude vector-magnitude-ref vector-scale vector-scale-ref
                 vector-cross vector-diff vector-copy vector-sum vector-pt-along-vector
                 vector-dir vector-combine dist-pt-pt dist-cpt-cpt-homo))

(defun vector-dot (a b)
  "A . B.  NLib: N_VectorDot."
  (declare (type point a b))
  (with-point (ax ay az) a (with-point (bx by bz) b (+ (+ (* ax bx) (* ay by)) (* az bz)))))

(defun vector-magnitude (a)
  "|A|.  NLib: N_VectorMagnitude."
  (declare (type point a))
  (sqrt (the (double-float 0d0) (vector-dot a a))))

(defun vector-magnitude-ref (a)
  "|A|.  NLib: N_VectorMagnitudeRef."
  (declare (type point a))
  (sqrt (the (double-float 0d0) (vector-dot a a))))

(defun vector-scale (a alpha b)
  "B = ALPHA*A.  NLib: N_VectorScale."
  (declare (type point a b) (type double-float alpha))
  (with-point (x y z) a (set-point b (* alpha x) (* alpha y) (* alpha z))))

(defun vector-scale-ref (a alpha b)
  "B = ALPHA*A.  NLib: N_VectorScaleRef."
  (declare (type point a b) (type double-float alpha))
  (with-point (x y z) a (set-point b (* alpha x) (* alpha y) (* alpha z))))

(defun vector-cross (a b c)
  "C = A x B.  NLib: N_VectorCross."
  (declare (type point a b c))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (set-point c (- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))))))

(defun vector-diff (a b c)
  "C = A - B.  NLib: N_VectorDiff."
  (declare (type point a b c))
  (with-point (ax ay az) a (with-point (bx by bz) b (set-point c (- ax bx) (- ay by) (- az bz)))))

(defun vector-copy (a b)
  "Copy A into B.  NLib: N_VectorCopy."
  (declare (type point a b))
  (with-point (x y z) a (set-point b x y z)))

(defun vector-sum (a b c)
  "C = A + B.  NLib: N_VectorSum."
  (declare (type point a b c))
  (with-point (ax ay az) a (with-point (bx by bz) b (set-point c (+ ax bx) (+ ay by) (+ az bz)))))

(defun vector-pt-along-vector (p alpha v q)
  "Q = P + ALPHA*V.  NLib: N_VectorPtAlongVector."
  (declare (type point p v q) (type double-float alpha))
  (with-point (px py pz) p
    (with-point (vx vy vz) v
      (set-point q (+ px (* alpha vx)) (+ py (* alpha vy)) (+ pz (* alpha vz))))))

(defun vector-dir (ps pe d)
  "D = PE - PS.  NLib: N_VectorDir."
  (declare (type point ps pe d))
  (with-point (sx sy sz) ps (with-point (ex ey ez) pe (set-point d (- ex sx) (- ey sy) (- ez sz)))))

(defun vector-combine (alf a bet b c)
  "C = ALF*A + BET*B.  NLib: N_VectorCombine."
  (declare (type double-float alf bet) (type point a b c))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (set-point c (+ (* alf ax) (* bet bx)) (+ (* alf ay) (* bet by)) (+ (* alf az) (* bet bz))))))

(defun dist-pt-pt (p q)
  "The distance between P and Q.  NLib: N_DistPtPt."
  (declare (type point p q))
  (with-point (px py pz) p
    (with-point (qx qy qz) q
      (let ((dx (- px qx)) (dy (- py qy)) (dz (- pz qz)))
        (sqrt (the (double-float 0d0) (+ (+ (* dx dx) (* dy dy)) (* dz dz))))))))

(defun dist-cpt-cpt-homo (pw qw)
  "The distance between PW and QW in homogeneous space, using z and w only
where both points have them.  NLib: N_DistCptCptHomo."
  (declare (type cpoint pw qw))
  (with-cpoint (px py pz pw*) pw
    (with-cpoint (qx qy qz qw*) qw
      (let ((disc (+ (* (- px qx) (- px qx)) (* (- py qy) (- py qy)))))
        (declare (type double-float disc))
        (when (and (/= pz +no-z+) (/= qz +no-z+))
          (incf disc (* (- pz qz) (- pz qz))))
        (when (and (/= pw* +no-w+) (/= qw* +no-w+))
          (incf disc (* (- pw* qw*) (- pw* qw*))))
        (sqrt (the (double-float 0d0) disc))))))

(defun vector-normalize-ref (a)
  "Scale A to unit length in place.  Signals :geo-err for a zero vector.
NLib: N_VectorNormalizeRef."
  (declare (type point a))
  (let ((d (vector-magnitude-ref a)))
    (declare (type double-float d))
    (when (float-op-is-bad 1d0 d +division+)
      (nurbs-error :geo-err 'vector-normalize-ref))
    (vector-scale-ref a (/ 1d0 d) a)))

(defun vector-normalize (a b)
  "Write A scaled to unit length into B, and return |A|.  For a zero vector B
and the magnitude are zeroed and :geo-err is signaled.  NLib: N_VectorNormalize."
  (declare (type point a b))
  (let ((d (vector-magnitude a)))
    (declare (type double-float d))
    (when (float-op-is-bad 1d0 d +division+)
      (set-point b 0d0 0d0 0d0)
      (nurbs-error :geo-err 'vector-normalize))
    (vector-scale a (/ 1d0 d) b)
    d))

;;; ------------------------------------------------------------ level 0, group 1

(defun vector-dot-ref (a b)
  "A . B, for vectors C passes by reference.
NLib: N_VectorDotRef."
  (declare (type point a b))
  (with-point (ax ay az) a (with-point (bx by bz) b (+ (+ (* ax bx) (* ay by)) (* az bz)))))

(defun vector-reverse (a b)
  "B = -A.
NLib: N_VectorReverse."
  (declare (type point a b))
  (with-point (x y z) a (set-point b (- x) (- y) (- z))))

(defun vector-reverse-in-place (a)
  "A = -A, in place.
NLib: N_VectorReverseInPlace."
  (declare (type point a))
  (with-point (x y z) a (set-point a (- x) (- y) (- z))))

(defun vector-create (x y z a)
  "Set vector A to (X Y Z).
NLib: N_VectorCreate."
  (declare (type double-float x y z) (type point a))
  (set-point a x y z))

(defun vector-perpendicular (a b)
  "B = the 2D vector A rotated 90 degrees counterclockwise, z zeroed.
NLib: N_VectorPerpendicular."
  (declare (type point a b))
  (with-point (x y z) a (set-point b (- y) x 0d0)))

(defun vectors-angle (a b)
  "The angle between vectors A and B in degrees, 0 to 180.  Signals :geo-err
when the cosine cannot be formed (a zero vector).
NLib: N_VectorsAngle."
  (declare (type point a b))
  (let ((dot (vector-dot a b))
        (ma (vector-magnitude a))
        (mb (vector-magnitude b)))
    (declare (type double-float dot ma mb))
    (when (float-op-is-bad dot (* ma mb) +division+)
      (nurbs-error :geo-err 'vectors-angle))
    (let ((arg (/ dot (* ma mb))))
      (declare (type double-float arg))
      (when (> arg 1d0) (setf arg 1d0))
      (when (< arg -1d0) (setf arg -1d0))
      (let ((rad (the double-float (acos arg))))
        (* 180d0 (/ rad +pi+))))))

(defun project-pt-line-param (p a b q)
  "Project point P onto the segment AB: write the projection into Q and return
the parameter t with Q = A + t*(B-A).  A degenerate segment projects to its
midpoint with t = 0.  Signals :num-err when the division fails.
NLib: N_ProjectPtLineParam."
  (declare (type point p a b q))
  (let ((ba (make-point)) (pa (make-point)) (tt (make-point)))
    (declare (dynamic-extent ba pa tt))
    (let ((dis (dist-pt-pt a b)))
      (declare (type double-float dis))
      (when (< dis *mtol*)
        (vector-sum a b tt)
        (vector-scale tt 0.5d0 q)
        (return-from project-pt-line-param 0d0))
      (setf dis (dist-pt-pt p a))
      (when (< dis *mtol*)
        (vector-copy a q)
        (return-from project-pt-line-param 0d0))
      (setf dis (dist-pt-pt p b))
      (when (< dis *mtol*)
        (vector-copy b q)
        (return-from project-pt-line-param 1d0))
      (vector-diff b a ba)
      (vector-diff p a pa)
      (let ((len (vector-dot ba ba))
            (dot (vector-dot pa ba)))
        (declare (type double-float len dot))
        (setf len (sqrt (the (double-float 0d0) len)))
        (when (float-op-is-bad dot len +division+)
          (nurbs-error :num-err 'project-pt-line-param))
        (setf dis (/ dot len))
        (cond ((< (abs dis) *mtol*)
               (vector-copy a q)
               0d0)
              ((< (abs (- dis len)) *mtol*)
               (vector-copy b q)
               1d0)
              (t
               (let ((u (/ dis len)))
                 (declare (type double-float u))
                 (vector-pt-along-vector a u ba q)
                 u)))))))

(defun dist-sq-pt-pt (p q)
  "The squared distance between P and Q.
NLib: N_DistSqPtPt."
  (declare (type point p q))
  (with-point (px py pz) p
    (with-point (qx qy qz) q
      (+ (+ (* (- px qx) (- px qx)) (* (- py qy) (- py qy))) (* (- pz qz) (- pz qz))))))

(defun dist-sq-pt-pt-2d (p q)
  "The squared distance between P and Q in x and y.
NLib: N_DistSqPtPt2d."
  (declare (type point p q))
  (with-point (px py pz) p
    (with-point (qx qy qz) q
      (+ (* (- px qx) (- px qx)) (* (- py qy) (- py qy))))))

(defun pts-are-equal (p q)
  "True when P and Q are exactly equal, coordinate by coordinate.
NLib: N_PtsAreEqual."
  (declare (type point p q))
  (with-point (px py pz) p
    (with-point (qx qy qz) q
      (and (= px qx) (= py qy) (= pz qz)))))

(defun dist-pt-line-seg (p a b)
  "The distance between point P and the segment AB: perpendicular when P
projects inside AB, otherwise to the nearer end.  A degenerate segment is its
midpoint.  Signals :num-err when the division fails.
NLib: N_DistPtLineSeg."
  (declare (type point p a b))
  (let ((ba (make-point)) (pa (make-point)) (bp (make-point)) (tt (make-point))
        (dotab (dist-pt-pt a b)))
    (declare (dynamic-extent ba pa bp tt) (type double-float dotab))
    (when (< dotab *mtol*)
      (vector-sum a b tt)
      (vector-scale tt 0.5d0 tt)
      (return-from dist-pt-line-seg (dist-pt-pt p tt)))
    (vector-diff b a ba)
    (vector-diff p a pa)
    (vector-diff b p bp)
    (let ((dotap (vector-dot pa ba))
          (dotbp (vector-dot bp ba)))
      (declare (type double-float dotap dotbp))
      (cond ((and (>= dotap 0d0) (>= dotbp 0d0))
             (vector-cross pa ba tt)
             (let ((num (vector-magnitude tt)) (den dotab))
               (declare (type double-float num den))
               (when (float-op-is-bad num den +division+)
                 (nurbs-error :num-err 'dist-pt-line-seg))
               (/ num den)))
            ((< dotap 0d0) (dist-pt-pt p a))
            ((< dotbp 0d0) (dist-pt-pt p b))
            ;; only a NaN reaches here; the C leaves *d unset
            (t 0d0)))))

(defun dist-perp-pt-line-seg (q a b)
  "The perpendicular distance between point Q and the segment AB, or -1 when Q
does not project inside AB.  A degenerate segment gives the distance to its
midpoint.  Signals :num-err when the division fails.
NLib: N_DistPerpPtLineSeg."
  (declare (type point q a b))
  (let ((ba (make-point)) (qa (make-point)) (bq (make-point)) (tt (make-point))
        (dotab (dist-pt-pt a b)))
    (declare (dynamic-extent ba qa bq tt) (type double-float dotab))
    (when (< dotab *mtol*)
      (vector-sum a b tt)
      (vector-scale tt 0.5d0 tt)
      (return-from dist-perp-pt-line-seg (dist-pt-pt q tt)))
    (vector-diff b a ba)
    (vector-diff q a qa)
    (vector-diff b q bq)
    (let ((dotaq (vector-dot qa ba))
          (dotbq (vector-dot bq ba)))
      (declare (type double-float dotaq dotbq))
      (when (or (< dotaq 0d0) (< dotbq 0d0))
        (return-from dist-perp-pt-line-seg -1d0))
      (setf dotab (vector-dot ba ba))
      (vector-cross qa ba tt)
      (let ((num (vector-magnitude tt)))
        (declare (type double-float num))
        (when (< dotab 0d0) (setf dotab 0d0))
        (let ((den (sqrt (the (double-float 0d0) dotab))))
          (declare (type double-float den))
          (when (float-op-is-bad num den +division+)
            (nurbs-error :num-err 'dist-perp-pt-line-seg))
          (/ num den))))))

(defun dist-polygon (p n)
  "The length of the polygon P[0..N].
NLib: N_DistPolygon."
  (declare (type point-array p) (type fixnum n))
  (let ((len 0d0) (a (make-point)) (b (make-point)))
    (declare (type double-float len) (dynamic-extent a b))
    (do ((i 1 (1+ i))) ((> i n) len)
      (declare (type fixnum i))
      (set-point a (aref p (1- i) 0) (aref p (1- i) 1) (aref p (1- i) 2))
      (set-point b (aref p i 0) (aref p i 1) (aref p i 2))
      (incf len (dist-pt-pt a b)))))

(defun dist-cpt-cpt (pw qw)
  "The Euclidean distance between control points PW and QW.
NLib: N_DistCptCpt."
  (declare (type cpoint pw qw))
  (let ((p (make-point)) (q (make-point)))
    (declare (dynamic-extent p q))
    (c-pt-to-pt-euclid pw p)
    (c-pt-to-pt-euclid qw q)
    (dist-pt-pt p q)))

(defun line-set-bounding-flag (lsg flg)
  "Set the bounding flag of line LSG to FLG (+bounded+ or +unbounded+).
NLib: N_LineSetBoundingFlag."
  (declare (type lineseg lsg) (type fixnum flg))
  (setf (lineseg-bounded lsg) flg)
  lsg)

(defun line-get-data (lsg p v)
  "Write the start point and direction of line LSG into P and V, and return
its bounding flag.
NLib: N_LineGetData."
  (declare (type lineseg lsg) (type point p v))
  (vector-copy (lineseg-p lsg) p)
  (vector-copy (lineseg-v lsg) v)
  (lineseg-bounded lsg))

(defun create-line-start-dir-vector (lsg p v flg)
  "Make LSG the line from P along V with bounding flag FLG.
NLib: N_CreateLineStartDirVector."
  (declare (type lineseg lsg) (type point p v) (type fixnum flg))
  (with-point (px py pz) p
    (with-point (vx vy vz) v
      (set-point (lineseg-p lsg) px py pz)
      (set-point (lineseg-v lsg) vx vy vz)))
  (setf (lineseg-bounded lsg) flg)
  lsg)

(defun create-plane-pt-normal (pln p n)
  "Make PLN the plane through P with normal N.
NLib: N_CreatePlanePtNormal."
  (declare (type plane pln) (type point p n))
  (with-point (px py pz) p
    (with-point (nx ny nz) n
      (set-point (plane-p pln) px py pz)
      (set-point (plane-n pln) nx ny nz)))
  pln)

(defun create-interval (i ul ur)
  "Set interval I (a real-vector of 2) to [UL, UR].
NLib: N_CreateInterval."
  (declare (type real-vector i) (type double-float ul ur))
  (setf (aref i 0) ul (aref i 1) ur)
  i)

(defun interval-get-data (i)
  "The ends of interval I as two values, ul and ur.
NLib: N_IntervalGetData."
  (declare (type real-vector i))
  (values (aref i 0) (aref i 1)))

(defun create-rectangle (r ul ur vb vt)
  "Set rectangle R (a real-vector of 4) to (UL UR VB VT).
NLib: N_CreateRectangle."
  (declare (type real-vector r) (type double-float ul ur vb vt))
  (setf (aref r 0) ul (aref r 1) ur (aref r 2) vb (aref r 3) vt)
  r)

(defun rectangle-get-data (r)
  "The sides of rectangle R as four values: ul ur vb vt.
NLib: N_RectangleGetData."
  (declare (type real-vector r))
  (values (aref r 0) (aref r 1) (aref r 2) (aref r 3)))

(defun get-b-box-data (box)
  "The extents of min-max box BOX as six values: xl xr yb yt zn zf.
NLib: N_GetBBoxData."
  (declare (type real-vector box))
  (values (aref box 0) (aref box 1) (aref box 2) (aref box 3) (aref box 4) (aref box 5)))

(defun b-box-define (box xl xr yb yt zn zf)
  "Set min-max box BOX (a real-vector of 6) to (XL XR YB YT ZN ZF).
NLib: N_BBoxDefine."
  (declare (type real-vector box) (type double-float xl xr yb yt zn zf))
  (setf (aref box 0) xl (aref box 1) xr (aref box 2) yb
        (aref box 3) yt (aref box 4) zn (aref box 5) zf)
  box)

(defun plane-get-data (pln p n)
  "Write the point and normal of plane PLN into P and N.
NLib: N_PlaneGetData."
  (declare (type plane pln) (type point p n))
  (vector-copy (plane-p pln) p)
  (vector-copy (plane-n pln) n)
  p)

(defun plane-get-normal (pln n)
  "Write the normal of plane PLN into N.
NLib: N_PlaneGetNormal."
  (declare (type plane pln) (type point n))
  (vector-copy (plane-n pln) n))

(defun calc-circ-arc-derivs (o x y r ang k d)
  "Write into D[0..K] (a point-array of K+1 rows, allocated by the caller) the
point and derivatives at sweep angle ANG degrees of the circle O + R*cos*X +
R*sin*Y, X and Y orthogonal unit axes.
NLib: N_CalcCircArcDerivs."
  (declare (type point o x y) (type double-float r ang) (type fixnum k) (type point-array d))
  (let* ((rad (/ (* +pi+ ang) 180d0))
         (gam (* r (cos rad)))
         (del (* r (sin rad)))
         (alf 0d0) (bet 0d0)
         (tmp (make-point)))
    (declare (type double-float rad gam del alf bet) (dynamic-extent tmp))
    (translate-sum-2-pts o gam x del y tmp)
    (setf (aref d 0 0) (aref tmp 0) (aref d 0 1) (aref tmp 1) (aref d 0 2) (aref tmp 2))
    (do ((i 1 (1+ i))) ((> i k) d)
      (declare (type fixnum i))
      (case (mod i 4)
        (1 (setf alf (- del) bet gam))
        (2 (setf alf (- gam) bet (- del)))
        (3 (setf alf del bet (- gam)))
        (0 (setf alf gam bet del)))
      (combine-2-pts alf x bet y tmp)
      (setf (aref d i 0) (aref tmp 0) (aref d i 1) (aref tmp 1) (aref d i 2) (aref tmp 2)))))

;;; ------------------------------------------------------------ level 1, group 4

(defun vector-cross-ref (a b c)
  "C = A x B, for vectors C passes by reference.
NLib: N_VectorCrossRef."
  (declare (type point a b c))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (set-point c (- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))))))

(defun vector-directed-angle (r v)
  "The directed angle in degrees, 0 to 360, from reference vector R to vector
V in 2D: counterclockwise angles beyond 180 come out as 360 minus the
unsigned angle.  Signals :geo-err when the cosine cannot be formed (a zero
vector).
NLib: N_VectorDirectedAngle."
  (declare (type point r v))
  (let ((dot (vector-dot r v))
        (mr (vector-magnitude r))
        (mv (vector-magnitude v)))
    (declare (type double-float dot mr mv))
    (when (float-op-is-bad dot (* mr mv) +division+)
      (nurbs-error :geo-err 'vector-directed-angle))
    (let ((arg (/ dot (* mr mv))))
      (declare (type double-float arg))
      (when (> arg 1d0) (setf arg 1d0))
      (when (< arg -1d0) (setf arg -1d0))
      (let* ((rad (the double-float (acos arg)))
             (ang (* 180d0 (/ rad +pi+)))
             (n (make-point)))
        (declare (type double-float rad ang) (dynamic-extent n))
        (vector-perpendicular r n)
        (setf dot (vector-dot v n))
        (when (< dot 0d0)
          (setf ang (- 360d0 ang)))
        ang))))

(defun isect-line-line (ls1 ls2 p)
  "Intersect lines LS1 and LS2 in the least-squares sense: write the midpoint of
the closest points into P.  Returns (values flg t1 t2), FLG +true+ when the
lines meet within their bounds and +false+ when they are parallel or a bounded
line misses; T1 and T2 (the parameters on each line) and P are written only
when FLG is +true+, and T1 and T2 are 0 otherwise.  Signals :geo-err for a
zero direction and :num-err when a division fails.
NLib: N_IsectLineLine."
  (declare (type lineseg ls1 ls2) (type point p))
  (let ((p21 (make-point)) (v1 (make-point)) (v2 (make-point))
        (p1 (make-point)) (p2 (make-point)))
    (declare (dynamic-extent p21 v1 v2 p1 p2))
    (let* ((bd1 (line-get-data ls1 p1 v1))
           (bd2 (line-get-data ls2 p2 v2))
           (m1 (vector-normalize v1 v1))
           (m2 (vector-normalize v2 v2)))
      (declare (type fixnum bd1 bd2) (type double-float m1 m2))
      (vector-dir p1 p2 p21)
      (let* ((a00 (vector-dot v1 v1))
             (a01 (vector-dot v1 v2))
             (a11 (vector-dot v2 v2))
             (b0 (vector-dot v1 p21))
             (b1 (vector-dot v2 p21))
             (den (- (* a00 a11) (* a01 a01))))
        (declare (type double-float a00 a01 a11 b0 b1 den))
        (when (< (abs den) (the double-float *zctl*))
          (return-from isect-line-line (values +false+ 0d0 0d0)))
        (let ((numa (- (* b0 a11) (* b1 a01)))
              (numb (- (* b0 a01) (* b1 a00))))
          (declare (type double-float numa numb))
          (when (or (float-op-is-bad numa den +division+) (float-op-is-bad numb den +division+))
            (nurbs-error :num-err 'isect-line-line))
          (let ((alpha (/ numa den))
                (beta (/ numb den)))
            (declare (type double-float alpha beta))
            (when (or (float-op-is-bad alpha m1 +division+) (float-op-is-bad beta m2 +division+))
              (nurbs-error :num-err 'isect-line-line))
            (let ((u1 (/ alpha m1))
                  (u2 (/ beta m2)))
              (declare (type double-float u1 u2))
              (when (and (= bd1 +true+) (or (< u1 0d0) (> u1 1d0)))
                (return-from isect-line-line (values +false+ 0d0 0d0)))
              (when (and (= bd2 +true+) (or (< u2 0d0) (> u2 1d0)))
                (return-from isect-line-line (values +false+ 0d0 0d0)))
              (vector-pt-along-vector p1 alpha v1 p1)
              (vector-pt-along-vector p2 beta v2 p2)
              (combine-2-pts 0.5d0 p1 0.5d0 p2 p)
              (values +true+ u1 u2))))))))

(defun isect-line-plane (lsg pln p)
  "Intersect line LSG with plane PLN, writing the intersection into P.  Returns
(values flg t): FLG +true+ when they meet, +false+ when they are parallel or a
bounded line misses; T is the line parameter, written (as the C writes it)
whenever they are not parallel and 0 when they are; P is written only when
FLG is +true+.  Signals :geo-err for a zero direction or normal and :num-err
when a division fails.
NLib: N_IsectLinePlane."
  (declare (type lineseg lsg) (type plane pln) (type point p))
  (let ((pl (make-point)) (v (make-point)) (n (make-point))
        (q (make-point)) (r (make-point)))
    (declare (dynamic-extent pl v n q r))
    (let ((bnd (line-get-data lsg q v)))
      (declare (type fixnum bnd))
      (plane-get-data pln r n)
      (let ((mv (vector-normalize v v)))
        (declare (type double-float mv))
        (vector-normalize n n)
        (vector-dir q r pl)
        (let ((den (vector-dot v n)))
          (declare (type double-float den))
          (when (< (abs den) (the double-float *ltol*))
            (return-from isect-line-plane (values +false+ 0d0)))
          (let ((num (vector-dot pl n)))
            (declare (type double-float num))
            (when (float-op-is-bad num den +division+)
              (nurbs-error :num-err 'isect-line-plane))
            (let ((alpha (/ num den)))
              (declare (type double-float alpha))
              (when (float-op-is-bad alpha mv +division+)
                (nurbs-error :num-err 'isect-line-plane))
              (let ((tt (/ alpha mv)))
                (declare (type double-float tt))
                (when (and (= bnd +true+) (or (< tt 0d0) (> tt 1d0)))
                  (return-from isect-line-plane (values +false+ tt)))
                (vector-pt-along-vector q alpha v p)
                (values +true+ tt)))))))))

(defun dist-c-polygon (pw n)
  "The length of the control polygon PW[0..N], summing the Euclidean distances
of consecutive control points.
NLib: N_DistCPolygon."
  (declare (type cpoint-array pw) (type fixnum n))
  (let ((len 0d0) (a (make-cpoint)) (b (make-cpoint)))
    (declare (type double-float len) (dynamic-extent a b))
    (do ((i 1 (1+ i))) ((> i n) len)
      (declare (type fixnum i))
      (%row-c-pt pw (1- i) a)
      (%row-c-pt pw i b)
      (incf len (dist-cpt-cpt a b)))))

(defun dist-signed-pt-plane (a b c d p)
  "The signed distance of point P from the plane a*x + b*y + c*z + d = 0, (A B
C) its unit normal.
NLib: N_DistSignedPtPlane."
  (declare (type double-float a b c d) (type point p))
  (multiple-value-bind (x y z) (pt-to-xyz p)
    (declare (type double-float x y z))
    (+ (+ (+ (* a x) (* b y)) (* c z)) d)))

(defun b-box-get-diagonal (box)
  "The length of the diagonal of min-max box BOX.
NLib: N_BBoxGetDiagonal."
  (declare (type real-vector box))
  (let ((pll (make-point)) (pur (make-point)))
    (declare (dynamic-extent pll pur))
    (multiple-value-bind (xl xr yb yt zn zf) (get-b-box-data box)
      (declare (type double-float xl xr yb yt zn zf))
      (pt-from-xyz xl yb zn pll)
      (pt-from-xyz xr yt zf pur)
      (dist-pt-pt pll pur))))

(defun create-line-pt-pt (lsg ps pe flg)
  "Make LSG the line from PS to PE (direction PE - PS) with bounding flag FLG
(+bounded+ or +unbounded+).  Returns LSG.
NLib: N_CreateLinePtPt."
  (declare (type lineseg lsg) (type point ps pe) (type fixnum flg))
  (let ((v (make-point)) (s (make-point)))
    (declare (dynamic-extent v s))
    (vector-copy ps s)
    (vector-dir s pe v)
    (create-line-start-dir-vector lsg s v flg)))

(defun plane-pt-normal-to-implicit (p n)
  "The normalized implicit equation a*x + b*y + c*z + d = 0 of the plane
through P with normal N, as (values a b c d).  Signals :geo-err for a zero
normal.
NLib: N_PlanePtNormalToImplicit."
  (declare (type point p n))
  (let ((nn (make-point)))
    (declare (dynamic-extent nn))
    (vector-copy n nn)
    (vector-normalize-ref nn)
    (multiple-value-bind (a b c) (pt-to-xyz nn)
      (declare (type double-float a b c))
      (let ((dot (vector-dot p nn)))
        (declare (type double-float dot))
        (values a b c (- dot))))))

(defun line-implicit-to-pt-vector (a b c p v)
  "Convert the 2D line a*x + b*y + c = 0 to a point P on it and its unit
direction V, both written.  Returns P.  Signals :geo-err when A and B are
both zero and :num-err when the division fails.
NLib: N_LineImplicitToPtVector."
  (declare (type double-float a b c) (type point p v))
  (pt-from-xyz b (- a) 0d0 v)
  (vector-normalize-ref v)
  (let ((x 0d0) (y 0d0))
    (declare (type double-float x y))
    (if (> (abs a) (abs b))
        (let ((num (- (+ (* b b) c))))
          (declare (type double-float num))
          (when (float-op-is-bad num a +division+)
            (nurbs-error :num-err 'line-implicit-to-pt-vector))
          (setf x (/ num a) y b))
        (let ((num (- (+ (* a a) c))))
          (declare (type double-float num))
          (when (float-op-is-bad num b +division+)
            (nurbs-error :num-err 'line-implicit-to-pt-vector))
          (setf x a y (/ num b))))
    (pt-from-xyz x y 0d0 p)))

(defun pts-are-colinear (a b c tol)
  "True when points A, B and C are collinear within TOL: two of them closer
than TOL, or the third point's perpendicular distance from the longest side
below TOL (a point projecting outside that side counts as collinear, as in
the C).  A failed distance computation counts as not collinear.
NLib: N_PtsAreColinear."
  (declare (type point a b c) (type double-float tol))
  (let ((dab (dist-sq-pt-pt a b))
        (dbc (dist-sq-pt-pt b c))
        (dca (dist-sq-pt-pt c a))
        (tol2 (* tol tol)))
    (declare (type double-float dab dbc dca tol2))
    (flet ((perp-below-tol (q s e)
             (let ((d (handler-case (dist-perp-pt-line-seg q s e)
                        (nurbs-error () (return-from pts-are-colinear nil)))))
               (declare (type double-float d))
               (< d tol))))
      (cond ((or (< dab tol2) (< dbc tol2) (< dca tol2)) t)
            ((and (> dab dbc) (> dab dca)) (perp-below-tol c a b))
            ((and (> dbc dca) (> dbc dab)) (perp-below-tol a b c))
            ((and (> dca dab) (> dca dbc)) (perp-below-tol b c a))
            (t nil)))))
