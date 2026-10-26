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
  (c-sqrt (vector-dot a a)))

(defun vector-magnitude-ref (a)
  "|A|.  NLib: N_VectorMagnitudeRef."
  (declare (type point a))
  (c-sqrt (vector-dot a a)))

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
        (c-sqrt (+ (+ (* dx dx) (* dy dy)) (* dz dz)))))))

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
        (c-sqrt disc)))))

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
        (setf len (c-sqrt len))
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
        (let ((den (c-sqrt dotab)))
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

;;; ------------------------------------------------------------ level 2, group 6

(defun isect-plane-plane (pl1 pl2 p v)
  "Intersect planes PL1 and PL2.  Returns +true+ and writes a point of the
intersection line into P and its direction (the cross product of the unit
normals, not normalized) into V; returns +false+, writing neither, when the
planes are parallel.  Signals :geo-err for a zero normal.
NLib: N_IsectPlanePlane."
  (declare (type plane pl1 pl2) (type point p v))
  (let ((p1 (make-point)) (p2 (make-point)) (n1 (make-point)) (n2 (make-point)))
    (declare (dynamic-extent p1 p2 n1 n2))
    (plane-get-data pl1 p1 n1)
    (plane-get-data pl2 p2 n2)
    (multiple-value-bind (a1 b1 c1 d1) (plane-pt-normal-to-implicit p1 n1)
      (declare (type double-float a1 b1 c1 d1))
      (multiple-value-bind (a2 b2 c2 d2) (plane-pt-normal-to-implicit p2 n2)
        (declare (type double-float a2 b2 c2 d2))
        (let* ((f (- (* b1 c2) (* b2 c1)))
               (g (- (* c1 a2) (* c2 a1)))
               (h (- (* a1 b2) (* a2 b1)))
               (det (+ (+ (* f f) (* g g)) (* h h))))
          (declare (type double-float f g h det))
          (when (< (abs det) (the double-float *zctl*))
            (return-from isect-plane-plane +false+))
          (let* ((fac (/ 1d0 det))
                 (dc (- (* d1 c2) (* c1 d2)))
                 (db (- (* d1 b2) (* b1 d2)))
                 (ad (- (* a1 d2) (* a2 d1)))
                 (x0 (* fac (- (* g dc) (* h db))))
                 (y0 (* (- fac) (+ (* f dc) (* h ad))))
                 (z0 (* fac (+ (* f db) (* g ad)))))
            (declare (type double-float fac dc db ad x0 y0 z0))
            (pt-from-xyz x0 y0 z0 p)
            (pt-from-xyz f g h v)
            +true+))))))

(defun project-pt-line (lsg p q)
  "Project point P onto line LSG, through the plane through P normal to the
line.  Writes the projection into Q and returns (values flg t), as
isect-line-plane does: FLG +false+ when the projection falls outside a bounded
line (Q then unwritten).  Signals what isect-line-plane signals.
NLib: N_ProjectPtLine."
  (declare (type lineseg lsg) (type point p q))
  (let ((nor (%make-plane (make-point) (make-point))))
    (declare (dynamic-extent nor))
    (create-plane-pt-normal nor p (lineseg-v lsg))
    (isect-line-plane lsg nor q)))

(defun project-pt-plane (pln p q)
  "Project point P onto plane PLN along its normal, writing the projection into
Q, which is left unwritten in the (numerically degenerate) case that the line
comes out parallel to the plane.  Returns Q.  Signals what isect-line-plane
signals (:geo-err for a zero normal).
NLib: N_ProjectPtPlane."
  (declare (type plane pln) (type point p q))
  (let ((n (make-point))
        (lsg (%make-lineseg (make-point) (make-point) +unbounded+)))
    (declare (dynamic-extent n lsg))
    (plane-get-normal pln n)
    (create-line-start-dir-vector lsg p n +unbounded+)
    (isect-line-plane lsg pln q)
    q))

(defun pts-1d-get-max-extent (p n)
  "The maximum extent of points P[0..N]: the diagonal of their min-max box.
NLib: N_Pts1dGetMaxExtent."
  (declare (type point-array p) (type fixnum n))
  (let ((box (make-real-vector 6)))
    (declare (dynamic-extent box))
    (polygon-get-b-box (e-polygon-from-pts n p) box)
    (b-box-get-diagonal box)))

(defun pts-2d-get-max-extent (p n m)
  "The maximum extent of points P[0..N][0..M]: the diagonal of their min-max
box.
NLib: N_Pts2dGetMaxExtent."
  (declare (type point-net p) (type fixnum n m))
  (let ((box (make-real-vector 6)))
    (declare (dynamic-extent box))
    (e-net-get-b-box (e-net-from-pts n m p) box)
    (b-box-get-diagonal box)))

(defun create-transform-matrix-from-axes (o1 x1 y1 z1 o2 x2 y2 z2 rma)
  "Write into RMA the 4x4 matrix (highest indexes 3, 3) mapping the orthonormal
frame O1, <X1 Y1 Z1> onto the frame O2, <X2 Y2 Z2>, applied from the left:
T(O2) * R2 * R1^T * T(-O1), multiplied in that order.  A NULL RMA
(make-rmatrix) gets storage; an existing one must be at least 3 x 3
(:sto-err).  Returns RMA.
NLib: N_CreateTransformMatrixFromAxes."
  (declare (type point o1 x1 y1 z1 o2 x2 y2 z2) (type rmatrix rma))
  (%check-mem-real-matrix rma 3 3 +mt-full+ 3 'create-transform-matrix-from-axes)
  (let ((rm1 (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3))
        (rmi (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3))
        (rt1 (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3))
        (rt2 (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3))
        (rmw (set-real-matrix (make-rmatrix) 3 3 +mt-full+ 3)))
    (let ((r1 (rmatrix-rm rm1)) (ri (rmatrix-rm rmi))
          (t1 (rmatrix-rm rt1)) (t2 (rmatrix-rm rt2)))
      (declare (type real-matrix r1 ri t1 t2))
      ;; X1, Y1, Z1 as the rows of R1 (the inverse rotation)
      (multiple-value-bind (x y z) (pt-to-xyz x1)
        (setf (aref r1 0 0) x (aref r1 0 1) y (aref r1 0 2) z))
      (setf (aref r1 0 3) 0d0)
      (multiple-value-bind (x y z) (pt-to-xyz y1)
        (setf (aref r1 1 0) x (aref r1 1 1) y (aref r1 1 2) z))
      (setf (aref r1 1 3) 0d0)
      (multiple-value-bind (x y z) (pt-to-xyz z1)
        (setf (aref r1 2 0) x (aref r1 2 1) y (aref r1 2 2) z))
      (setf (aref r1 2 3) 0d0
            (aref r1 3 0) 0d0 (aref r1 3 1) 0d0 (aref r1 3 2) 0d0 (aref r1 3 3) 1d0)
      ;; X2, Y2, Z2 as the columns of RI (the forward rotation)
      (multiple-value-bind (x y z) (pt-to-xyz x2)
        (setf (aref ri 0 0) x (aref ri 1 0) y (aref ri 2 0) z))
      (setf (aref ri 0 3) 0d0)
      (multiple-value-bind (x y z) (pt-to-xyz y2)
        (setf (aref ri 0 1) x (aref ri 1 1) y (aref ri 2 1) z))
      (setf (aref ri 1 3) 0d0)
      (multiple-value-bind (x y z) (pt-to-xyz z2)
        (setf (aref ri 0 2) x (aref ri 1 2) y (aref ri 2 2) z))
      (setf (aref ri 2 3) 0d0
            (aref ri 3 0) 0d0 (aref ri 3 1) 0d0 (aref ri 3 2) 0d0 (aref ri 3 3) 1d0)
      ;; -O1 into T1 (the inverse translation)
      (multiple-value-bind (x y z) (pt-to-xyz o1)
        (declare (type double-float x y z))
        (setf (aref t1 0 0) 1d0 (aref t1 0 1) 0d0 (aref t1 0 2) 0d0 (aref t1 0 3) (- x)
              (aref t1 1 0) 0d0 (aref t1 1 1) 1d0 (aref t1 1 2) 0d0 (aref t1 1 3) (- y)
              (aref t1 2 0) 0d0 (aref t1 2 1) 0d0 (aref t1 2 2) 1d0 (aref t1 2 3) (- z)
              (aref t1 3 0) 0d0 (aref t1 3 1) 0d0 (aref t1 3 2) 0d0 (aref t1 3 3) 1d0))
      ;; +O2 into T2 (the forward translation)
      (multiple-value-bind (x y z) (pt-to-xyz o2)
        (declare (type double-float x y z))
        (setf (aref t2 0 0) 1d0 (aref t2 0 1) 0d0 (aref t2 0 2) 0d0 (aref t2 0 3) x
              (aref t2 1 0) 0d0 (aref t2 1 1) 1d0 (aref t2 1 2) 0d0 (aref t2 1 3) y
              (aref t2 2 0) 0d0 (aref t2 2 1) 0d0 (aref t2 2 2) 1d0 (aref t2 2 3) z
              (aref t2 3 0) 0d0 (aref t2 3 1) 0d0 (aref t2 3 2) 0d0 (aref t2 3 3) 1d0)))
    ;; rmw = rt2 * rmi; rt2 = rmw * rm1; rma = rt2 * rt1
    (real-matrix-multiply rt2 rmi rmw)
    (real-matrix-multiply rmw rm1 rt2)
    (real-matrix-multiply rt2 rt1 rma)
    rma))

(defun line-fit-2d-pts (p n)
  "The least-squares line through the 2D points P[0..N] (z ignored), as the
coefficients (values a b c) of its normalized equation a*x + b*y + c = 0, (a
b) the unit normal: the line passes through the centroid along the singular
vector of the scatter matrix with the smallest singular value.  As in the C,
a failure of the decomposition is ignored.  Signals :num-err for N below 0.
NLib: N_LineFit2dPts."
  (declare (type point-array p) (type fixnum n))
  (let* ((ld 1)
         (y (make-point-array (1+ n)))
         (mw (make-real-vector (1+ ld)))
         (mu (make-real-matrix (1+ ld) (1+ ld)))
         (mv (make-real-matrix (1+ ld) (1+ ld)))
         (qx 0d0) (qy 0d0)
         (cnt (+ (float n 1d0) 1d0)))
    (declare (type fixnum ld) (type point-array y) (type real-vector mw)
             (type real-matrix mu mv) (type double-float qx qy cnt))
    (when (float-op-is-bad 1d0 cnt +division+)
      (nurbs-error :num-err 'line-fit-2d-pts))
    ;; Q, the centroid
    (fast
      (dotimes (i (1+ n))
        (setf qx (+ qx (aref p i 0))
              qy (+ qy (aref p i 1)))))
    (setf qx (/ qx cnt)
          qy (/ qy cnt))
    ;; Y[i] = P[i] - Q
    (fast
      (dotimes (i (1+ n))
        (setf (aref y i 0) (- (aref p i 0) qx)
              (aref y i 1) (- (aref p i 1) qy))))
    ;; M = Sum_i (Yi.Yi I - Yi Yi^T)
    (fast
      (dotimes (i (1+ n))
        (let* ((yx (aref y i 0)) (yy (aref y i 1))
               (ydoty (+ (* yx yx) (* yy yy))))
          (declare (type double-float yx yy ydoty))
          (setf (aref mu 0 0) (+ (aref mu 0 0) (- ydoty (* yx yx)))
                (aref mu 0 1) (+ (aref mu 0 1) (- (* yx yy)))
                (aref mu 1 0) (+ (aref mu 1 0) (- (* yy yx)))
                (aref mu 1 1) (+ (aref mu 1 1) (- ydoty (* yy yy)))))))
    ;; the C ignores the decomposition's error
    (handler-case (real-matrix-single-value-decompose mu mw mv ld ld)
      (nurbs-error () nil))
    ;; the singular vector of the smallest singular value
    (let ((imin 0) (sval (aref mw 0)))
      (declare (type fixnum imin) (type double-float sval))
      (do ((i 1 (1+ i))) ((> i ld))
        (declare (type fixnum i))
        (when (> sval (aref mw i))
          (setf imin i sval (aref mw i))))
      (let* ((a (- (aref mu 1 imin)))
             (b (aref mu 0 imin))
             (c (- (+ (* a qx) (* b qy)))))
        (declare (type double-float a b c))
        (values a b c)))))

;;; ------------------------------------------------------------ level 3, group 4

(defun dist-pt-inf-line (lsg p)
  "The distance of point P from the infinite line through LSG (its bounding
flag is ignored): the distance to the projection of P on the line, or to the
origin when the line comes out parallel to its normal plane.  Signals what
project-pt-line signals.
NLib: N_DistPtInfLine."
  (declare (type lineseg lsg) (type point p))
  (let ((q (make-point))
        (lin (copy-structure lsg)))
    (declare (dynamic-extent q))
    (line-set-bounding-flag lin +unbounded+)
    (project-pt-line lin p q)
    (dist-pt-pt p q)))

(defconstant +geom--cpt-reflect-noz-tol+ 1.0d-6
  "NOZ_TOL in N_CptReflect.")

(defun cpt-reflect (pw pln qw)
  "Write into QW the reflection of control point PW through plane PLN: the
Euclidean point is reflected and weighted with PW's weight.  A 2D PW stays 2D
when the plane's unit normal has (nearly) no z component.  Signals :geo-err
for a zero normal.  Returns QW.
NLib: N_CptReflect."
  (declare (type cpoint pw qw) (type plane pln))
  (let ((p (make-point)) (q (make-point)) (r (make-point)))
    (declare (dynamic-extent p q r))
    (c-pt-to-pt-euclid pw p)
    (project-pt-plane pln p q)
    (vector-dir p q r)
    (sum-2-pts q r q)
    (multiple-value-bind (x y z) (pt-to-xyz q)
      (declare (type double-float x y z))
      (multiple-value-bind (dd1 dd2 wz w) (c-pt-to-wx-wy-wz pw)
        (declare (type double-float wz w) (ignore dd1 dd2))
        (when (/= w +no-w+)
          (setf x (* x w)
                y (* y w)
                z (* z w)))
        (plane-get-normal pln r)
        (vector-normalize r r)
        (let ((vz (aref r 2)))
          (declare (type double-float vz))
          (when (and (= wz +no-z+) (<= (abs vz) +geom--cpt-reflect-noz-tol+))
            (setf z +no-z+)))
        (c-pt-from-wx-wy-wz x y z w qw)))))

(defun line-fit-3d-pts (p n q v)
  "The least-squares line through the points P[0..N]: writes the centroid into
Q and into V the unit direction minimizing the sum of the squared distances,
the singular vector of the smallest singular value of the matrix
Sum_i (Yi.Yi I - Yi Yi^T), Yi = P[i] - Q.  As in the C, a failure of the
decomposition is ignored.  Signals :num-err for N below 0 (Q zeroed first).
Returns Q.
NLib: N_LineFit3dPts."
  (declare (type point-array p) (type fixnum n) (type point q v))
  (let* ((ld 2)
         (y (make-point-array (1+ (max n 0))))
         (mw (make-real-vector (1+ ld)))
         (mu (make-real-matrix (1+ ld) (1+ ld)))
         (mv (make-real-matrix (1+ ld) (1+ ld)))
         (qx 0d0) (qy 0d0) (qz 0d0)
         (cnt (+ (float n 1d0) 1d0)))
    (declare (type fixnum ld) (type point-array y) (type real-vector mw)
             (type real-matrix mu mv) (type double-float qx qy qz cnt))
    ;; Q, the centroid
    (set-point q 0d0 0d0 0d0)
    (when (float-op-is-bad 1d0 cnt +division+)
      (nurbs-error :num-err 'line-fit-3d-pts))
    (fast
      (dotimes (i (1+ n))
        (setf qx (+ (aref p i 0) qx)
              qy (+ (aref p i 1) qy)
              qz (+ (aref p i 2) qz))))
    (let ((fac (/ 1d0 cnt)))
      (declare (type double-float fac))
      (setf qx (* fac qx)
            qy (* fac qy)
            qz (* fac qz)))
    (set-point q qx qy qz)
    ;; Y[i] = P[i] - Q
    (fast
      (dotimes (i (1+ n))
        (setf (aref y i 0) (- (aref p i 0) qx)
              (aref y i 1) (- (aref p i 1) qy)
              (aref y i 2) (- (aref p i 2) qz))))
    ;; M = Sum_i (Yi.Yi I - Yi Yi^T)
    (fast
      (dotimes (i (1+ n))
        (let* ((yx (aref y i 0)) (yy (aref y i 1)) (yz (aref y i 2))
               (ydoty (+ (+ (* yx yx) (* yy yy)) (* yz yz))))
          (declare (type double-float yx yy yz ydoty))
          (setf (aref mu 0 0) (+ (aref mu 0 0) (- ydoty (* yx yx)))
                (aref mu 0 1) (+ (aref mu 0 1) (- (* yx yy)))
                (aref mu 0 2) (+ (aref mu 0 2) (- (* yx yz)))
                (aref mu 1 0) (+ (aref mu 1 0) (- (* yy yx)))
                (aref mu 1 1) (+ (aref mu 1 1) (- ydoty (* yy yy)))
                (aref mu 1 2) (+ (aref mu 1 2) (- (* yy yz)))
                (aref mu 2 0) (+ (aref mu 2 0) (- (* yz yx)))
                (aref mu 2 1) (+ (aref mu 2 1) (- (* yz yy)))
                (aref mu 2 2) (+ (aref mu 2 2) (- ydoty (* yz yz)))))))
    ;; the C ignores the decomposition's error
    (handler-case (real-matrix-single-value-decompose mu mw mv ld ld)
      (nurbs-error () nil))
    ;; the singular vector of the smallest singular value
    (let ((imin 0) (sval (aref mw 0)))
      (declare (type fixnum imin) (type double-float sval))
      (do ((i 1 (1+ i))) ((> i ld))
        (declare (type fixnum i))
        (when (> sval (aref mw i))
          (setf imin i sval (aref mw i))))
      (set-point v (aref mu 0 imin) (aref mu 1 imin) (aref mu 2 imin)))
    q))

(defun calc-circ-center-and-radius (a b c d)
  "The circle through the points A, B and C, in 2D or 3D: writes its center
into D and returns (values r cfl), R the radius and CFL +true+; CFL is +false+
(D unwritten, R 0) when the points are collinear or two coincide, or when a
pair of perpendicular bisectors does not meet.  The center is the average of
the pairwise intersections of the bisectors, the radius the average distance.
Signals what isect-plane-plane and isect-line-line signal.
NLib: N_CalcCircCenterAndRadius."
  (declare (type point a b c d))
  (when (pts-are-colinear a b c *mtol*)
    (return-from calc-circ-center-and-radius (values 0d0 +false+)))
  (let ((za (dist-pt-pt a b))
        (zb (dist-pt-pt b c))
        (zc (dist-pt-pt c a)))
    (declare (type double-float za zb zc))
    (when (or (< za *mtol*) (< zb *mtol*) (< zc *mtol*))
      (return-from calc-circ-center-and-radius (values 0d0 +false+)))
    (let ((fac (/ 1d0 3d0))
          (mab (make-point)) (mbc (make-point)) (mca (make-point))
          (vab (make-point)) (vbc (make-point)) (vca (make-point))
          (lab (%make-lineseg (make-point) (make-point) +unbounded+))
          (lbc (%make-lineseg (make-point) (make-point) +unbounded+))
          (lca (%make-lineseg (make-point) (make-point) +unbounded+))
          (pa (make-point)) (pb (make-point)) (pc (make-point)) (cc (make-point)))
      (declare (type double-float fac) (dynamic-extent mab mbc mca vab vbc vca pa pb pc cc))
      ;; some points and vectors
      (vector-combine 0.5d0 a 0.5d0 b mab)
      (vector-diff b a vab)
      (vector-combine 0.5d0 b 0.5d0 c mbc)
      (vector-diff c b vbc)
      (vector-combine 0.5d0 c 0.5d0 a mca)
      (vector-diff a c vca)
      ;; perpendicular bisectors
      (setf za (aref a 2) zb (aref b 2) zc (aref c 2))
      (if (and (< (abs za) *mtol*) (< (abs zb) *mtol*) (< (abs zc) *mtol*))
          ;; the points are in the <x,y> plane
          (let ((wab (make-point)) (wbc (make-point)) (wca (make-point)))
            (declare (dynamic-extent wab wbc wca))
            (vector-perpendicular vab wab)
            (vector-perpendicular vbc wbc)
            (vector-perpendicular vca wca)
            (create-line-start-dir-vector lab mab wab +unbounded+)
            (create-line-start-dir-vector lbc mbc wbc +unbounded+)
            (create-line-start-dir-vector lca mca wca +unbounded+))
          ;; the points are not in the <x,y> plane
          (let ((na (make-point)) (nb (make-point)) (nc (make-point))
                (nv (make-point)) (g (make-point)) (pp (make-point)) (vv (make-point))
                (pln (%make-plane (make-point) (make-point)))
                (pab (%make-plane (make-point) (make-point)))
                (pbc (%make-plane (make-point) (make-point)))
                (pca (%make-plane (make-point) (make-point))))
            (declare (dynamic-extent na nb nc nv g pp vv))
            (vector-cross vab vbc nb)
            (vector-cross vbc vca nc)
            (vector-cross vca vab na)
            (vector-create 0d0 0d0 0d0 nv)
            (vector-create 0d0 0d0 0d0 g)
            (vector-sum a b g)
            (vector-sum c g g)
            (vector-sum na nb nv)
            (vector-sum nc nv nv)
            (vector-scale g fac g)
            (vector-scale nv fac nv)
            (create-plane-pt-normal pln g nv)
            (create-plane-pt-normal pab mab vab)
            (create-plane-pt-normal pbc mbc vbc)
            (create-plane-pt-normal pca mca vca)
            (flet ((bisector (pl lsg)
                     (when (= (isect-plane-plane pln pl pp vv) +false+)
                       (return-from calc-circ-center-and-radius (values 0d0 +false+)))
                     (create-line-start-dir-vector lsg pp vv +unbounded+)))
              (bisector pab lab)
              (bisector pbc lbc)
              (bisector pca lca))))
      ;; intersect the three bisectors
      (flet ((meet (l1 l2 pt)
               (when (= (the fixnum (isect-line-line l1 l2 pt)) +false+)
                 (return-from calc-circ-center-and-radius (values 0d0 +false+)))))
        (meet lab lbc pb)
        (meet lbc lca pc)
        (meet lca lab pa))
      ;; the center and the radius
      (vector-create 0d0 0d0 0d0 cc)
      (vector-sum pa pb cc)
      (vector-sum pc cc cc)
      (vector-scale cc fac cc)
      (let* ((ra (dist-pt-pt cc a))
             (rb (dist-pt-pt cc b))
             (rc (dist-pt-pt cc c))
             (rr (* fac (+ (+ ra rb) rc))))
        (declare (type double-float ra rb rc rr))
        (vector-copy cc d)
        (values rr +true+)))))

;;; ------------------------------------------------------------- L4-g1

(defun line-fit-pts (p n c v)
  "The best fitting line to the points P[0..N], 2D or 3D: writes the centroid
into C and the unit direction into V (line-fit-3d-pts) and returns the
average and the maximum distance of the points from the line.  When P[0] is
2D the z of every point is zeroed for the fit and set to +no-z+ afterwards,
in P itself, and C and V come out 2D.  As in the C, a failure of the fit is
ignored unless there are no points to measure (N below 0), when it is
signaled after all.  Signals what dist-pt-inf-line signals.
NLib: N_LineFitPts."
  (declare (type point-array p) (type fixnum n) (type point c v))
  (let ((is-2d nil)
        (fit-error nil)
        (lsg (%make-lineseg (make-point) (make-point) +unbounded+))
        (pi* (make-point))
        (da 0d0) (dm 0d0))
    (declare (type double-float da dm) (type point pi*) (dynamic-extent pi*))
    ;; 2D points: z temporarily 0
    (when (= (aref p 0 2) +no-z+)
      (setf is-2d t)
      (dotimes (ii (1+ n))
        (setf (aref p ii 2) 0d0)))
    ;; the fit; the C keeps its error flag only when no distance follows
    (handler-case (line-fit-3d-pts p n c v)
      (nurbs-error (e) (setf fit-error e)))
    (when is-2d
      (dotimes (ii (1+ n))
        (setf (aref p ii 2) +no-z+))
      (setf (aref c 2) +no-z+
            (aref v 2) +no-z+))
    ;; the errors
    (create-line-start-dir-vector lsg c v +unbounded+)
    (setf da 0d0
          dm -1d0)
    (dotimes (ii (1+ n))
      (%row-pt p ii pi*)
      (let ((d (dist-pt-inf-line lsg pi*)))
        (declare (type double-float d))
        (when (> d dm)
          (setf dm d))
        (setf da (+ da d))))
    (setf da (/ da (+ (float n 1d0) 1d0)))
    (when (and fit-error (< n 0))
      (error fit-error))
    (values da dm)))

;;; ------------------------------------------------------------- L5-g1

(defun create-circ-from-2d-pts (pp n tol etp)
  "The least-squares circle through the points PP[0..N] (a point-array; x and
y are fitted, the points' z only enters the collinearity check), minimizing
Sum (r^2 - (x-cx)^2 - (y-cy)^2)^2 for the first guess and improved by up to 31
Gauss-Newton steps, stopped when the relative average error changes or the
center moves by less than TOL.  Returns (values cx cy rr era erm cfl): the
center, the radius, the average and maximum absolute errors and CFL +true+.
CFL is +false+ when the points are collinear -- the average (ETP +average+) or
maximum (ETP +maximum+) distance of LINE-FIT-PTS below TOL, ERA and ERM then
being the line's -- or, for N = 2, when CALC-CIRC-CENTER-AND-RADIUS finds no
circle; cx, cy and rr, which the C leaves unwritten then, are returned as 0.
A failing least-squares step ends the iteration with the best circle so far.
Signals :cal-err for another ETP, :num-err for a bad division, and what
line-fit-pts, calc-circ-center-and-radius and the 3 x 3 solver signal.
NLib: N_CreateCircFrom2dPts."
  (declare (type point-array pp) (type fixnum n etp) (type double-float tol))
  (let ((cc (make-point))
        (vv (make-point)))
    (declare (type point cc vv) (dynamic-extent cc vv))
    ;; see if the points are collinear
    (multiple-value-bind (era erm) (line-fit-pts pp n cc vv)
      (declare (type double-float era erm))
      (cond ((= etp +maximum+)
             (when (< erm tol)
               (return-from create-circ-from-2d-pts (values 0d0 0d0 0d0 era erm +false+))))
            ((= etp +average+)
             (when (< era tol)
               (return-from create-circ-from-2d-pts (values 0d0 0d0 0d0 era erm +false+))))
            (t (nurbs-error :cal-err 'create-circ-from-2d-pts)))
      ;; the special case n = 2
      (when (= n 2)
        (let ((p0 (make-point)) (p1 (make-point)) (p2 (make-point)))
          (declare (dynamic-extent p0 p1 p2))
          (%row-pt pp 0 p0)
          (%row-pt pp 1 p1)
          (%row-pt pp 2 p2)
          (multiple-value-bind (rr cfl) (calc-circ-center-and-radius p0 p1 p2 cc)
            (declare (type double-float rr) (type fixnum cfl))
            (when (= cfl +false+)
              (return-from create-circ-from-2d-pts (values 0d0 0d0 0d0 era erm +false+)))
            (multiple-value-bind (cx cy) (pt-to-xyz cc)
              (declare (type double-float cx cy))
              (return-from create-circ-from-2d-pts (values cx cy rr 0d0 0d0 +true+))))))))
  ;; various entities
  (let ((a 0d0) (b 0d0) (c 0d0) (d 0d0) (e 0d0) (f 0d0) (g 0d0) (h 0d0)
        (fac (+ (float n 1d0) 1d0))
        (sol (make-real-vector 3))
        (ind (make-index-vector 3))
        (rma (make-rmatrix))
        (itlm 30))
    (declare (type double-float a b c d e f g h fac) (type real-vector sol)
             (type index-vector ind) (type fixnum itlm))
    (dotimes (i (1+ n))
      (let* ((x (aref pp i 0))
             (y (aref pp i 1))
             (x2 (* x x))
             (y2 (* y y))
             (xy (* x y))
             (x2py2 (+ x2 y2)))
        (declare (type double-float x y x2 y2 xy x2py2))
        (setf a (+ a x)
              b (+ b y)
              c (+ c x2)
              d (+ d y2)
              e (+ e xy)
              g (+ g (* x x2py2))
              h (+ h (* y x2py2)))))
    (setf f (+ c d))
    ;; solve the system of equations
    (set-real-matrix rma 2 2 +mt-full+ 2)
    (let ((rm (rmatrix-rm rma)))
      (declare (type real-matrix rm))
      (setf (aref rm 0 0) (* 2.0d0 c)
            (aref rm 0 1) (* 2.0d0 e)
            (aref rm 0 2) a
            (aref sol 0) g
            (aref rm 1 0) (* 2.0d0 e)
            (aref rm 1 1) (* 2.0d0 d)
            (aref rm 1 2) b
            (aref sol 1) h
            (aref rm 2 0) (* 2.0d0 a)
            (aref rm 2 1) (* 2.0d0 b)
            (aref rm 2 2) fac
            (aref sol 2) f))
    (real-matrix-lu-decompose-pivot rma ind)
    (real-matrix-for-back-pivot rma ind sol)
    ;; the initial center and radius and the error measures
    (let* ((xc (aref sol 0))
           (yc (aref sol 1))
           (rc (c-sqrt (+ (+ (* xc xc) (* yc yc)) (aref sol 2))))
           (eavp 0d0) (eavn 0d0) (eavo 0d0)
           (dm -1.0d0)
           (xo 0d0) (yo 0d0) (ro 0d0) (dao 0d0) (dmo 0d0))
      (declare (type double-float xc yc rc eavp eavn eavo dm xo yo ro dao dmo))
      (dotimes (i (1+ n))
        (let* ((x (aref pp i 0))
               (y (aref pp i 1))
               (ri (c-sqrt (+ (* (- x xc) (- x xc)) (* (- y yc) (- y yc))))))
          (declare (type double-float x y ri))
          (setf eavp (+ eavp (abs (- rc ri))))
          (when (> (abs (- rc ri)) dm)
            (setf dm (abs (- rc ri))))))
      (when (float-op-is-bad eavp fac +division+)
        (nurbs-error :num-err 'create-circ-from-2d-pts))
      (setf eavp (/ eavp fac))
      (setf xo xc
            yo yc
            ro rc
            dao eavp
            dmo dm)
      (when (< eavp (the double-float *mtol*))
        (return-from create-circ-from-2d-pts (values xo yo ro eavp dm +true+)))
      (when (float-op-is-bad eavp rc +division+)
        (nurbs-error :num-err 'create-circ-from-2d-pts))
      (setf eavp (/ eavp rc))
      (setf eavo eavp)
      ;; improve the fit by Gauss-Newton iteration: the least-squares problem
      (let ((rmj (make-rmatrix))
            (rhs (make-real-vector (1+ n)))
            (itc 0))
        (declare (type real-vector rhs) (type fixnum itc))
        (set-real-matrix rmj n 2 +mt-full+ 2)
        (let ((jm (rmatrix-rm rmj)))
          (declare (type real-matrix jm))
          (do () ((> itc itlm))
            ;; the Jacobian and the right hand side
            (dotimes (i (1+ n))
              (let* ((x (aref pp i 0))
                     (y (aref pp i 1))
                     (ri (c-sqrt (+ (* (- x xc) (- x xc)) (* (- y yc) (- y yc))))))
                (declare (type double-float x y ri))
                (when (float-op-is-bad 1.0d0 ri +division+)
                  (nurbs-error :num-err 'create-circ-from-2d-pts))
                (let ((inv (/ 1.0d0 ri)))
                  (declare (type double-float inv))
                  (setf (aref jm i 0) (* (- xc x) inv)
                        (aref jm i 1) (* (- yc y) inv)
                        (aref jm i 2) -1.0d0
                        (aref rhs i) (- rc ri)))))
            ;; solve by least squares; a failure ends the iteration
            (unless (handler-case (progn (real-matrix-lst-sq-solve rmj rhs sol) t)
                      (nurbs-error () (err-clear) nil))
              (return))
            (setf xc (+ xc (aref sol 0))
                  yc (+ yc (aref sol 1))
                  rc (+ rc (aref sol 2)))
            ;; recompute the average error
            (setf eavn 0d0
                  dm -1.0d0)
            (dotimes (i (1+ n))
              (let* ((x (aref pp i 0))
                     (y (aref pp i 1))
                     (ri (c-sqrt (+ (* (- x xc) (- x xc)) (* (- y yc) (- y yc))))))
                (declare (type double-float x y ri))
                (setf eavn (+ eavn (abs (- rc ri))))
                (when (> (abs (- rc ri)) dm)
                  (setf dm (abs (- rc ri))))))
            (when (float-op-is-bad eavn fac +division+)
              (nurbs-error :num-err 'create-circ-from-2d-pts))
            (setf eavn (/ eavn fac))
            (let ((da eavn)
                  (dis 0d0))
              (declare (type double-float da dis))
              (when (float-op-is-bad eavn rc +division+)
                (nurbs-error :num-err 'create-circ-from-2d-pts))
              (setf eavn (/ eavn rc))
              ;; check for convergence and save the best result
              (setf dis (c-sqrt (+ (+ (* (aref sol 0) (aref sol 0)) (* (aref sol 1) (aref sol 1)))
                                         (* (aref sol 2) (aref sol 2)))))
              (when (float-op-is-bad dis rc +division+)
                (nurbs-error :num-err 'create-circ-from-2d-pts))
              (setf dis (/ dis rc))
              ;; the best result has the smallest average error
              (when (< eavn eavo)
                (setf xo xc
                      yo yc
                      ro rc
                      dao da
                      dmo dm
                      eavo eavn))
              ;; out when the average error does not change, or the center
              ;; does not move
              (when (< (abs (- eavp eavn)) tol)
                (return))
              (when (< dis tol)
                (return)))
            (setf eavp eavn)
            (incf itc))))
      ;; the output
      (values xo yo ro dao dmo +true+))))

;;; ------------------------------------------------------------- L10-g1

(defun dist-pt-arc (p c r alf bet ps pe)
  "The distance between point P and the circular arc of center C, radius R,
start and end angles ALF and BET (degrees, the arc sweeping less than 180
degrees) and end points PS and PE: |dist(C, P) - R| when P lies within
*mtol* of an end point or its direction from C lies within the arc's sweep,
else +undefined+.  R = +infinite+ marks a straight segment PS-PE: the
perpendicular distance of P from it, +undefined+ where P does not project
inside it (signals what DIST-PERP-PT-LINE-SEG signals).  When P is C the C
reads an uninitialized angle (it ignores N_VectorDirectedAngle's error);
the Lisp takes that angle as 0.
NLib: N_DistPtArc."
  (declare (type point p c ps pe) (type double-float r alf bet))
  ;; see if the segment is a line
  (when (= r +infinite+)
    (let ((dis (dist-perp-pt-line-seg p ps pe)))
      (declare (type double-float dis))
      (return-from dist-pt-arc (if (>= dis 0d0) dis +undefined+))))
  ;; check the special case
  (let* ((ds (dist-pt-pt p ps))
         (de (dist-pt-pt p pe))
         (dis (dist-pt-pt c p))
         (dif (abs (- dis r)))
         (mtol *mtol*))
    (declare (type double-float ds de dis dif mtol))
    (when (or (<= ds mtol) (<= de mtol))
      (return-from dist-pt-arc dif))
    ;; process the circle case now
    (let ((rr (make-point)) (v (make-point)) (ang 0d0))
      (declare (dynamic-extent rr v) (type double-float ang))
      (vector-create 1d0 0d0 0d0 rr)
      (vector-diff p c v)
      ;; the C ignores N_VectorDirectedAngle's error
      (handler-case (setf ang (vector-directed-angle rr v))
        (nurbs-error () nil))
      (cond ((<= (abs (- alf bet)) 180d0)
             (cond ((< alf bet)
                    (if (and (<= alf ang) (<= ang bet)) dif +undefined+))
                   ((> alf bet)
                    (if (and (>= alf ang) (>= ang bet)) dif +undefined+))
                   (t +undefined+)))
            ((> (abs (- alf bet)) 180d0)
             (cond ((and (< alf 180d0) (> bet 180d0))
                    (if (or (and (<= 0d0 ang) (<= ang alf)) (and (<= bet ang) (<= ang 360d0)))
                        dif +undefined+))
                   ((and (> alf 180d0) (< bet 180d0))
                    (if (or (and (<= 0d0 ang) (<= ang bet)) (and (<= alf ang) (<= ang 360d0)))
                        dif +undefined+))
                   (t +undefined+)))
            (t +undefined+)))))

;;; ------------------------------------------------------------- L11-g0

(defun dist-max-line-arc (c r alf bet ps pe p q)
  "The maximum radial distance between the segment P Q and the circular arc
of center C, radius R, start and end angles ALF and BET (degrees, sweeping
less than 180) and end points PS and PE: the largest |dist(C, X) - R| over
the part of the segment inside the arc's wedge (its ends, clipped to the
wedge's radial lines, and C's projection onto it).  +undefined+ when the
segment lies outside the wedge.  R = +infinite+ marks a straight segment
PS PE: the larger perpendicular distance of P and Q from it, +undefined+
where either does not project inside it (signals what dist-perp-pt-line-seg
signals).  Signals :num-err as project-pt-line-param does.  Where a callee
fails on the way the C leaves through EXIT returning no error and *d unset;
the Lisp returns 0 there (dist-pt-arc's and isect-line-line's errors, e.g.
C on PS or PE).
NLib: N_DistMaxLineArc."
  (declare (type point c ps pe p q) (type double-float r alf bet))
  ;; the segment a line?
  (when (= r +infinite+)
    (let ((d1 (dist-perp-pt-line-seg p ps pe))
          (d2 (dist-perp-pt-line-seg q ps pe)))
      (declare (type double-float d1 d2))
      (return-from dist-max-line-arc
        (if (or (< d1 0d0) (< d2 0d0)) +undefined+ (nl-max d1 d2)))))
  (let ((pp (make-point)) (qq (make-point)) (rr (make-point)) (m (make-point))
        (lpq (%make-lineseg (make-point) (make-point) 0))
        (lcs (%make-lineseg (make-point) (make-point) 0))
        (lce (%make-lineseg (make-point) (make-point) 0))
        (d1 0d0) (d2 0d0) (d3 0d0) (d 0d0))
    (declare (type point pp qq rr m) (dynamic-extent pp qq rr m)
             (type double-float d1 d2 d3 d))
    ;; C passes P and Q by value and may overwrite its copies
    (vector-copy p pp)
    (vector-copy q qq)
    (handler-case
        (progn
          ;; classify the points with respect to the wedge
          (setf d1 (dist-pt-arc pp c r alf bet ps pe)
                d2 (dist-pt-arc qq c r alf bet ps pe))
          (when (and (= d1 +undefined+) (= d2 +undefined+))
            (return-from dist-max-line-arc +undefined+))
          ;; split the line segment if necessary
          (create-line-pt-pt lpq pp qq +bounded+)
          (create-line-pt-pt lcs c ps +unbounded+)
          (create-line-pt-pt lce c pe +unbounded+)
          (flet ((clip (end)
                   ;; R: where the segment crosses a radial line of the wedge
                   (let ((ifl (isect-line-line lpq lcs rr)))
                     (declare (type fixnum ifl))
                     (setf d3 (if (= ifl +true+) (dist-pt-pt rr end) +bigd+))
                     (when (or (= ifl +false+) (< d3 *mtol*))
                       (setf ifl (isect-line-line lpq lce rr))
                       (when (= ifl +false+)
                         (return-from dist-max-line-arc +undefined+))))))
            (when (and (/= d1 +undefined+) (= d2 +undefined+))
              (clip pp)
              (vector-copy rr qq))
            (when (and (= d1 +undefined+) (/= d2 +undefined+))
              (clip qq)
              (vector-copy rr pp))))
      (nurbs-error () (return-from dist-max-line-arc 0d0)))
    ;; the circle case
    (setf d1 (abs (- (dist-pt-pt c pp) r))
          d2 (abs (- (dist-pt-pt c qq) r))
          d (nl-max d1 d2))
    (let ((tt (project-pt-line-param c pp qq m)))
      (declare (type double-float tt))
      (if (and (>= tt 0d0) (<= tt 1d0))
          (setf d3 (abs (- (dist-pt-pt c m) r)))
          (setf d3 -1d0)))
    (nl-max d d3)))

;;; ------------------------------------------------------------- L14-g0

(defun vector-diff-c-pts (pw qw v)
  "Write into V the vector from the Euclidean point of control point QW to
that of PW (PW - QW).  Returns V.
NLib: N_VectorDiffCPts."
  (declare (type cpoint pw qw) (type point v))
  (let ((p (make-point)) (q (make-point)))
    (declare (type point p q) (dynamic-extent p q))
    (c-pt-to-pt-euclid pw p)
    (c-pt-to-pt-euclid qw q)
    (vector-diff p q v)))

(defun line-segs-2d-are-intersecting (x1 y1 x2 y2 x3 y3 x4 y4 tol)
  "True when the 2D segments (X1,Y1)-(X2,Y2) and (X3,Y3)-(X4,Y4) intersect, to
parameter tolerance TOL: a bounding box check, then the signs of the two
segments' parameters scaled by the denominator.  Segments within *ltol* of
parallel (collinear ones included) do not intersect.
NLib: N_LineSegs2dAreIntersecting."
  (declare (type double-float x1 y1 x2 y2 x3 y3 x4 y4 tol))
  ;; bounding box check
  (when (< (+ (nl-max x3 x4) tol) (- (nl-min x1 x2) tol))
    (return-from line-segs-2d-are-intersecting nil))
  (when (< (+ (nl-max x1 x2) tol) (- (nl-min x3 x4) tol))
    (return-from line-segs-2d-are-intersecting nil))
  (when (< (+ (nl-max y3 y4) tol) (- (nl-min y1 y2) tol))
    (return-from line-segs-2d-are-intersecting nil))
  (when (< (+ (nl-max y1 y2) tol) (- (nl-min y3 y4) tol))
    (return-from line-segs-2d-are-intersecting nil))
  ;; intersection
  (let* ((ax (- x2 x1))
         (ay (- y2 y1))
         (bx (- x3 x4))
         (by (- y3 y4))
         (cx (- x1 x3))
         (cy (- y1 y3))
         (den (- (* ay bx) (* ax by))))
    (declare (type double-float ax ay bx by cx cy den))
    (when (< (abs den) (the double-float *ltol*))
      (return-from line-segs-2d-are-intersecting nil))
    (let ((numa (- (* by cx) (* bx cy)))
          (denp (+ den tol))
          (denm (- den tol)))
      (declare (type double-float numa denp denm))
      (if (> den 0d0)
          (when (or (< numa (- tol)) (> numa denp))
            (return-from line-segs-2d-are-intersecting nil))
          (when (or (> numa tol) (< numa denm))
            (return-from line-segs-2d-are-intersecting nil)))
      (let ((numb (- (* ax cy) (* ay cx))))
        (declare (type double-float numb))
        (if (> den 0d0)
            (when (or (< numb (- tol)) (> numb denp))
              (return-from line-segs-2d-are-intersecting nil))
            (when (or (> numb tol) (< numb denm))
              (return-from line-segs-2d-are-intersecting nil)))
        t))))

(defun b-box-2d-calc-overlap (xl xr yb yt ul ur vb vt tol)
  "The area of overlap of the 2D boxes [XL,XR]x[YB,YT] and [UL,UR]x[VB,VT]:
-1 when they are apart by more than TOL in x or y, 0 when they overlap by
less than TOL in x or y, else the area of the overlap.
NLib: N_BBox2dCalcOverlap."
  (declare (type double-float xl xr yb yt ul ur vb vt tol))
  (let ((w (- (nl-min xr ur) (nl-max xl ul)))
        (h (- (nl-min yt vt) (nl-max yb vb))))
    (declare (type double-float w h))
    (cond ((or (< w (- tol)) (< h (- tol))) -1d0)
          ((or (< w tol) (< h tol)) 0d0)
          (t (* w h)))))
