;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Math.cpp (usd-brep source/SMLib/src).  Modified:
;;;; translated to Common Lisp.  NL_POINT and NL_CPOINT arguments, which C
;;;; passes by value, are arrays here: every routine reads all of its input
;;;; coordinates before writing any output, so an output may be an input.

(in-package :batten)

(defvar *bypass-float-checks* nil
  "When true, FLOAT-OP-IS-BAD answers false without checking, as NLib does
when its NLIBCHKBYPASS is set.  NLib's default, and this one, is to check.")

(declaim (inline float-op-is-bad))
(defun float-op-is-bad (x y op)
  "True when X OP Y would overflow or underflow a double.  OP is +addition+,
+subtraction+, +multiplication+ or +division+.  NLib: N_FloatOpIsBad."
  (declare (type double-float x y) (type fixnum op))
  (cond
    (*bypass-float-checks* nil)
    ((= op +division+)
     (cond ((= y 0d0) t)
           ((= x 0d0) nil)
           (t (let ((x (abs x)) (y (abs y)))
                (cond ((and (<= y 1d0) (< x (* +bigd+ y))) nil)
                      ((and (>= y 1d0) (> x (* +smad+ y))) nil)
                      (t t))))))
    ((= op +multiplication+)
     (cond ((= x 0d0) nil)
           ((= y 0d0) nil)
           (t (let ((x (abs x)) (y (abs y)))
                (cond ((and (>= x 1d0) (<= y 1d0)) nil)
                      ((and (<= x 1d0) (>= y 1d0)) nil)
                      ((>= x 1d0) (not (or (<= y (/ +bigd+ x)) (<= x (/ +bigd+ y)))))
                      (t (not (> x (/ +smad+ y)))))))))
    ((= op +addition+)
     (if (/= (nl-sign x) (nl-sign y))
         nil
         (not (< (abs y) (- +bigd+ (abs x))))))
    ((= op +subtraction+)
     (if (> (* (nl-sign x) (nl-sign y)) 0d0)
         nil
         (not (> (- (abs y)) (- (abs x) +bigd+)))))
    (t t)))

;;; ------------------------------------------------- homogeneous to Euclidean

(declaim (inline %homogeneous-to-euclid %cpoint-row-to-euclid))
(defun %homogeneous-to-euclid (wx wy wz w)
  "The Euclidean point of homogeneous (WX WY WZ W) as three values, honoring
+no-z+ and +no-w+, multiplying by the reciprocal of w as N_CPtToPtEuclid does
(which rounds differently from dividing by w)."
  (declare (type double-float wx wy wz w))
  (fast
    (if (/= w +no-w+)
        (let ((r (/ 1d0 w)))
          (values (* wx r) (* wy r) (if (/= wz +no-z+) (* wz r) 0d0)))
        (values wx wy (if (/= wz +no-z+) wz 0d0)))))

(defun %cpoint-row-to-euclid (cpoints i)
  "The Euclidean point of row I of a cpoint-array, as three values."
  (declare (type cpoint-array cpoints) (type index i))
  (fast (%homogeneous-to-euclid (aref cpoints i 0) (aref cpoints i 1)
                                (aref cpoints i 2) (aref cpoints i 3))))

(declaim (inline c-pt-to-pt-euclid))
(defun c-pt-to-pt-euclid (pw p)
  "Write the Euclidean point of control point PW into P.  NLib: N_CPtToPtEuclid."
  (declare (type cpoint pw) (type point p))
  (multiple-value-bind (x y z) (%homogeneous-to-euclid (aref pw 0) (aref pw 1) (aref pw 2) (aref pw 3))
    (setf (aref p 0) x (aref p 1) y (aref p 2) z)
    p))

;;; ---------------------------------------------------- control point arithmetic

(defmacro with-cpoint ((x y z w) cp &body body)
  (let ((c (gensym)))
    `(let* ((,c ,cp) (,x (aref ,c 0)) (,y (aref ,c 1)) (,z (aref ,c 2)) (,w (aref ,c 3)))
       (declare (type double-float ,x ,y ,z ,w) (ignorable ,x ,y ,z ,w))
       ,@body)))

(defmacro with-point ((x y z) pt &body body)
  (let ((c (gensym)))
    `(let* ((,c ,pt) (,x (aref ,c 0)) (,y (aref ,c 1)) (,z (aref ,c 2)))
       (declare (type double-float ,x ,y ,z) (ignorable ,x ,y ,z))
       ,@body)))

(defmacro set-cpoint (cp x y z w)
  (let ((c (gensym)))
    `(let ((,c ,cp)) (setf (aref ,c 0) ,x (aref ,c 1) ,y (aref ,c 2) ,z (aref ,c 3) ,w) ,c)))

(defmacro set-point (pt x y z)
  (let ((c (gensym)))
    `(let ((,c ,pt)) (setf (aref ,c 0) ,x (aref ,c 1) ,y (aref ,c 2) ,z) ,c)))

;;; ------------------------------------------------- rows of arrays and nets

;;; C passes Pw[i] (a row of a cpoint array) or Pw[i][j] (an entry of a net)
;;; straight to N_CopyCPt, N_Combine2CPts and the like.  In Batten those take
;;; cpoints, so a port moves the row into a scratch cpoint, calls the routine
;;; and moves the result back: a copy changes no bits, so the arithmetic stays
;;; NLib's.  Use these rather than private copies of them.

(declaim (inline %row-c-pt %c-pt-row %row-pt %pt-row %net-c-pt %c-pt-net %net-pt %pt-net
                 %copy-c-pt-row %copy-pt-row %copy-c-pt-net))

(defun %row-c-pt (arr i cp)
  "Copy row I of cpoint-array ARR into cpoint CP; return CP."
  (declare (type cpoint-array arr) (type index i) (type cpoint cp))
  (fast (set-cpoint cp (aref arr i 0) (aref arr i 1) (aref arr i 2) (aref arr i 3))))

(defun %c-pt-row (cp arr i)
  "Copy cpoint CP into row I of cpoint-array ARR; return ARR."
  (declare (type cpoint-array arr) (type index i) (type cpoint cp))
  (fast (setf (aref arr i 0) (aref cp 0) (aref arr i 1) (aref cp 1)
              (aref arr i 2) (aref cp 2) (aref arr i 3) (aref cp 3)))
  arr)

(defun %row-pt (arr i p)
  "Copy row I of point-array ARR into point P; return P."
  (declare (type point-array arr) (type index i) (type point p))
  (fast (set-point p (aref arr i 0) (aref arr i 1) (aref arr i 2))))

(defun %pt-row (p arr i)
  "Copy point P into row I of point-array ARR; return ARR."
  (declare (type point-array arr) (type index i) (type point p))
  (fast (setf (aref arr i 0) (aref p 0) (aref arr i 1) (aref p 1) (aref arr i 2) (aref p 2)))
  arr)

(defun %net-c-pt (net i j cp)
  "Copy entry I,J of cpoint-net NET into cpoint CP; return CP."
  (declare (type cpoint-net net) (type index i j) (type cpoint cp))
  (fast (set-cpoint cp (aref net i j 0) (aref net i j 1) (aref net i j 2) (aref net i j 3))))

(defun %c-pt-net (cp net i j)
  "Copy cpoint CP into entry I,J of cpoint-net NET; return NET."
  (declare (type cpoint-net net) (type index i j) (type cpoint cp))
  (fast (setf (aref net i j 0) (aref cp 0) (aref net i j 1) (aref cp 1)
              (aref net i j 2) (aref cp 2) (aref net i j 3) (aref cp 3)))
  net)

(defun %net-pt (net i j p)
  "Copy entry I,J of point-net NET into point P; return P."
  (declare (type point-net net) (type index i j) (type point p))
  (fast (set-point p (aref net i j 0) (aref net i j 1) (aref net i j 2))))

(defun %pt-net (p net i j)
  "Copy point P into entry I,J of point-net NET; return NET."
  (declare (type point-net net) (type index i j) (type point p))
  (fast (setf (aref net i j 0) (aref p 0) (aref net i j 1) (aref p 1) (aref net i j 2) (aref p 2)))
  net)

(defun %copy-c-pt-row (src i dst j)
  "N_CopyCPt( src[i], &dst[j] ) on rows of cpoint-arrays; return DST."
  (declare (type cpoint-array src dst) (type index i j))
  (fast (let ((a (aref src i 0)) (b (aref src i 1)) (c (aref src i 2)) (d (aref src i 3)))
          (setf (aref dst j 0) a (aref dst j 1) b (aref dst j 2) c (aref dst j 3) d)))
  dst)

(defun %copy-pt-row (src i dst j)
  "N_CopyPt( src[i], &dst[j] ) on rows of point-arrays; return DST."
  (declare (type point-array src dst) (type index i j))
  (fast (let ((a (aref src i 0)) (b (aref src i 1)) (c (aref src i 2)))
          (setf (aref dst j 0) a (aref dst j 1) b (aref dst j 2) c)))
  dst)

(defun %copy-c-pt-net (src i j dst k l)
  "N_CopyCPt( src[i][j], &dst[k][l] ) on entries of cpoint-nets; return DST."
  (declare (type cpoint-net src dst) (type index i j k l))
  (fast (let ((a (aref src i j 0)) (b (aref src i j 1)) (c (aref src i j 2)) (d (aref src i j 3)))
          (setf (aref dst k l 0) a (aref dst k l 1) b (aref dst k l 2) c (aref dst k l 3) d)))
  dst)

(declaim (inline %combine-2-c-pt-rows %vector-blend-c-pt-row))

(defun %combine-2-c-pt-rows (alpha a ai beta b bi c ci)
  "N_Combine2CPts( alpha, a[ai], beta, b[bi], &c[ci] ) on rows of
cpoint-arrays: row CI of C = ALPHA * row AI of A + BETA * row BI of B,
honoring +no-z+ and +no-w+ as N_Combine2CPts does.  Reads before writing."
  (declare (type double-float alpha beta) (type cpoint-array a b c) (type fixnum ai bi ci))
  (let ((ax (aref a ai 0)) (ay (aref a ai 1)) (az (aref a ai 2)) (aw (aref a ai 3))
        (bx (aref b bi 0)) (by (aref b bi 1)) (bz (aref b bi 2)) (bw (aref b bi 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref c ci 0) (+ (* alpha ax) (* beta bx))
          (aref c ci 1) (+ (* alpha ay) (* beta by))
          (aref c ci 2) (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ (* alpha az) (* beta bz)))
                              ((/= az +no-z+) (* alpha az))
                              ((/= bz +no-z+) (* beta bz))
                              (t +no-z+))
          (aref c ci 3) (if (and (/= aw +no-w+) (/= bw +no-w+)) (+ (* alpha aw) (* beta bw)) +no-w+))
    c))

(defun %vector-blend-c-pt-row (alpha a ai b bi)
  "N_VectorBlendCPt( alpha, a[ai], &b[bi] ) on rows of cpoint-arrays: row BI
of B = row BI of B + ALPHA * row AI of A, honoring +no-z+ and +no-w+ of A as
N_VectorBlendCPt does."
  (declare (type double-float alpha) (type cpoint-array a b) (type fixnum ai bi))
  (let ((ax (aref a ai 0)) (ay (aref a ai 1)) (az (aref a ai 2)) (aw (aref a ai 3))
        (bx (aref b bi 0)) (by (aref b bi 1)) (bz (aref b bi 2)) (bw (aref b bi 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref b bi 0) (+ bx (* alpha ax))
          (aref b bi 1) (+ by (* alpha ay))
          (aref b bi 2) (if (/= az +no-z+) (+ bz (* alpha az)) +no-z+)
          (aref b bi 3) (if (/= aw +no-w+) (+ bw (* alpha aw)) +no-w+))
    b))

(declaim (inline copy-c-pt c-pt-get-w c-pt-from-wx-wy-wz c-pt-to-wx-wy-wz pt-to-c-pt weight
                 combine-2-c-pts sum-2-c-pts diff-2-c-pts scale-c-pt combine-4-c-pts
                 translate-c-pt translate-sum-2-c-pts vector-blend-c-pt))

(defun copy-c-pt (aw bw)
  "Copy control point AW into BW.  NLib: N_CopyCPt."
  (declare (type cpoint aw bw))
  (with-cpoint (x y z w) aw (set-cpoint bw x y z w)))

(defun c-pt-get-w (pw)
  "The weight of PW.  NLib: N_CPtGetW."
  (declare (type cpoint pw))
  (aref pw 3))

(defun c-pt-from-wx-wy-wz (wx wy wz w pw)
  "Set PW to (WX WY WZ W).  NLib: N_CPtFromWxWyWz."
  (declare (type double-float wx wy wz w) (type cpoint pw))
  (set-cpoint pw wx wy wz w))

(defun c-pt-to-wx-wy-wz (pw)
  "The coordinates of PW as four values.  NLib: N_CPtToWxWyWz."
  (declare (type cpoint pw))
  (values (aref pw 0) (aref pw 1) (aref pw 2) (aref pw 3)))

(defun pt-to-c-pt (p pw)
  "Set PW to point P, non-rational.  NLib: N_PtToCPt."
  (declare (type point p) (type cpoint pw))
  (with-point (x y z) p (set-cpoint pw x y z +no-w+)))

(defun weight (p w pw)
  "Set PW to point P with weight W: (w*x w*y w*z w), or P unweighted when W
is +no-w+.  NLib: N_Weight."
  (declare (type point p) (type double-float w) (type cpoint pw))
  (with-point (x y z) p
    (if (= w +no-w+)
        (set-cpoint pw x y z +no-w+)
        (set-cpoint pw (* w x) (* w y) (* w z) w))))

(defun combine-2-c-pts (alpha aw beta bw cw)
  "CW = ALPHA*AW + BETA*BW, honoring +no-z+ and +no-w+.  NLib: N_Combine2CPts."
  (declare (type double-float alpha beta) (type cpoint aw bw cw))
  (with-cpoint (ax ay az aw*) aw
    (with-cpoint (bx by bz bw*) bw
      (set-cpoint cw
                  (+ (* alpha ax) (* beta bx))
                  (+ (* alpha ay) (* beta by))
                  (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ (* alpha az) (* beta bz)))
                        ((/= az +no-z+) (* alpha az))
                        ((/= bz +no-z+) (* beta bz))
                        (t +no-z+))
                  (if (and (/= aw* +no-w+) (/= bw* +no-w+)) (+ (* alpha aw*) (* beta bw*)) +no-w+)))))

(defun sum-2-c-pts (aw bw cw)
  "CW = AW + BW, honoring +no-z+ and +no-w+.  NLib: N_Sum2CPts."
  (declare (type cpoint aw bw cw))
  (with-cpoint (ax ay az aw*) aw
    (with-cpoint (bx by bz bw*) bw
      (set-cpoint cw (+ ax bx) (+ ay by)
                  (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ az bz))
                        ((/= az +no-z+) az)
                        ((/= bz +no-z+) bz)
                        (t +no-z+))
                  (if (and (/= aw* +no-w+) (/= bw* +no-w+)) (+ aw* bw*) +no-w+)))))

(defun diff-2-c-pts (aw bw cw)
  "CW = AW - BW, honoring +no-z+ and +no-w+.  NLib: N_Diff2CPts."
  (declare (type cpoint aw bw cw))
  (with-cpoint (ax ay az aw*) aw
    (with-cpoint (bx by bz bw*) bw
      (set-cpoint cw (- ax bx) (- ay by)
                  (cond ((and (/= az +no-z+) (/= bz +no-z+)) (- az bz))
                        ((/= az +no-z+) az)
                        ((/= bz +no-z+) (- bz))
                        (t +no-z+))
                  (if (and (/= aw* +no-w+) (/= bw* +no-w+)) (- aw* bw*) +no-w+)))))

(defun scale-c-pt (alpha aw bw)
  "BW = ALPHA*AW, honoring +no-z+ and +no-w+.  NLib: N_ScaleCPt."
  (declare (type double-float alpha) (type cpoint aw bw))
  (with-cpoint (x y z w) aw
    (set-cpoint bw (* alpha x) (* alpha y)
                (if (/= z +no-z+) (* alpha z) +no-z+)
                (if (/= w +no-w+) (* alpha w) +no-w+))))

(defun combine-4-c-pts (alf aw bet bw gam cw del dw ew)
  "EW = ALF*AW + BET*BW + GAM*CW + DEL*DW.  NLib: N_Combine4CPts."
  (declare (type double-float alf bet gam del) (type cpoint aw bw cw dw ew))
  (with-cpoint (ax ay az aw*) aw
    (with-cpoint (bx by bz bw*) bw
      (with-cpoint (cx cy cz cw*) cw
        (with-cpoint (dx dy dz dw*) dw
          (set-cpoint ew
                      (+ (+ (+ (* alf ax) (* bet bx)) (* gam cx)) (* del dx))
                      (+ (+ (+ (* alf ay) (* bet by)) (* gam cy)) (* del dy))
                      (if (and (/= az +no-z+) (/= bz +no-z+) (/= cz +no-z+) (/= dz +no-z+))
                          (+ (+ (+ (* alf az) (* bet bz)) (* gam cz)) (* del dz))
                          +no-z+)
                      (if (and (/= aw* +no-w+) (/= bw* +no-w+) (/= cw* +no-w+) (/= dw* +no-w+))
                          (+ (+ (+ (* alf aw*) (* bet bw*)) (* gam cw*)) (* del dw*))
                          +no-w+)))))))

(defun translate-c-pt (pw tv qw)
  "QW = PW translated by vector TV in Euclidean space, keeping PW's weight.
NLib: N_TranslateCPt."
  (declare (type cpoint pw qw) (type point tv))
  (with-cpoint (x y z w) pw
    (with-point (tx ty tz) tv
      (when (/= w +no-w+)
        (setf x (/ x w) y (/ y w))
        (when (/= z +no-z+) (setf z (/ z w))))
      (let ((qx (+ x tx)) (qy (+ y ty)) (qz (if (/= z +no-z+) (+ z tz) +no-z+)))
        (declare (type double-float qx qy qz))
        (if (/= w +no-w+)
            (set-cpoint qw (* qx w) (* qy w) (if (/= z +no-z+) (* qz w) qz) w)
            (set-cpoint qw qx qy qz +no-w+))))))

(defun translate-sum-2-c-pts (cw alpha aw beta bw dw)
  "DW = CW + ALPHA*AW + BETA*BW.  NLib: N_TranslateSum2CPts."
  (declare (type double-float alpha beta) (type cpoint cw aw bw dw))
  (with-cpoint (cx cy cz cw*) cw
    (with-cpoint (ax ay az aw*) aw
      (with-cpoint (bx by bz bw*) bw
        (set-cpoint dw
                    (+ (+ cx (* alpha ax)) (* beta bx))
                    (+ (+ cy (* alpha ay)) (* beta by))
                    (if (and (/= az +no-z+) (/= bz +no-z+) (/= cz +no-z+))
                        (+ (+ cz (* alpha az)) (* beta bz))
                        +no-z+)
                    (if (and (/= aw* +no-w+) (/= bw* +no-w+) (/= cw* +no-w+))
                        (+ (+ cw* (* alpha aw*)) (* beta bw*))
                        +no-w+))))))

(defun vector-blend-c-pt (alpha aw bw)
  "BW = BW + ALPHA*AW, honoring +no-z+ and +no-w+ of AW.  NLib: N_VectorBlendCPt."
  (declare (type double-float alpha) (type cpoint aw bw))
  (with-cpoint (ax ay az aw*) aw
    (with-cpoint (bx by bz bw*) bw
      (set-cpoint bw (+ bx (* alpha ax)) (+ by (* alpha ay))
                  (if (/= az +no-z+) (+ bz (* alpha az)) +no-z+)
                  (if (/= aw* +no-w+) (+ bw* (* alpha aw*)) +no-w+)))))

;;; ------------------------------------------------------------ point arithmetic

(declaim (inline copy-pt combine-2-pts scale-pt vector-blend-pt translate-sum-2-pts
                 diff-2-pts sum-2-pts combine-4-pts dot-2-pts pt-magnitude))

(defun copy-pt (a b)
  "Copy point A into B.  NLib: N_CopyPt."
  (declare (type point a b))
  (with-point (x y z) a (set-point b x y z)))

(defun combine-2-pts (alpha a beta b c)
  "C = ALPHA*A + BETA*B.  NLib: N_Combine2Pts."
  (declare (type double-float alpha beta) (type point a b c))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (set-point c (+ (* alpha ax) (* beta bx)) (+ (* alpha ay) (* beta by)) (+ (* alpha az) (* beta bz))))))

(defun scale-pt (alpha a b)
  "B = ALPHA*A.  NLib: N_ScalePt."
  (declare (type double-float alpha) (type point a b))
  (with-point (x y z) a (set-point b (* alpha x) (* alpha y) (* alpha z))))

(defun vector-blend-pt (alpha a b)
  "B = B + ALPHA*A.  NLib: N_VectorBlendPt."
  (declare (type double-float alpha) (type point a b))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (set-point b (+ bx (* alpha ax)) (+ by (* alpha ay)) (+ bz (* alpha az))))))

(defun translate-sum-2-pts (c alpha a beta b d)
  "D = C + ALPHA*A + BETA*B.  NLib: N_TranslateSum2Pts."
  (declare (type double-float alpha beta) (type point c a b d))
  (with-point (cx cy cz) c
    (with-point (ax ay az) a
      (with-point (bx by bz) b
        (set-point d (+ (+ cx (* alpha ax)) (* beta bx))
                   (+ (+ cy (* alpha ay)) (* beta by))
                   (+ (+ cz (* alpha az)) (* beta bz)))))))

(defun diff-2-pts (a b c)
  "C = A - B.  NLib: N_Diff2Pts."
  (declare (type point a b c))
  (with-point (ax ay az) a (with-point (bx by bz) b (set-point c (- ax bx) (- ay by) (- az bz)))))

(defun sum-2-pts (a b c)
  "C = A + B.  NLib: N_Sum2Pts."
  (declare (type point a b c))
  (with-point (ax ay az) a (with-point (bx by bz) b (set-point c (+ ax bx) (+ ay by) (+ az bz)))))

(defun combine-4-pts (alf a bet b gam c del d e)
  "E = ALF*A + BET*B + GAM*C + DEL*D.  NLib: N_Combine4Pts."
  (declare (type double-float alf bet gam del) (type point a b c d e))
  (with-point (ax ay az) a
    (with-point (bx by bz) b
      (with-point (cx cy cz) c
        (with-point (dx dy dz) d
          (set-point e
                     (+ (+ (+ (* alf ax) (* bet bx)) (* gam cx)) (* del dx))
                     (+ (+ (+ (* alf ay) (* bet by)) (* gam cy)) (* del dy))
                     (+ (+ (+ (* alf az) (* bet bz)) (* gam cz)) (* del dz))))))))

(defun dot-2-pts (a b)
  "A . B.  NLib: N_Dot2Pts."
  (declare (type point a b))
  (with-point (ax ay az) a (with-point (bx by bz) b (+ (+ (* ax bx) (* ay by)) (* az bz)))))

(defun pt-magnitude (a)
  "|A|.  NLib: N_PtMagnitude."
  (declare (type point a))
  (c-sqrt (dot-2-pts a a)))

;;; ------------------------------------------------- more control point access

(defun c-pt-to-pt-and-w (pw p)
  "Write the Euclidean point of PW into P and return the weight of PW.
NLib: N_CPtToPtAndW."
  (declare (type cpoint pw) (type point p))
  (let ((w (c-pt-get-w pw)))
    (declare (type double-float w))
    (c-pt-to-pt-euclid pw p)
    w))

(defun scale-c-pt-xyz (alpha aw bw)
  "BW = ALPHA*AW in the first three coordinates (z only when not +no-z+),
keeping AW's weight.  NLib: N_ScaleCPtXYZ."
  (declare (type double-float alpha) (type cpoint aw bw))
  (with-cpoint (x y z w) aw
    (set-cpoint bw (* alpha x) (* alpha y)
                (if (/= z +no-z+) (* alpha z) +no-z+)
                w)))

(defun scale-c-pt-with-pt-and-vector (pw c f qw)
  "QW = PW scaled about center C by the vector-valued factor F, keeping PW's
weight.  NLib: N_ScaleCPtWithPtAndVector."
  (declare (type cpoint pw qw) (type point c f))
  (with-cpoint (px py pz pw*) pw
    (with-point (cx cy cz) c
      (with-point (fx fy fz) f
        (let* ((gx (- 1d0 fx)) (gy (- 1d0 fy)) (gz (- 1d0 fz))
               (w (if (/= pw* +no-w+) pw* 1d0))
               (cwx (* w cx)) (cwy (* w cy)) (cwz (* w cz)))
          (declare (type double-float gx gy gz w cwx cwy cwz))
          (set-cpoint qw
                      (+ (* fx px) (* gx cwx))
                      (+ (* fy py) (* gy cwy))
                      (if (/= pz +no-z+) (+ (* fz pz) (* gz cwz)) +no-z+)
                      (if (/= pw* +no-w+) pw* +no-w+)))))))

(defun pt-to-xyz (p)
  "The coordinates of P as three values.  NLib: N_PtToXYZ."
  (declare (type point p))
  (values (aref p 0) (aref p 1) (aref p 2)))

(defun pt-from-xyz (x y z p)
  "Set P to (X Y Z).  NLib: N_PtFromXYZ."
  (declare (type double-float x y z) (type point p))
  (set-point p x y z))

(defun c-pt-set-x (wx pw)
  "Set the x coordinate of PW to WX.  NLib: N_CPtSetX."
  (declare (type double-float wx) (type cpoint pw))
  (setf (aref pw 0) wx)
  pw)

(defun c-pt-set-y (wy pw)
  "Set the y coordinate of PW to WY.  NLib: N_CPtSetY."
  (declare (type double-float wy) (type cpoint pw))
  (setf (aref pw 1) wy)
  pw)

(defun c-pt-set-z (wz pw)
  "Set the z coordinate of PW to WZ.  NLib: N_CPtSetZ."
  (declare (type double-float wz) (type cpoint pw))
  (setf (aref pw 2) wz)
  pw)

(defun c-pt-set-w (w pw)
  "Set the weight of PW to W.  NLib: N_CPtSetW."
  (declare (type double-float w) (type cpoint pw))
  (setf (aref pw 3) w)
  pw)

(defun c-pt-rat-deriv (pw d wd dw)
  "DW = (wd*x + dx*w, wd*y + dy*w, wd*z + dz*w, wd), the homogeneous
derivative from control point PW = (wx wy wz w), the Euclidean derivative D and
the weight derivative WD.  NLib: N_CPtRatDeriv."
  (declare (type cpoint pw dw) (type point d) (type double-float wd))
  (let ((p (make-point)) (w (c-pt-get-w pw)) (pz (aref pw 2)))
    (declare (dynamic-extent p) (type double-float w pz))
    (c-pt-to-pt-euclid pw p)
    (with-point (x y z) p
      (with-point (dx dy dz) d
        (set-cpoint dw
                    (+ (* wd x) (* dx w))
                    (+ (* wd y) (* dy w))
                    (if (/= pz +no-z+) (+ (* wd z) (* dz w)) +no-z+)
                    wd)))))

(defun c-pt-to-pt (pw p)
  "Write PW's first three coordinates into P (z = 0 for +no-z+), dropping the
weight without dividing by it.  NLib: N_CPtToPt."
  (declare (type cpoint pw) (type point p))
  (with-cpoint (x y z w) pw
    (set-point p x y (if (/= z +no-z+) z 0d0))))

(defun swap-c-pts (pw qw)
  "Exchange control points PW and QW in place.  NLib: N_SwapCPts."
  (declare (type cpoint pw qw))
  (let ((aw (make-cpoint)))
    (declare (dynamic-extent aw))
    (copy-c-pt pw aw)
    (copy-c-pt qw pw)
    (copy-c-pt aw qw)
    pw))

;;; ------------------------------------------------------- reals and integers

(declaim (inline swap-reals))
(defun swap-reals (a b)
  "The two reals exchanged: (values B A).  NLib swaps through pointers; Lisp
returns the in-out pair.  NLib: N_SwapReals."
  (declare (type double-float a b))
  (values b a))

(defun shell-sort-real (a n)
  "Sort A[0..N] ascending in place by Shellsort, with NLib's increments.
NLib: N_ShellSortReal."
  (declare (type real-vector a) (type fixnum n))
  (fast
    (let ((k (1+ n)))
      (declare (type fixnum k))
      (do () ((<= k 1))
        (if (>= k 5)
            (setf k (truncate (- (* 5 k) 1) 11))
            (setf k 1))
        (do ((i (- n k) (1- i))) ((< i 0))
          (declare (type fixnum i))
          (let ((b (aref a i)) (j (+ i k)))
            (declare (type double-float b) (type fixnum j))
            (do () ((not (and (<= j n) (> b (aref a j)))))
              (setf (aref a (- j k)) (aref a j))
              (incf j k))
            (setf (aref a (- j k)) b)))))
    a))

(defun pascal-tri-row (tri row)
  "Fill TRI[0..ROW] with the rows of Pascal's triangle.  TRI, a fixnum matrix
of at least ROW+1 rows and columns, is allocated by the caller.
NLib: N_PascalTriRow."
  (declare (type (simple-array fixnum (* *)) tri) (type fixnum row))
  (fast
    (setf (aref tri 0 0) 1)
    (when (= row 0) (return-from pascal-tri-row tri))
    (setf (aref tri 1 0) 1 (aref tri 1 1) 1)
    (when (= row 1) (return-from pascal-tri-row tri))
    (do ((k 2 (1+ k))) ((> k row))
      (declare (type fixnum k))
      (setf (aref tri k 0) 1)
      (let ((r (truncate k 2)) (tmp2 1))
        (declare (type fixnum r tmp2))
        (do ((j 1 (1+ j))) ((> j r))
          (declare (type fixnum j))
          (let ((tmp1 (aref tri (1- k) j)))
            (declare (type fixnum tmp1))
            (setf (aref tri k j) (+ (aref tri (1- k) j) tmp2))
            (setf (aref tri k (- k j)) (aref tri k j))
            (setf tmp2 tmp1))))
      (setf (aref tri k k) 1)))
  tri)

;;; -------------------------------------------------------------- real matrices

;;; NL_RMATRIX is the rmatrix structure in types.lisp.

(defun create-real-matrix (rma n m rm mtp bw)
  "Set every member of RMA: highest indexes N and M, storage RM, type MTP and
bandwidth BW.  NLib: N_CreateRealMatrix."
  (declare (type rmatrix rma) (type fixnum n m mtp bw) (type (or null real-matrix) rm))
  (setf (rmatrix-n rma) n
        (rmatrix-m rma) m
        (rmatrix-rm rma) rm
        (rmatrix-mtp rma) mtp
        (rmatrix-bw rma) bw)
  rma)

(defun real-matrix-define (rma n m mtp bw)
  "Set the highest indexes, type and bandwidth of RMA, leaving its storage.
NLib: N_RealMatrixDefine."
  (declare (type rmatrix rma) (type fixnum n m mtp bw))
  (setf (rmatrix-n rma) n
        (rmatrix-m rma) m
        (rmatrix-mtp rma) mtp
        (rmatrix-bw rma) bw)
  rma)

(defun get-max-index-real-matrix (rma)
  "The highest row and column indexes of RMA, as two values.
NLib: N_GetMaxIndexRealMatrix."
  (declare (type rmatrix rma))
  (values (rmatrix-n rma) (rmatrix-m rma)))

(defun real-matrix-is-null (rma)
  "True when RMA is the NULL matrix (needs storage): a negative index, no
storage or a bandwidth below 1.  NLib: N_RealMatrixIsNULL."
  (declare (type rmatrix rma))
  (or (< (rmatrix-n rma) 0) (< (rmatrix-m rma) 0) (null (rmatrix-rm rma)) (<= (rmatrix-bw rma) 0)))

(defun %check-mem-real-matrix (rma n m mtp bw rname)
  "N_CheckMemRealMatrix (which vanishes): give a NULL RMA storage for an
N x M matrix of type MTP and bandwidth BW, as N_SetRealMatrix does; otherwise
check that RMA is defined at least N x M (:sto-err for RNAME, as
N_CheckRealMatStorage does) and redefine it."
  (declare (type rmatrix rma) (type fixnum n m mtp bw))
  (if (real-matrix-is-null rma)
      (create-real-matrix rma n m
                          (if (= mtp +mt-banded+)
                              (make-real-matrix (1+ n) (if (= bw 1) (1+ bw) bw))
                              (make-real-matrix (1+ n) (1+ m)))
                          mtp bw)
      (progn
        (when (or (< (rmatrix-n rma) n) (< (rmatrix-m rma) m))
          (nurbs-error :sto-err rname))
        (real-matrix-define rma n m mtp bw))))

(defun transform-c-pt (pw rma qw)
  "QW = PW transformed by the 4x4 matrix RMA (rotation and scale upper left,
translation in the last column, projection in the last row), applied to the
Euclidean point and re-weighted.  NLib: N_TransformCPt."
  (declare (type cpoint pw qw) (type rmatrix rma))
  (let ((rm (rmatrix-rm rma)))
    (declare (type real-matrix rm))
    (with-cpoint (x y z w) pw
      (when (/= w +no-w+)
        (setf x (/ x w) y (/ y w))
        (when (/= z +no-z+) (setf z (/ z w))))
      (let ((rx 0d0) (ry 0d0) (rz 0d0))
        (declare (type double-float rx ry rz))
        (if (/= z +no-z+)
            (let ((f (/ 1d0 (+ (+ (+ (* x (aref rm 3 0)) (* y (aref rm 3 1))) (* z (aref rm 3 2))) (aref rm 3 3)))))
              (declare (type double-float f))
              (setf rx (* f (+ (+ (+ (* x (aref rm 0 0)) (* y (aref rm 0 1))) (* z (aref rm 0 2))) (aref rm 0 3)))
                    ry (* f (+ (+ (+ (* x (aref rm 1 0)) (* y (aref rm 1 1))) (* z (aref rm 1 2))) (aref rm 1 3)))
                    rz (* f (+ (+ (+ (* x (aref rm 2 0)) (* y (aref rm 2 1))) (* z (aref rm 2 2))) (aref rm 2 3)))))
            (let ((f (/ 1d0 (+ (+ (* x (aref rm 3 0)) (* y (aref rm 3 1))) (aref rm 3 3)))))
              (declare (type double-float f))
              (setf rx (* f (+ (+ (* x (aref rm 0 0)) (* y (aref rm 0 1))) (aref rm 0 3)))
                    ry (* f (+ (+ (* x (aref rm 1 0)) (* y (aref rm 1 1))) (aref rm 1 3)))
                    rz +no-z+)))
        (if (/= w +no-w+)
            (set-cpoint qw (* rx w) (* ry w) (if (/= z +no-z+) (* rz w) rz) w)
            (set-cpoint qw rx ry rz +no-w+))))))

(defun real-matrix-multiply (rma rmb rmc)
  "RMC = RMA * RMB.  RMC is lower-left when both factors are, upper-right when
both are, full otherwise; a NULL RMC (make-rmatrix) gets storage, an existing
one must be defined large enough (:sto-err).  Signals :ind-err when the inner
dimensions differ and :inp-err for a banded factor.  Returns RMC.
NLib: N_RealMatrixMultiply."
  (declare (type rmatrix rma rmb rmc))
  (let ((na (rmatrix-n rma)) (ma (rmatrix-m rma)) (a (rmatrix-rm rma)) (mta (rmatrix-mtp rma))
        (nb (rmatrix-n rmb)) (mb (rmatrix-m rmb)) (b (rmatrix-rm rmb)) (mtb (rmatrix-mtp rmb)))
    (declare (type fixnum na ma nb mb mta mtb))
    (when (/= ma nb)
      (nurbs-error :ind-err 'real-matrix-multiply))
    (when (or (= mta +mt-banded+) (= mtb +mt-banded+))
      (nurbs-error :inp-err 'real-matrix-multiply))
    (let ((mtc (cond ((and (= mta +mt-lowerleft+) (= mtb +mt-lowerleft+)) +mt-lowerleft+)
                     ((and (= mta +mt-upperright+) (= mtb +mt-upperright+)) +mt-upperright+)
                     (t +mt-full+))))
      (%check-mem-real-matrix rmc na mb mtc na 'real-matrix-multiply))
    (let ((c (rmatrix-rm rmc)))
      (declare (type real-matrix a b c))
      (macrolet ((product (jlo jhi klo khi)
                   ;; C[i][j] = 0, then C[i][j] += A[i][k]*B[k][j], in memory
                   `(fast
                      (do ((i 0 (1+ i))) ((> i na))
                        (declare (type fixnum i))
                        (do ((j ,jlo (1+ j))) ((> j ,jhi))
                          (declare (type fixnum j))
                          (setf (aref c i j) 0d0)
                          (do ((k ,klo (1+ k)) (kh ,khi)) ((> k kh))
                            (declare (type fixnum k kh))
                            (setf (aref c i j) (+ (aref c i j) (* (aref a i k) (aref b k j))))))))))
        (cond ((and (= mta +mt-full+) (= mtb +mt-full+)) (product 0 mb 0 ma))
              ((and (= mta +mt-full+) (= mtb +mt-lowerleft+)) (product 0 mb j ma))
              ((and (= mta +mt-full+) (= mtb +mt-upperright+)) (product 0 mb 0 j))
              ((and (= mta +mt-lowerleft+) (= mtb +mt-full+)) (product 0 mb 0 i))
              ((and (= mta +mt-lowerleft+) (= mtb +mt-lowerleft+)) (product 0 i j i))
              ((and (= mta +mt-lowerleft+) (= mtb +mt-upperright+)) (product 0 mb 0 (nl-min i j)))
              ((and (= mta +mt-upperright+) (= mtb +mt-full+)) (product 0 mb i ma))
              ((and (= mta +mt-upperright+) (= mtb +mt-lowerleft+)) (product 0 mb (nl-max i j) ma))
              ((and (= mta +mt-upperright+) (= mtb +mt-upperright+)) (product i mb i j)))))
    rmc))

(defmacro %with-matrix-band ((n m a mtp bw sbw) rma &body body)
  "Bind the members of real matrix RMA, and SBW = BW/2."
  (let ((r (gensym)))
    `(let* ((,r ,rma) (,n (rmatrix-n ,r)) (,m (rmatrix-m ,r)) (,a (rmatrix-rm ,r))
            (,mtp (rmatrix-mtp ,r)) (,bw (rmatrix-bw ,r)) (,sbw (truncate ,bw 2)))
       (declare (type fixnum ,n ,m ,mtp ,bw ,sbw) (ignorable ,n ,m ,a ,mtp ,bw ,sbw))
       ,@body)))

(defun real-matrix-multiply-pt-array (rma p q)
  "Q = RMA * P for a point array P; Q is allocated by the caller.  A banded
RMA is indexed by its band.  Returns Q.  NLib: N_RealMatrixMultiplyPtArray."
  (declare (type rmatrix rma) (type point-array p q))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (let ((s (make-point)) (pj (make-point)))
      (declare (dynamic-extent s pj) (type real-matrix a))
      (macrolet ((rows (jlo jhi index)
                   ;; Q[i] = 0, then Q[i] += A[i][index] * P[j]
                   `(fast
                      (do ((i 0 (1+ i))) ((> i n))
                        (declare (type fixnum i))
                        (set-point s 0d0 0d0 0d0)
                        (do ((j ,jlo (1+ j)) (jh ,jhi)) ((> j jh))
                          (declare (type fixnum j jh))
                          (set-point pj (aref p j 0) (aref p j 1) (aref p j 2))
                          (vector-blend-pt (aref a i ,index) pj s))
                        (setf (aref q i 0) (aref s 0) (aref q i 1) (aref s 1) (aref q i 2) (aref s 2))))))
        (cond ((= mtp +mt-lowerleft+) (rows 0 i j))
              ((= mtp +mt-upperright+) (rows i m j))
              ((= mtp +mt-full+) (rows 0 m j))
              ((= mtp +mt-banded+) (rows (nl-max 0 (- i sbw)) (nl-min n (+ i sbw)) (+ (- j i) sbw)))))))
  q)

(defun real-matrix-multiply-c-pt-array (rma pw qw)
  "QW = RMA * PW for a control point array PW; QW is allocated by the caller.
A banded RMA is indexed by its band.  Returns QW.
NLib: N_RealMatrixMultiplyCPtArray."
  (declare (type rmatrix rma) (type cpoint-array pw qw))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (let ((s (make-cpoint)) (pj (make-cpoint)) (czero (make-cpoint 0d0 0d0 0d0 0d0)))
      (declare (dynamic-extent s pj czero) (type real-matrix a))
      (macrolet ((rows (jlo jhi index)
                   ;; Qw[i] = NL_CZERO, then Qw[i] += A[i][index] * Pw[j]
                   `(fast
                      (do ((i 0 (1+ i))) ((> i n))
                        (declare (type fixnum i))
                        (copy-c-pt czero s)
                        (do ((j ,jlo (1+ j)) (jh ,jhi)) ((> j jh))
                          (declare (type fixnum j jh))
                          (set-cpoint pj (aref pw j 0) (aref pw j 1) (aref pw j 2) (aref pw j 3))
                          (vector-blend-c-pt (aref a i ,index) pj s))
                        (setf (aref qw i 0) (aref s 0) (aref qw i 1) (aref s 1)
                              (aref qw i 2) (aref s 2) (aref qw i 3) (aref s 3))))))
        (cond ((= mtp +mt-lowerleft+) (rows 0 i j))
              ((= mtp +mt-upperright+) (rows i m j))
              ((= mtp +mt-full+) (rows 0 m j))
              ((= mtp +mt-banded+) (rows (nl-max 0 (- i sbw)) (nl-min n (+ i sbw)) (+ (- j i) sbw)))))))
  qw)

(defun real-matrix-lu-decompose (rma)
  "Decompose the square full or banded matrix RMA into LU in place, without
pivoting.  Signals :ind-err for a non-square matrix, :inp-err for another
type and :num-err for a pivot below *ludt*.  Returns RMA.
NLib: N_RealMatrixLuDecompose."
  (declare (type rmatrix rma))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (when (/= m n)
      (nurbs-error :ind-err 'real-matrix-lu-decompose))
    (when (and (/= mtp +mt-full+) (/= mtp +mt-banded+))
      (nurbs-error :inp-err 'real-matrix-lu-decompose))
    (let ((ludt *ludt*))
      (declare (type double-float ludt) (type real-matrix a))
      (if (= mtp +mt-full+)
          (fast
            (do ((i 0 (1+ i))) ((> i n))
              (declare (type fixnum i))
              (do ((j i (1+ j))) ((> j n))
                (declare (type fixnum j))
                (let ((sum 0d0))
                  (declare (type double-float sum))
                  (dotimes (k i)
                    (setf sum (+ sum (* (aref a i k) (aref a k j)))))
                  (setf (aref a i j) (- (aref a i j) sum))))
              (do ((j (1+ i) (1+ j))) ((> j n))
                (declare (type fixnum j))
                (let ((sum 0d0))
                  (declare (type double-float sum))
                  (dotimes (k i)
                    (setf sum (+ sum (* (aref a j k) (aref a k i)))))
                  (when (< (abs (aref a i i)) ludt)
                    (nurbs-error :num-err 'real-matrix-lu-decompose))
                  (setf (aref a j i) (/ (- (aref a j i) sum) (aref a i i)))))))
          (fast
            (do ((i 0 (1+ i))) ((> i n))
              (declare (type fixnum i))
              (let ((jh (nl-min n (+ i sbw))))
                (declare (type fixnum jh))
                (do ((j i (1+ j))) ((> j jh))
                  (declare (type fixnum j))
                  (let ((sum 0d0))
                    (declare (type double-float sum))
                    (do ((k (nl-max 0 (- j sbw)) (1+ k))) ((>= k i))
                      (declare (type fixnum k))
                      (setf sum (+ sum (* (aref a i (+ (- k i) sbw)) (aref a k (+ (- j k) sbw))))))
                    (setf (aref a i (+ (- j i) sbw)) (- (aref a i (+ (- j i) sbw)) sum))))
                (do ((j (1+ i) (1+ j))) ((> j jh))
                  (declare (type fixnum j))
                  (let ((sum 0d0))
                    (declare (type double-float sum))
                    (do ((k (nl-max 0 (- j sbw)) (1+ k))) ((>= k i))
                      (declare (type fixnum k))
                      (setf sum (+ sum (* (aref a j (+ (- k j) sbw)) (aref a k (+ (- i k) sbw))))))
                    (when (< (abs (aref a i sbw)) ludt)
                      (nurbs-error :num-err 'real-matrix-lu-decompose))
                    (setf (aref a j (+ (- i j) sbw))
                          (/ (- (aref a j (+ (- i j) sbw)) sum) (aref a i sbw)))))))))))
  rma)

(defun real-matrix-for-back (rma a ptp qw)
  "Solve LU-decomposed RMA * QW = A by forward elimination and back
substitution.  A is a point-array when PTP is +epoint+ (the solution is then
non-rational) or a cpoint-array when PTP is +hpoint+; QW, a cpoint-array, is
allocated by the caller.  Signals :ind-err for a non-square matrix, :inp-err
for a type other than full or banded, :num-err for a bad pivot and :cal-err
for another PTP.  Returns QW.  NLib: N_RealMatrixForBack."
  (declare (type rmatrix rma) (type fixnum ptp) (type cpoint-array qw))
  (%with-matrix-band (n m b mtp bw sbw) rma
    (when (/= m n)
      (nurbs-error :ind-err 'real-matrix-for-back))
    (when (and (/= mtp +mt-full+) (/= mtp +mt-banded+))
      (nurbs-error :inp-err 'real-matrix-for-back))
    (locally (declare (type real-matrix b))
      (macrolet ((diag (i) `(if (= mtp +mt-full+) (aref b ,i ,i) (aref b ,i sbw)))
                 (elt* (i j) `(if (= mtp +mt-full+) (aref b ,i ,j) (aref b ,i (+ (- ,j ,i) sbw))))
                 (jlo (i) `(if (= mtp +mt-full+) 0 (nl-max 0 (- ,i sbw))))
                 (jhi (i) `(if (= mtp +mt-full+) n (nl-min n (+ ,i sbw)))))
        (cond
          ((= ptp +epoint+)
           (let ((p a) (z (make-point-array (1+ n))) (q (make-point-array (1+ n)))
                 (s (make-point)) (t1 (make-point)) (t2 (make-point)) (qc (make-cpoint)))
             (declare (type point-array p z q) (dynamic-extent s t1 t2 qc))
             (macrolet ((row->pt (arr i pt) `(set-point ,pt (aref ,arr ,i 0) (aref ,arr ,i 1) (aref ,arr ,i 2)))
                        (pt->row (pt arr i) `(setf (aref ,arr ,i 0) (aref ,pt 0) (aref ,arr ,i 1) (aref ,pt 1)
                                                   (aref ,arr ,i 2) (aref ,pt 2))))
               (fast
                 (do ((i 0 (1+ i))) ((> i n))
                   (declare (type fixnum i))
                   (set-point s 0d0 0d0 0d0)
                   (do ((j (jlo i) (1+ j))) ((>= j i))
                     (declare (type fixnum j))
                     (row->pt z j t1)
                     (vector-blend-pt (elt* i j) t1 s))
                   (row->pt p i t1)
                   (diff-2-pts t1 s t2)
                   (pt->row t2 z i))
                 (do ((i n (1- i))) ((< i 0))
                   (declare (type fixnum i))
                   (when (float-op-is-bad 1d0 (diag i) +division+)
                     (nurbs-error :num-err 'real-matrix-for-back))
                   (let ((fact (/ 1d0 (diag i))))
                     (declare (type double-float fact))
                     (set-point s 0d0 0d0 0d0)
                     (do ((j (1+ i) (1+ j)) (jh (jhi i))) ((> j jh))
                       (declare (type fixnum j jh))
                       (row->pt q j t1)
                       (vector-blend-pt (elt* i j) t1 s))
                     (row->pt z i t1)
                     (combine-2-pts fact t1 (- fact) s t2)
                     (pt->row t2 q i)
                     (pt-to-c-pt t2 qc)
                     (setf (aref qw i 0) (aref qc 0) (aref qw i 1) (aref qc 1)
                           (aref qw i 2) (aref qc 2) (aref qw i 3) (aref qc 3))))))))
          ((= ptp +hpoint+)
           (let ((pw a) (zw (make-cpoint-array (1+ n)))
                 (sw (make-cpoint)) (t1 (make-cpoint)) (t2 (make-cpoint)) (czero (make-cpoint 0d0 0d0 0d0 0d0)))
             (declare (type cpoint-array pw zw) (dynamic-extent sw t1 t2 czero))
             (macrolet ((row->cpt (arr i cp) `(set-cpoint ,cp (aref ,arr ,i 0) (aref ,arr ,i 1)
                                                          (aref ,arr ,i 2) (aref ,arr ,i 3)))
                        (cpt->row (cp arr i) `(setf (aref ,arr ,i 0) (aref ,cp 0) (aref ,arr ,i 1) (aref ,cp 1)
                                                    (aref ,arr ,i 2) (aref ,cp 2) (aref ,arr ,i 3) (aref ,cp 3))))
               (fast
                 (do ((i 0 (1+ i))) ((> i n))
                   (declare (type fixnum i))
                   (copy-c-pt czero sw)
                   (do ((j (jlo i) (1+ j))) ((>= j i))
                     (declare (type fixnum j))
                     (row->cpt zw j t1)
                     (vector-blend-c-pt (elt* i j) t1 sw))
                   (row->cpt pw i t1)
                   (diff-2-c-pts t1 sw t2)
                   (cpt->row t2 zw i))
                 (do ((i n (1- i))) ((< i 0))
                   (declare (type fixnum i))
                   (when (float-op-is-bad 1d0 (diag i) +division+)
                     (nurbs-error :num-err 'real-matrix-for-back))
                   (let ((fact (/ 1d0 (diag i))))
                     (declare (type double-float fact))
                     (copy-c-pt czero sw)
                     (do ((j (1+ i) (1+ j)) (jh (jhi i))) ((> j jh))
                       (declare (type fixnum j jh))
                       (row->cpt qw j t1)
                       (vector-blend-c-pt (elt* i j) t1 sw))
                     (row->cpt zw i t1)
                     (combine-2-c-pts fact t1 (- fact) sw t2)
                     (cpt->row t2 qw i)))))))
          (t (nurbs-error :cal-err 'real-matrix-for-back))))))
  qw)

(defun real-matrix-right-for-back (rma y x)
  "Solve LU-decomposed RMA * X = Y for real vectors by forward elimination
and back substitution; X is allocated by the caller.  Signals :ind-err for a
non-square matrix, :inp-err for a type other than full or banded and :num-err
for a bad division.  Returns X.  NLib: N_RealMatrixRightForBack."
  (declare (type rmatrix rma) (type real-vector y x))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (when (/= m n)
      (nurbs-error :ind-err 'real-matrix-right-for-back))
    (when (and (/= mtp +mt-full+) (/= mtp +mt-banded+))
      (nurbs-error :inp-err 'real-matrix-right-for-back))
    (let ((z (make-real-vector (1+ n))))
      (declare (type real-matrix a) (type real-vector z))
      (macrolet ((solve (jlo jhi el diag)
                   `(fast
                      (do ((i 0 (1+ i))) ((> i n))
                        (declare (type fixnum i))
                        (let ((sum 0d0))
                          (declare (type double-float sum))
                          (do ((j ,jlo (1+ j))) ((>= j i))
                            (declare (type fixnum j))
                            (setf sum (+ sum (* ,el (aref z j)))))
                          (setf (aref z i) (- (aref y i) sum))))
                      (do ((i n (1- i))) ((< i 0))
                        (declare (type fixnum i))
                        (let ((sum 0d0))
                          (declare (type double-float sum))
                          (do ((j (1+ i) (1+ j)) (jh ,jhi)) ((> j jh))
                            (declare (type fixnum j jh))
                            (setf sum (+ sum (* ,el (aref x j)))))
                          (when (float-op-is-bad (- (aref z i) sum) ,diag +division+)
                            (nurbs-error :num-err 'real-matrix-right-for-back))
                          (setf (aref x i) (/ (- (aref z i) sum) ,diag)))))))
        (if (= mtp +mt-full+)
            (solve 0 n (aref a i j) (aref a i i))
            (solve (nl-max 0 (- i sbw)) (nl-min n (+ i sbw)) (aref a i (+ (- j i) sbw)) (aref a i sbw))))))
  x)

(defun real-matrix-for-back-pivot (rma ind b)
  "Solve RMA * X = B in place in B, RMA full and LU-decomposed with partial
pivoting, IND its permutation.  Signals :ind-err for a non-square matrix,
:inp-err for a type other than full and :num-err for a bad division.
Returns B.  NLib: N_RealMatrixForBackPivot."
  (declare (type rmatrix rma) (type index-vector ind) (type real-vector b))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (when (/= m n)
      (nurbs-error :ind-err 'real-matrix-for-back-pivot))
    (when (/= mtp +mt-full+)
      (nurbs-error :inp-err 'real-matrix-for-back-pivot))
    (let ((ii -1) (ludt *ludt*))
      (declare (type fixnum ii) (type double-float ludt) (type real-matrix a))
      (fast
        (do ((i 0 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (let* ((ip (aref ind i)) (sum (aref b ip)))
            (declare (type fixnum ip) (type double-float sum))
            (setf (aref b ip) (aref b i))
            (if (>= ii 0)
                (do ((j ii (1+ j))) ((>= j i))
                  (declare (type fixnum j))
                  (setf sum (- sum (* (aref a i j) (aref b j)))))
                (when (> (abs sum) ludt)
                  (setf ii i)))
            (setf (aref b i) sum)))
        (do ((i n (1- i))) ((< i 0))
          (declare (type fixnum i))
          (let ((sum (aref b i)))
            (declare (type double-float sum))
            (do ((j (1+ i) (1+ j))) ((> j n))
              (declare (type fixnum j))
              (setf sum (- sum (* (aref a i j) (aref b j)))))
            (when (float-op-is-bad sum (aref a i i) +division+)
              (nurbs-error :num-err 'real-matrix-for-back-pivot))
            (setf (aref b i) (/ sum (aref a i i))))))))
  b)

(defun set-real-matrix (rma n m mtp bw)
  "Give real matrix RMA fresh storage and define it N x M (highest indexes) of
type MTP and bandwidth BW: N+1 rows of M+1 elements, or of BW elements (2 when
BW is 1) for a banded matrix.  Returns RMA.  NLib: N_SetRealMatrix."
  (declare (type rmatrix rma) (type fixnum n m mtp bw))
  (let ((rm (if (= mtp +mt-banded+)
                (if (= bw 1)
                    (make-real-matrix (1+ n) (1+ bw))
                    (make-real-matrix (1+ n) bw))
                (make-real-matrix (1+ n) (1+ m)))))
    (create-real-matrix rma n m rm mtp bw)))

(defun real-matrix-transpose (rma rmb)
  "RMB = the transpose of RMA.  Lower-left becomes upper-right and the other
way about, banded stays banded (stored by its band), anything else is full.
When RMB is RMA the transposition is in place: a triangular matrix keeps its
storage, any other gets fresh (m+1) x (n+1) storage.  Otherwise a NULL RMB
gets storage and an existing one must be defined large enough (:sto-err).
Returns RMB.  NLib: N_RealMatrixTranspose."
  (declare (type rmatrix rma rmb))
  (%with-matrix-band (n m a mta bw sbw) rma
    (let ((mtb (cond ((= mta +mt-lowerleft+) +mt-upperright+)
                     ((= mta +mt-upperright+) +mt-lowerleft+)
                     ((= mta +mt-banded+) +mt-banded+)
                     (t +mt-full+)))
          (b nil))
      (declare (type fixnum mtb))
      (if (eq rma rmb)
          (progn
            (setf b (if (or (= mta +mt-lowerleft+) (= mta +mt-upperright+))
                        a
                        (make-real-matrix (1+ m) (1+ n))))
            (create-real-matrix rma m n b mtb bw))
          (progn
            (%check-mem-real-matrix rmb m n mtb bw 'real-matrix-transpose)
            (setf b (rmatrix-rm rmb))))
      (let ((a a) (b b))
        (declare (type real-matrix a b))
        (fast
          (cond
            ((= mta +mt-lowerleft+)
             (do ((i 0 (1+ i))) ((> i n))
               (declare (type fixnum i))
               (do ((j i (1+ j))) ((> j m))
                 (declare (type fixnum j))
                 (setf (aref b i j) (aref a j i)))))
            ((= mta +mt-upperright+)
             (do ((j 0 (1+ j))) ((> j m))
               (declare (type fixnum j))
               (do ((i j (1+ i))) ((> i n))
                 (declare (type fixnum i))
                 (setf (aref b i j) (aref a j i)))))
            ((= mta +mt-banded+)
             (do ((i 0 (1+ i))) ((> i n))
               (declare (type fixnum i))
               (let ((jl (nl-max 0 (- i sbw)))
                     (jh (nl-min n (+ i sbw))))
                 (declare (type fixnum jl jh))
                 (do ((j jl (1+ j))) ((> j jh))
                   (declare (type fixnum j))
                   (setf (aref b j (+ (- i j) sbw)) (aref a i (+ (- j i) sbw)))))))
            ((= mta +mt-full+)
             (do ((i 0 (1+ i))) ((> i n))
               (declare (type fixnum i))
               (do ((j 0 (1+ j))) ((> j m))
                 (declare (type fixnum j))
                 (setf (aref b j i) (aref a i j))))))))))
  rmb)

(defun real-matrix-lu-decompose-pivot (rma ind)
  "Decompose the square full matrix RMA into LU in place by Crout's method
with partial pivoting (the LU of a row permutation of RMA); IND, an
index-vector of n+1 allocated by the caller, receives the row interchanges.
Signals :ind-err for a non-square matrix, :inp-err for another type and
:num-err for a row or pivot below *ludt*.  Returns RMA.
NLib: N_RealMatrixLuDecomposePivot."
  (declare (type rmatrix rma) (type index-vector ind))
  (%with-matrix-band (n m a mtp bw sbw) rma
    (when (/= m n)
      (nurbs-error :ind-err 'real-matrix-lu-decompose-pivot))
    (when (/= mtp +mt-full+)
      (nurbs-error :inp-err 'real-matrix-lu-decompose-pivot))
    (let ((v (make-real-vector (1+ n)))
          (ludt *ludt*))
      (declare (type real-vector v) (type double-float ludt) (type real-matrix a))
      (fast
        ;; scaling info
        (do ((i 0 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (let ((big 0d0))
            (declare (type double-float big))
            (do ((j 0 (1+ j))) ((> j n))
              (declare (type fixnum j))
              (let ((temp (abs (aref a i j))))
                (declare (type double-float temp))
                (when (> temp big)
                  (setf big temp))))
            (when (< big ludt)
              (nurbs-error :num-err 'real-matrix-lu-decompose-pivot))
            (setf (aref v i) (/ 1d0 big))))
        ;; for each column
        (do ((j 0 (1+ j))) ((> j n))
          (declare (type fixnum j))
          (let ((imax j) (big 0d0))
            (declare (type fixnum imax) (type double-float big))
            (do ((i 0 (1+ i))) ((>= i j))
              (declare (type fixnum i))
              (let ((sum (aref a i j)))
                (declare (type double-float sum))
                (do ((k 0 (1+ k))) ((>= k i))
                  (declare (type fixnum k))
                  (setf sum (- sum (* (aref a i k) (aref a k j)))))
                (setf (aref a i j) sum)))
            (do ((i j (1+ i))) ((> i n))
              (declare (type fixnum i))
              (let ((sum (aref a i j)))
                (declare (type double-float sum))
                (do ((k 0 (1+ k))) ((>= k j))
                  (declare (type fixnum k))
                  (setf sum (- sum (* (aref a i k) (aref a k j)))))
                (setf (aref a i j) sum)
                (let ((dum (* (aref v i) (abs sum))))
                  (declare (type double-float dum))
                  (when (>= dum big)
                    (setf big dum
                          imax i)))))
            ;; exchange rows if necessary
            (when (/= j imax)
              (do ((k 0 (1+ k))) ((> k n))
                (declare (type fixnum k))
                (multiple-value-bind (x y) (swap-reals (aref a imax k) (aref a j k))
                  (declare (type double-float x y))
                  (setf (aref a imax k) x
                        (aref a j k) y)))
              (setf (aref v imax) (aref v j)))
            (setf (aref ind j) imax)
            ;; divide by the pivot
            (when (< (abs (aref a j j)) ludt)
              (nurbs-error :num-err 'real-matrix-lu-decompose-pivot))
            (when (< j n)
              (let ((dum (/ 1d0 (aref a j j))))
                (declare (type double-float dum))
                (do ((i (1+ j) (1+ i))) ((> i n))
                  (declare (type fixnum i))
                  (setf (aref a i j) (* (aref a i j) dum))))))))))
  rma)

(defun real-matrix-inverse (rma lud rmb)
  "RMB = the inverse of the square full or banded matrix RMA, without
pivoting, column by column through RMA's LU decomposition.  LUD +yes+ says
RMA is already LU-decomposed; +no+ decomposes it, in place.  A NULL RMB gets
full storage, an existing one must be defined large enough (:sto-err).
Signals :ind-err for a non-square RMA, :inp-err for another type, and the
errors of the decomposition and the solver.  Returns RMB.
NLib: N_RealMatrixInverse."
  (declare (type rmatrix rma rmb) (type fixnum lud))
  (let ((n (rmatrix-n rma)) (m (rmatrix-m rma)) (mta (rmatrix-mtp rma)))
    (declare (type fixnum n m mta))
    (when (/= n m)
      (nurbs-error :ind-err 'real-matrix-inverse))
    (when (and (/= mta +mt-full+) (/= mta +mt-banded+))
      (nurbs-error :inp-err 'real-matrix-inverse))
    (%check-mem-real-matrix rmb n m +mt-full+ n 'real-matrix-inverse)
    (let ((b (rmatrix-rm rmb)))
      (declare (type real-matrix b))
      (when (= lud +no+)
        (real-matrix-lu-decompose rma))
      (let ((x (make-real-vector (1+ n)))
            (y (make-real-vector (1+ n))))
        (declare (type real-vector x y))
        (setf (aref y 0) 1d0)
        (do ((i 1 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (setf (aref y i) 0d0))
        (real-matrix-right-for-back rma y x)
        (dotimes (i (1+ n))
          (setf (aref b i 0) (aref x i)))
        (do ((i 1 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (setf (aref y (1- i)) 0d0
                (aref y i) 1d0)
          (real-matrix-right-for-back rma y x)
          (dotimes (j (1+ n))
            (setf (aref b j i) (aref x j)))))))
  rmb)

(defconstant +real-matrix-single-value-decompose--mxits+ 30
  "MXITS in N_RealMatrixSingleValueDecompose.")

(defun real-matrix-single-value-decompose (a w v m n)
  "Singular value decomposition A = U W V^T of the (M+1) x (N+1) matrix A,
M >= N (else :inp-err), after Numerical Recipes: A is overwritten by U, the
singular values go to W (N+1) and the orthogonal V to V ((N+1) x (N+1)), both
allocated by the caller.  Signals :con-err when an iteration does not
converge.  Returns A.  NLib: N_RealMatrixSingleValueDecompose."
  (declare (type real-matrix a v) (type real-vector w) (type fixnum m n))
  (when (> n m)
    (nurbs-error :inp-err 'real-matrix-single-value-decompose))
  (let ((rv1 (make-real-vector (1+ n)))
        (i1 0) (j1 0) (flag 0) (kk 0) (ll 0) (nm 0)
        (c 0d0) (f 0d0) (h 0d0) (s 0d0) (x 0d0) (y 0d0) (z 0d0)
        (anorm 0d0) (g 0d0) (scale 0d0))
    (declare (type real-vector rv1) (type fixnum i1 j1 flag kk ll nm)
             (type double-float c f h s x y z anorm g scale))
    (fast
      ;; Householder reduction to bidiagonal form
      (setf i1 0)
      (do () ((> i1 n))
        (setf ll (1+ i1))
        (setf (aref rv1 i1) (* scale g))
        (setf g 0d0 s 0d0 scale 0d0)
        (when (<= i1 m)
          (setf kk i1)
          (do () ((> kk m)) (setf scale (+ scale (abs (aref a kk i1)))) (incf kk))
          (when (/= scale 0d0)
            (setf kk i1)
            (do () ((> kk m))
              (setf (aref a kk i1) (/ (aref a kk i1) scale))
              (setf s (+ s (* (aref a kk i1) (aref a kk i1))))
              (incf kk))
            (setf f (aref a i1 i1))
            (setf g (- (nl-siab (c-sqrt s) f)))
            (setf h (- (* f g) s))
            (setf (aref a i1 i1) (- f g))
            (when (/= i1 n)
              (setf j1 ll)
              (do () ((> j1 n))
                (setf s 0d0)
                (setf kk i1)
                (do () ((> kk m)) (setf s (+ s (* (aref a kk i1) (aref a kk j1)))) (incf kk))
                (setf f (/ s h))
                (setf kk i1)
                (do () ((> kk m)) (setf (aref a kk j1) (+ (aref a kk j1) (* f (aref a kk i1)))) (incf kk))
                (incf j1)))
            (setf kk i1)
            (do () ((> kk m)) (setf (aref a kk i1) (* (aref a kk i1) scale)) (incf kk))))
        (setf (aref w i1) (* scale g))
        (setf g 0d0 s 0d0 scale 0d0)
        (when (and (<= i1 m) (/= i1 n))
          (setf kk ll)
          (do () ((> kk n)) (setf scale (+ scale (abs (aref a i1 kk)))) (incf kk))
          (when (/= scale 0d0)
            (setf kk ll)
            (do () ((> kk n))
              (setf (aref a i1 kk) (/ (aref a i1 kk) scale))
              (setf s (+ s (* (aref a i1 kk) (aref a i1 kk))))
              (incf kk))
            (setf f (aref a i1 ll))
            (setf g (- (nl-siab (c-sqrt s) f)))
            (setf h (- (* f g) s))
            (setf (aref a i1 ll) (- f g))
            (setf kk ll)
            (do () ((> kk n)) (setf (aref rv1 kk) (/ (aref a i1 kk) h)) (incf kk))
            (when (/= i1 m)
              (setf j1 ll)
              (do () ((> j1 m))
                (setf s 0d0)
                (setf kk ll)
                (do () ((> kk n)) (setf s (+ s (* (aref a j1 kk) (aref a i1 kk)))) (incf kk))
                (setf kk ll)
                (do () ((> kk n)) (setf (aref a j1 kk) (+ (aref a j1 kk) (* s (aref rv1 kk)))) (incf kk))
                (incf j1)))
            (setf kk ll)
            (do () ((> kk n)) (setf (aref a i1 kk) (* (aref a i1 kk) scale)) (incf kk))))
        (setf anorm (nl-max anorm (+ (abs (aref w i1)) (abs (aref rv1 i1)))))
        (incf i1))
      ;; Accumulation of right hand transformations
      (setf i1 n ll i1)
      (do () ((< i1 0))
        (when (< i1 n)
          (when (/= g 0d0)
            (setf j1 ll)
            (do () ((> j1 n))
              (setf (aref v j1 i1) (/ (/ (aref a i1 j1) (aref a i1 ll)) g))
              (incf j1))
            (setf j1 ll)
            (do () ((> j1 n))
              (setf s 0d0)
              (setf kk ll)
              (do () ((> kk n)) (setf s (+ s (* (aref a i1 kk) (aref v kk j1)))) (incf kk))
              (setf kk ll)
              (do () ((> kk n)) (setf (aref v kk j1) (+ (aref v kk j1) (* s (aref v kk i1)))) (incf kk))
              (incf j1)))
          (setf j1 ll)
          (do () ((> j1 n))
            (setf (aref v j1 i1) 0d0)
            (setf (aref v i1 j1) 0d0)
            (incf j1)))
        (setf (aref v i1 i1) 1d0)
        (setf g (aref rv1 i1))
        (setf ll i1)
        (decf i1))
      ;; Accumulation of left hand transformations
      (setf i1 n)
      (do () ((< i1 0))
        (setf ll (1+ i1))
        (setf g (aref w i1))
        (when (< i1 n)
          (setf j1 ll)
          (do () ((> j1 n)) (setf (aref a i1 j1) 0d0) (incf j1)))
        (if (/= g 0d0)
            (progn
              (setf g (/ 1d0 g))
              (when (/= i1 n)
                (setf j1 ll)
                (do () ((> j1 n))
                  (setf s 0d0)
                  (setf kk ll)
                  (do () ((> kk m)) (setf s (+ s (* (aref a kk i1) (aref a kk j1)))) (incf kk))
                  (setf f (* (/ s (aref a i1 i1)) g))
                  (setf kk i1)
                  (do () ((> kk m)) (setf (aref a kk j1) (+ (aref a kk j1) (* f (aref a kk i1)))) (incf kk))
                  (incf j1)))
              (setf j1 i1)
              (do () ((> j1 m)) (setf (aref a j1 i1) (* (aref a j1 i1) g)) (incf j1)))
            (progn
              (setf j1 i1)
              (do () ((> j1 m)) (setf (aref a j1 i1) 0d0) (incf j1))))
        (setf (aref a i1 i1) (+ (aref a i1 i1) 1d0))
        (decf i1))
      ;; Diagonalization of the bidiagonal form
      (setf kk n)
      (do () ((< kk 0))
        (block iterations
          (do ((its 1 (1+ its))) ((> its +real-matrix-single-value-decompose--mxits+))
            (declare (type fixnum its))
            (setf flag 1)
            ;; test for splitting; rv1[0] is always 0
            (setf ll kk)
            (block split
              (do () ((< ll 0))
                (setf nm (1- ll))
                (when (= (+ (abs (aref rv1 ll)) anorm) anorm)
                  (setf flag 0)
                  (return-from split))
                (when (= (+ (abs (aref w nm)) anorm) anorm)
                  (return-from split))
                (decf ll)))
            (when (/= flag 0)
              (setf c 0d0 s 1d0)
              (setf i1 ll)
              (do () ((> i1 kk))
                (setf f (* s (aref rv1 i1)))
                (when (/= (+ (abs f) anorm) anorm)
                  (setf g (aref w i1))
                  (setf h (c-sqrt (+ (* f f) (* g g))))
                  (setf (aref w i1) h)
                  (setf h (/ 1d0 h))
                  (setf c (* g h))
                  (setf s (* (- f) h))
                  (setf j1 0)
                  (do () ((> j1 m))
                    (setf y (aref a j1 nm))
                    (setf z (aref a j1 i1))
                    (setf (aref a j1 nm) (+ (* y c) (* z s)))
                    (setf (aref a j1 i1) (- (* z c) (* y s)))
                    (incf j1)))
                (incf i1)))
            (setf z (aref w kk))
            (when (= ll kk)
              (when (< z 0d0)
                (setf (aref w kk) (- z))
                (setf j1 0)
                (do () ((> j1 n)) (setf (aref v j1 kk) (- (aref v j1 kk))) (incf j1)))
              (return-from iterations))
            (when (= its +real-matrix-single-value-decompose--mxits+)
              (nurbs-error :con-err 'real-matrix-single-value-decompose))
            (setf x (aref w ll))
            (setf nm (1- kk))
            (setf y (aref w nm))
            (setf g (aref rv1 nm))
            (setf h (aref rv1 kk))
            (setf f (/ (+ (* (- y z) (+ y z)) (* (- g h) (+ g h))) (* (* 2d0 h) y)))
            (setf g (c-sqrt (+ (* f f) 1d0)))
            (setf f (/ (+ (* (- x z) (+ x z)) (* h (- (/ y (+ f (nl-siab g f))) h))) x))
            (setf c 1d0 s 1d0)
            (setf j1 ll)
            (do () ((> j1 nm))
              (setf i1 (1+ j1))
              (setf g (aref rv1 i1))
              (setf y (aref w i1))
              (setf h (* s g))
              (setf g (* c g))
              (setf z (c-sqrt (+ (* f f) (* h h))))
              (setf (aref rv1 j1) z)
              (setf c (/ f z))
              (setf s (/ h z))
              (setf f (+ (* x c) (* g s)))
              (setf g (- (* g c) (* x s)))
              (setf h (* y s))
              (setf y (* y c))
              (dotimes (jj (1+ n))
                (setf x (aref v jj j1))
                (setf z (aref v jj i1))
                (setf (aref v jj j1) (+ (* x c) (* z s)))
                (setf (aref v jj i1) (- (* z c) (* x s))))
              (setf z (c-sqrt (+ (* f f) (* h h))))
              (setf (aref w j1) z)
              (when (/= z 0d0)
                (setf z (/ 1d0 z))
                (setf c (* f z))
                (setf s (* h z)))
              (setf f (+ (* c g) (* s y)))
              (setf x (- (* c y) (* s g)))
              (dotimes (jj (1+ m))
                (setf y (aref a jj j1))
                (setf z (aref a jj i1))
                (setf (aref a jj j1) (+ (* y c) (* z s)))
                (setf (aref a jj i1) (- (* z c) (* y s))))
              (incf j1))
            (setf (aref rv1 ll) 0d0)
            (setf (aref rv1 kk) f)
            (setf (aref w kk) x)))
        (decf kk))))
  a)

;;; ------------------------------------------------------------- minimization

(defconstant +func-find-minima--mxits+ 100 "MXITS in N_FuncFindMinima.")

(defun func-find-minima (xa xb xc fb f xtol ftol)
  "Minimize the function F (double-float to double-float) by Brent's method,
the minimum bracketed by XA, XB, XC with FB = F(XB).  Stops when the bracket
is within the fractional tolerance XTOL or F falls to FTOL or below.  Returns
(values xmin fmin).  NLib: N_FuncFindMinima."
  (declare (type double-float xa xb xc fb xtol ftol) (type function f))
  (let* ((gold 0.3819660d0) (zeps 1.0d-14) (e 0d0) (d 0d0)
         (a (if (< xa xc) xa xc)) (b (if (< xa xc) xc xa))
         (x xb) (w xb) (v xb) (fx fb) (fw fb) (fv fb)
         (u 0d0) (fu 0d0) (xm 0d0) (tol1 0d0) (tol2 0d0)
         (p 0d0) (q 0d0) (r 0d0) (etemp 0d0))
    (declare (type double-float gold zeps e d a b x w v fx fw fv u fu xm tol1 tol2 p q r etemp))
    (do ((its 1 (1+ its))) ((> its +func-find-minima--mxits+))
      (declare (type fixnum its))
      (setf xm (* 0.5d0 (+ a b)))
      (setf tol1 (+ (* xtol (abs x)) zeps))
      (setf tol2 (* 2d0 tol1))
      (when (<= (abs (- x xm)) (- tol2 (* 0.5d0 (- b a))))
        (return-from func-find-minima (values x fx)))
      (if (> (abs e) tol1)
          (progn
            (setf r (* (- x w) (- fx fv)))
            (setf q (* (- x v) (- fx fw)))
            (setf p (- (* q (- x v)) (* r (- x w))))
            (setf q (* 2d0 (- q r)))
            (when (> q 0d0) (setf p (- p)))
            (setf q (abs q))
            (setf etemp e)
            (setf e d)
            (if (or (>= (abs p) (abs (* (* 0.5d0 q) etemp)))
                    (<= p (* q (- a x)))
                    (>= p (* q (- b x))))
                (progn
                  (setf e (if (>= x xm) (- a x) (- b x)))
                  (setf d (* gold e)))
                (progn
                  (setf d (/ p q))
                  (setf u (+ x d))
                  (when (or (< (- u a) tol2) (< (- b u) tol2))
                    (setf d (nl-siab tol1 (- xm x)))))))
          (progn
            (setf e (if (>= x xm) (- a x) (- b x)))
            (setf d (* gold e))))
      (setf u (if (>= (abs d) tol1) (+ x d) (+ x (nl-siab tol1 d))))
      (setf fu (the double-float (funcall f u)))
      (when (<= fu ftol)
        (return-from func-find-minima (values u fu)))
      (if (<= fu fx)
          (progn
            (if (>= u x) (setf a x) (setf b x))
            (setf v w w x x u)
            (setf fv fw fw fx fx fu))
          (progn
            (if (< u x) (setf a u) (setf b u))
            (cond ((or (<= fu fw) (= w x))
                   (setf v w w u fv fw fw fu))
                  ((or (<= fu fv) (= v x) (= v w))
                   (setf v u fv fu))))))
    (values x fx)))

;;; ------------------------------------------------------------------- L2-g3

(defun real-matrix-lst-sq-solve (rma rhs sol)
  "The least squares solution of RMA x = RHS, RMA having m+1 rows and n+1
columns, m > n (its storage read as a full matrix), through the normal
equations RMA^T RMA x = RMA^T RHS, solved by LU decomposition with partial
pivoting.  SOL, a real-vector of at least n+1 allocated by the caller,
receives x.  Signals :ind-err unless m > n, and the errors of the transpose,
product, decomposition and solver.  Returns SOL.
NLib: N_RealMatrixLstSqSolve."
  (declare (type rmatrix rma) (type real-vector rhs sol))
  (let ((m (rmatrix-n rma))
        (n (rmatrix-m rma))
        (rm (rmatrix-rm rma)))
    (declare (type fixnum m n))
    (when (<= m n)
      (nurbs-error :ind-err 'real-matrix-lst-sq-solve))
    ;; the transpose of RMA times RMA: the square (n+1)x(n+1) system
    (let ((rmt (make-rmatrix))
          (rmta (make-rmatrix)))
      (real-matrix-transpose rma rmt)
      (real-matrix-multiply rmt rma rmta)
      ;; the transpose times RHS: the new right hand side
      (let ((rm rm))
        (declare (type real-matrix rm))
        (fast
          (do ((ii 0 (1+ ii))) ((> ii n))
            (declare (type fixnum ii))
            (setf (aref sol ii) 0d0)
            (do ((kk 0 (1+ kk))) ((> kk m))
              (declare (type fixnum kk))
              (setf (aref sol ii) (+ (aref sol ii) (* (aref rm kk ii) (aref rhs kk))))))))
      ;; decompose, then forward/backward substitution
      (let ((ind (make-index-vector (1+ n))))
        (declare (type index-vector ind))
        (real-matrix-lu-decompose-pivot rmta ind)
        (real-matrix-for-back-pivot rmta ind sol))))
  sol)
