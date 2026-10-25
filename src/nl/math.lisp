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
  (sqrt (the (double-float 0d0) (dot-2-pts a a))))
