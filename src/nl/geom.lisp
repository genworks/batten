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
