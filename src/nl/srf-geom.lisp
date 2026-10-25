;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfGeom.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun srf-get-average-len (sur)
  "Estimate the average lengths of surface SUR in u and v from its control
net, and a bound on its area as the product of the longest control polygon
in each direction.  Returns (values a lu lv).  NLib: N_SrfGetAverageLen."
  (declare (type surface sur))
  (let* ((pw (surface-cpoints sur))
         (n (surface-n sur))
         (m (surface-m sur))
         (cp (make-cpoint))
         (prev (make-point))
         (cur (make-point))
         (maxu 0d0) (sumu 0d0) (maxv 0d0) (sumv 0d0))
    (declare (type cpoint-net pw) (type fixnum n m) (type cpoint cp) (type point prev cur)
             (type double-float maxu sumu maxv sumv)
             (dynamic-extent cp prev cur))
    (flet ((euclid (i j p)
             (declare (type fixnum i j) (type point p))
             (set-cpoint cp (aref pw i j 0) (aref pw i j 1) (aref pw i j 2) (aref pw i j 3))
             (c-pt-to-pt-euclid cp p)))
      (declare (inline euclid))
      (dotimes (j (1+ m))
        (euclid 0 j prev)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((i 1 (1+ i))) ((> i n))
            (declare (type fixnum i))
            (euclid i j cur)
            (incf len (dist-pt-pt cur prev))
            (rotatef prev cur))
          (incf sumu len)
          (when (> len maxu) (setf maxu len))))
      (dotimes (i (1+ n))
        (euclid i 0 prev)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((j 1 (1+ j))) ((> j m))
            (declare (type fixnum j))
            (euclid i j cur)
            (incf len (dist-pt-pt cur prev))
            (rotatef prev cur))
          (incf sumv len)
          (when (> len maxv) (setf maxv len)))))
    (values (* maxu maxv)
            (/ sumu (+ (float m 1d0) 1d0))
            (/ sumv (+ (float n 1d0) 1d0)))))
