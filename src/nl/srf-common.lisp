;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfCommon.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun create-srf-corner-pts (p00 p10 p01 p11)
  "A fresh bilinear, non-rational surface through corner points P00, P10, P01
and P11 (Pij is the corner at u = i, v = j), on [0,1] x [0,1].
NLib: N_CreateSrfCornerPts."
  (declare (type point p00 p10 p01 p11))
  (let ((pw (make-cpoint-net 2 2))
        (u (make-real-vector 4))
        (v (make-real-vector 4))
        (cp (make-cpoint)))
    (declare (type cpoint-net pw) (type real-vector u v) (type cpoint cp))
    (flet ((put (p i j)
             (pt-to-c-pt p cp)
             (dotimes (k 4) (setf (aref pw i j k) (aref cp k)))))
      (put p00 0 0)
      (put p10 1 0)
      (put p01 0 1)
      (put p11 1 1))
    (setf (aref u 0) 0d0 (aref u 1) 0d0 (aref u 2) 1d0 (aref u 3) 1d0
          (aref v 0) 0d0 (aref v 1) 0d0 (aref v 2) 1d0 (aref v 3) 1d0)
    (%make-surface 1 1 u v pw)))
