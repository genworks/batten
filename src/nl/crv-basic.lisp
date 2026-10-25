;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(declaim (inline is-crv-rat))
(defun is-crv-rat (cur)
  "True when curve CUR is rational (its first control point has a weight).
NLib: N_IsCrvRat."
  (declare (type curve cur))
  (/= (aref (curve-cpoints cur) 0 3) +no-w+))

(defun crv-get-param-bounds (cur)
  "The parameter domain of CUR, (values U[0] U[m]).  NLib: N_CrvGetParamBounds."
  (declare (type curve cur))
  (let ((k (curve-knots cur)))
    (values (aref k 0) (aref k (1- (length k))))))
