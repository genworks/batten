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
