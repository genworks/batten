;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_Util.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun xyz-to-2d-c-pt-array (wx wy wz w n m)
  "A fresh (N+1) x (M+1) cpoint-net whose [i][j] control point is
(WX[i][j] WY[i][j] WZ[i][j] W[i][j]), from four real-matrices.  W may hold
+no-w+ for a non-rational net.
NLib: N_XYZTo2dCPtArray."
  (declare (type real-matrix wx wy wz w) (type fixnum n m))
  (let ((pw (make-cpoint-net (1+ n) (1+ m))))
    (declare (type cpoint-net pw))
    (fast
      (dotimes (i (1+ n))
        (dotimes (j (1+ m))
          (setf (aref pw i j 0) (aref wx i j)
                (aref pw i j 1) (aref wy i j)
                (aref pw i j 2) (aref wz i j)
                (aref pw i j 3) (aref w i j)))))
    pw))

(defun check-real-mat-storage (rma nr mc rname)
  "Signal :sto-err in routine RNAME unless real matrix RMA is defined with
highest indexes at least NR (rows) and MC (columns).  Returns NIL otherwise.
NLib: N_CheckRealMatStorage."
  (declare (type rmatrix rma) (type fixnum nr mc))
  (multiple-value-bind (n m) (get-max-index-real-matrix rma)
    (declare (type fixnum n m))
    (when (or (< n nr) (< m mc))
      (nurbs-error :sto-err rname))
    nil))

(defun file-close (fptr)
  "Close the stream FPTR, as NLib wraps fclose for streams it opened.
NLib: N_FileClose."
  (declare (type stream fptr))
  (close fptr)
  nil)

(defun err-clear ()
  "Reset NLib's error record to \"no error\".  Batten reports errors as
NURBS-ERROR conditions and keeps no error record, so there is nothing to
reset; callers that clear the C record after handling an error have
already handled the condition.
NLib: N_ErrClear."
  nil)
