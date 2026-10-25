;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_DataStruct.h, NL_Geometry.h and NL_Nurbsdef.h
;;;; (usd-brep source/SMLib/inc).  Modified: Lisp structures and conditions
;;;; in place of C structs, error flags and memory stacks.

(in-package :batten)

(deftype real-vector () '(simple-array double-float (*)))
(deftype real-matrix () '(simple-array double-float (* *)))
(deftype index () '(integer 0 #.(1- array-dimension-limit)))

(defmacro fast (&body body)
  "BODY compiled for speed with safety off.  Used inside routines whose
arguments are already declared or checked, as NLib's C trusts its inputs."
  `(locally (declare (optimize (speed 3) (safety 0) (debug 0))) ,@body))

;;; NL_MAXDEG and NL_MAXDER.
(defconstant +max-degree+ 17)
(defconstant +max-derivative+ 4)

;;; NL_NOZ and NL_NOW: a control point whose z is +no-z+ is 2D, and one whose
;;; w is +no-w+ is non-rational.  Kept as NLib has them, inside the
;;; coordinate arrays, so that every ported routine reads them as its C
;;; source does.
(defconstant +no-z+ most-positive-double-float)
(defconstant +no-w+ most-negative-double-float)

;;; NL_SMAD and NL_BIGD.
(defconstant +smallest-double+ least-positive-normalized-double-float)
(defconstant +biggest-double+ most-positive-double-float)

(define-condition nurbs-error (error)
  ((code :initarg :code :reader nurbs-error-code)
   (routine :initarg :routine :reader nurbs-error-routine))
  (:report (lambda (c stream)
             (format stream "NURBS error ~s in ~s."
                     (nurbs-error-code c) (nurbs-error-routine c))))
  (:documentation "Signaled where NLib sets NL_ERROR and returns NL_YES.  The
code names the NLib error: :parameter (NL_PAR_ERR), :numeric (NL_NUM_ERR),
:call (NL_CAL_ERR), :input (NL_INP_ERR), :max-derivative (NL_MXD_ERR)."))

(defun nurbs-error (code routine)
  (error 'nurbs-error :code code :routine routine))

;;; A NURBS curve: NL_CURVE with its control polygon and knot vector flattened.
;;; CPOINTS holds the homogeneous control points Pw[0..n] four doubles apiece,
;;; (wx wy wz w), with +no-z+ and +no-w+ as above.  KNOTS is U[0..m].
(defstruct (curve (:constructor %make-curve (degree knots cpoints)))
  (degree 0 :type fixnum)
  (knots (make-array 0 :element-type 'double-float) :type real-vector)
  (cpoints (make-array 0 :element-type 'double-float) :type real-vector))

(declaim (inline curve-n curve-m))
(defun curve-n (curve)
  "Highest control point index, NLib's n."
  (1- (floor (length (curve-cpoints curve)) 4)))
(defun curve-m (curve)
  "Highest knot index, NLib's m."
  (1- (length (curve-knots curve))))

(defun curve-rational-p (curve)
  "True when the curve's control points carry weights.  NLib: N_IsCrvRat."
  (/= (aref (curve-cpoints curve) 3) +no-w+))

(defun make-curve (degree knots points &key weights two-d)
  "A NURBS curve from Euclidean POINTS (a sequence of 3-sequences), KNOTS and
DEGREE.  With WEIGHTS the curve is rational and its homogeneous control points
are (w*x w*y w*z w); without, w is +no-w+.  With TWO-D, z is +no-z+."
  (let* ((n+1 (length points))
         (cpoints (make-array (* 4 n+1) :element-type 'double-float))
         (knots (map 'real-vector (lambda (k) (float k 1d0)) knots)))
    (unless (= (length knots) (+ n+1 degree 1))
      (nurbs-error :input 'make-curve))
    (let ((i 0))
      (map nil (lambda (point weight)
                 (let ((w (if weights (float weight 1d0) 1d0))
                       (base (* 4 i)))
                   (setf (aref cpoints base) (* w (float (elt point 0) 1d0))
                         (aref cpoints (+ base 1)) (* w (float (elt point 1) 1d0))
                         (aref cpoints (+ base 2)) (if two-d +no-z+ (* w (float (elt point 2) 1d0)))
                         (aref cpoints (+ base 3)) (if weights w +no-w+))
                   (incf i)))
           points (or weights (make-list n+1 :initial-element 1))))
    (%make-curve degree knots cpoints)))
