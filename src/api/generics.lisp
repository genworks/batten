;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The operations curves and surfaces share, as generic functions: one name
;;;; for a curve's domain and a surface's, a curve's knots and a surface's.
;;;; The methods live with the rest of each section (curves.lisp,
;;;; surfaces.lisp ...).  A call costs a dispatch; the ! functions and the
;;;; port routines are there for inner loops.

(in-package :batten)

;;; What every API function accepts: any real for a parameter or a length,
;;; any sequence of three reals for a point or a vector (a Batten point is
;;; passed through, not copied).

(declaim (inline %real))
(defun %real (x)
  "X as a double-float."
  (float x 1d0))

(defun %point (p)
  "P, a point or any sequence of three reals, as a Batten point."
  (if (typep p 'point)
      p
      (make-point (float (elt p 0) 1d0) (float (elt p 1) 1d0) (float (elt p 2) 1d0))))

(defun %point-array (points)
  "POINTS, a sequence of points (as %point takes them), as a point-array."
  (let* ((n (length points)) (array (make-point-array n)) (i 0))
    (map nil #'(lambda (p)
                 (let ((p (%point p)))
                   (setf (aref array i 0) (aref p 0) (aref array i 1) (aref p 1) (aref array i 2) (aref p 2)))
                 (incf i))
         points)
    array))

(defgeneric domain (object)
  (:documentation "The parameter domain: (values u0 u1) for a curve, (values u0 u1 v0 v1)
for a surface."))

(defgeneric bounding-box (object)
  (:documentation "The axis-aligned box around the control points (which contains the
object): (values min-point max-point)."))

(defgeneric copy (object)
  (:documentation "A copy sharing nothing with OBJECT."))

(defgeneric transform (object matrix)
  (:documentation "A fresh copy of OBJECT with its control points transformed by MATRIX, a
4x4 real-matrix applied to homogeneous points."))

(defgeneric translate (object vector)
  (:documentation "A fresh copy of OBJECT moved by VECTOR."))

(defgeneric insert-knots (object knots &key direction)
  (:documentation "A fresh copy of OBJECT with KNOTS (a sequence, repeats allowed) inserted;
the shape is unchanged.  For a surface, DIRECTION is :u or :v."))

(defgeneric elevate-degree (object by &key direction)
  (:documentation "A fresh copy of OBJECT with its degree raised BY; the shape is unchanged.
For a surface, DIRECTION is :u or :v."))

(defgeneric reparameterize (object &rest bounds)
  (:documentation "A fresh copy of OBJECT with its domain mapped linearly onto BOUNDS:
u0 u1 for a curve, u0 u1 v0 v1 for a surface."))

(defgeneric split (object at &key direction)
  (:documentation "OBJECT cut at parameter AT into two fresh pieces, (values before after).
For a surface, DIRECTION is :u or :v."))
