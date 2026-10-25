;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvGeom.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------ level 2, group 6

(defun crv-closest-pt (cur p u0 top toc q &optional ucell)
  "Find the point of curve CUR closest to point P (inversion when P is on the
curve, projection when it is off it) by Newton iteration from the guess U0,
wrapping the parameter on a closed curve and clamping it on an open one.  TOP
is the point coincidence tolerance, TOC the zero cosine tolerance.  Writes the
best point found into Q and returns its parameter u.

Signals :con-err when the iteration fails or does not converge within *itlim*
steps, and what crv-derivs signals.  The C returns the best u and Q found even
then: Q is written as it improves, and so is UCELL, an optional real-vector of
1 receiving u, for callers that carry on after the error.
NLib: N_CrvClosestPt."
  (declare (type curve cur) (type point p q) (type double-float u0 top toc))
  (let* ((knots (curve-knots cur))
         (m (1- (length knots)))
         (closed (crv-is-closed cur))
         (k 0) (ucount 0)
         (itlim *itlim*)
         (u 0d0)
         (unew u0) (uold 0d0)
         (best-dis +bigd+)
         (num 0d0) (den 0d0) (der 0d0) (dot 0d0) (dis 0d0) (zco 0d0)
         (pp (make-point))
         (d (make-point-array 3))
         (d0 (make-point)) (d1 (make-point)) (d2 (make-point))
         (v (make-point)))
    (declare (type real-vector knots) (type fixnum m k ucount itlim)
             (type double-float u unew uold best-dis num den der dot dis zco)
             (dynamic-extent pp d d0 d1 d2 v))
    ;; P is passed by value in C; Q may be the same array here
    (copy-pt p pp)
    (flet ((con-err () (nurbs-error :con-err 'crv-closest-pt)))
      (do () ((>= k itlim))
        ;; derivatives
        (crv-derivs cur unew +left+ 2 d)
        (%row-pt d 0 d0)
        (%row-pt d 1 d1)
        (%row-pt d 2 d2)
        ;; prepare for the convergence test
        (vector-diff d0 pp v)
        (setf num (vector-dot d1 v))
        (setf der (vector-magnitude d1))
        ;; check 1: point coincidence
        (setf dis (vector-magnitude v))
        (when (< dis best-dis)
          (setf u unew)
          (when ucell (setf (aref (the real-vector ucell) 0) u))
          (copy-pt d0 q)
          (setf best-dis dis))
        (when (<= dis top) (return))
        ;; check 2: zero cosine
        (when (float-op-is-bad num (* der dis) +division+) (con-err))
        (setf zco (/ (abs num) (* der dis)))
        (when (and (<= zco toc) (= dis best-dis)) (return))
        ;; no convergence: a new parameter
        (setf dot (vector-dot d2 v))
        (setf den (+ dot (* der der)))
        (when (float-op-is-bad num den +division+) (con-err))
        (setf uold unew
              unew (- uold (/ num den)))
        ;; check 3: parameter range
        (let ((ua (aref knots 0)) (ub (aref knots m)))
          (declare (type double-float ua ub))
          (if closed
              (cond ((< unew ua)
                     (cond ((< (- ua unew) (* 100d0 (- ub ua)))
                            (do () ((not (< unew ua)))
                              (setf unew (+ (- ub ua) unew)))
                            (setf ucount 0))
                           (t
                            (when (>= ucount 2) (con-err))
                            (setf unew uold)
                            (incf ucount))))
                    ((> unew ub)
                     (cond ((< (- unew ub) (* 100d0 (- ub ua)))
                            (do () ((not (> unew ub)))
                              (setf unew (+ (- ua ub) unew)))
                            (setf ucount 0))
                           (t
                            (when (>= ucount 2) (con-err))
                            (setf unew uold)
                            (incf ucount))))
                    (t (setf ucount 0)))
              (progn
                (when (< unew ua) (setf unew ua))
                (when (> unew ub) (setf unew ub)))))
        ;; check 4: parameter change
        (when (and (<= (abs (* (- unew uold) der)) top) (= dis best-dis))
          (return))
        (incf k))
      (when (>= k itlim) (con-err))
      u)))
