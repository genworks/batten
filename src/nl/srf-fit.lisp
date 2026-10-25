;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfFit.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(declaim (inline %srf-fit-net-cpt))
(defun %srf-fit-net-cpt (pw row col dst)
  "Copy control point PW[ROW][COL] of a cpoint-net into the cpoint DST; return DST."
  (declare (type cpoint-net pw) (type fixnum row col) (type cpoint dst))
  (set-cpoint dst (aref pw row col 0) (aref pw row col 1) (aref pw row col 2) (aref pw row col 3)))

(defun fit-srf-removal-boundary (sur r s f l dir br mr)
  "Update the knot removal error bounds BR (a real-matrix indexed by control
point row and column) and the maximum error MR for removing once the knot of
index R and multiplicity S of surface SUR in direction DIR (+udir+ or
+vdir+), over rows/columns F..L of the other direction.  T[r] must be an
interior knot with T[r] /= T[r+1] and S > 0; these are not checked.  MR equal
to +bigd+ is first reset to -1.  BR is updated in place; returns the new MR.
As in the C, which returns NL_NO even then, a DIR other than +udir+ or +vdir+
updates nothing and signals nothing.  NLib: N_FitSrfRemovalBoundary."
  (declare (type surface sur) (type fixnum r s f l dir) (type real-matrix br)
           (type double-float mr))
  (let* ((p (surface-degree-u sur))
         (q (surface-degree-v sur))
         (pw (surface-cpoints sur))
         (udir (= dir +udir+))
         (deg (if udir p q))
         (knots (if udir (surface-knots-u sur) (surface-knots-v sur)))
         (size (1+ (* 2 (nl-max p q))))
         (alf (make-real-vector size))
         (oma (make-real-vector size))
         (bet (make-real-vector size))
         (omb (make-real-vector size))
         (rw (make-array size))
         (a (make-cpoint))
         (pt (make-cpoint)))
    (declare (type fixnum p q deg size) (type cpoint-net pw) (type real-vector knots alf oma bet omb)
             (type simple-vector rw) (type cpoint a pt))
    (dotimes (k size) (setf (svref rw k) (make-cpoint)))
    (when (= mr +bigd+)
      (setf mr -1d0))
    (unless (or (= dir +udir+) (= dir +vdir+))
      ;; NL_ERROR(NL_CAL_ERR), but the routine returns NL_NO
      (return-from fit-srf-removal-boundary mr))
    (flet ((net (k other)
             ;; Pw[k][col] removing in u, Pw[row][k] removing in v
             (if udir
                 (%srf-fit-net-cpt pw k other pt)
                 (%srf-fit-net-cpt pw other k pt)))
           (rw (k) (the cpoint (svref rw k))))
      (declare (inline net rw))
      (let* ((first (- r deg))
             (last (- r s))
             (off (1- first))
             (i first)
             (j last)
             (del 0d0)
             (omd 0d0))
        (declare (type fixnum first last off i j) (type double-float del omd))
        ;; save some parameters
        (do () ((not (> (- j i) 0)))
          (setf (aref alf (- i first)) (/ (- (aref knots (+ i deg 1)) (aref knots i))
                                          (- (aref knots r) (aref knots i)))
                (aref oma (- i first)) (- 1d0 (aref alf (- i first)))
                (aref bet (- j first)) (/ (- (aref knots (+ j deg 1)) (aref knots j))
                                          (- (aref knots (+ j deg 1)) (aref knots r)))
                (aref omb (- j first)) (- 1d0 (aref bet (- j first))))
          (incf i)
          (decf j))
        (setf del (/ (- (aref knots r) (aref knots i))
                     (- (aref knots (+ i deg 1)) (aref knots i)))
              omd (- 1d0 del))
        ;; update the maximum error for the requested rows/columns
        (do ((other f (1+ other))) ((> other l))
          (declare (type fixnum other))
          (let ((ii 1)
                (jj (- last off))
                (dw 0d0))
            (declare (type fixnum ii jj) (type double-float dw))
            (setf i first
                  j last)
            (copy-c-pt (net off other) (rw 0))
            (copy-c-pt (net (1+ last) other) (rw (- (1+ last) off)))
            ;; new control points for one removal step
            (do () ((not (> (- j i) 0)))
              (combine-2-c-pts (aref alf (- i first)) (net i other)
                               (aref oma (- i first)) (rw (1- ii)) (rw ii))
              (combine-2-c-pts (aref bet (- j first)) (net j other)
                               (aref omb (- j first)) (rw (1+ jj)) (rw jj))
              (incf i)
              (decf j)
              (incf ii)
              (decf jj))
            ;; the error
            (if (< (- j i) 0)
                (setf dw (dist-cpt-cpt-homo (rw (1- ii)) (rw (1+ jj))))
                (progn
                  (combine-2-c-pts del (rw (1+ jj)) omd (rw (1- ii)) a)
                  (setf dw (dist-cpt-cpt-homo (net i other) a))))
            (if udir
                (incf (aref br r other) dw)
                (incf (aref br other r) dw))
            (when (> dw mr)
              (setf mr dw))))))
    mr))
