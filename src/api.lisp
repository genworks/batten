;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: keyword arguments and fresh results over the ported
;;;; routines.

(in-package :batten)

(defun side-flag (side)
  (ecase side (:left +left+) (:right +right+)))

(defun find-span (knots p u &key (side :left))
  "The knot span of U in KNOTS for degree P (see BASIS-FIND-SPAN)."
  (basis-find-span knots p (float u 1d0) (side-flag side)))

(defun curve-rational-p (curve) (is-crv-rat curve))
(defun surface-rational-p (surface) (is-srf-rat surface))

(defun curve-point! (curve u result &key (side :left))
  "Write the point of CURVE at U into RESULT, a point, and return it."
  (crv-eval curve u (side-flag side) result))

(defun curve-point (curve u &key (side :left))
  "The point of CURVE at U, fresh."
  (crv-eval curve (float u 1d0) (side-flag side) (make-point)))

(defun curve-derivatives! (curve u der result &key (side :left))
  "Write the point of CURVE at U and its derivatives to order DER into
RESULT, a point-array of (1+ DER) rows, and return it."
  (crv-derivs curve u (side-flag side) der result))

(defun curve-derivatives (curve u der &key (side :left))
  "The point of CURVE at U and its derivatives to order DER, a fresh
point-array of (1+ DER) rows."
  (crv-derivs curve (float u 1d0) (side-flag side) der (make-point-array (1+ der))))

(defun surface-point! (surface u v result &key (u-side :left) (v-side :left))
  "Write the point of SURFACE at (U, V) into RESULT and return it."
  (srf-eval-pt surface u v (side-flag u-side) (side-flag v-side) result))

(defun surface-point (surface u v &key (u-side :left) (v-side :left))
  "The point of SURFACE at (U, V), fresh."
  (srf-eval-pt surface (float u 1d0) (float v 1d0) (side-flag u-side) (side-flag v-side)
               (make-point)))

(defun surface-derivatives! (surface u v udr vdr result &key (u-side :left) (v-side :left))
  "Write the partial derivatives of SURFACE at (U, V) to orders UDR, VDR into
RESULT, a point-net, and return it: RESULT[k][l] is d^(k+l)S/du^k dv^l."
  (srf-derivs surface u v (side-flag u-side) (side-flag v-side) +false+ udr vdr result))

(defun surface-derivatives (surface u v udr vdr &key (u-side :left) (v-side :left))
  "The partial derivatives of SURFACE at (U, V) to orders UDR, VDR, a fresh
point-net."
  (srf-derivs surface (float u 1d0) (float v 1d0) (side-flag u-side) (side-flag v-side)
              +false+ udr vdr (make-point-net (1+ udr) (1+ vdr))))
