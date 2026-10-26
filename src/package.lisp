;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;
;;; The exported API is small and readable.  Every ported NLib routine also
;;; exists inside the package under its mechanical name (N_CrvEval is
;;; BATTEN::CRV-EVAL; PORTING.md gives the rule), with the C signature mapped
;;; as PORTING.md describes.

(defpackage :batten
  (:use :common-lisp)
  (:export
   ;; types
   #:real-vector #:real-matrix #:point #:cpoint #:point-array #:cpoint-array
   #:point-net #:cpoint-net
   #:make-point #:make-cpoint #:make-real-vector #:make-real-matrix
   #:make-point-array #:make-cpoint-array #:make-point-net #:make-cpoint-net
   #:+no-z+ #:+no-w+ #:+maxdeg+ #:+maxder+ #:+left+ #:+right+ #:+true+ #:+false+
   ;; conditions
   #:nurbs-error #:nurbs-error-code #:nurbs-error-routine #:*bypass-float-checks*
   ;; curves
   #:curve #:curve-p #:make-curve #:%make-curve
   #:curve-degree #:curve-knots #:curve-cpoints #:curve-n #:curve-m #:curve-rational-p
   #:curve-point #:curve-point! #:curve-derivatives #:curve-derivatives!
   ;; surfaces
   #:surface #:surface-p #:%make-surface
   #:surface-degree-u #:surface-degree-v #:surface-knots-u #:surface-knots-v
   #:surface-cpoints #:surface-n #:surface-m #:surface-rational-p
   #:surface-point #:surface-point! #:surface-derivatives #:surface-derivatives!
   ;; scalar functions
   #:cfun #:%make-cfun #:cfun-degree #:cfun-knots #:cfun-vals
   #:sfun #:%make-sfun #:sfun-degree-u #:sfun-degree-v #:sfun-knots-u #:sfun-knots-v #:sfun-vals
   ;; basis functions
   #:find-span
   ;; shared by curves and surfaces (src/api/generics.lisp)
   #:domain #:bounding-box #:copy #:transform #:translate #:insert-knots #:elevate-degree
   #:reparameterize #:split
   ;; tessellation
   #:tessellate-surface))
