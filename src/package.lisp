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
   ;; construction (src/api/construction.lisp)
   #:make-surface #:line-segment #:circular-arc #:elliptical-arc #:conic-arc
   #:interpolate-curve #:approximate-curve #:fit-arc #:fit-line
   #:interpolate-surface #:approximate-surface #:ruled-surface #:extruded-surface
   #:revolved-surface #:skinned-surface #:gordon-surface #:coons-surface #:sphere
   ;; curve queries (src/api/curve-queries.lisp)
   #:curve-tangent #:curve-curvature #:curve-frenet-frame #:curve-length #:curve-point-at-length
   #:curve-closest-point #:curve-project-points #:curve-polyline
   #:curve-closed-p #:curve-planar-p #:curve-line-p #:curve-degenerate-p #:curves-coincident-p
   #:curve-type
   ;; curve operations (src/api/curve-operations.lisp)
   #:remove-knots #:subcurve #:bezier-segments #:reverse-curve #:extend-curve #:offset-curve
   #:make-compatible #:curve-to-2d #:curve-to-3d #:curve-to-rational #:curve-to-polynomial
   ;; surfaces (src/api/surfaces.lisp)
   #:surface-normal #:surface-grid #:surface-closest-point #:isocurve #:subsurface
   #:bezier-patches #:offset-surface #:surface-swap-uv #:surface-closed-p #:surface-flat-p
   #:surface-type
   ;; geometry and scalar functions (src/api/geometry.lisp)
   #:distance #:dot #:cross #:normalize #:angle-between
   #:project-point-to-line #:project-point-to-plane
   #:intersect-lines #:intersect-line-plane #:intersect-planes
   #:make-cfun #:cfun-value #:cfun-derivatives #:cfun-derivative #:cfun* #:cfun+ #:cfun-
   #:make-sfun
   ;; tessellation
   #:tessellate-surface))
