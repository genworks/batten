;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0

(defpackage :batten
  (:use :common-lisp)
  (:export
   ;; types and limits
   #:real-vector #:+max-degree+ #:+max-derivative+ #:+no-z+ #:+no-w+
   ;; conditions
   #:nurbs-error #:nurbs-error-code #:nurbs-error-routine
   #:*bypass-float-checks*
   ;; curves
   #:curve #:curve-p #:make-curve #:%make-curve
   #:curve-degree #:curve-knots #:curve-cpoints #:curve-n #:curve-m
   #:curve-rational-p
   ;; basis functions
   #:find-span #:basis-functions! #:basis-derivatives!
   ;; curve evaluation
   #:curve-point! #:curve-point #:curve-derivatives! #:curve-derivatives))
