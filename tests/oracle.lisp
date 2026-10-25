;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The oracle: NLib's own C functions, from an SMLib 10 build, through CFFI.
;;;; Every ported routine is checked against the routine it came from.

(defpackage :batten-oracle
  (:use :common-lisp :batten)
  (:export #:load-oracle #:with-c-curve
           #:c-find-span #:c-basis-functions #:c-basis-derivatives
           #:c-curve-point #:c-curve-derivatives))

(in-package :batten-oracle)

(defun load-oracle (&optional (path (or (uiop:getenv "BATTEN_ORACLE_LIB")
                                        "/projects/tmp-shots/smlib-port/out/smlib64.so")))
  "Load the SMLib 10 library whose NLib functions are the oracle."
  (cffi:load-foreign-library path))

;;; NL_INDEX is long, NL_FLAG and NL_DEGREE are short.
(cffi:defcstruct nl-cpolygon (n :long) (pw :pointer))
(cffi:defcstruct nl-knotvector (m :long) (u :pointer))
(cffi:defcstruct nl-curve (pol :pointer) (p :short) (knt :pointer))

(defconstant +nl-left+ 1)
(defconstant +nl-right+ 2)

(defun side-flag (side) (ecase side (:left +nl-left+) (:right +nl-right+)))

(cffi:defcfun ("N_BasisFindSpan" %find-span) :short
  (knt :pointer) (p :short) (u :double) (flg :short) (spn :pointer))
(cffi:defcfun ("N_BasisEval" %basis-eval) :short
  (knt :pointer) (p :short) (u :double) (flg :short) (n :pointer) (spn :pointer))
(cffi:defcfun ("N_BasisDerivs" %basis-derivs) :short
  (knt :pointer) (p :short) (u :double) (flg :short) (der :long) (nd :pointer) (spn :pointer))
(cffi:defcfun ("N_CrvEval" %crv-eval) :short
  (cur :pointer) (u :double) (flg :short) (c :pointer))
(cffi:defcfun ("N_CrvDerivs" %crv-derivs) :short
  (cur :pointer) (u :double) (flg :short) (der :long) (cd :pointer))

(defun copy-to-foreign (vector)
  (let* ((n (length vector))
         (ptr (cffi:foreign-alloc :double :count (max n 1))))
    (dotimes (i n ptr)
      (setf (cffi:mem-aref ptr :double i) (aref vector i)))))

(defmacro with-c-knots ((ptr knots) &body body)
  (let ((u (gensym "U")))
    `(let* ((,u (copy-to-foreign ,knots))
            (,ptr (cffi:foreign-alloc '(:struct nl-knotvector))))
       (unwind-protect
            (progn
              (setf (cffi:foreign-slot-value ,ptr '(:struct nl-knotvector) 'm) (1- (length ,knots))
                    (cffi:foreign-slot-value ,ptr '(:struct nl-knotvector) 'u) ,u)
              ,@body)
         (cffi:foreign-free ,ptr)
         (cffi:foreign-free ,u)))))

(defmacro with-c-curve ((ptr curve) &body body)
  "Bind PTR to a freshly built NL_CURVE copy of CURVE for BODY."
  (let ((c (gensym "C")) (pw (gensym "PW")) (pol (gensym "POL")) (knt (gensym "KNT")))
    `(let* ((,c ,curve)
            (,pw (copy-to-foreign (curve-cpoints ,c)))
            (,pol (cffi:foreign-alloc '(:struct nl-cpolygon)))
            (,ptr (cffi:foreign-alloc '(:struct nl-curve))))
       (with-c-knots (,knt (curve-knots ,c))
         (unwind-protect
              (progn
                (setf (cffi:foreign-slot-value ,pol '(:struct nl-cpolygon) 'n) (curve-n ,c)
                      (cffi:foreign-slot-value ,pol '(:struct nl-cpolygon) 'pw) ,pw
                      (cffi:foreign-slot-value ,ptr '(:struct nl-curve) 'pol) ,pol
                      (cffi:foreign-slot-value ,ptr '(:struct nl-curve) 'p) (curve-degree ,c)
                      (cffi:foreign-slot-value ,ptr '(:struct nl-curve) 'knt) ,knt)
                ,@body)
           (cffi:foreign-free ,ptr)
           (cffi:foreign-free ,pol)
           (cffi:foreign-free ,pw))))))

;;; Each c- function answers (values result error-p) where error-p is true
;;; when the C routine returned NL_YES.

(defun c-find-span (knots p u side)
  (with-c-knots (knt knots)
    (cffi:with-foreign-object (spn :long)
      (let ((err (%find-span knt p u (side-flag side) spn)))
        (values (cffi:mem-ref spn :long) (/= err 0))))))

(defun c-basis-functions (knots p u side)
  "(values n-vector span error-p)"
  (with-c-knots (knt knots)
    (cffi:with-foreign-objects ((n :double (1+ p)) (spn :long))
      (let ((err (%basis-eval knt p u (side-flag side) n spn))
            (out (make-array (1+ p) :element-type 'double-float)))
        (dotimes (i (1+ p)) (setf (aref out i) (cffi:mem-aref n :double i)))
        (values out (cffi:mem-ref spn :long) (/= err 0))))))

(defun c-basis-derivatives (knots p u side der)
  "(values nd-matrix span error-p)"
  (with-c-knots (knt knots)
    (let* ((rows (1+ der)) (cols (1+ p))
           (data (cffi:foreign-alloc :double :count (* rows cols)))
           (rowptrs (cffi:foreign-alloc :pointer :count rows)))
      (unwind-protect
           (cffi:with-foreign-object (spn :long)
             (dotimes (k rows)
               (setf (cffi:mem-aref rowptrs :pointer k)
                     (cffi:inc-pointer data (* 8 k cols))))
             (let ((err (%basis-derivs knt p u (side-flag side) der rowptrs spn))
                   (out (make-array (list rows cols) :element-type 'double-float)))
               (dotimes (k rows)
                 (dotimes (j cols)
                   (setf (aref out k j) (cffi:mem-aref data :double (+ (* k cols) j)))))
               (values out (cffi:mem-ref spn :long) (/= err 0))))
        (cffi:foreign-free rowptrs)
        (cffi:foreign-free data)))))

(defun c-curve-point (curve u side)
  (with-c-curve (cur curve)
    (cffi:with-foreign-object (pt :double 3)
      (let ((err (%crv-eval cur u (side-flag side) pt))
            (out (make-array 3 :element-type 'double-float)))
        (dotimes (i 3) (setf (aref out i) (cffi:mem-aref pt :double i)))
        (values out (/= err 0))))))

(defun c-curve-derivatives (curve u der side)
  (with-c-curve (cur curve)
    (cffi:with-foreign-object (cd :double (* 3 (1+ der)))
      (let ((err (%crv-derivs cur u (side-flag side) der cd))
            (out (make-array (list (1+ der) 3) :element-type 'double-float)))
        (dotimes (k (1+ der))
          (dotimes (i 3) (setf (aref out k i) (cffi:mem-aref cd :double (+ (* 3 k) i)))))
        (values out (/= err 0))))))
