;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The oracle: NLib's own C, built by tools/build-oracle.sh, through CFFI.
;;;; oracle-bindings.lisp declares C-<lisp name> for every routine; this file
;;;; loads the libraries and defines NLib's structs; marshal.lisp moves Batten
;;;; objects in and out of them.

(defpackage :batten-oracle
  (:use :common-lisp :batten)
  (:export #:load-oracle #:*oracle-dir*
           #:with-stacks #:with-foreign #:with-foreign-out #:with-scalars
           #:to-foreign #:from-foreign #:foreign-free-object
           #:c-flag-error-p))

(in-package :batten-oracle)

(defvar *oracle-dir*
  (or (uiop:getenv "BATTEN_ORACLE_DIR")
      (namestring (asdf:system-relative-pathname "batten" "oracle/")))
  "Where build-oracle.sh put libnlib-oracle.so and libnlib-oracle-shims.so:
$BATTEN_ORACLE_DIR, else oracle/ in the Batten checkout.")

(defun load-oracle (&optional (dir *oracle-dir*))
  "Load NLib and the shims from DIR."
  (let ((dir (uiop:ensure-directory-pathname dir)))
    (cffi:load-foreign-library (merge-pathnames "libnlib-oracle.so" dir))
    (cffi:load-foreign-library (merge-pathnames "libnlib-oracle-shims.so" dir))))

(defun c-flag-error-p (flag)
  "True when an NL_FLAG return value says the C routine failed (NL_YES)."
  (and (integerp flag) (/= flag 0)))

;;; NLib's structs, as laid out by the C compiler: NL_INDEX long, NL_FLAG and
;;; NL_DEGREE short, NL_REAL double.
(cffi:defcstruct nl-point (x :double) (y :double) (z :double))
(cffi:defcstruct nl-cpoint (x :double) (y :double) (z :double) (w :double))
(cffi:defcstruct nl-knotvector (m :long) (u :pointer))
(cffi:defcstruct nl-cpolygon (n :long) (pw :pointer))
(cffi:defcstruct nl-cnet (n :long) (m :long) (pw :pointer))
(cffi:defcstruct nl-curve (pol :pointer) (p :short) (knt :pointer))
(cffi:defcstruct nl-surface (net :pointer) (p :short) (q :short) (knu :pointer) (knv :pointer))
(cffi:defcstruct nl-cvalue (n :long) (fu :pointer))
(cffi:defcstruct nl-cfun (cvl :pointer) (p :short) (knt :pointer))
(cffi:defcstruct nl-svalue (n :long) (m :long) (fuv :pointer))
(cffi:defcstruct nl-sfun (svl :pointer) (p :short) (q :short) (knu :pointer) (knv :pointer))
(cffi:defcstruct nl-epolygon (n :long) (pts :pointer))
(cffi:defcstruct nl-enet (n :long) (m :long) (pts :pointer))
(cffi:defcstruct nl-interval (ul :double) (ur :double))
(cffi:defcstruct nl-rectangle (ul :double) (ur :double) (vb :double) (vt :double))
(cffi:defcstruct nl-minmaxbox (xl :double) (xr :double) (yb :double) (yt :double) (zn :double) (zf :double))
(cffi:defcstruct nl-lineseg (p (:struct nl-point)) (v (:struct nl-point)) (bounded :short))
(cffi:defcstruct nl-plane (p (:struct nl-point)) (n (:struct nl-point)))
(cffi:defcstruct nl-rmatrix (n :long) (m :long) (rm :pointer) (mtp :int) (bw :long))

;;; NL_STACKS is opaque here: 352 bytes in SMLib 10, allocated generously.
(defconstant +stacks-bytes+ 1024)

(defmacro with-stacks ((var) &body body)
  "Bind VAR to an initialized NL_STACKS for BODY; everything C allocates from it
is freed afterwards, so read results out inside BODY."
  `(let ((,var (cffi:foreign-alloc :uint8 :count +stacks-bytes+ :initial-element 0)))
     (unwind-protect
          (progn (c-init-nurbs ,var) ,@body)
       (c-end-nurbs ,var)
       (cffi:foreign-free ,var))))

(defmacro with-scalars (bindings &body body)
  "BINDINGS: ((var type) ...), each a fresh foreign cell (zeroed) for a C scalar
output argument; read with (cffi:mem-ref var type).  TYPE :pointer gives a cell
for an output pointer (NL_REAL **, NL_POINT ** ...), initially NULL."
  `(cffi:with-foreign-objects ,(mapcar (lambda (b) (list (first b) (second b))) bindings)
     ,@(mapcar (lambda (b) `(setf (cffi:mem-ref ,(first b) ,(second b)) ,(ecase (second b)
                                                                            ((:double) 0d0)
                                                                            ((:long :short :int :unsigned-long) 0)
                                                                            ((:pointer) '(cffi:null-pointer)))))
               bindings)
     ,@body))
