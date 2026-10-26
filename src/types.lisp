;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_DataStruct.h, NL_Geometry.h, NL_Nurbsdef.h and
;;;; NL_Globals.cpp (usd-brep source/SMLib).  Modified: Lisp structures,
;;;; arrays and conditions in place of C structs, error flags and memory
;;;; stacks.  PORTING.md describes how each C type maps.

(in-package :batten)

;;; ------------------------------------------------------------------ types

(deftype real-vector () '(simple-array double-float (*)))
(deftype real-matrix () '(simple-array double-float (* *)))
(deftype point () '(simple-array double-float (3)))        ; NL_POINT, NL_VECTOR
(deftype cpoint () '(simple-array double-float (4)))       ; NL_CPOINT (wx wy wz w)
(deftype point-array () '(simple-array double-float (* 3)))   ; NL_POINT[], NL_EPOLYGON
(deftype cpoint-array () '(simple-array double-float (* 4)))  ; NL_CPOINT[], NL_CPOLYGON
(deftype point-net () '(simple-array double-float (* * 3)))   ; NL_POINT**, NL_ENET
(deftype cpoint-net () '(simple-array double-float (* * 4)))  ; NL_CPOINT**, NL_CNET
(deftype index-vector () '(simple-array fixnum (*)))         ; NL_INDEX[]
(deftype index-matrix () '(simple-array fixnum (* *)))       ; NL_INTEGER**, NL_INDEX**
(deftype real-3d-array () '(simple-array double-float (* * *))) ; NL_REAL***
(deftype index () '(integer 0 #.(1- array-dimension-limit)))

(defmacro fast (&body body)
  "BODY compiled for speed with safety off.  Used inside routines whose
arguments are already declared or checked, as NLib's C trusts its inputs."
  `(locally (declare (optimize (speed 3) (safety 0) (debug 0))) ,@body))

(defmacro define-inline (name lambda-list &body body)
  "DEFUN NAME, declaimed inline.  Allegro CL ignores INLINE declarations on
user functions -- each call stays a full call, boxing its double arguments
and results -- but honors compiler macros, so there a compiler macro expands
every call into the body as a lambda application."
  `(progn
     (declaim (inline ,name))
     (defun ,name ,lambda-list ,@body)
     #+allegro
     (define-compiler-macro ,name (&rest args)
       (list* '(lambda ,lambda-list ,@body) args))))

;;; Scratch-array sizes as named constants: Allegro stack-allocates a
;;; DYNAMIC-EXTENT array only when its size is a literal or a constant symbol,
;;; not a constant expression such as (1+ +maxdeg+) (measured: 176 bytes a
;;; call for that form, none for the others).
(defconstant +deg+1+ (1+ +maxdeg+))
(defconstant +deg+2+ (+ +maxdeg+ 2))
(defconstant +2deg+2+ (* 2 (1+ +maxdeg+)))
(defconstant +deg+1-squared+ (* (1+ +maxdeg+) (1+ +maxdeg+)))
(defconstant +der+1+ (1+ +maxder+))
(defconstant +der+1-squared+ (* (1+ +maxder+) (1+ +maxder+)))
(defconstant +der+1-deg+1+ (* (1+ +maxder+) (1+ +maxdeg+)))

(defmacro with-scratch-doubles ((var size bound) &body body)
  "Bind VAR to a scratch vector of at least SIZE doubles for BODY.  When SIZE
is at most BOUND (a constant), VAR is a stack-allocated vector of BOUND
elements: implementations allocate a DYNAMIC-EXTENT array on the stack only
when its size is constant (CCL and Allegro heap-allocate a variable-size one).
Larger requests fall back to the heap.  BODY is expanded twice; keep it a
call.  The contents start undefined, as C's stack arrays do."
  (let ((n (gensym "N")))
    `(let ((,n ,size))
       (declare (type fixnum ,n))
       (if (<= ,n ,bound)
           (let ((,var (make-array ,bound :element-type 'double-float)))
             (declare (dynamic-extent ,var))
             ,@body)
           (let ((,var (make-array ,n :element-type 'double-float)))
             ,@body)))))

(defmacro with-c-float-environment (&body body)
  "Run BODY with floating-point traps masked, as C runs: a division by zero
gives an infinity and an invalid operation a NaN instead of signaling, so a
ported routine computes the same bits as its C source on every input,
including the ones where NLib computes an infinity and tests for it
afterwards.  Batten's routines assume this environment; the exported API
establishes it.  Allegro CL, EGCL and CLISP offer no way to mask the traps
(CLISP has no infinities or NaNs at all): there an input on which NLib
divides by zero signals instead (every other result is the same)."
  #+sbcl `(sb-int:with-float-traps-masked (:overflow :invalid :divide-by-zero :inexact :underflow)
            ,@body)
  #+ccl (let ((mode (gensym)))
          `(let ((,mode (ccl:get-fpu-mode)))
             (unwind-protect
                  (progn (ccl:set-fpu-mode :overflow nil :underflow nil :division-by-zero nil
                                           :invalid nil :inexact nil)
                         ,@body)
               (apply #'ccl:set-fpu-mode ,mode))))
  #+ecl (let ((mask (gensym)))
          `(let ((,mask (ext:trap-fpe 'last t)))
             (unwind-protect
                  (progn (ext:trap-fpe t nil) ,@body)
               (ext:trap-fpe ,mask t))))
  #-(or sbcl ccl ecl) `(progn ,@body))

;;; ---------------------------------------------- constants not in constants.lisp

;;; NL_BIGD, NL_SMAD; NL_NOZ marks a 2D control point (z = +no-z+), NL_NOW a
;;; non-rational one (w = +no-w+).  They stay inside the coordinate arrays,
;;; as in NLib, so every ported routine reads them as its C source does.
(defconstant +bigd+ most-positive-double-float)
(defconstant +smad+ least-positive-normalized-double-float)
(defconstant +no-z+ most-positive-double-float)
(defconstant +no-w+ most-negative-double-float)
(defconstant +infinite+ most-positive-double-float)
(defconstant +undefined+ most-negative-double-float)

;;; NL_PascalTri: binomial coefficients through degree 17.
(defconstant +pascal-tri+
  (if (boundp '+pascal-tri+)
      (symbol-value '+pascal-tri+)
      (let ((a (make-array '(18 18) :element-type 'fixnum :initial-element 0)))
        (dotimes (k 18 a)
          (setf (aref a k 0) 1)
          (do ((j 1 (1+ j))) ((> j k))
            (setf (aref a k j) (+ (aref a (1- k) (1- j))
                                  (if (< j k) (aref a (1- k) j) 0)))))))
  "NLib's NL_PascalTri[18][18].")

;;; ----------------------------------------------------- NLib's global tolerances

;;; NL_PUBLIC variables in NL_Globals.cpp: ported routines read these.
(defvar *mtol* 1d-7 "NL_MTOL, model space tolerance.")
(defvar *ptol* 1d-12 "NL_PTOL, parameter space tolerance.")
(defvar *wtol* 1d-3 "NL_WTOL, weight tolerance.")
(defvar *ltol* 1d-9 "NL_LTOL, parallel line tolerance (radians).")
(defvar *ludt* 1d-9 "NL_LUDT, LU-decomposition tolerance.")
(defvar *zdtl* 1d-8 "NL_ZDTL, zero discriminant tolerance.")
(defvar *zctl* 1d-12 "NL_ZCTL, zero coefficient tolerance.")
(defvar *wmin* 1d-3 "NL_WMIN, minimum weight.")
(defvar *wmax* 100d0 "NL_WMAX, maximum weight.")
(defvar *dmax* 17 "NL_DMAX, maximum degree.")
(defvar *itlim* 50 "NL_ITLIM, maximum number of iterations.")
(defvar *recurselim* 20 "NL_RECURSELIM, maximum recursion levels.")
(defvar *ccplim* 1000 "NL_CCPLIM, maximum curve control points.")

;;; ------------------------------------------------------------ NLib's macros

(defmacro nl-min (x y) "NL_MIN: (x < y) ? x : y" (let ((a (gensym)) (b (gensym))) `(let ((,a ,x) (,b ,y)) (if (< ,a ,b) ,a ,b))))
(defmacro nl-max (x y) "NL_MAX: (x > y) ? x : y" (let ((a (gensym)) (b (gensym))) `(let ((,a ,x) (,b ,y)) (if (> ,a ,b) ,a ,b))))
(defmacro nl-min3 (x y z)
  "NL_MIN3" (let ((a (gensym)) (b (gensym)) (c (gensym)))
              `(let ((,a ,x) (,b ,y) (,c ,z)) (if (and (< ,a ,b) (< ,a ,c)) ,a (if (< ,b ,c) ,b ,c)))))
(defmacro nl-max3 (x y z)
  "NL_MAX3" (let ((a (gensym)) (b (gensym)) (c (gensym)))
              `(let ((,a ,x) (,b ,y) (,c ,z)) (if (and (> ,a ,b) (> ,a ,c)) ,a (if (> ,b ,c) ,b ,c)))))

(define-inline nl-sign (x)
  "NL_SIGN: 1.0 for positive X, else -1.0 -- zero included."
  (declare (type double-float x))
  (if (> x 0d0) 1d0 -1d0))
(define-inline nl-siab (a b)
  "NL_SIAB: |a| with the sign of b (b >= 0 counts as positive)."
  (declare (type double-float a b))
  (if (>= b 0d0) (abs a) (- (abs a))))

;;; ------------------------------------------------ C's libm on the real line

;;; CL's SQRT and EXPT go complex where C's sqrt and pow give a NaN, and a
;;; declaration such as (the (double-float 0d0) x) signals on a NaN that C
;;; carries on with.  Ported routines call these instead.

(defvar *c-nan*
  ;; Allegro has no trap masking (see with-c-float-environment): its bits,
  ;; the x86-64 default NaN, are written out instead.
  #+allegro (excl:shorts-to-double-float #xfff8 0 0 0)
  ;; nil where the implementation can neither mask the trap nor represent a
  ;; NaN (CLISP, EGCL): there %C-NAN signals instead
  #-allegro (ignore-errors
             (with-c-float-environment
               ;; read at load time so that no compiler folds (and traps on) 0/0
               (let ((zero (read-from-string "0d0"))) (/ zero zero))))
  "The NaN an invalid operation produces (0/0), as C gets from sqrt(-1), or
nil on an implementation without one.")

(defun %c-nan (operation &rest operands)
  "*C-NAN*, or where the implementation has none, the floating-point error C
would have masked."
  (or *c-nan*
      (error 'floating-point-invalid-operation :operation operation :operands operands)))

(define-inline c-sqrt (x)
  "C's sqrt of the double X: X itself for a NaN, a NaN for a negative X."
  (declare (type double-float x))
  (cond ((>= x 0d0) (sqrt (the (double-float 0d0) x)))
        ((/= x x) x)
        (t (%c-nan 'sqrt x))))

(define-inline c-pow (x y)
  "C's pow of the doubles X and Y, Y not a whole number: a NaN for a negative
X, where CL's EXPT would return a complex."
  (declare (type double-float x y))
  (if (< x 0d0)
      (%c-nan 'expt x y)
      (the double-float (expt x y))))

;;; ------------------------------------------------------------------ errors

(define-condition nurbs-error (error)
  ((code :initarg :code :reader nurbs-error-code)
   (routine :initarg :routine :reader nurbs-error-routine))
  (:report (lambda (c stream)
             (format stream "NURBS error ~s in ~s." (nurbs-error-code c) (nurbs-error-routine c))))
  (:documentation "Signaled where NLib sets NL_ERROR and returns NL_YES.  The code is
the NLib error constant's name as a keyword: NL_PAR_ERR is :par-err, NL_NUM_ERR
is :num-err, NL_GEO_ERR is :geo-err, and so on."))

(defun nurbs-error (code routine)
  "Signal NURBS-ERROR.  CODE is the NLib error constant's name as a keyword;
where the C returns NL_YES from an NL_QUIT without setting an error code, the
code is :quit."
  (error 'nurbs-error :code code :routine routine))

;;; --------------------------------------------------------- points and arrays

(define-inline make-point (&optional (x 0d0) (y 0d0) (z 0d0))
  (let ((p (make-array 3 :element-type 'double-float)))
    (setf (aref p 0) (float x 1d0) (aref p 1) (float y 1d0) (aref p 2) (float z 1d0))
    p))
(define-inline make-cpoint (&optional (wx 0d0) (wy 0d0) (wz 0d0) (w +no-w+))
  (let ((p (make-array 4 :element-type 'double-float)))
    (setf (aref p 0) (float wx 1d0) (aref p 1) (float wy 1d0)
          (aref p 2) (float wz 1d0) (aref p 3) (float w 1d0))
    p))

(defun make-real-vector (n) (make-array n :element-type 'double-float :initial-element 0d0))
(defun make-real-matrix (rows cols) (make-array (list rows cols) :element-type 'double-float :initial-element 0d0))
(defun make-point-array (n) (make-array (list n 3) :element-type 'double-float :initial-element 0d0))
(defun make-cpoint-array (n) (make-array (list n 4) :element-type 'double-float :initial-element 0d0))
(defun make-point-net (n m) (make-array (list n m 3) :element-type 'double-float :initial-element 0d0))
(defun make-cpoint-net (n m) (make-array (list n m 4) :element-type 'double-float :initial-element 0d0))
(defun make-index-vector (n) (make-array n :element-type 'fixnum :initial-element 0))

;;; ------------------------------------------------------------------ curves

;;; NL_CURVE: degree p, knot vector U[0..m], control polygon Pw[0..n] as an
;;; (n+1) x 4 array of homogeneous control points.
(defstruct (curve (:constructor %make-curve (degree knots cpoints)))
  (degree 0 :type fixnum)
  (knots (make-real-vector 0) :type real-vector)
  (cpoints (make-cpoint-array 0) :type cpoint-array))

(declaim (inline curve-n curve-m))
(defun curve-n (curve) "NLib's n: highest control point index." (1- (array-dimension (curve-cpoints curve) 0)))
(defun curve-m (curve) "NLib's m: highest knot index." (1- (length (curve-knots curve))))

;;; ----------------------------------------------------------------- surfaces

;;; NL_SURFACE: degrees p (u) and q (v), knot vectors U and V, control net
;;; Pw[0..n][0..m] as an (n+1) x (m+1) x 4 array.
(defstruct (surface (:constructor %make-surface (degree-u degree-v knots-u knots-v cpoints)))
  (degree-u 0 :type fixnum)
  (degree-v 0 :type fixnum)
  (knots-u (make-real-vector 0) :type real-vector)
  (knots-v (make-real-vector 0) :type real-vector)
  (cpoints (make-cpoint-net 0 0) :type cpoint-net))

(declaim (inline surface-n surface-m))
(defun surface-n (s) "NLib's n: highest control point index in u." (1- (array-dimension (surface-cpoints s) 0)))
(defun surface-m (s) "NLib's m: highest control point index in v." (1- (array-dimension (surface-cpoints s) 1)))

;;; ------------------------------------------------ curve and surface functions

;;; NL_CFUN: a scalar B-spline function -- degree p, knots U, values fu[0..n].
(defstruct (cfun (:constructor %make-cfun (degree knots vals)))
  (degree 0 :type fixnum)
  (knots (make-real-vector 0) :type real-vector)
  (vals (make-real-vector 0) :type real-vector))

;;; NL_SFUN: a scalar tensor-product function -- degrees p, q, knots U, V,
;;; values fuv[0..n][0..m].
(defstruct (sfun (:constructor %make-sfun (degree-u degree-v knots-u knots-v vals)))
  (degree-u 0 :type fixnum)
  (degree-v 0 :type fixnum)
  (knots-u (make-real-vector 0) :type real-vector)
  (knots-v (make-real-vector 0) :type real-vector)
  (vals (make-real-matrix 0 0) :type real-matrix))

;;; ------------------------------------------------------- small geometric objects

;;; NL_LINESEG and NL_PLANE.  NL_INTERVAL, NL_RECTANGLE and NL_MINMAXBOX are
;;; real-vectors of 2, 4 and 6 in the C field order (ul ur), (ul ur vb vt),
;;; (xl xr yb yt zn zf).
(defstruct (lineseg (:constructor %make-lineseg (p v bounded)))
  (p (make-point) :type point)
  (v (make-point) :type point)
  (bounded 0 :type fixnum))

(defstruct (plane (:constructor %make-plane (p n)))
  (p (make-point) :type point)
  (n (make-point) :type point))

;;; ------------------------------------------------------------ real matrices

;;; NL_MATRIXTYPE, an enum in NL_DataStruct.h.
(defconstant +mt-full+ 0)
(defconstant +mt-lowerleft+ 1)
(defconstant +mt-upperright+ 2)
(defconstant +mt-banded+ 3)

;;; NL_RMATRIX: highest row and column indexes N and M, the element storage RM
;;; (a real-matrix, or nil for C's NULL), the matrix type MTP and the
;;; bandwidth BW.  As in C, N and M are recorded apart from RM's dimensions: a
;;; matrix may be defined smaller than its storage.  A banded matrix stores
;;; rows of BW elements, the diagonal in column BW/2.  (make-rmatrix) is
;;; N_InitRealMatrix's NULL matrix.  An optional caller-allocated NL_RMATRIX
;;; output is an rmatrix argument; a NULL one is given storage as
;;; N_CheckMemRealMatrix does (%check-mem-real-matrix in nl/math.lisp).
(defstruct (rmatrix (:constructor make-rmatrix (&optional (n -1) (m -1) rm (mtp +mt-full+) (bw 0))))
  (n -1 :type fixnum)
  (m -1 :type fixnum)
  (rm nil :type (or null real-matrix))
  (mtp +mt-full+ :type fixnum)
  (bw 0 :type fixnum))

;;; ------------------------------------------------------ convenient construction

(defun make-curve (degree knots points &key weights two-d)
  "A NURBS curve from Euclidean POINTS (a sequence of 3-sequences), KNOTS and
DEGREE.  With WEIGHTS the curve is rational: control points (w*x w*y w*z w);
without, w is +no-w+.  With TWO-D, z is +no-z+."
  (let* ((n+1 (length points))
         (cpoints (make-cpoint-array n+1))
         (knots (map 'real-vector (lambda (k) (float k 1d0)) knots)))
    (unless (= (length knots) (+ n+1 degree 1))
      (nurbs-error :inp-err 'make-curve))
    (let ((i 0))
      (map nil (lambda (point weight)
                 (let ((w (if weights (float weight 1d0) 1d0)))
                   (setf (aref cpoints i 0) (* w (float (elt point 0) 1d0))
                         (aref cpoints i 1) (* w (float (elt point 1) 1d0))
                         (aref cpoints i 2) (if two-d +no-z+ (* w (float (elt point 2) 1d0)))
                         (aref cpoints i 3) (if weights w +no-w+))
                   (incf i)))
           points (or weights (make-list n+1 :initial-element 1))))
    (%make-curve degree knots cpoints)))
