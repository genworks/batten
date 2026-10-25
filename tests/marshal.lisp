;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Moving Batten objects into NLib's C structs and back.
;;;;
;;;;   (with-foreign ((cur curve :curve) (pts points :point-array)) ...)
;;;;     binds foreign copies; every allocation is freed on exit.
;;;;   (with-foreign-out ((out :curve) (pt :point) (nd :real-matrix 3 4)) ...)
;;;;     binds empty foreign objects for C to fill.
;;;;   (from-foreign ptr :curve)  (from-foreign ptr :point-array n) ...
;;;;     copies a C object back into a fresh Batten object.
;;;;
;;;; Types and their C forms:
;;;;   :real-vector NL_REAL*          :index-vector NL_INDEX*
;;;;   :point NL_POINT*               :cpoint NL_CPOINT*
;;;;   :point-array NL_POINT*         :cpoint-array NL_CPOINT*
;;;;   :real-matrix NL_REAL**         :point-net NL_POINT**     :cpoint-net NL_CPOINT**
;;;;   :knots NL_KNOTVECTOR*          :curve NL_CURVE*          :surface NL_SURFACE*
;;;;   :cfun NL_CFUN*                 :sfun NL_SFUN*
;;;;   :epolygon NL_EPOLYGON*         :enet NL_ENET*            :cpolygon NL_CPOLYGON*
;;;;   :interval NL_INTERVAL*         :rectangle NL_RECTANGLE*  :box NL_MINMAXBOX*
;;;;   :lineseg NL_LINESEG*           :plane NL_PLANE*
;;;;   :basis-derivatives NL_BASISDERIVATIVES (fixed [5][5][18][18])
;;;;   :rmatrix NL_RMATRIX* (batten::rmatrix; nil storage is RM NULL)
;;;;   :cnet NL_CNET* (a cpoint-net)  :index-matrix NL_INTEGER**  :real-3d NL_REAL***
;;;; Dimensions after the type give sizes where C does not record them:
;;;;   (from-foreign p :point-array n) (from-foreign p :real-matrix rows cols)
;;;;   (from-foreign p :basis-derivatives udr+1 vdr+1 p+1 q+1)
;;;;   (from-foreign p :index-matrix rows cols) (from-foreign p :real-3d a b c)
;;;;   (from-foreign p :rmatrix [rows cols]) -- storage size; by default n+1 rows
;;;;     and m+1 columns, or bw columns (2 when bw is 1) for a banded matrix
;;;;   (empty-foreign :rmatrix) is N_InitRealMatrix's NULL matrix;
;;;;   (empty-foreign :rmatrix rows cols) has full storage, defined rows-1 x cols-1.

(in-package :batten-oracle)

(defvar *allocations* nil "Foreign pointers to free when the enclosing WITH-FOREIGN exits.")

(defun falloc (type count)
  (let ((p (cffi:foreign-alloc type :count (max count 1))))
    ;; zero the block: NLib's allocators hand out cleared memory too
    (dotimes (i (* (max count 1) (cffi:foreign-type-size type)))
      (setf (cffi:mem-aref p :uint8 i) 0))
    (push p *allocations*)
    p))

(defun put-doubles (ptr offset &rest values)
  (let ((i offset))
    (dolist (v values) (setf (cffi:mem-aref ptr :double i) v) (incf i))))

(defun fill-doubles (array)
  "A foreign copy of every element of a double ARRAY, in row-major order."
  (let* ((n (array-total-size array)) (p (falloc :double n)))
    (dotimes (i n p) (setf (cffi:mem-aref p :double i) (row-major-aref array i)))))

(defun row-pointers (data rows stride-bytes)
  (let ((rp (falloc :pointer rows)))
    (dotimes (r rows rp)
      (setf (cffi:mem-aref rp :pointer r) (cffi:inc-pointer data (* r stride-bytes))))))

(defun knots-to-foreign (knots)
  (let ((k (falloc '(:struct nl-knotvector) 1)))
    (setf (cffi:foreign-slot-value k '(:struct nl-knotvector) 'm) (1- (length knots))
          (cffi:foreign-slot-value k '(:struct nl-knotvector) 'u) (fill-doubles knots))
    k))

(defun net-to-foreign (net)
  "NL_CPOINT** / NL_POINT** over a (* * k) array: rows of contiguous points."
  (let* ((rows (array-dimension net 0)) (cols (array-dimension net 1)) (k (array-dimension net 2)))
    (row-pointers (fill-doubles net) rows (* cols k 8))))

(defun index-matrix-to-foreign (object)
  (let* ((rows (array-dimension object 0)) (cols (array-dimension object 1))
         (data (falloc :long (* rows cols))))
    (dotimes (i (* rows cols)) (setf (cffi:mem-aref data :long i) (row-major-aref object i)))
    (row-pointers data rows (* 8 cols))))

(defun real-3d-to-foreign (object)
  "NL_REAL***: a pointer per first index to a NL_REAL** over the rest."
  (destructuring-bind (a b c) (array-dimensions object)
    (let ((data (fill-doubles object)) (top (falloc :pointer a)))
      (dotimes (i a top)
        (setf (cffi:mem-aref top :pointer i)
              (row-pointers (cffi:inc-pointer data (* i b c 8)) b (* c 8)))))))

(defun cnet-to-foreign (net)
  (let ((s (falloc '(:struct nl-cnet) 1)))
    (setf (cffi:foreign-slot-value s '(:struct nl-cnet) 'n) (1- (array-dimension net 0))
          (cffi:foreign-slot-value s '(:struct nl-cnet) 'm) (1- (array-dimension net 1))
          (cffi:foreign-slot-value s '(:struct nl-cnet) 'pw) (net-to-foreign net))
    s))

(defun rmatrix-to-foreign (object)
  (let ((s (falloc '(:struct nl-rmatrix) 1)) (rm (batten::rmatrix-rm object)))
    (setf (cffi:foreign-slot-value s '(:struct nl-rmatrix) 'n) (batten::rmatrix-n object)
          (cffi:foreign-slot-value s '(:struct nl-rmatrix) 'm) (batten::rmatrix-m object)
          (cffi:foreign-slot-value s '(:struct nl-rmatrix) 'rm)
          (if rm (to-foreign rm :real-matrix) (cffi:null-pointer))
          (cffi:foreign-slot-value s '(:struct nl-rmatrix) 'mtp) (batten::rmatrix-mtp object)
          (cffi:foreign-slot-value s '(:struct nl-rmatrix) 'bw) (batten::rmatrix-bw object))
    s))

(defun to-foreign (object type)
  "A foreign copy of OBJECT as C type TYPE (see the file header).  The memory is
released by the enclosing WITH-FOREIGN."
  (ecase type
    (:rmatrix (rmatrix-to-foreign object))
    (:cnet (cnet-to-foreign object))
    (:index-matrix (index-matrix-to-foreign object))
    (:real-3d (real-3d-to-foreign object))
    ((:real-vector :point :cpoint :point-array :cpoint-array) (fill-doubles object))
    (:index-vector (let ((p (falloc :long (length object))))
                     (dotimes (i (length object) p) (setf (cffi:mem-aref p :long i) (aref object i)))))
    (:real-matrix (row-pointers (fill-doubles object) (array-dimension object 0)
                                (* 8 (array-dimension object 1))))
    ((:point-net :cpoint-net) (net-to-foreign object))
    (:knots (knots-to-foreign object))
    (:curve
     (let ((pol (falloc '(:struct nl-cpolygon) 1)) (c (falloc '(:struct nl-curve) 1)))
       (setf (cffi:foreign-slot-value pol '(:struct nl-cpolygon) 'n) (curve-n object)
             (cffi:foreign-slot-value pol '(:struct nl-cpolygon) 'pw) (fill-doubles (curve-cpoints object))
             (cffi:foreign-slot-value c '(:struct nl-curve) 'pol) pol
             (cffi:foreign-slot-value c '(:struct nl-curve) 'p) (curve-degree object)
             (cffi:foreign-slot-value c '(:struct nl-curve) 'knt) (knots-to-foreign (curve-knots object)))
       c))
    (:surface
     (let ((net (falloc '(:struct nl-cnet) 1)) (s (falloc '(:struct nl-surface) 1)))
       (setf (cffi:foreign-slot-value net '(:struct nl-cnet) 'n) (surface-n object)
             (cffi:foreign-slot-value net '(:struct nl-cnet) 'm) (surface-m object)
             (cffi:foreign-slot-value net '(:struct nl-cnet) 'pw) (net-to-foreign (surface-cpoints object))
             (cffi:foreign-slot-value s '(:struct nl-surface) 'net) net
             (cffi:foreign-slot-value s '(:struct nl-surface) 'p) (surface-degree-u object)
             (cffi:foreign-slot-value s '(:struct nl-surface) 'q) (surface-degree-v object)
             (cffi:foreign-slot-value s '(:struct nl-surface) 'knu) (knots-to-foreign (surface-knots-u object))
             (cffi:foreign-slot-value s '(:struct nl-surface) 'knv) (knots-to-foreign (surface-knots-v object)))
       s))
    (:cfun
     (let ((cv (falloc '(:struct nl-cvalue) 1)) (f (falloc '(:struct nl-cfun) 1)))
       (setf (cffi:foreign-slot-value cv '(:struct nl-cvalue) 'n) (1- (length (cfun-vals object)))
             (cffi:foreign-slot-value cv '(:struct nl-cvalue) 'fu) (fill-doubles (cfun-vals object))
             (cffi:foreign-slot-value f '(:struct nl-cfun) 'cvl) cv
             (cffi:foreign-slot-value f '(:struct nl-cfun) 'p) (cfun-degree object)
             (cffi:foreign-slot-value f '(:struct nl-cfun) 'knt) (knots-to-foreign (cfun-knots object)))
       f))
    (:sfun
     (let ((sv (falloc '(:struct nl-svalue) 1)) (f (falloc '(:struct nl-sfun) 1))
           (vals (sfun-vals object)))
       (setf (cffi:foreign-slot-value sv '(:struct nl-svalue) 'n) (1- (array-dimension vals 0))
             (cffi:foreign-slot-value sv '(:struct nl-svalue) 'm) (1- (array-dimension vals 1))
             (cffi:foreign-slot-value sv '(:struct nl-svalue) 'fuv) (to-foreign vals :real-matrix)
             (cffi:foreign-slot-value f '(:struct nl-sfun) 'svl) sv
             (cffi:foreign-slot-value f '(:struct nl-sfun) 'p) (sfun-degree-u object)
             (cffi:foreign-slot-value f '(:struct nl-sfun) 'q) (sfun-degree-v object)
             (cffi:foreign-slot-value f '(:struct nl-sfun) 'knu) (knots-to-foreign (sfun-knots-u object))
             (cffi:foreign-slot-value f '(:struct nl-sfun) 'knv) (knots-to-foreign (sfun-knots-v object)))
       f))
    (:epolygon
     (let ((e (falloc '(:struct nl-epolygon) 1)))
       (setf (cffi:foreign-slot-value e '(:struct nl-epolygon) 'n) (1- (array-dimension object 0))
             (cffi:foreign-slot-value e '(:struct nl-epolygon) 'pts) (fill-doubles object))
       e))
    (:enet
     (let ((e (falloc '(:struct nl-enet) 1)))
       (setf (cffi:foreign-slot-value e '(:struct nl-enet) 'n) (1- (array-dimension object 0))
             (cffi:foreign-slot-value e '(:struct nl-enet) 'm) (1- (array-dimension object 1))
             (cffi:foreign-slot-value e '(:struct nl-enet) 'pts) (net-to-foreign object))
       e))
    (:cpolygon
     (let ((e (falloc '(:struct nl-cpolygon) 1)))
       (setf (cffi:foreign-slot-value e '(:struct nl-cpolygon) 'n) (1- (array-dimension object 0))
             (cffi:foreign-slot-value e '(:struct nl-cpolygon) 'pw) (fill-doubles object))
       e))
    ((:interval :rectangle :box) (fill-doubles object))
    (:lineseg
     (let ((l (falloc '(:struct nl-lineseg) 1)) (p (batten::lineseg-p object)) (v (batten::lineseg-v object)))
       (put-doubles l 0 (aref p 0) (aref p 1) (aref p 2) (aref v 0) (aref v 1) (aref v 2))
       (setf (cffi:foreign-slot-value l '(:struct nl-lineseg) 'bounded) (batten::lineseg-bounded object))
       l))
    (:plane
     (let ((l (falloc '(:struct nl-plane) 1)) (p (batten::plane-p object)) (n (batten::plane-n object)))
       (put-doubles l 0 (aref p 0) (aref p 1) (aref p 2) (aref n 0) (aref n 1) (aref n 2))
       l))
    (:basis-derivatives
     (let ((bd (falloc :double (* 5 5 18 18))))
       (destructuring-bind (a b c d) (array-dimensions object)
         (dotimes (k a bd)
           (dotimes (l b)
             (dotimes (i c)
               (dotimes (j d)
                 (setf (cffi:mem-aref bd :double (+ (* (+ (* (+ (* k 5) l) 18) i) 18) j))
                       (aref object k l i j)))))))))))

(defun empty-foreign (type &rest dims)
  "An empty, zeroed foreign object of TYPE for C to fill (see WITH-FOREIGN-OUT)."
  (ecase type
    (:rmatrix (if dims
                  (rmatrix-to-foreign (batten::make-rmatrix (1- (first dims)) (1- (second dims))
                                                            (make-real-matrix (first dims) (second dims))
                                                            batten::+mt-full+ (first dims)))
                  (rmatrix-to-foreign (batten::make-rmatrix))))
    (:cnet (if dims
               (cnet-to-foreign (make-cpoint-net (first dims) (second dims)))
               (falloc '(:struct nl-cnet) 1)))
    (:index-matrix (row-pointers (falloc :long (* (first dims) (second dims))) (first dims) (* 8 (second dims))))
    (:real-3d (real-3d-to-foreign (make-array dims :element-type 'double-float :initial-element 0d0)))
    (:point (falloc :double 3))
    (:cpoint (falloc :double 4))
    (:point-array (falloc :double (* 3 (first dims))))
    (:cpoint-array (falloc :double (* 4 (first dims))))
    (:real-vector (falloc :double (first dims)))
    (:index-vector (falloc :long (first dims)))
    (:real-matrix (row-pointers (falloc :double (* (first dims) (second dims))) (first dims) (* 8 (second dims))))
    (:point-net (row-pointers (falloc :double (* 3 (first dims) (second dims))) (first dims) (* 24 (second dims))))
    (:cpoint-net (row-pointers (falloc :double (* 4 (first dims) (second dims))) (first dims) (* 32 (second dims))))
    (:knots (falloc '(:struct nl-knotvector) 1))
    (:curve (falloc '(:struct nl-curve) 1))
    (:surface (falloc '(:struct nl-surface) 1))
    (:cfun (falloc '(:struct nl-cfun) 1))
    (:sfun (falloc '(:struct nl-sfun) 1))
    (:epolygon (falloc '(:struct nl-epolygon) 1))
    (:enet (falloc '(:struct nl-enet) 1))
    (:cpolygon (falloc '(:struct nl-cpolygon) 1))
    (:interval (falloc :double 2))
    (:rectangle (falloc :double 4))
    (:box (falloc :double 6))
    (:lineseg (falloc '(:struct nl-lineseg) 1))
    (:plane (falloc '(:struct nl-plane) 1))
    (:basis-derivatives (falloc :double (* 5 5 18 18)))))

(defun read-doubles (ptr dims)
  (let ((a (make-array dims :element-type 'double-float)))
    (dotimes (i (array-total-size a) a)
      (setf (row-major-aref a i) (cffi:mem-aref ptr :double i)))))

(defun read-net (rows-ptr rows cols k)
  (let ((a (make-array (list rows cols k) :element-type 'double-float)))
    (dotimes (r rows a)
      (let ((row (cffi:mem-aref rows-ptr :pointer r)))
        (dotimes (c cols)
          (dotimes (j k)
            (setf (aref a r c j) (cffi:mem-aref row :double (+ (* c k) j)))))))))

(defun read-knots (k)
  (read-doubles (cffi:foreign-slot-value k '(:struct nl-knotvector) 'u)
                (list (1+ (cffi:foreign-slot-value k '(:struct nl-knotvector) 'm)))))

(defun from-foreign (ptr type &rest dims)
  "A fresh Batten object copied from C object PTR of TYPE (see the file header)."
  (ecase type
    (:rmatrix
     (let* ((n (cffi:foreign-slot-value ptr '(:struct nl-rmatrix) 'n))
            (m (cffi:foreign-slot-value ptr '(:struct nl-rmatrix) 'm))
            (rm (cffi:foreign-slot-value ptr '(:struct nl-rmatrix) 'rm))
            (mtp (cffi:foreign-slot-value ptr '(:struct nl-rmatrix) 'mtp))
            (bw (cffi:foreign-slot-value ptr '(:struct nl-rmatrix) 'bw))
            (rows (or (first dims) (1+ n)))
            (cols (or (second dims) (if (= mtp batten::+mt-banded+) (if (= bw 1) 2 bw) (1+ m)))))
       (batten::make-rmatrix n m (if (cffi:null-pointer-p rm) nil (from-foreign rm :real-matrix rows cols))
                             mtp bw)))
    (:cnet
     (let ((n (cffi:foreign-slot-value ptr '(:struct nl-cnet) 'n))
           (m (cffi:foreign-slot-value ptr '(:struct nl-cnet) 'm)))
       (read-net (cffi:foreign-slot-value ptr '(:struct nl-cnet) 'pw) (1+ n) (1+ m) 4)))
    (:index-matrix
     (let ((a (make-array dims :element-type 'fixnum)))
       (dotimes (r (first dims) a)
         (let ((row (cffi:mem-aref ptr :pointer r)))
           (dotimes (c (second dims))
             (setf (aref a r c) (cffi:mem-aref row :long c)))))))
    (:real-3d
     (destructuring-bind (a b c) dims
       (let ((out (make-array dims :element-type 'double-float)))
         (dotimes (i a out)
           (let ((rows (cffi:mem-aref ptr :pointer i)))
             (dotimes (j b)
               (let ((row (cffi:mem-aref rows :pointer j)))
                 (dotimes (k c)
                   (setf (aref out i j k) (cffi:mem-aref row :double k))))))))))
    (:point (read-doubles ptr '(3)))
    (:cpoint (read-doubles ptr '(4)))
    (:point-array (read-doubles ptr (list (first dims) 3)))
    (:cpoint-array (read-doubles ptr (list (first dims) 4)))
    (:real-vector (read-doubles ptr (list (first dims))))
    ((:interval) (read-doubles ptr '(2)))
    ((:rectangle) (read-doubles ptr '(4)))
    ((:box) (read-doubles ptr '(6)))
    (:index-vector (let ((a (make-array (first dims) :element-type 'fixnum)))
                     (dotimes (i (first dims) a) (setf (aref a i) (cffi:mem-aref ptr :long i)))))
    (:real-matrix (let ((a (make-array dims :element-type 'double-float)))
                    (dotimes (r (first dims) a)
                      (let ((row (cffi:mem-aref ptr :pointer r)))
                        (dotimes (c (second dims))
                          (setf (aref a r c) (cffi:mem-aref row :double c)))))))
    (:point-net (read-net ptr (first dims) (second dims) 3))
    (:cpoint-net (read-net ptr (first dims) (second dims) 4))
    (:knots (read-knots ptr))
    (:curve
     (let* ((pol (cffi:foreign-slot-value ptr '(:struct nl-curve) 'pol))
            (n (cffi:foreign-slot-value pol '(:struct nl-cpolygon) 'n)))
       (%make-curve (cffi:foreign-slot-value ptr '(:struct nl-curve) 'p)
                    (read-knots (cffi:foreign-slot-value ptr '(:struct nl-curve) 'knt))
                    (read-doubles (cffi:foreign-slot-value pol '(:struct nl-cpolygon) 'pw) (list (1+ n) 4)))))
    (:surface
     (let* ((net (cffi:foreign-slot-value ptr '(:struct nl-surface) 'net))
            (n (cffi:foreign-slot-value net '(:struct nl-cnet) 'n))
            (m (cffi:foreign-slot-value net '(:struct nl-cnet) 'm)))
       (%make-surface (cffi:foreign-slot-value ptr '(:struct nl-surface) 'p)
                      (cffi:foreign-slot-value ptr '(:struct nl-surface) 'q)
                      (read-knots (cffi:foreign-slot-value ptr '(:struct nl-surface) 'knu))
                      (read-knots (cffi:foreign-slot-value ptr '(:struct nl-surface) 'knv))
                      (read-net (cffi:foreign-slot-value net '(:struct nl-cnet) 'pw) (1+ n) (1+ m) 4))))
    (:cfun
     (let* ((cv (cffi:foreign-slot-value ptr '(:struct nl-cfun) 'cvl))
            (n (cffi:foreign-slot-value cv '(:struct nl-cvalue) 'n)))
       (%make-cfun (cffi:foreign-slot-value ptr '(:struct nl-cfun) 'p)
                   (read-knots (cffi:foreign-slot-value ptr '(:struct nl-cfun) 'knt))
                   (read-doubles (cffi:foreign-slot-value cv '(:struct nl-cvalue) 'fu) (list (1+ n))))))
    (:sfun
     (let* ((sv (cffi:foreign-slot-value ptr '(:struct nl-sfun) 'svl))
            (n (cffi:foreign-slot-value sv '(:struct nl-svalue) 'n))
            (m (cffi:foreign-slot-value sv '(:struct nl-svalue) 'm)))
       (%make-sfun (cffi:foreign-slot-value ptr '(:struct nl-sfun) 'p)
                   (cffi:foreign-slot-value ptr '(:struct nl-sfun) 'q)
                   (read-knots (cffi:foreign-slot-value ptr '(:struct nl-sfun) 'knu))
                   (read-knots (cffi:foreign-slot-value ptr '(:struct nl-sfun) 'knv))
                   (from-foreign (cffi:foreign-slot-value sv '(:struct nl-svalue) 'fuv)
                                 :real-matrix (1+ n) (1+ m)))))
    (:epolygon
     (let ((n (cffi:foreign-slot-value ptr '(:struct nl-epolygon) 'n)))
       (read-doubles (cffi:foreign-slot-value ptr '(:struct nl-epolygon) 'pts) (list (1+ n) 3))))
    (:enet
     (let ((n (cffi:foreign-slot-value ptr '(:struct nl-enet) 'n))
           (m (cffi:foreign-slot-value ptr '(:struct nl-enet) 'm)))
       (read-net (cffi:foreign-slot-value ptr '(:struct nl-enet) 'pts) (1+ n) (1+ m) 3)))
    (:cpolygon
     (let ((n (cffi:foreign-slot-value ptr '(:struct nl-cpolygon) 'n)))
       (read-doubles (cffi:foreign-slot-value ptr '(:struct nl-cpolygon) 'pw) (list (1+ n) 4))))
    (:lineseg (batten::%make-lineseg (read-doubles ptr '(3))
                                     (read-doubles (cffi:inc-pointer ptr 24) '(3))
                                     (cffi:foreign-slot-value ptr '(:struct nl-lineseg) 'bounded)))
    (:plane (batten::%make-plane (read-doubles ptr '(3)) (read-doubles (cffi:inc-pointer ptr 24) '(3))))
    (:basis-derivatives
     (destructuring-bind (a b c d) dims
       (let ((out (make-array dims :element-type 'double-float)))
         (dotimes (k a out)
           (dotimes (l b)
             (dotimes (i c)
               (dotimes (j d)
                 (setf (aref out k l i j)
                       (cffi:mem-aref ptr :double (+ (* (+ (* (+ (* k 5) l) 18) i) 18) j))))))))))))

(defmacro with-foreign (bindings &body body)
  "BINDINGS: ((var object type) ...).  Bind each VAR to a foreign copy of OBJECT
as C type TYPE; free every allocation when BODY exits."
  `(let ((*allocations* nil))
     (unwind-protect
          (let* ,(mapcar (lambda (b) `(,(first b) (to-foreign ,(second b) ,(third b)))) bindings)
            ,@body)
       (mapc #'cffi:foreign-free *allocations*))))

(defmacro with-foreign-out (bindings &body body)
  "BINDINGS: ((var type &rest dims) ...).  Bind each VAR to an empty, zeroed
foreign object of TYPE for C to fill; read it back with FROM-FOREIGN inside
BODY.  Freed on exit (memory C allocated from NL_STACKS is freed by
WITH-STACKS)."
  `(let ((*allocations* nil))
     (unwind-protect
          (let* ,(mapcar (lambda (b) `(,(first b) (empty-foreign ,@(rest b)))) bindings)
            ,@body)
       (mapc #'cffi:foreign-free *allocations*))))
