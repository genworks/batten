;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_FrameBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.
;;;;
;;;; NLib's polygon and net structs (NL_CPOLYGON, NL_EPOLYGON, NL_CNET,
;;;; NL_ENET) are, in Batten, the point arrays themselves: the highest indexes
;;;; n and m are the array dimensions less one.  The C routines that only
;;;; assign a struct's fields therefore take no struct here: the ...-FROM-...
;;;; routines return the polygon or net, and the ...-GET-PTS routines return
;;;; the highest indexes and the array.

(in-package :batten)

(defun %frame-basic-rows (array rows)
  "ARRAY (a point-array or cpoint-array) as a polygon of ROWS points: ARRAY
itself when it has exactly ROWS rows (shared, as the C shares the pointer),
else a fresh copy of its first ROWS rows."
  (declare (type (simple-array double-float (* *)) array) (type fixnum rows))
  (if (= rows (array-dimension array 0))
      array
      (let* ((k (array-dimension array 1))
             (out (make-array (list rows k) :element-type 'double-float)))
        (declare (type fixnum k))
        (fast (dotimes (i rows out)
                (dotimes (c k)
                  (setf (aref out i c) (aref array i c))))))))

(defun %frame-basic-net (array rows cols)
  "ARRAY (a point-net or cpoint-net) as a net of ROWS by COLS points: ARRAY
itself when it has exactly those dimensions (shared, as the C shares the
pointer), else a fresh copy of that corner of it."
  (declare (type (simple-array double-float (* * *)) array) (type fixnum rows cols))
  (if (and (= rows (array-dimension array 0)) (= cols (array-dimension array 1)))
      array
      (let* ((k (array-dimension array 2))
             (out (make-array (list rows cols k) :element-type 'double-float)))
        (declare (type fixnum k))
        (fast (dotimes (i rows out)
                (dotimes (j cols)
                  (dotimes (c k)
                    (setf (aref out i j c) (aref array i j c)))))))))

(defun c-polygon-from-c-pts (pw n)
  "The control polygon of control points PW[0..N]: PW itself when N is its
highest index, else its first N+1 points.  (The C fills the struct POL the
caller passes; here the polygon is the array.)
NLib: N_CPolygonFromCPts."
  (declare (type cpoint-array pw) (type fixnum n))
  (the cpoint-array (%frame-basic-rows pw (1+ n))))

(defun e-polygon-from-pts (n p)
  "The polygon with vertices P[0..N]: P itself when N is its highest index,
else its first N+1 points.  (The C fills the struct PPL the caller passes;
here the polygon is the array.)
NLib: N_EPolygonFromPts."
  (declare (type fixnum n) (type point-array p))
  (the point-array (%frame-basic-rows p (1+ n))))

(defun e-polygon-get-pts (ppl)
  "The highest vertex index of polygon PPL and its vertex array (PPL itself),
as two values.
NLib: N_EPolygonGetPts."
  (declare (type point-array ppl))
  (values (1- (array-dimension ppl 0)) ppl))

(defun c-net-from-c-pts (pw n m)
  "The control net of control points PW[0..N][0..M]: PW itself when N and M
are its highest indexes, else that corner of it.  (The C fills the struct NET
the caller passes; here the net is the array.)
NLib: N_CNetFromCPts."
  (declare (type cpoint-net pw) (type fixnum n m))
  (the cpoint-net (%frame-basic-net pw (1+ n) (1+ m))))

(defun e-net-from-pts (n m p)
  "The point net with vertices P[0..N][0..M]: P itself when N and M are its
highest indexes, else that corner of it.  (The C fills the struct NTL the
caller passes; here the net is the array.)
NLib: N_ENetFromPts."
  (declare (type fixnum n m) (type point-net p))
  (the point-net (%frame-basic-net p (1+ n) (1+ m))))

(defun e-net-get-pts (ntl)
  "The highest indexes of point net NTL and its vertex array (NTL itself), as
three values.
NLib: N_ENetGetPts."
  (declare (type point-net ntl))
  (values (1- (array-dimension ntl 0)) (1- (array-dimension ntl 1)) ntl))

;;; ------------------------------------------------------------ level 1, group 1

(defun polygon-get-b-box (ppl box)
  "Write the min-max box of polygon PPL (a point-array) into BOX (a
real-vector of 6: xl xr yb yt zn zf) and return BOX.
NLib: N_PolygonGetBBox."
  (declare (type point-array ppl) (type real-vector box))
  (multiple-value-bind (n p) (e-polygon-get-pts ppl)
    (declare (type fixnum n) (type point-array p))
    (let* ((x (aref p 0 0)) (y (aref p 0 1)) (z (aref p 0 2))
           (xl x) (xr x) (yb y) (yt y) (zn z) (zf z))
      (declare (type double-float x y z xl xr yb yt zn zf))
      (fast
        (do ((i 1 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (setf x (aref p i 0) y (aref p i 1) z (aref p i 2))
          (when (< x xl) (setf xl x))
          (when (> x xr) (setf xr x))
          (when (< y yb) (setf yb y))
          (when (> y yt) (setf yt y))
          (when (< z zn) (setf zn z))
          (when (> z zf) (setf zf z))))
      (b-box-define box xl xr yb yt zn zf))))

(defun e-net-get-b-box (ntl box)
  "Write the min-max box of point net NTL (a point-net) into BOX (a
real-vector of 6: xl xr yb yt zn zf) and return BOX.
NLib: N_ENetGetBBox."
  (declare (type point-net ntl) (type real-vector box))
  (multiple-value-bind (n m p) (e-net-get-pts ntl)
    (declare (type fixnum n m) (type point-net p))
    (let* ((x (aref p 0 0 0)) (y (aref p 0 0 1)) (z (aref p 0 0 2))
           (xl x) (xr x) (yb y) (yt y) (zn z) (zf z))
      (declare (type double-float x y z xl xr yb yt zn zf))
      (fast
        (dotimes (i (1+ n))
          (dotimes (j (1+ m))
            (setf x (aref p i j 0) y (aref p i j 1) z (aref p i j 2))
            (when (< x xl) (setf xl x))
            (when (> x xr) (setf xr x))
            (when (< y yb) (setf yb y))
            (when (> y yt) (setf yt y))
            (when (< z zn) (setf zn z))
            (when (> z zf) (setf zf z)))))
      (b-box-define box xl xr yb yt zn zf))))

(defun e-net-is-closed (ntl dir)
  "True when point net NTL is closed in direction DIR: for +udir+ every
column's first and last points, for +vdir+ every row's, lie within *mtol*
of each other.  Any other DIR answers false.
NLib: N_ENetIsClosed."
  (declare (type point-net ntl) (type fixnum dir))
  (multiple-value-bind (n m p) (e-net-get-pts ntl)
    (declare (type fixnum n m) (type point-net p))
    (let ((a (make-point)) (b (make-point)) (mtol *mtol*))
      (declare (dynamic-extent a b) (type double-float mtol))
      (cond ((= dir +udir+)
             (dotimes (j (1+ m) t)
               (%net-pt p 0 j a)
               (%net-pt p n j b)
               (when (> (dist-pt-pt a b) mtol)
                 (return nil))))
            ((= dir +vdir+)
             (dotimes (i (1+ n) t)
               (%net-pt p i 0 a)
               (%net-pt p i m b)
               (when (> (dist-pt-pt a b) mtol)
                 (return nil))))
            (t nil)))))
