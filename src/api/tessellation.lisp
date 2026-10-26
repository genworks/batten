;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The readable API: surfaces as triangles, by NLib's tessellator
;;;; (tess-trimmed-srf), trimmed or not, holes included.

(in-package :batten)

(defun %parameter-line (u0 v0 u1 v1)
  "The straight 2D curve from (U0, V0) to (U1, V1) in a surface's parameter
space."
  (make-curve 1 '(0 0 1 1) (list (list u0 v0 0) (list u1 v1 0)) :two-d t))

(defun %domain-loop (surface)
  "SURFACE's parameter domain as one counterclockwise loop of four lines."
  (multiple-value-bind (u0 u1 v0 v1) (srf-get-parameter-bounds surface)
    (vector (%parameter-line u0 v0 u1 v0) (%parameter-line u1 v0 u1 v1)
            (%parameter-line u1 v1 u0 v1) (%parameter-line u0 v1 u0 v0))))

(defun %index-vector (list)
  (make-array (length list) :element-type 'fixnum :initial-contents list))

(defun tessellate-surface (surface &key tolerance trim-tolerance (parameter-tolerance 1d-7)
                                        outer-loops holes normals)
  "SURFACE as triangles, within TOLERANCE of it.  Values: POINTS, a point-array
of the vertices; TRIANGLES, an array of N x 3 vertex indexes, each triangle
counterclockwise in the parameter space (so facing the way the surface's
u-cross-v normal faces); PARAMETERS, an array of the vertices' (u v); and,
with NORMALS true, a point-array of the surface's unit normals at the
vertices (a zero vector where the surface has none, at a pole).

Untrimmed, the whole parameter domain is meshed.  To trim, give OUTER-LOOPS,
a list of loops, each a list of 2D curves in SURFACE's parameter space joined
end to start, counterclockwise; and HOLES, inside them, clockwise: for one
outer loop a list of loops, for several a list with one such list per outer
loop.  The loops must be closed, simple and disjoint.

TOLERANCE (default a thousandth of the surface's size) bounds how far a
triangle strays from the surface; TRIM-TOLERANCE (default TOLERANCE) how far
the polygons standing for the trimming curves stray from them;
PARAMETER-TOLERANCE is the distance at which two parameters are one.
SURFACE is not changed.
Port: tess-trimmed-srf, tess-get-triangles, srf-eval-pt, srf-eval-pt-pt-deriv-normal."
  (with-c-float-environment
    (let* ((s (srf-copy surface))   ; tess-trimmed-srf rescales the knots
           (eps (if tolerance
                    (%real tolerance)
                    (* 1d-3 (sqrt (max 1d-12 (values (tess-srf-area s)))))))
           (epc (if trim-tolerance (%real trim-tolerance) eps))
           (outer (if outer-loops
                      (mapcar #'(lambda (loop) (coerce loop 'simple-vector)) outer-loops)
                      (list (%domain-loop s))))
           (inner (mapcar #'(lambda (loops)
                              (coerce (mapcar #'(lambda (loop) (coerce loop 'simple-vector)) loops)
                                      'simple-vector))
                          (if (= (length outer) 1)
                              (list holes)
                              (or holes (make-list (length outer))))))
           (m (1- (length outer)))
           (cuo (coerce outer 'simple-vector))
           (cui (coerce inner 'simple-vector))
           (ho (%index-vector (mapcar #'(lambda (loop) (1- (length loop))) outer)))
           (hi (%index-vector (mapcar #'(lambda (loops) (1- (length loops))) inner)))
           (hs (coerce (mapcar #'(lambda (loops) (%index-vector (map 'list #'(lambda (loop) (1- (length loop))) loops)))
                               inner)
                       'simple-vector)))
      (multiple-value-bind (u v n dt hd)
          (tess-trimmed-srf s cuo cui m ho hi hs epc eps (%real parameter-tolerance))
        (declare (type real-vector u v) (type fixnum n))
        (multiple-value-bind (alf bet gam k) (tess-get-triangles dt hd n)
          (declare (type index-vector alf bet gam) (type fixnum k))
          (let ((points (make-point-array (1+ n)))
                (parameters (make-array (list (1+ n) 2) :element-type 'double-float))
                (triangles (make-array (list (1+ k) 3) :element-type 'fixnum))
                (point (make-point)))
            (dotimes (i (1+ n))
              (srf-eval-pt s (aref u i) (aref v i) +left+ +left+ point)
              (setf (aref points i 0) (aref point 0) (aref points i 1) (aref point 1)
                    (aref points i 2) (aref point 2)
                    (aref parameters i 0) (aref u i) (aref parameters i 1) (aref v i)))
            (dotimes (j (1+ k))
              (setf (aref triangles j 0) (aref alf j) (aref triangles j 1) (aref bet j)
                    (aref triangles j 2) (aref gam j)))
            (if normals
                (let ((normal-array (make-point-array (1+ n)))
                      (p (make-point)) (su (make-point)) (sv (make-point)) (nv (make-point)))
                  (dotimes (i (1+ n))
                    (when (handler-case
                              (progn (srf-eval-pt-pt-deriv-normal s (aref u i) (aref v i) +left+ +left+ p su sv nv) t)
                            (nurbs-error () nil))
                      (setf (aref normal-array i 0) (aref nv 0) (aref normal-array i 1) (aref nv 1)
                            (aref normal-array i 2) (aref nv 2))))
                  (values points triangles parameters normal-array))
                (values points triangles parameters))))))))
