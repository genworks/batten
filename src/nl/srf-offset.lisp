;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfOffset.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------- L9-g0

(defconstant +srf-offset--slm+ 1000.0d0
  "slm in NL_SrfOffset.cpp: the cap on the sampling count of one Bezier patch.")

(defconstant +srf-offset--totalmax+ 10000
  "totalmax in NL_SrfOffset.cpp: the cap on the number of sampled points.")

(defun %srf-offset-c-nan ()
  "The NaN x86-64 produces for an invalid operation (the sign bit set)."
  (sb-kernel:make-double-float -524288 0))

(defun %srf-offset-c-sqrt (x)
  "C's sqrt of the double X: X itself for a NaN, a NaN for a negative X, where
CL's SQRT would return a complex."
  (declare (type double-float x))
  (cond ((>= x 0d0) (sqrt (the (double-float 0d0) x)))
        ((/= x x) x)
        (t (%srf-offset-c-nan))))

(defun %srf-offset-c-pow (x y)
  "C's pow of the doubles X and Y (Y not a whole number): a NaN for a negative
X, where CL's EXPT would return a complex."
  (declare (type double-float x y))
  (if (< x 0d0)
      (%srf-offset-c-nan)
      (expt x y)))

(defun %srf-offset-c-ceil (x)
  "C's ceil of the double X, passing a NaN, an infinity or a double too large
to have a fraction through unchanged (where FCEILING would not)."
  (declare (type double-float x))
  (if (and (= x x) (< (abs x) 4.503599627370496d15))
      (values (fceiling x))
      x))

(defun %srf-offset-c-index (x)
  "(NL_INDEX)X: C's truncating conversion of the double X to a 64-bit integer,
x86-64's integer indefinite (-2^63) for a NaN or a value out of range."
  (declare (type double-float x))
  (if (and (>= x -9.223372036854775808d18) (< x 9.223372036854775808d18))
      (values (truncate x))
      (- (expt 2 63))))

(defun %srf-offset-orient (sur-q np sd)
  "The tail of N_SrfOffset's revolved and extruded cases: SUR-Q with u and v
swapped when its normal at its start corner points against NP.  The C ignores
an evaluation failure there when it swaps (the swap's success is returned) and
reports it otherwise; the normal is used in the state the failed evaluation
left it in."
  (declare (type surface sur-q) (type point np) (type point-net sd))
  (let ((pnt (make-point)) (b (make-point)) (c (make-point)) (nq (make-point))
        (failure nil))
    (multiple-value-bind (u0 u1 v0 v1) (srf-get-parameter-bounds sur-q)
      (declare (type double-float u0 v0) (ignore u1 v1))
      (handler-case (srf-eval-pt-pt-deriv-normal-fast sur-q u0 v0 +right+ +right+ pnt b c nq sd)
        (nurbs-error (e) (setf failure e))))
    (cond ((< (vector-dot np nq) 0d0) (swap-uv sur-q))
          (failure (error failure))
          (t sur-q))))

(defun srf-offset (sur-p d pp qq par eps tol)
  "An approximation of the offset at distance D of the G1 surface SUR-P, as a
fresh surface.  The surface type is found to within EPS (srf-type).  A plane
is offset exactly: translated (PAR +inherited+) or rebuilt bilinear on the
unit square from its moved corners.  A sphere, torus, cone, cylinder or
surface of revolution with a planar profile revolves the profile's offset
(crv-offset-pt-sampling, of degree 1 for a cone or cylinder, else QQ or PP
across the revolution); an extrusion of a planar profile perpendicular to it
extrudes the profile's offset; either is turned (u and v swapped) to face
SUR-P's normal at the start corner.  Anything else is sampled: per Bezier
patch a count from the offset's second-derivative bound and TOL (at most
1000, then the area estimate of the patch's offset control net), at most
10000 points in all; the offset points are interpolated with degrees PP, QQ
at SUR-P's parameters (PAR +inherited+) or approximated by least squares,
bicubic, at parameters by PAR (+uniform+, +chordlength+, +centripetal+), and
knots are removed within TOL.  Signals :num-err when 1/TOL is bad, :quit
for a negative sampling count, and what the classification, normals,
evaluations, curve offset, revolution, extrusion, fitting and knot removal
signal (a failed normal at the turned surface's start corner only when it is
not turned).
NLib: N_SrfOffset."
  (declare (type surface sur-p) (type double-float d eps tol) (type fixnum pp qq par))
  (let ((aa (make-point)) (bb (make-point)) (cc (make-point)) (dd (make-point))
        (v1 (make-point)) (v2 (make-point)) (v3 (make-point))
        (nn (make-point)) (cross (make-point)) (ta (make-point))
        (np (make-point)) (pnt (make-point)) (b (make-point)) (c (make-point))
        (zero (make-point 0d0 0d0 0d0))
        (sd (make-point-net 2 2)))
    (declare (type point aa bb cc dd v1 v2 v3 nn cross ta np pnt b c zero) (type point-net sd))
    ;; the surface type
    (multiple-value-bind (cur-a cur-b r1 r2 ac al h stp dir) (srf-type sur-p eps aa bb cc dd v1 v2)
      (declare (ignore cur-b r2 ac h) (type double-float r1 al) (type fixnum stp dir))
      (let ((flt (if (pts-are-equal cc zero) +no+ +yes+)))
        (declare (type fixnum flt))
        ;; offset a planar surface
        (when (= stp +nplane+)
          (flet ((normalize-ok (vec)
                   (handler-case (progn (vector-normalize-ref vec) t)
                     (nurbs-error () nil))))
            (vector-diff bb aa v1)
            (vector-diff dd aa v2)
            (vector-cross v1 v2 nn)
            (unless (normalize-ok nn)
              ;; assumed due to coincident BB, AA or DD, AA
              (vector-diff cc aa v3)
              (vector-cross v3 v2 nn)
              (unless (normalize-ok nn)
                (vector-cross v1 v3 nn)
                (vector-normalize-ref nn))))
          (return-from srf-offset
            (if (= par +inherited+)
                (progn
                  (vector-scale-ref nn d nn)
                  (srf-translate (srf-copy sur-p) nn))
                (progn
                  (vector-pt-along-vector aa d nn aa)
                  (vector-pt-along-vector bb d nn bb)
                  (vector-pt-along-vector cc d nn cc)
                  (vector-pt-along-vector dd d nn dd)
                  (create-srf-corner-pts aa bb dd cc)))))
        ;; offset a surface of revolution
        (let ((rev (or (= stp +nsphere+) (= stp +ntorus+) (= stp +ncone+)
                       (= stp +ncylinder+) (= stp +nrevolution+))))
          (when (and rev (= flt +yes+))
            (let ((pq (if (or (= stp +ncone+) (= stp +ncylinder+))
                          1
                          (if (= dir +udir+) qq pp))))
              (declare (type fixnum pq))
              (multiple-value-bind (u0 u1 v0 v1*) (srf-get-parameter-bounds sur-p)
                (declare (type double-float u0 v0) (ignore u1 v1*))
                (srf-eval-pt-pt-deriv-normal-fast sur-p u0 v0 +right+ +right+ pnt b c np sd))
              (let ((u0 (crv-get-param-bounds cur-a)))
                (declare (type double-float u0))
                (handler-case (crv-eval-tangent cur-a u0 +right+ pnt ta)
                  (nurbs-error () nil)))
              (vector-cross cc ta cross)
              (when (< (vector-dot np cross) 0d0)
                (vector-scale cc -1d0 cc))
              (let* ((cur-p (crv-offset-pt-sampling cur-a cc d pq +no+ par tol))
                     (sur-q (create-revolved-srf cur-p aa bb al +quadratic+)))
                (return-from srf-offset (%srf-offset-orient sur-q np sd))))))
        ;; offset an extruded surface
        (when (= stp +nextrusion+)
          (let* ((pflt (crv-is-planar cur-a eps dd nn))
                 (f1 (progn (vector-cross aa nn v1) (vector-magnitude v1))))
            (declare (type fixnum pflt) (type double-float f1))
            (when (and (= pflt +true+) (< f1 *mtol*))
              (let ((pq (if (= dir +udir+) qq pp)))
                (declare (type fixnum pq))
                ;; the C ignores this evaluation's failure
                (multiple-value-bind (u0 u1 v0 v1*) (srf-get-parameter-bounds sur-p)
                  (declare (type double-float u0 v0) (ignore u1 v1*))
                  (handler-case (srf-eval-pt-pt-deriv-normal-fast sur-p u0 v0 +right+ +right+
                                                                  pnt b c np sd)
                    (nurbs-error () nil)))
                (let ((u0 (crv-get-param-bounds cur-a)))
                  (declare (type double-float u0))
                  (handler-case (crv-eval-tangent cur-a u0 +right+ pnt ta)
                    (nurbs-error () nil)))
                (vector-cross nn ta cross)
                (when (< (vector-dot np cross) 0d0)
                  (vector-scale nn -1d0 nn))
                (let* ((cur-p (crv-offset-pt-sampling cur-a nn d pq +no+ par tol))
                       (sur-q (create-srf-extrude-crv cur-p aa r1 dir)))
                  (return-from srf-offset (%srf-offset-orient sur-q np sd)))))))
        ;; offset a general surface
        (%srf-offset-general sur-p d pp qq par tol)))))

(defun %srf-offset-general (sur-p d pp qq par tol)
  "N_SrfOffset's free-form case: sampling, interpolation or least-squares
approximation, knot removal (see srf-offset)."
  (declare (type surface sur-p) (type double-float d tol) (type fixnum pp qq par))
  (let* ((up (surface-knots-u sur-p)) (vp (surface-knots-v sur-p))
         (rp (1- (length up))) (sp (1- (length vp)))
         (p (surface-degree-u sur-p)) (q (surface-degree-v sur-p))
         ;; output control point counts, for the least-squares fit only
         (nnn (1+ rp)) (mmm (1+ sp))
         (ex (if (or (= pp 1) (= qq 1)) 0.5d0 0.34d0))
         (f1 0d0) (f2 0d0) (du 0d0) (dv 0d0)
         (pq 0)
         (mu 0) (mv 0))
    (declare (type real-vector up vp) (type fixnum rp sp p q nnn mmm pq)
             (type double-float ex f1 f2 du dv) (type integer mu mv))
    (when (float-op-is-bad 1.0d0 tol +division+)
      (nurbs-error :num-err 'srf-offset))
    (setf f1 (/ 1.0d0 tol))
    (setf f1 (%srf-offset-c-pow f1 ex))
    (setf pq (nl-max pp qq))
    (setf f2 (* 2.0d0 (float pq 1d0)))
    (setf du (* 0.00005d0 (- (aref up rp) (aref up 0))))
    (setf dv (* 0.00005d0 (- (aref vp sp) (aref vp 0))))
    ;; the Bezier patches
    (multiple-value-bind (bez k l) (srf-decompose-to-bez sur-p)
      (declare (type (simple-array t (* *)) bez) (type fixnum k l))
      (let ((nps (make-array (list (1+ k) (1+ l)) :initial-element 0))
            (nu (make-array (1+ k) :initial-element 0))
            (nv (make-array (1+ l) :initial-element 0))
            (nnet (make-point-net (1+ p) (1+ q)))
            (sd (make-point-net 2 2))
            (unit (make-real-vector 4))
            (aa (make-point)) (bb (make-point)) (cc (make-point)) (nij (make-point))
            (cp (make-cpoint)))
        (declare (type (simple-array t (* *)) nps) (type simple-vector nu nv)
                 (type point-net nnet sd) (type real-vector unit) (type point aa bb cc nij)
                 (type cpoint cp))
        (create-rectangle unit 0d0 1d0 0d0 1d0)
        ;; the number of sampling points of each patch
        (dotimes (i (1+ k))
          (dotimes (j (1+ l))
            (let ((patch (aref bez i j)))
              (declare (type surface patch))
              (multiple-value-bind (ul ur vl vr) (srf-get-parameter-bounds patch)
                (declare (type double-float ul ur vl vr))
                (if (or (<= (- ur ul) du) (<= (- vr vl) dv))
                    (setf (aref nps i j) 1)
                    (progn
                      (srf-reparam-to-interval patch unit +uvdir+)
                      (multiple-value-bind (muu muv mvv) (srf-offset-get-max-deriv patch d)
                        (declare (type double-float muu muv mvv))
                        (let* ((num (+ (* 0.125d0 (+ muu mvv)) (* 0.25d0 muv)))
                               (gro (nl-max f2 (%srf-offset-c-sqrt num)))
                               (del (%srf-offset-c-ceil (* f1 gro))))
                          (declare (type double-float num gro del))
                          (if (>= del +srf-offset--slm+)
                              (let ((pw (surface-cpoints patch))
                                    (uinc (/ 1.0d0 (float p 1d0)))
                                    (vinc (/ 1.0d0 (float q 1d0))))
                                (declare (type cpoint-net pw) (type double-float uinc vinc))
                                (setf (aref nps i j) (%srf-offset-c-index +srf-offset--slm+))
                                ;; normals at the patch's nodes
                                (dotimes (ii (1+ p))
                                  (let ((u (* (float ii 1d0) uinc)))
                                    (declare (type double-float u))
                                    (when (> u 1.0d0) (setf u 1.0d0))
                                    (dotimes (jj (1+ q))
                                      (let ((v (* (float jj 1d0) vinc)))
                                        (declare (type double-float v))
                                        (when (> v 1.0d0) (setf v 1.0d0))
                                        (handler-case
                                            (srf-eval-pt-pt-deriv-normal-fast
                                             patch u v +left+ +left+ aa bb cc nij sd)
                                          (nurbs-error ()
                                            (srf-eval-pt-normal-at-pole
                                             patch u v +left+ +left+ aa bb cc nij sd)))
                                        (%pt-net nij nnet ii jj)))))
                                ;; the patch's offset control net
                                (dotimes (ii (1+ p))
                                  (dotimes (jj (1+ q))
                                    (%net-c-pt pw ii jj cp)
                                    (let ((w (c-pt-to-pt-and-w cp aa)))
                                      (declare (type double-float w))
                                      (%net-pt nnet ii jj nij)
                                      (vector-pt-along-vector aa d nij aa)
                                      (weight aa w cp)
                                      (%c-pt-net cp pw ii jj))))
                                (let* ((a (values (srf-get-average-len patch)))
                                       (del2 (%srf-offset-c-ceil (%srf-offset-c-sqrt (/ a tol))))
                                       (cnt (%srf-offset-c-index del2)))
                                  (declare (type double-float a del2))
                                  (when (< cnt (aref nps i j))
                                    (setf (aref nps i j) cnt))))
                              (setf (aref nps i j) (%srf-offset-c-index del)))))))))))
        ;; the most points to be sampled
        (setf mu 0)
        (dotimes (i (1+ k))
          (setf (svref nu i) (aref nps i 0))
          (do ((j 1 (1+ j))) ((> j l))
            (declare (type fixnum j))
            (when (> (aref nps i j) (svref nu i))
              (setf (svref nu i) (aref nps i j))))
          (setf mu (+ mu (svref nu i))))
        (setf mv 0)
        (dotimes (j (1+ l))
          (setf (svref nv j) (aref nps 0 j))
          (do ((i 1 (1+ i))) ((> i k))
            (declare (type fixnum i))
            (when (> (aref nps i j) (svref nv j))
              (setf (svref nv j) (aref nps i j))))
          (setf mv (+ mv (svref nv j))))
        (when (> (* (1+ mu) (1+ mv)) +srf-offset--totalmax+)
          (let* ((ddd (* (+ (float mu 1d0) 1.0d0) (+ (float mv 1d0) 1.0d0)))
                 (per (%srf-offset-c-sqrt (/ (float +srf-offset--totalmax+ 1d0) ddd))))
            (declare (type double-float ddd per))
            (setf mu 0)
            (dotimes (i (1+ k))
              (setf ddd (float (svref nu i) 1d0))
              (setf (svref nu i) (%srf-offset-c-index (* ddd per)))
              (when (< (svref nu i) p)
                (setf (svref nu i) p))
              (setf mu (+ mu (svref nu i))))
            (setf mv 0)
            (dotimes (j (1+ l))
              (setf ddd (float (svref nv j) 1d0))
              (setf (svref nv j) (%srf-offset-c-index (* ddd per)))
              (when (< (svref nv j) q)
                (setf (svref nv j) q))
              (setf mv (+ mv (svref nv j))))))
        (when (or (< mu 0) (< mv 0))
          (nurbs-error :quit 'srf-offset))
        ;; the output control point counts
        (when (> nnn mu) (setf nnn mu))
        (when (> mmm mv) (setf mmm mv))
        ;; the parameters
        (let ((u (make-real-vector (1+ mu)))
              (v (make-real-vector (1+ mv))))
          (declare (type real-vector u v))
          (let ((ul (aref up p)) (ur 0d0) (i (1+ p)) (j 0) (ll 0))
            (declare (type double-float ul ur) (type fixnum i j ll))
            (setf (aref u 0) ul)
            (do () ((>= i rp))
              (do () ((not (and (< i rp) (= (aref up i) (aref up (1+ i))))))
                (incf i))
              (setf ur (aref up i))
              (let ((uinc (/ (- ur ul) (float (svref nu j) 1d0))))
                (declare (type double-float uinc))
                (do ((kk 1 (1+ kk))) ((>= kk (svref nu j)))
                  (declare (type fixnum kk))
                  (setf (aref u (incf ll)) (+ ul (* (float kk 1d0) uinc)))))
              (setf (aref u (incf ll)) ur)
              (setf ul ur)
              (incf i)
              (incf j)))
          (let ((vl (aref vp q)) (vr 0d0) (i (1+ q)) (j 0) (ll 0))
            (declare (type double-float vl vr) (type fixnum i j ll))
            (setf (aref v 0) vl)
            (do () ((>= i sp))
              (do () ((not (and (< i sp) (= (aref vp i) (aref vp (1+ i))))))
                (incf i))
              (setf vr (aref vp i))
              (let ((vinc (/ (- vr vl) (float (svref nv j) 1d0))))
                (declare (type double-float vinc))
                (do ((kk 1 (1+ kk))) ((>= kk (svref nv j)))
                  (declare (type fixnum kk))
                  (setf (aref v (incf ll)) (+ vl (* (float kk 1d0) vinc)))))
              (setf (aref v (incf ll)) vr)
              (setf vl vr)
              (incf i)
              (incf j)))
          ;; sample the surface
          (let ((pts (make-point-net (1+ mu) (1+ mv)))
                (nrm (make-point-net (1+ mu) (1+ mv)))
                (pt (make-point)) (nv* (make-point)))
            (declare (type point-net pts nrm) (type point pt nv*))
            (srf-get-grid-pts-normals sur-p u v mu mv +left+ +left+ pts nrm)
            (dotimes (i (1+ mu))
              (dotimes (j (1+ mv))
                (%net-pt pts i j pt)
                (%net-pt nrm i j nv*)
                (vector-pt-along-vector pt d nv* pt)
                (%pt-net pt pts i j)))
            ;; the knot vectors
            (let ((knu (make-real-vector (+ mu pp 2)))
                  (knv (make-real-vector (+ mv qq 2)))
                  (sur-q nil))
              (declare (type real-vector knu knv))
              (when (/= par +inherited+)
                (fit-calc-srf-param-values pts mu mv +epoint+ par u v))
              (fit-crv-calc-knot-vector u mu pp knu)
              (fit-crv-calc-knot-vector v mv qq knv)
              ;; interpolate or approximate the points
              (setf sur-q (if (= par +inherited+)
                              (fit-srf-to-pts-knots pts mu mv +epoint+ u v knu knv pp qq)
                              (fit-srf-lst-sq-approx pts mu mv nnn mmm 3 3 par)))
              ;; remove knots
              (setf sur-q (srf-remove-all-knots sur-q tol +uvdir+))
              ;; compact the output surface
              (when (or (< (surface-n sur-q) mu) (< (surface-m sur-q) mv))
                (srf-compress sur-q))
              ;; rescale the knot vectors if necessary
              (when (/= par +inherited+)
                (let* ((uq (surface-knots-u sur-q)) (vq (surface-knots-v sur-q))
                       (rq (1- (length uq))) (sq (1- (length vq))))
                  (declare (type real-vector uq vq) (type fixnum rq sq))
                  (when (or (/= (aref up 0) (aref uq 0)) (/= (aref up rp) (aref uq rq))
                            (/= (aref vp 0) (aref vq 0)) (/= (aref vp sp) (aref vq sq)))
                    (let ((rect (make-real-vector 4)))
                      (declare (type real-vector rect))
                      (create-rectangle rect (aref up 0) (aref up rp) (aref vp 0) (aref vp sp))
                      (cond ((and (= (aref up 0) (aref uq 0)) (= (aref up rp) (aref uq rq)))
                             (srf-reparam-to-interval sur-q rect +vdir+))
                            ((and (= (aref vp 0) (aref vq 0)) (= (aref vp sp) (aref vq sq)))
                             (srf-reparam-to-interval sur-q rect +udir+))
                            (t (srf-reparam-to-interval sur-q rect +uvdir+)))))))
              sur-q)))))))
