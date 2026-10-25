;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvGeom.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------ level 2, group 6

(defun crv-closest-pt (cur p u0 top toc q &optional ucell)
  "Find the point of curve CUR closest to point P (inversion when P is on the
curve, projection when it is off it) by Newton iteration from the guess U0,
wrapping the parameter on a closed curve and clamping it on an open one.  TOP
is the point coincidence tolerance, TOC the zero cosine tolerance.  Writes the
best point found into Q and returns its parameter u.

Signals :con-err when the iteration fails or does not converge within *itlim*
steps, and what crv-derivs signals.  The C returns the best u and Q found even
then: Q is written as it improves, and so is UCELL, an optional real-vector of
1 receiving u, for callers that carry on after the error.
NLib: N_CrvClosestPt."
  (declare (type curve cur) (type point p q) (type double-float u0 top toc))
  (let* ((knots (curve-knots cur))
         (m (1- (length knots)))
         (closed (crv-is-closed cur))
         (k 0) (ucount 0)
         (itlim *itlim*)
         (u 0d0)
         (unew u0) (uold 0d0)
         (best-dis +bigd+)
         (num 0d0) (den 0d0) (der 0d0) (dot 0d0) (dis 0d0) (zco 0d0)
         (pp (make-point))
         (d (make-point-array 3))
         (d0 (make-point)) (d1 (make-point)) (d2 (make-point))
         (v (make-point)))
    (declare (type real-vector knots) (type fixnum m k ucount itlim)
             (type double-float u unew uold best-dis num den der dot dis zco)
             (dynamic-extent pp d d0 d1 d2 v))
    ;; P is passed by value in C; Q may be the same array here
    (copy-pt p pp)
    (flet ((con-err () (nurbs-error :con-err 'crv-closest-pt)))
      (do () ((>= k itlim))
        ;; derivatives
        (crv-derivs cur unew +left+ 2 d)
        (%row-pt d 0 d0)
        (%row-pt d 1 d1)
        (%row-pt d 2 d2)
        ;; prepare for the convergence test
        (vector-diff d0 pp v)
        (setf num (vector-dot d1 v))
        (setf der (vector-magnitude d1))
        ;; check 1: point coincidence
        (setf dis (vector-magnitude v))
        (when (< dis best-dis)
          (setf u unew)
          (when ucell (setf (aref (the real-vector ucell) 0) u))
          (copy-pt d0 q)
          (setf best-dis dis))
        (when (<= dis top) (return))
        ;; check 2: zero cosine
        (when (float-op-is-bad num (* der dis) +division+) (con-err))
        (setf zco (/ (abs num) (* der dis)))
        (when (and (<= zco toc) (= dis best-dis)) (return))
        ;; no convergence: a new parameter
        (setf dot (vector-dot d2 v))
        (setf den (+ dot (* der der)))
        (when (float-op-is-bad num den +division+) (con-err))
        (setf uold unew
              unew (- uold (/ num den)))
        ;; check 3: parameter range
        (let ((ua (aref knots 0)) (ub (aref knots m)))
          (declare (type double-float ua ub))
          (if closed
              (cond ((< unew ua)
                     (cond ((< (- ua unew) (* 100d0 (- ub ua)))
                            (do () ((not (< unew ua)))
                              (setf unew (+ (- ub ua) unew)))
                            (setf ucount 0))
                           (t
                            (when (>= ucount 2) (con-err))
                            (setf unew uold)
                            (incf ucount))))
                    ((> unew ub)
                     (cond ((< (- unew ub) (* 100d0 (- ub ua)))
                            (do () ((not (> unew ub)))
                              (setf unew (+ (- ua ub) unew)))
                            (setf ucount 0))
                           (t
                            (when (>= ucount 2) (con-err))
                            (setf unew uold)
                            (incf ucount))))
                    (t (setf ucount 0)))
              (progn
                (when (< unew ua) (setf unew ua))
                (when (> unew ub) (setf unew ub)))))
        ;; check 4: parameter change
        (when (and (<= (abs (* (- unew uold) der)) top) (= dis best-dis))
          (return))
        (incf k))
      (when (>= k itlim) (con-err))
      u)))

;;; ------------------------------------------------------------- L5-g0

(defconstant +crv-geom--arc-length-ctl+ 1d-9
  "N_CrvArcLength's ctl: the line/circle tolerance, relative to the extent.")

(defun crv-arc-length (cur ul ur tol flg)
  "The arc length of curve CUR between parameters UL and UR (0 when they are
within *ptol*), to tolerance TOL, absolute (FLG +absolute+) or relative to the
segment's control-polygon extent (+relative+).  The segment is extracted; a
line or circular arc (to 1e-9 of the extent) is measured exactly, anything
else Bezier piece by Bezier piece, decomposed on the fly, each to TOL over the
number of spans.  Signals :cal-err for another FLG (only for a free-form
segment), :quit for an interior knot of multiplicity above the degree, and
what crv-extract-crv-seg, crv-get-type and bez-calc-length signal.
NLib: N_CrvArcLength."
  (declare (type curve cur) (type double-float ul ur tol) (type fixnum flg))
  ;; parameter degeneracy
  (when (< (abs (- ur ul)) (the double-float *ptol*))
    (return-from crv-arc-length 0d0))
  (let* ((cura (crv-extract-crv-seg cur ul ur))
         (d (crv-get-b-box-max-diag-dist cura))
         (lctl (* d +crv-geom--arc-length-ctl+))
         (aa (make-point)) (bb (make-point)) (xx (make-point)) (yy (make-point)))
    (declare (type double-float d lctl) (dynamic-extent aa bb xx yy))
    ;; line and circle
    (multiple-value-bind (rad al typ) (crv-get-type cura lctl aa bb xx yy)
      (declare (type double-float rad al) (type fixnum typ))
      (when (= typ +nline+)
        (vector-sum aa bb bb)
        (return-from crv-arc-length (dist-pt-pt aa bb)))
      (when (= typ +ncircle+)
        (return-from crv-arc-length (* (* (/ al 180d0) +pi+) rad))))
    (let* ((ltol (cond ((= flg +relative+) (* d tol))
                       ((= flg +absolute+) tol)
                       (t (nurbs-error :cal-err 'crv-arc-length))))
           (pw (curve-cpoints cura))
           (p (curve-degree cura))
           (m (curve-m cura))
           (u (curve-knots cura))
           (nsp (basis-get-span-count u p))
           (eps (/ ltol (float nsp 1d0)))
           (len 0d0)
           (bw (make-cpoint-array (1+ p)))
           (nw (make-cpoint-array (1+ p)))
           (alfs (make-real-vector (1+ p)))
           (omas (make-real-vector (1+ p)))
           (a p) (b (1+ p)) (r -1))
      (declare (type double-float ltol eps len) (type cpoint-array pw bw nw)
               (type real-vector u alfs omas) (type fixnum p m nsp a b r))
      (dotimes (i (1+ p))
        (%copy-c-pt-row pw i bw i))
      ;; decompose the curve and measure each piece
      (do () ((not (< b m)))
        (let ((i b))
          (declare (type fixnum i))
          (do () ((not (and (< b m) (= (aref u b) (aref u (1+ b))))))
            (incf b))
          (let ((mlt (+ (- b i) 1)))
            (declare (type fixnum mlt))
            (setf r (- p mlt))
            (when (> r 0)
              (let ((num (- (aref u b) (aref u a))))
                (declare (type double-float num))
                (do ((k p (1- k))) ((<= k mlt))
                  (declare (type fixnum k))
                  (setf (aref alfs (- k mlt 1)) (/ num (- (aref u (+ a k)) (aref u a))))
                  (setf (aref omas (- k mlt 1)) (- 1d0 (aref alfs (- k mlt 1))))))
              (do ((j 1 (1+ j))) ((> j r))
                (declare (type fixnum j))
                (let ((save (- r j)) (s (+ mlt j)))
                  (declare (type fixnum save s))
                  (do ((k p (1- k))) ((< k s))
                    (declare (type fixnum k))
                    (%combine-2-c-pt-rows (aref alfs (- k s)) bw k (aref omas (- k s)) bw (1- k) bw k))
                  (%copy-c-pt-row bw p nw save))))))
        ;; the length of the Bezier segment
        (incf len (the double-float (bez-calc-length bw p eps)))
        ;; initialize for the next pass
        (when (< b m)
          (when (< r 0)
            ;; an interior knot of full multiplicity plus one: the C returns
            ;; NL_YES with no error code
            (nurbs-error :quit 'crv-arc-length))
          (dotimes (i r)
            (%copy-c-pt-row nw i bw i))
          (do ((i r (1+ i))) ((> i p))
            (declare (type fixnum i))
            (%copy-c-pt-row pw (+ (- b p) i) bw i))
          (setf a b)
          (incf b)))
      len)))

(defun %crv-geom-grow-rows (arr rows)
  "A copy of 2D array ARR (double-float) with ROWS rows, the old rows copied."
  (declare (type (simple-array double-float (* *)) arr) (type fixnum rows))
  (let* ((cols (array-dimension arr 1))
         (new (make-array (list rows cols) :element-type 'double-float :initial-element 0d0)))
    (dotimes (i (min rows (array-dimension arr 0)) new)
      (dotimes (j cols)
        (setf (aref new i j) (aref arr i j))))))

(defun approx-crv-with-polyline (cur tol tflg rflg)
  "A polyline through points of curve CUR deviating from it by at most TOL
(TFLG +absolute+), or TOL times the curve's length (+relative+): the curve is
decomposed into Bezier pieces, each subdivided until its control polygon is
within TOL of its chord.  A degree-1 curve gives its control points and the
parameters U[i+1].  Returns (values pts pars n): PTS a point-array of the
N+1 points (nil when RFLG is +parameters+), PARS a real-vector of their
parameters (nil when RFLG is +points+).  Signals what crv-decompose-bez and
bez-calc-length signal.
NLib: N_ApproxCrvWithPolyline."
  (declare (type curve cur) (type double-float tol) (type fixnum tflg rflg))
  (let ((p (curve-degree cur))
        (want-pts (/= rflg +parameters+))
        (want-pars (/= rflg +points+)))
    (declare (type fixnum p))
    ;; p = 1 as a special case
    (when (= p 1)
      (let* ((pw (curve-cpoints cur))
             (u (curve-knots cur))
             (nn (curve-n cur))
             (pnts (when want-pts (make-point-array (1+ nn))))
             (parms (when want-pars (make-real-vector (1+ nn))))
             (cw (make-cpoint)) (pt (make-point)))
        (declare (type cpoint-array pw) (type real-vector u) (type fixnum nn)
                 (dynamic-extent cw pt))
        (dotimes (ii (1+ nn))
          (when want-pts
            (c-pt-to-pt-euclid (%row-c-pt pw ii cw) pt)
            (%pt-row pt pnts ii))
          (when want-pars
            (setf (aref (the real-vector parms) ii) (aref u (1+ ii)))))
        (return-from approx-crv-with-polyline (values pnts parms nn))))
    ;; decompose the curve into Bezier segments
    (multiple-value-bind (bezs nn) (crv-decompose-bez cur)
      (declare (type simple-vector bezs) (type fixnum nn))
      ;; adjust the tolerance
      (when (= tflg +relative+)
        (let ((length 0d0))
          (declare (type double-float length))
          (dotimes (ii (1+ nn))
            (let ((bz (svref bezs ii)))
              (setf p (curve-degree bz))
              (incf length (the double-float (bez-calc-length (curve-cpoints bz) p 0.01d0)))))
          (setf tol (* tol length))))
      ;; the stacks of Bezier segments
      (let* ((s-size (cond ((< nn 6) 100) ((< nn 15) (* 10 nn)) (t (* 5 nn))))
             (cp-stack (make-array (1+ s-size)))
             (pm-stack (when want-pars (make-real-matrix (1+ s-size) 2)))
             (top -1)
             (p-size 250)
             (pnts (when want-pts (make-point-array (1+ p-size))))
             (parms (when want-pars (make-real-vector (1+ p-size))))
             (p1 (make-point)) (p2 (make-point)) (pp (make-point)) (qq (make-point))
             (last-pt (make-point)) (v12 (make-point)) (v1p (make-point))
             (cw (make-cpoint))
             (tolsq 0d0) (dd 0d0) (dot 0d0)
             (ii 0))
        (declare (type fixnum s-size top p-size ii) (type simple-vector cp-stack)
                 (type double-float tolsq dd dot)
                 (dynamic-extent p1 p2 pp qq last-pt v12 v1p cw))
        (dotimes (k (1+ s-size))
          (setf (svref cp-stack k) (make-cpoint-array (1+ p))))
        (do ((i nn (1- i))) ((< i 0))
          (declare (type fixnum i))
          (let* ((bz (svref bezs i))
                 (pw (curve-cpoints bz))
                 (u (curve-knots bz))
                 (kk (curve-m bz)))
            (declare (type cpoint-array pw) (type real-vector u) (type fixnum kk))
            (setf p (curve-degree bz))
            (incf top)
            (when want-pars
              (setf (aref (the real-matrix pm-stack) top 0) (aref u 0)
                    (aref (the real-matrix pm-stack) top 1) (aref u kk)))
            (dotimes (jj (1+ p))
              (%copy-c-pt-row pw jj (svref cp-stack top) jj))))
        ;; the points and parameters
        (flet ((top-pt (k pt)
                 (c-pt-to-pt-euclid (%row-c-pt (svref cp-stack top) k cw) pt)))
          (when want-pts
            (top-pt 0 pp)
            (%pt-row pp pnts 0))
          (when want-pars
            (setf (aref (the real-vector parms) 0) (aref (the real-matrix pm-stack) top 0)))
          (top-pt 0 last-pt)
          (setf nn 0)
          ;; subdivide the segments until the stack is empty
          (setf tolsq (* tol tol))
          (do () ((not (>= top 0)))
            ;; has this segment converged?
            (top-pt 0 p1)
            (top-pt p p2)
            (vector-dir p1 p2 v12)
            (if (handler-case (progn (vector-normalize v12 v12) t)
                  (nurbs-error () nil))
                (progn
                  (setf ii 1)
                  (do () ((not (< ii p)))
                    (top-pt ii pp)
                    (vector-dir p1 pp v1p)
                    (setf dot (vector-dot v12 v1p))
                    (vector-pt-along-vector p1 dot v12 qq)
                    (setf dd (dist-sq-pt-pt pp qq))
                    (when (> dd tolsq)
                      (return))
                    (incf ii)))
                ;; the segment is closed or degenerate
                (progn
                  (top-pt (truncate (1+ p) 2) pp)
                  (setf dd (dist-sq-pt-pt pp p1))
                  (if (< dd tolsq)
                      (setf ii p)       ; degenerate
                      (setf ii -1))))   ; closed
            (cond
              ((>= ii p)
               ;; converged
               (when (>= nn p-size)
                 (incf p-size 150)
                 (when want-pts
                   (setf pnts (%crv-geom-grow-rows pnts (1+ p-size))))
                 (when want-pars
                   (let ((new (make-real-vector (1+ p-size))))
                     (replace new (the real-vector parms))
                     (setf parms new))))
               ;; load the point and parameter
               (setf dd (dist-pt-pt p2 last-pt))
               (when (= top 0)
                 (cond ((= nn 0)
                        (setf dd tol))
                       ((<= dd (* 0.01d0 tol))
                        (decf nn)
                        (setf dd tol))))
               (when (> dd (* 0.01d0 tol))
                 (incf nn)
                 (when want-pts
                   (%pt-row p2 pnts nn))
                 (when want-pars
                   (setf (aref (the real-vector parms) nn) (aref (the real-matrix pm-stack) top 1)))
                 (copy-pt p2 last-pt))
               (decf top))
              (t
               ;; split the segment
               (when (>= top s-size)
                 (let ((new (make-array (+ s-size 75 1))))
                   (replace new cp-stack)
                   (do ((k (1+ s-size) (1+ k))) ((> k (+ s-size 75)))
                     (declare (type fixnum k))
                     (setf (svref new k) (make-cpoint-array (1+ p))))
                   (setf cp-stack new))
                 (incf s-size 75)
                 (when want-pars
                   (setf pm-stack (%crv-geom-grow-rows pm-stack (1+ s-size)))))
               (let ((a (svref cp-stack top)) (b (svref cp-stack (1+ top))))
                 (declare (type cpoint-array a b))
                 (%copy-c-pt-row a 0 b 0)
                 (do ((kk 1 (1+ kk))) ((> kk p))
                   (declare (type fixnum kk))
                   (dotimes (i (1+ (- p kk)))
                     (%combine-2-c-pt-rows 0.5d0 a i 0.5d0 a (1+ i) a i))
                   (%copy-c-pt-row a 0 b kk)))
               (when want-pars
                 (let ((pm pm-stack))
                   (declare (type real-matrix pm))
                   (setf (aref pm (1+ top) 0) (aref pm top 0))
                   (setf (aref pm top 0) (* 0.5d0 (+ (aref pm top 0) (aref pm top 1))))
                   (setf (aref pm (1+ top) 1) (aref pm top 0))))
               (incf top)))))
        ;; the output
        (values (when want-pts
                  (let ((lpts (make-point-array (1+ nn))))
                    (dotimes (i (1+ nn) lpts)
                      (%copy-pt-row pnts i lpts i))))
                (when want-pars
                  (subseq (the real-vector parms) 0 (1+ nn)))
                nn)))))
