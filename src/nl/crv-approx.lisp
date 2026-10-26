;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvApprox.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; ------------------------------------------------------------- L6-g1

(defun approx-nurbs-with-non-rat-crv (curp tol deg par tn)
  "A fresh non-rational curve of degree DEG approximating curve CURP to
within TOL: each Bezier segment of CURP is sampled at a number of evenly
spaced parameters set by TOL and its second derivative bound, the samples are
interpolated (TN +yes+: with CURP's end tangents, scaled to its arc length
unless PAR is +inherited+), parametrized by PAR (+uniform+, +chordlength+,
+centripetal+, or +inherited+: CURP's own sample parameters), and then as
many knots are removed as TOL allows (keeping the end derivative knots for a
tangent-constrained quadratic).  Unless PAR is +inherited+, the result is
reparametrized onto CURP's domain.  A 2D CURP gives a 2D curve.  TOL must be
positive.  Signals :num-err when 1/TOL is out of range, and what the
decomposition, the sampling, the parametrization (:cal-err for another PAR),
the interpolation (:inp-err for TN +yes+ and DEG 1), the normalization of a
zero end tangent and the knot removal signal.
NLib: N_ApproxNurbsWithNonRatCrv."
  (declare (type curve curp) (type double-float tol) (type fixnum deg par tn))
  (let* ((mp (curve-m curp))
         (up (curve-knots curp))
         (p (curve-degree curp))
         (twod (not (crv-is-3d curp)))
         (ex (if (= deg 1) 0.5d0 0.34d0))
         (pct 0.01d0)
         (f1 0d0) (f2 0d0))
    (declare (type fixnum mp p) (type real-vector up) (type double-float ex pct f1 f2))
    (when (float-op-is-bad 1d0 tol +division+)
      (nurbs-error :num-err 'approx-nurbs-with-non-rat-crv))
    (setf f1 (/ 1d0 tol))
    (setf f1 (c-pow f1 ex))
    (setf f2 (float deg 1d0))
    ;; the Bezier segments
    (multiple-value-bind (bez k) (crv-decompose-bez curp)
      (declare (type simple-vector bez) (type fixnum k))
      ;; the number of sampling points per segment
      (let ((nu (make-index-vector (1+ k)))
            (mu 0)
            (unitspan (make-real-vector 2)))
        (declare (type index-vector nu) (type fixnum mu) (type real-vector unitspan))
        (create-interval unitspan 0d0 1d0)
        (dotimes (i (1+ k))
          (crv-reparam-to-interval (svref bez i) unitspan)
          (let* ((muu (crv-get-max-2nd-deriv (svref bez i)))
                 (num (* 0.125d0 muu))
                 (gro (nl-max f2 (c-sqrt num)))
                 (del (* f1 gro)))
            (declare (type double-float muu num gro del))
            (setf del (fceiling del))
            (setf (aref nu i) (truncate del))
            (setf mu (+ mu (aref nu i)))))
        ;; the parameters
        (let ((u (make-real-vector (1+ mu)))
              (ul (aref up p))
              (ur 0d0)
              (i (1+ p))
              (j 0)
              (l 0))
          (declare (type real-vector u) (type double-float ul ur) (type fixnum i j l))
          (setf (aref u 0) ul)
          (do () ((not (< i mp)))
            (do () ((not (and (< i mp) (= (aref up i) (aref up (1+ i))))))
              (incf i))
            (setf ur (aref up i))
            (let ((uinc (/ (- ur ul) (float (aref nu j) 1d0))))
              (declare (type double-float uinc))
              (do ((kk 1 (1+ kk))) ((not (< kk (aref nu j))))
                (declare (type fixnum kk))
                (setf (aref u (incf l)) (+ ul (* (float kk 1d0) uinc))))
              (setf (aref u (incf l)) ur))
            (setf ul ur)
            (incf i)
            (incf j))
          ;; sample the curve
          (let ((pts (make-point-array (1+ mu)))
                (pt (make-point))
                (knu (make-real-vector (if (= tn +yes+) (+ mu deg 4) (+ mu deg 2))))
                (curq nil))
            (declare (type point-array pts) (type point pt) (type real-vector knu))
            (dotimes (ii (1+ mu))
              (crv-eval curp (aref u ii) +left+ pt)
              (%pt-row pt pts ii))
            ;; the knot vector
            (when (/= par +inherited+)
              (fit-calc-crv-param-values pts mu +epoint+ par u))
            (if (= tn +yes+)
                (fit-calc-knot-vector-end-derivs u mu deg knu)
                (fit-crv-calc-knot-vector u mu deg knu))
            ;; interpolate the points
            (cond
              ((= tn +yes+)
               ;; the end tangents
               (let ((cd (make-point-array 2))
                     (ts (make-point))
                     (te (make-point)))
                 (declare (type point-array cd) (type point ts te))
                 (crv-derivs curp (aref up 0) +left+ 1 cd)
                 (%row-pt cd 1 ts)
                 (crv-derivs curp (aref up mp) +left+ 1 cd)
                 (%row-pt cd 1 te)
                 (when (/= par +inherited+)
                   ;; adjust the magnitudes
                   (let ((del (crv-arc-length curp (aref up 0) (aref up mp) pct +relative+)))
                     (declare (type double-float del))
                     (vector-normalize-ref ts)
                     (vector-normalize-ref te)
                     (vector-scale ts del ts)
                     (vector-scale te del te)))
                 (setf curq (fit-crv-knots-and-derivs pts mu +epoint+ u knu deg ts te))))
              (t
               (setf curq (fit-crv-interp-given-params pts mu +epoint+ u knu deg))))
            (when twod
              (crv-3d-to-2d curq))
            ;; remove knots
            (let ((der (if (and (= tn +yes+) (= deg 2)) 1 0)))
              (declare (type fixnum der))
              (setf curq (crv-remove-knots-deriv-constraints curq tol +both+ der)))
            ;; compact the output curve
            (when (< (curve-n curq) mu)
              (crv-compress curq))
            ;; rescale the knot vector if necessary
            (when (/= par +inherited+)
              (let* ((uq (curve-knots curq))
                     (mq (1- (length uq))))
                (declare (type real-vector uq) (type fixnum mq))
                (when (or (/= (aref up 0) (aref uq 0)) (/= (aref up mp) (aref uq mq)))
                  (let ((iv (make-real-vector 2)))
                    (create-interval iv (aref up 0) (aref up mp))
                    (crv-reparam-to-interval curq iv)))))
            curq))))))

;;; ------------------------------------------------------------- L12-g0

(defun crv-is-polygon-within-tol (pw n p k eps)
  "Whether the polygon P[0..K] (a point-array) lies within EPS of the biarc
curve whose control points are PW[0..N] (N even; each triple PW[2i],
PW[2i+1], PW[2i+2] a circular arc, or a straight segment when collinear):
every vertex must be within EPS of its closest arc, and every polygon leg
within EPS of the arcs closest to its two ends.  Returns +yes+ or +no+.
Signals what bez-center-radius-from-3-c-pts, dist-pt-arc and
dist-max-line-arc signal.  Where dist-max-line-arc leaves its distance unset
the C reads a stale local; the Lisp gets the 0 that dist-max-line-arc
returns there.
NLib: N_CrvIsPolygonWithinTol."
  (declare (type cpoint-array pw) (type fixnum n k) (type point-array p)
           (type double-float eps))
  (let* ((n2 (truncate n 2))
         (cs (make-array (max n2 0)))
         (pss (make-array (max n2 0)))
         (pes (make-array (max n2 0)))
         (r (make-real-vector (max n2 0)))
         (alf (make-real-vector (max n2 0)))
         (bet (make-real-vector (max n2 0)))
         (ind (make-index-vector (1+ k)))
         (p0w (make-cpoint)) (p1w (make-cpoint)) (p2w (make-cpoint))
         (pt (make-point)) (qt (make-point)))
    (declare (type fixnum n2) (type simple-vector cs pss pes)
             (type real-vector r alf bet) (type index-vector ind)
             (type cpoint p0w p1w p2w) (type point pt qt)
             (dynamic-extent p0w p1w p2w pt qt))
    ;; center and radius of every arc; a straight segment has r = +infinite+
    (dotimes (i n2)
      (let ((j (* 2 i))
            (c (make-point)) (ps (make-point)) (pe (make-point)))
        (declare (type fixnum j) (type point c ps pe))
        (setf (svref cs i) c (svref pss i) ps (svref pes i) pe)
        (multiple-value-bind (ri ai bi)
            (bez-center-radius-from-3-c-pts (%row-c-pt pw j p0w) (%row-c-pt pw (+ j 1) p1w)
                                            (%row-c-pt pw (+ j 2) p2w) c ps pe)
          (declare (type double-float ri ai bi))
          (setf (aref r i) ri (aref alf i) ai (aref bet i) bi))))
    ;; for each polygon vertex, the closest arc
    (dotimes (i (1+ k))
      (let ((dmin +bigd+))
        (declare (type double-float dmin))
        (setf (aref ind i) -1)
        (%row-pt p i pt)
        (dotimes (j n2)
          (let ((d (dist-pt-arc pt (svref cs j) (aref r j) (aref alf j) (aref bet j)
                                (svref pss j) (svref pes j))))
            (declare (type double-float d))
            (when (/= d +undefined+)
              (when (< d dmin)
                (setf dmin d
                      (aref ind i) j)))))
        (when (or (< (aref ind i) 0) (> dmin eps))
          (return-from crv-is-polygon-within-tol +no+))))
    ;; for each polygon leg, its distance from the arcs at its ends
    (dotimes (i k)
      (let ((dmax -1d0))
        (declare (type double-float dmax))
        (%row-pt p i pt)
        (%row-pt p (1+ i) qt)
        (dotimes (j n2)
          (unless (and (/= j (aref ind i)) (/= j (aref ind (1+ i))))
            (let ((d (dist-max-line-arc (svref cs j) (aref r j) (aref alf j) (aref bet j)
                                        (svref pss j) (svref pes j) pt qt)))
              (declare (type double-float d))
              (when (/= d +undefined+)
                (when (> d dmax)
                  (setf dmax d))))))
        (when (or (< dmax 0d0) (> dmax eps))
          (return-from crv-is-polygon-within-tol +no+))))
    +yes+))

;;; ------------------------------------------------------------- L13-g0

(defun %crv-approx-grow-reals (r new-size)
  "A copy of real-vector R, NEW-SIZE long (N_Realloc1dRealArray)."
  (declare (type real-vector r) (type fixnum new-size))
  (let ((s (make-real-vector new-size)))
    (replace s r)
    s))

(defun %crv-approx-biarcs (curp pp u k eps)
  "The biarc fitting of N_ApproxContinuousCrvWithArcs for a curve CURP of
degree above 1, whose polyline is PP[0..K] at parameters U[0..K]: returns
(values curq n), CURQ a working curve of degree 2 whose control points
0..N hold the biarcs (its arrays may be longer; the knots are unset)."
  (declare (type curve curp) (type point-array pp) (type real-vector u)
           (type fixnum k) (type double-float eps))
  (let* ((tt (make-point-array (1+ k)))
         (cl (make-point-array 2))
         (cr (make-point-array 2))
         (q (make-point-array 2))
         (bw (make-cpoint-array 18))
         (ps (make-point)) (ts (make-point)) (pe (make-point)) (te (make-point))
         (std 20)
         (ul (make-real-vector (1+ std)))
         (ur (make-real-vector (1+ std)))
         ;; initial size: two biarcs per polyline leg
         (n (* 8 k))
         (curq (%make-curve 2 (make-real-vector (+ n 4)) (make-cpoint-array (1+ n))))
         (closed nil) (frt nil) (sbd nil)
         (l 0) (le 0) (ks 0) (ke 0) (kl 0) (kr 0) (kd k) (top 0)
         (um 0d0) (u1 0d0) (u2 0d0))
    (declare (type point-array tt cl cr q) (type cpoint-array bw) (type point ps ts pe te)
             (type fixnum std n l le ks ke kl kr kd top) (type real-vector ul ur)
             (type double-float um u1 u2) (dynamic-extent ps ts pe te))
    ;; the tangents at the polyline points
    (dotimes (i (1+ k))
      (crv-derivs curp (aref u i) +left+ 1 cl)
      (%copy-pt-row cl 1 tt i))
    ;; open or closed?
    (setf closed (< (dist-pt-pt (%row-pt pp 0 ps) (%row-pt pp k pe)) *mtol*))
    (setf n 0)
    (flet ((emit ()
             ;; Bw[0..l] into Pw[n..n+l], expanding the arrays as the C does
             (let ((cur-n (curve-n curq)))
               (declare (type fixnum cur-n))
               (when (> (+ n l) cur-n)
                 ;; the C expands to 2*cur_n and would write past it if that
                 ;; were still short
                 (let ((new-n (max (* 2 cur-n) (+ n l))))
                   (declare (type fixnum new-n))
                   (crv-expand curq new-n (+ new-n 3)))))
             (let ((pw (curve-cpoints curq)))
               (declare (type cpoint-array pw))
               (dotimes (j (1+ l))
                 (%copy-c-pt-row bw j pw (+ n j))))))
      (do () ((not (< ks k)))
        (setf ke (if closed (truncate k 2) (nl-min (+ ks kd) k)))
        (setf kl ks kr ke frt t sbd nil)
        (do () ((not (/= kl ke)))
          ;; a biarc from P[ks] to P[ke]
          (setf l (fit-arc-to-end-pts-and-tangents (%row-pt pp ks ps) (%row-pt tt ks ts)
                                                   (%row-pt pp ke pe) (%row-pt tt ke te) bw))
          (if (= (let ((poly (make-point-array (1+ (- ke ks)))))
                   (dotimes (j (1+ (- ke ks)))
                     (%copy-pt-row pp (+ ks j) poly j))
                   (crv-is-polygon-within-tol bw l poly (- ke ks) eps))
                 +yes+)
              ;; success: keep it, and try further unless bisecting
              (progn
                (emit)
                (setf kl ke le l)
                (when (and frt (< ke k) (not closed))
                  (setf ke (nl-min (+ ke kd) k)
                        kr ke)))
              (if (= (- ke ks) 1)
                  ;; a single leg fails: subdivide the curve over it
                  (progn
                    (setf um (* 0.5d0 (+ (aref u ks) (aref u ke)))
                          sbd t)
                    (setf (aref ul 0) um (aref ur 0) (aref u ke)
                          (aref ul 1) (aref u ks) (aref ur 1) um
                          top 1)
                    (do () ((not (>= top 0)))
                      (setf u1 (aref ul top) u2 (aref ur top))
                      (decf top)
                      (crv-derivs curp u1 +left+ 1 cl)
                      (crv-derivs curp u2 +left+ 1 cr)
                      (setf l (fit-arc-to-end-pts-and-tangents (%row-pt cl 0 ps) (%row-pt cl 1 ts)
                                                               (%row-pt cr 0 pe) (%row-pt cr 1 te) bw))
                      (%copy-pt-row cl 0 q 0)
                      (%copy-pt-row cr 0 q 1)
                      (if (= (crv-is-polygon-within-tol bw l q 1 eps) +yes+)
                          (progn
                            (emit)
                            (incf n l))
                          (progn
                            (setf um (* 0.5d0 (+ u1 u2)))
                            (when (> (+ top 2) std)
                              (setf ul (%crv-approx-grow-reals ul (1+ (+ std std)))
                                    ur (%crv-approx-grow-reals ur (1+ (+ std std))))
                              (incf std std))
                            (setf (aref ul (+ top 1)) um (aref ur (+ top 1)) u2
                                  (aref ul (+ top 2)) u1 (aref ur (+ top 2)) um)
                            (incf top 2))))
                    (setf kl ke))
                  ;; several legs fail: bisect
                  (setf kr ke
                        frt nil)))
          (unless frt
            (setf ke (truncate (+ kl kr) 2))))
        (unless sbd
          (incf n le))
        (setf kd (nl-max 1 (- ke ks)))
        (setf ks ke
              closed nil)))
    (values curq n)))

(defun approx-continuous-crv-with-arcs (curp tol)
  "A fresh piecewise circular (biarc) curve approximating curve CURP, which
must be xy-planar and C1 (interior knot multiplicities below its degree), to
within the absolute tolerance TOL: a degree-2 curve whose control point
triples are rational Bezier arcs, with double interior knots at chord-length
parameters of its arc ends.  CURP is first replaced by a polyline within
TOL/2; a degree-1 CURP gives the quadratic through the polyline points and
leg midpoints, otherwise biarcs are fitted from point to point (with the
polyline tangents), over as many legs as stay within TOL/2 of the polyline,
subdividing the curve itself where a single leg fails.  Signals what
approx-crv-with-polyline, crv-derivs, fit-arc-to-end-pts-and-tangents
(:geo-err for coincident end points), crv-is-polygon-within-tol and
fit-calc-crv-param-values signal.
NLib: N_ApproxContinuousCrvWithArcs."
  (declare (type curve curp) (type double-float tol))
  (let ((p (curve-degree curp))
        (eps (* 0.5d0 tol))
        (curq nil)
        (n 0))
    (declare (type fixnum p n) (type double-float eps))
    ;; the polyline
    (multiple-value-bind (pp u k) (approx-crv-with-polyline curp eps +absolute+ +both+)
      (declare (type point-array pp) (type real-vector u) (type fixnum k))
      (if (= p 1)
          ;; degree 1: the polyline points and the leg midpoints
          (let ((a (make-point)) (b (make-point)) (pm (make-point)) (cw (make-cpoint))
                (j 0))
            (declare (type point a b pm) (type cpoint cw) (type fixnum j)
                     (dynamic-extent a b pm cw))
            (setf n (* 2 k))
            (setf curq (%make-curve 2 (make-real-vector (+ n 4)) (make-cpoint-array (1+ n))))
            (let ((pw (curve-cpoints curq)))
              (declare (type cpoint-array pw))
              (%c-pt-row (weight (%row-pt pp 0 a) 1d0 cw) pw 0)
              (dotimes (i k)
                (%row-pt pp i a)
                (%row-pt pp (1+ i) b)
                (combine-2-pts 0.5d0 a 0.5d0 b pm)
                (%c-pt-row (weight pm 1d0 cw) pw (incf j))
                (%c-pt-row (weight b 1d0 cw) pw (incf j)))))
          (setf (values curq n) (%crv-approx-biarcs curp pp u k eps))))
    ;; the knot vector: chord-length parameters of the arc ends
    (let* ((l (truncate n 2))
           (tv (make-real-vector (1+ l)))
           (r (make-point-array (1+ l)))
           (pw (curve-cpoints curq))
           ;; the C writes 2l+4 knots (6 for l = 0)
           (uq (make-real-vector (max (+ n 4) (+ (* 2 l) 4))))
           (cw (make-cpoint)) (pt (make-point))
           (m -1))
      (declare (type fixnum l m) (type real-vector tv uq) (type point-array r)
               (type cpoint-array pw) (type cpoint cw) (type point pt)
               (dynamic-extent cw pt))
      (dotimes (i (1+ l))
        (c-pt-to-pt-euclid (%row-c-pt pw (* 2 i) cw) pt)
        (%pt-row pt r i))
      (fit-calc-crv-param-values r l +epoint+ +chordlength+ tv)
      (dotimes (i 3)
        (setf (aref uq (incf m)) (aref tv 0)))
      (do ((i 1 (1+ i))) ((not (< i l)))
        (declare (type fixnum i))
        (setf (aref uq (incf m)) (aref tv i))
        (setf (aref uq (incf m)) (aref tv i)))
      (dotimes (i 3)
        (setf (aref uq (incf m)) (aref tv l)))
      ;; the output arrays at their exact sizes (N_CrvSetSizeIndices and
      ;; N_CrvCompress)
      (let ((qw (make-cpoint-array (1+ n)))
            (v (make-real-vector (+ n 4))))
        (declare (type cpoint-array qw) (type real-vector v))
        (dotimes (i (1+ n))
          (%copy-c-pt-row pw i qw i))
        (replace v uq)
        (%make-curve 2 v qw)))))
