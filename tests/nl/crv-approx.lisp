;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-approx.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L6-g1

(defun %crv-approx-curve ()
  "A random curve for N_ApproxNurbsWithNonRatCrv: degree up to 5, a few
spans, rational or not, 2D or not, its control points within 10 of the
origin; now and then a single Bezier segment, a straight line (control
points on one), or a degenerate end (two equal end control points, a zero
end tangent)."
  (let* ((c (random-curve :max-degree 5 :max-extra 5))
         (pw (curve-cpoints c))
         (n+1 (array-dimension pw 0))
         (p (curve-degree c)))
    ;; scale down (homogeneous coordinates scale alike)
    (dotimes (i n+1)
      (dotimes (k 3)
        (unless (and (= k 2) (= (aref pw i 2) +no-z+))
          (setf (aref pw i k) (* 0.1d0 (aref pw i k))))))
    (case (rnd 12)
      (0 ;; a Bezier segment
       (let ((knots (curve-knots c)) (q (make-cpoint-array (1+ p))) (u (make-real-vector (* 2 (1+ p)))))
         (dotimes (i (1+ p))
           (setf (aref u i) (aref knots 0)
                 (aref u (+ p 1 i)) (aref knots (1- (length knots))))
           (dotimes (k 4) (setf (aref q i k) (aref pw i k))))
         (setf c (%make-curve p u q))))
      (1 ;; a straight line
       ;; from the origin to the first control point's Euclidean position
       (let* ((w0 (if (= (aref pw 0 3) +no-w+) 1d0 (aref pw 0 3)))
              (e (list (/ (aref pw 0 0) w0) (/ (aref pw 0 1) w0)
                       (if (= (aref pw 0 2) +no-z+) 0d0 (/ (aref pw 0 2) w0)))))
         (dotimes (i n+1)
           (let ((w (if (= (aref pw i 3) +no-w+) 1d0 (aref pw i 3)))
                 (s (/ (float i 1d0) (float (1- n+1) 1d0))))
             (dotimes (k 3)
               (unless (and (= k 2) (= (aref pw i 2) +no-z+))
                 (setf (aref pw i k) (* w (* s (elt e k))))))))))
      (2 ;; a zero end tangent; not for degree 1, whose first span would
         ;; then be a single point: its coincident samples give chord-length
         ;; parameters, so knots, of multiplicity above the degree, on which
         ;; N_CrvRemoveKnotMaxErr indexes outside its scratch array (Rw[p-s+2])
       (when (> p 1)
         (dotimes (k 4) (setf (aref pw 1 k) (aref pw 0 k))))))
    c))

(define-oracle-test approx-nurbs-with-non-rat-crv (:file "crv-approx" :trials 3000)
  ;; TOL positive (a negative one makes pow() a NaN the C casts to an index),
  ;; or now and then 0, which the division check rejects
  (let* ((c (%crv-approx-curve))
         (tol (case (rnd 50)
                (0 0d0)
                (t (case (rnd 6) (0 1d-3) (1 1d-2) (2 0.1d0) (3 1d0) (4 5d0) (t (rnd-range 1d-2 1d0))))))
         (deg (case (rnd 10) (0 1) (1 2) (2 2) (3 3) (4 3) (5 4) (6 5) (t (1+ (rnd 3)))))
         (par (case (rnd 25)
                (0 0)
                (1 7)
                (t (elt (list +uniform+ +chordlength+ +centripetal+ +inherited+) (rnd 4)))))
         (tn (if (zerop (rnd 2)) +yes+ +no+)))
    (agree (list :curve c :tol tol :deg deg :par par :tan tn)
           (list (batten::approx-nurbs-with-non-rat-crv c tol deg par tn))
           (with-foreign ((cp c :curve))
             (with-stacks (sg)
               (with-foreign-out ((q :curve))
                 (let ((failed (c-flag-error-p (c-approx-nurbs-with-non-rat-crv cp tol deg par tn q sg))))
                   (values failed (unless failed (list (from-foreign q :curve)))))))))))

;;; ------------------------------------------------------------- L12-g0

(defun %crv-approx-arc-point (x0 y0 x1 y1 w x2 y2 tt)
  "The point at TT of the rational quadratic with control points (X0 Y0),
(X1 Y1) of weight W and (X2 Y2), as (values x y)."
  (let* ((b0 (* (- 1d0 tt) (- 1d0 tt)))
         (b1 (* 2d0 tt (- 1d0 tt) w))
         (b2 (* tt tt))
         (s (+ b0 b1 b2)))
    (values (/ (+ (* b0 x0) (* b1 x1) (* b2 x2)) s)
            (/ (+ (* b0 y0) (* b1 y1) (* b2 y2)) s))))

(defun %crv-approx-biarcs ()
  "A chain of N2 circular arcs in the xy plane, as (values pw n arcs): PW the
2*N2+1 control points (z +no-z+ or 0; each arc's middle weight cos phi for an
arc turning 2 phi from its chord, rational, or a straight segment, collinear
and now and then non-rational), N = 2*N2, ARCS a list of (x0 y0 x1 y1 w x2
y2) per arc.  The tangent carries over from arc to arc most of the time
(tangent continuous) and turns at random otherwise; now and then N2 is 0."
  (let* ((n2 (if (zerop (rnd 30)) 0 (1+ (rnd 4))))
         (n (* 2 n2))
         (pw (make-cpoint-array (1+ n)))
         (no-z (zerop (rnd 2)))
         (rational (not (zerop (rnd 6))))
         (x (rnd-range -10d0 10d0)) (y (rnd-range -10d0 10d0))
         (th (rnd (* 2 pi)))
         (arcs nil))
    (flet ((put (i px py w)
             (setf (aref pw i 0) (* (if rational w 1d0) px)
                   (aref pw i 1) (* (if rational w 1d0) py)
                   (aref pw i 2) (if no-z +no-z+ 0d0)
                   (aref pw i 3) (if rational w +no-w+))))
      (put 0 x y 1d0)
      (dotimes (a n2)
        ;; the chord leaves the tangent direction TH at angle PHI
        (let* ((phi (case (rnd 8)
                      (0 0d0)
                      (t (* (if (zerop (rnd 2)) 1d0 -1d0) (rnd-range 0.05d0 1.45d0)))))
               (phi (if rational phi 0d0))
               (len (rnd-range 0.5d0 8d0))
               (ch (+ th phi))
               (x2 (+ x (* len (cos ch))))
               (y2 (+ y (* len (sin ch))))
               (hl (/ (* 0.5d0 len) (cos phi)))
               (x1 (+ x (* hl (cos th))))
               (y1 (+ y (* hl (sin th))))
               (w (cos phi)))
          (put (+ (* 2 a) 1) x1 y1 w)
          (put (+ (* 2 a) 2) x2 y2 1d0)
          (push (list x y x1 y1 w x2 y2) arcs)
          (setf x x2 y y2
                th (if (zerop (rnd 4)) (rnd (* 2 pi)) (+ th (* 2d0 phi)))))))
    (values pw n (nreverse arcs))))

(defun %crv-approx-polygon (arcs eps)
  "Polygon vertices along the arcs ARCS for tolerance EPS, as (values p k):
points on the arcs in order, densely or sparsely, moved off them by up to
0.6 EPS for most polygons and by up to 1.2 EPS (or 3 EPS, or far away) for
some, now and then exactly an arc end; a single vertex now and then."
  (let* ((na (max 1 (length arcs)))
         (k (case (rnd 10)
              (0 0)
              (1 (1+ (rnd 4)))
              (t (+ (* na (+ 1 (ceiling (exp (rnd (log 150d0)))))) (rnd 3)))))
         (p (make-point-array (1+ k)))
         (bad (zerop (rnd 3))))
    (dotimes (i (1+ k))
      (let* ((s (if (= k 0) (rnd 1d0) (/ (float i 1d0) (float k 1d0))))
             (ai (min (1- na) (floor (* s na))))
             (tt (- (* s na) ai))
             (arc (nth ai arcs)))
        (when (and arc (zerop (rnd 12)))
          ;; exactly the arc's start
          (setf tt 0d0))
        (multiple-value-bind (px py)
            (cond ((null arc) (values (rnd-range -10d0 10d0) (rnd-range -10d0 10d0)))
                  ((= tt 0d0) (values (first arc) (second arc)))
                  (t (apply #'%crv-approx-arc-point (append arc (list tt)))))
          (let* ((off (cond ((zerop (rnd 5)) 0d0)
                            ;; at an arc end a point off the arc sits on the
                            ;; wedge's edge, where rounding decides
                            ((and (not bad) (or (< tt 1d-3) (> tt 0.999d0))) 0d0)
                            ((not bad) (rnd (* 0.6d0 eps)))
                            ((zerop (rnd 6)) (* 3d0 eps))
                            (t (rnd (* 1.2d0 eps)))))
                 (far (and bad (zerop (rnd 40)))))
            ;; off along the arc's normal mostly (staying inside its wedge),
            ;; in any direction for the bad polygons now and then
            (multiple-value-bind (nx ny)
                (if (or (null arc) (and bad (zerop (rnd 3))))
                    (let ((ang (rnd (* 2 pi)))) (values (cos ang) (sin ang)))
                    (multiple-value-bind (ax ay)
                        (apply #'%crv-approx-arc-point (append arc (list (max 0d0 (- tt 1d-6)))))
                      (multiple-value-bind (bx by)
                          (apply #'%crv-approx-arc-point (append arc (list (min 1d0 (+ tt 1d-6)))))
                        (let ((l (sqrt (+ (expt (- bx ax) 2) (expt (- by ay) 2))))
                              (sg (if (zerop (rnd 2)) 1d0 -1d0)))
                          (values (* sg (/ (- ay by) l)) (* sg (/ (- bx ax) l)))))))
              (setf (aref p i 0) (+ px (* off nx) (if far 5d0 0d0))
                    (aref p i 1) (+ py (* off ny))
                    (aref p i 2) 0d0))))))
    (values p k)))

(define-oracle-test crv-is-polygon-within-tol (:file "crv-approx" :trials 4000)
  (multiple-value-bind (pw n arcs) (%crv-approx-biarcs)
    (when (and (> n 0) (zerop (rnd 30)))
      ;; a degenerate arc: its middle control point on its start (a straight
      ;; segment), or its end on its start (no circle)
      (let ((a (* 2 (rnd (floor n 2)))))
        (dotimes (q 4)
          (setf (aref pw (+ a (if (zerop (rnd 2)) 1 2)) q) (aref pw a q)))))
    (let ((eps (case (rnd 5) (0 1d-3) (1 1d-2) (2 0.1d0) (3 1d0) (t (rnd-range 1d-3 0.5d0)))))
      (multiple-value-bind (p k) (%crv-approx-polygon arcs eps)
        (agree (list :pw pw :n n :p p :k k :eps eps)
               (list (batten::crv-is-polygon-within-tol pw n p k eps))
               (with-foreign ((cpw pw :cpoint-array) (cp p :point-array))
                 (with-scalars ((aef :short))
                   (let ((failed (c-flag-error-p (c-crv-is-polygon-within-tol cpw n cp k eps aef))))
                     (values failed (unless failed (list (cffi:mem-ref aef :short))))))))))))
