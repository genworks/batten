;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_SrfGeom.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

(defun srf-get-average-len (sur)
  "Estimate the average lengths of surface SUR in u and v from its control
net, and a bound on its area as the product of the longest control polygon
in each direction.  Returns (values a lu lv).  NLib: N_SrfGetAverageLen."
  (declare (type surface sur))
  (let* ((pw (surface-cpoints sur))
         (n (surface-n sur))
         (m (surface-m sur))
         (cp (make-cpoint))
         (prev (make-point))
         (cur (make-point))
         (maxu 0d0) (sumu 0d0) (maxv 0d0) (sumv 0d0))
    (declare (type cpoint-net pw) (type fixnum n m) (type cpoint cp) (type point prev cur)
             (type double-float maxu sumu maxv sumv)
             (dynamic-extent cp prev cur))
    (flet ((euclid (i j p)
             (declare (type fixnum i j) (type point p))
             (set-cpoint cp (aref pw i j 0) (aref pw i j 1) (aref pw i j 2) (aref pw i j 3))
             (c-pt-to-pt-euclid cp p)))
      (declare (inline euclid))
      (dotimes (j (1+ m))
        (euclid 0 j prev)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((i 1 (1+ i))) ((> i n))
            (declare (type fixnum i))
            (euclid i j cur)
            (incf len (dist-pt-pt cur prev))
            (rotatef prev cur))
          (incf sumu len)
          (when (> len maxu) (setf maxu len))))
      (dotimes (i (1+ n))
        (euclid i 0 prev)
        (let ((len 0d0))
          (declare (type double-float len))
          (do ((j 1 (1+ j))) ((> j m))
            (declare (type fixnum j))
            (euclid i j cur)
            (incf len (dist-pt-pt cur prev))
            (rotatef prev cur))
          (incf sumv len)
          (when (> len maxv) (setf maxv len)))))
    (values (* maxu maxv)
            (/ sumu (+ (float m 1d0) 1d0))
            (/ sumv (+ (float n 1d0) 1d0)))))

;;; ------------------------------------------------------------ level 1, group 4

(defun get-closest-pt-on-srf (sur p u0 v0 top toc q d &optional uv)
  "Find the point of surface SUR closest to point P by Newton iteration from
the guess (U0, V0), on objective functions built from the unit first
derivatives, with the second-derivative terms limited and the parameters
wrapped (closed directions) or clamped (open ones).  TOP is the point
coincidence tolerance, TOC the zero cosine tolerance.  Writes the best point
found into Q and the last derivatives into D, a point-net of at least 3 x 3
(D[k][l] differentiated k times in u, l times in v, k + l <= 2).  Returns
(values u v), the parameters of Q.

Signals :con-err when the iteration fails or does not converge within
*itlim* steps, and what srf-derivs and srf-eval-pt signal.  The C returns the
best u, v and Q found even then: Q is written as it improves, and so is UV,
an optional real-vector of 2 receiving u and v, for callers that carry on
after the error.
NLib: N_GetClosestPtOnSrf."
  (declare (type surface sur) (type point p q) (type double-float u0 v0 top toc)
           (type point-net d))
  (let* ((knu (surface-knots-u sur)) (knv (surface-knots-v sur))
         (r (1- (length knu))) (s (1- (length knv)))
         (uclosed (srf-is-closed sur +udir+))
         (vclosed (srf-is-closed sur +vdir+))
         (k 0) (obndy 0) (ucount 0) (vcount 0)
         (itlim *itlim*)
         (u 0d0) (v 0d0)
         (unew u0) (vnew v0) (uold 0d0) (vold 0d0) (ub 0d0) (vb 0d0)
         (best-dis +bigd+)
         (a (make-point)) (b (make-point)) (gap (make-point))
         (ss (make-point))
         (ssu (make-point)) (ssv (make-point)) (ssuu (make-point)) (ssuv (make-point)) (ssvv (make-point))
         (ssun (make-point)) (ssvn (make-point))
         (ssunu (make-point)) (ssunv (make-point)) (ssvnu (make-point)) (ssvnv (make-point))
         (tmp1 (make-point)) (tmp2 (make-point)))
    (declare (type real-vector knu knv) (type fixnum r s k obndy ucount vcount itlim)
             (type double-float u v unew vnew uold vold ub vb best-dis)
             (dynamic-extent a b gap ss ssu ssv ssuu ssuv ssvv ssun ssvn ssunu ssunv ssvnu ssvnv
                             tmp1 tmp2))
    (flet ((con-err () (nurbs-error :con-err 'get-closest-pt-on-srf)))
      (do () ((>= k itlim))
        ;; derivatives
        (srf-derivs sur unew vnew +left+ +left+ +true+ 2 2 d)
        (%net-pt d 0 0 ss)
        (%net-pt d 1 0 ssu)
        (%net-pt d 0 1 ssv)
        (%net-pt d 2 0 ssuu)
        (%net-pt d 1 1 ssuv)
        (%net-pt d 0 2 ssvv)
        (vector-diff ss p gap)
        ;; convergence test 1: point coincidence
        (let ((dis (vector-magnitude gap)))
          (declare (type double-float dis))
          (when (< dis best-dis)
            (setf u unew v vnew)
            (when uv (setf (aref (the real-vector uv) 0) u (aref (the real-vector uv) 1) v))
            (copy-pt ss q)
            (setf best-dis dis))
          (when (<= dis top) (return))
          (let ((dmu (vector-magnitude ssu))
                (dmv (vector-magnitude ssv)))
            (declare (type double-float dmu dmv))
            (when (float-op-is-bad 1d0 dmu +division+) (con-err))
            (when (float-op-is-bad 1d0 dmv +division+) (con-err))
            (vector-scale ssu (/ 1d0 dmu) ssun)
            (vector-scale ssv (/ 1d0 dmv) ssvn)
            (let ((fval (vector-dot gap ssun))
                  (gval (vector-dot gap ssvn)))
              (declare (type double-float fval gval))
              ;; convergence test 2: zero cosine
              (when (float-op-is-bad fval (* dmu dis) +division+) (con-err))
              (when (float-op-is-bad gval (* dmv dis) +division+) (con-err))
              (let ((zcu (/ (abs fval) (* dmu dis)))
                    (zcv (/ (abs gval) (* dmv dis))))
                (declare (type double-float zcu zcv))
                (when (and (<= zcu toc) (<= zcv toc) (= dis best-dis)) (return)))
              ;; derivatives of the objective functions
              (let ((fu 0d0) (fv 0d0) (gu 0d0) (gv 0d0))
                (declare (type double-float fu fv gu gv))
                (let* ((dmuu (vector-dot ssun ssuu)))
                  (declare (type double-float dmuu))
                  (vector-scale ssun dmuu tmp1)
                  (vector-diff ssuu tmp1 tmp2)
                  (vector-scale tmp2 (/ 1d0 dmu) ssunu)
                  (let ((fu1 (vector-dot ssu ssun))
                        (fu2 (vector-dot gap ssunu)))
                    (declare (type double-float fu1 fu2))
                    (when (< fu2 (* (- fu1) 0.5d0))
                      (setf fu2 (* (- fu1) 0.5d0)))
                    (setf fu (+ fu1 fu2))))
                (let ((dmvv (vector-dot ssvn ssvv)))
                  (declare (type double-float dmvv))
                  (vector-scale ssvn dmvv tmp1)
                  (vector-diff ssvv tmp1 tmp2)
                  (vector-scale tmp2 (/ 1d0 dmv) ssvnv)
                  (let ((gv1 (vector-dot ssv ssvn))
                        (gv2 (vector-dot gap ssvnv)))
                    (declare (type double-float gv1 gv2))
                    (when (< gv2 (* (- gv1) 0.5d0))
                      (setf gv2 (* (- gv1) 0.5d0)))
                    (setf gv (+ gv1 gv2))))
                (let ((dmuv (vector-dot ssun ssuv)))
                  (declare (type double-float dmuv))
                  (vector-scale ssun dmuv tmp1)
                  (vector-diff ssuv tmp1 tmp2)
                  (vector-scale tmp2 (/ 1d0 dmu) ssunv)
                  (setf fv (+ (vector-dot ssv ssun) (vector-dot gap ssunv))))
                (let ((dmvu (vector-dot ssvn ssuv)))
                  (declare (type double-float dmvu))
                  (vector-scale ssvn dmvu tmp1)
                  (vector-diff ssuv tmp1 tmp2)
                  (vector-scale tmp2 (/ 1d0 dmv) ssvnu)
                  (setf gu (+ (vector-dot ssu ssvn) (vector-dot gap ssvnu))))
                ;; Cramer's rule for the shifts
                (let ((den (- (* fu gv) (* fv gu)))
                      (num (- (* fval gv) (* fv gval)))
                      (delta-u 0d0) (delta-v 0d0))
                  (declare (type double-float den num delta-u delta-v))
                  (when (float-op-is-bad num den +division+) (con-err))
                  (setf delta-u (/ (- num) den))
                  (setf num (- (* fu gval) (* fval gu)))
                  (when (float-op-is-bad num den +division+) (con-err))
                  (setf delta-v (/ (- num) den))
                  (setf uold unew vold vnew
                        unew (+ uold delta-u)
                        vnew (+ vold delta-v)))
                ;; check 3: parameter range
                (setf obndy -1)
                (let ((u0* (aref knu 0)) (ur (aref knu r)))
                  (declare (type double-float u0* ur))
                  (cond
                    (uclosed
                     (cond ((< unew u0*)
                            (cond ((< (- u0* unew) (* 100d0 (- ur u0*)))
                                   (do () ((not (< unew u0*)))
                                     (setf unew (+ (- ur u0*) unew)))
                                   (setf ucount 0))
                                  (t
                                   (when (>= ucount 2) (con-err))
                                   (setf unew uold)
                                   (incf ucount))))
                           ((> unew ur)
                            (cond ((< (- unew ur) (* 100d0 (- ur u0*)))
                                   (do () ((not (> unew ur)))
                                     (setf unew (+ (- u0* ur) unew)))
                                   (setf ucount 0))
                                  (t
                                   (when (>= ucount 2) (con-err))
                                   (setf unew uold)
                                   (incf ucount))))
                           (t (setf ucount 0))))
                    ((or (< unew u0*) (> unew ur))
                     (when (< unew u0*) (setf unew u0*))
                     (when (> unew ur) (setf unew ur))
                     (unless (float-op-is-bad gval gv +division+)
                       (setf vb (- vold (/ gval gv)))
                       (when (and (>= vb (aref knv 0)) (<= vb (aref knv s)))
                         (setf obndy 1))))))
                (let ((v0* (aref knv 0)) (vs (aref knv s)))
                  (declare (type double-float v0* vs))
                  (cond
                    (vclosed
                     (cond ((< vnew v0*)
                            (cond ((< (- v0* vnew) (* 100d0 (- vs v0*)))
                                   (do () ((not (< vnew v0*)))
                                     (setf vnew (+ (- vs v0*) vnew)))
                                   (setf vcount 0))
                                  (t
                                   (when (>= vcount 2) (con-err))
                                   (setf vnew vold)
                                   (incf vcount))))
                           ((> vnew vs)
                            (cond ((< (- vnew vs) (* 100d0 (- vs v0*)))
                                   (do () ((not (> vnew vs)))
                                     (setf vnew (+ (- v0* vs) vnew)))
                                   (setf vcount 0))
                                  (t
                                   (when (>= vcount 2) (con-err))
                                   (setf vnew vold)
                                   (incf vcount))))
                           (t (setf vcount 0))))
                    ((or (< vnew v0*) (> vnew vs))
                     (when (< vnew v0*) (setf vnew v0*))
                     (when (> vnew vs) (setf vnew vs))
                     (unless (float-op-is-bad fval fu +division+)
                       (setf ub (- uold (/ fval fu)))
                       (when (and (>= ub (aref knu 0)) (<= ub (aref knu r)))
                         (setf obndy (if (> obndy 0) 3 2)))))))
                ;; which of the two possibilities is better
                (when (> obndy 0)
                  (when (< obndy 2) (setf ub unew))
                  (when (= obndy 2) (setf vb vnew))
                  (srf-eval-pt sur unew vnew +left+ +left+ a)
                  (srf-eval-pt sur ub vb +left+ +left+ b)
                  (let ((d1 (dist-sq-pt-pt a p))
                        (d2 (dist-sq-pt-pt b p)))
                    (declare (type double-float d1 d2))
                    (when (< d2 d1)
                      (setf unew ub vnew vb))))
                ;; check 4: parameter change
                (let ((delta-u (- unew uold))
                      (delta-v (- vnew vold)))
                  (declare (type double-float delta-u delta-v))
                  (combine-2-pts delta-u ssu delta-v ssv a)
                  (when (and (<= (pt-magnitude a) top) (= dis best-dis))
                    (return)))))))
        (incf k))
      (when (>= k itlim) (con-err))
      (values u v))))

;;; ------------------------------------------------------------- L4-g1

(declaim (inline %srf-geom-grow))
(defun %srf-geom-grow (a need)
  "A real-vector holding A's elements with room for index NEED: A itself when
it is long enough, else a copy twice as long (or longer)."
  (declare (type real-vector a) (type fixnum need))
  (if (< need (length a))
      a
      (let ((b (make-real-vector (max (* 2 (length a)) (1+ need)))))
        (replace b a)
        b)))

(defun approx-srf-with-quad-srf (sur tol tfl rfl)
  "A piecewise quadrilateral approximation of surface SUR: its parameter
rectangle is split in half, in u or v as srf-is-flat-cheap advises, until
every patch is flat to within TOL (TFL +absolute+; +relative+ scales TOL by
the diagonal of SUR's bounding box).  Returns (values P U V N M): the grid of
surface points P[0..N][0..M] (a point-net, when RFL is +points+ or +both+,
else nil), the grid parameters U[0..N] and V[0..M] (when RFL is +parameters+
or +both+, else nil) and the highest indexes N and M.  Signals what
srf-extract-patch, srf-is-flat-cheap and srf-eval-pt signal.  (The C's
growing work arrays, sized by function-static counts, only hold the stack
and the grid; their sizes do not affect the result.)
NLib: N_ApproxSrfWithQuadSrf."
  (declare (type surface sur) (type double-float tol) (type fixnum tfl rfl))
  (let* ((us (surface-knots-u sur))
         (vs (surface-knots-v sur))
         (rs (1- (length us)))
         (ss (1- (length vs)))
         (s-ul (make-real-vector 41))
         (s-ur (make-real-vector 41))
         (s-vb (make-real-vector 41))
         (s-vt (make-real-vector 41))
         (uu (make-real-vector 101))
         (vv (make-real-vector 101))
         (ts 0) (tu 1) (tv 1)
         (pts nil) (uo nil) (vo nil))
    (declare (type real-vector us vs s-ul s-ur s-vb s-vt uu vv) (type fixnum rs ss ts tu tv))
    ;; adjust the tolerance
    (when (= tfl +relative+)
      (setf tol (* (srf-max-diag-dist-b-box sur) tol)))
    ;; initialize
    (setf (aref s-ul 0) (aref us 0)
          (aref s-ur 0) (aref us rs)
          (aref s-vb 0) (aref vs 0)
          (aref s-vt 0) (aref vs ss))
    (setf (aref uu 0) (aref us 0)
          (aref uu 1) (aref us rs)
          (aref vv 0) (aref vs 0)
          (aref vv 1) (aref vs ss))
    ;; subdivide the surface
    (do () ((< ts 0))
      (let ((ul (aref s-ul ts))
            (ur (aref s-ur ts))
            (vb (aref s-vb ts))
            (vt (aref s-vt ts)))
        (declare (type double-float ul ur vb vt))
        (decf ts)
        ;; extract the patch and check its flatness
        (multiple-value-bind (flt dir) (srf-is-flat-cheap (srf-extract-patch sur ul ur vb vt) tol)
          (declare (type fixnum flt dir))
          (when (= flt +no+)
            ;; room for the new rectangles and the new grid value
            (setf s-ul (%srf-geom-grow s-ul (+ ts 2))
                  s-ur (%srf-geom-grow s-ur (+ ts 2))
                  s-vb (%srf-geom-grow s-vb (+ ts 2))
                  s-vt (%srf-geom-grow s-vt (+ ts 2)))
            (if (= dir +udir+)
                (setf vv (%srf-geom-grow vv (1+ tv)))
                (setf uu (%srf-geom-grow uu (1+ tu))))
            ;; subdivision points, and the new rectangles on the stack
            (if (= dir +udir+)
                (let ((vm (* 0.5d0 (+ vb vt)))
                      (first 0)
                      (last tv))
                  (declare (type double-float vm) (type fixnum first last))
                  (do () ((> first last))
                    (let ((mid (truncate (+ first last) 2)))
                      (declare (type fixnum mid))
                      (when (< (abs (- (aref vv mid) vm)) *ptol*)
                        (return))
                      (if (< vm (aref vv mid))
                          (setf last (1- mid))
                          (setf first (1+ mid)))))
                  (when (> first last)
                    (do ((j tv (1- j))) ((< j (1+ last)))
                      (declare (type fixnum j))
                      (setf (aref vv (1+ j)) (aref vv j)))
                    (setf (aref vv (1+ last)) vm)
                    (incf tv))
                  (setf (aref s-ul (+ ts 1)) ul
                        (aref s-ur (+ ts 1)) ur
                        (aref s-vb (+ ts 1)) vb
                        (aref s-vt (+ ts 1)) vm
                        (aref s-ul (+ ts 2)) ul
                        (aref s-ur (+ ts 2)) ur
                        (aref s-vb (+ ts 2)) vm
                        (aref s-vt (+ ts 2)) vt))
                (let ((um (* 0.5d0 (+ ul ur)))
                      (first 0)
                      (last tu))
                  (declare (type double-float um) (type fixnum first last))
                  (do () ((> first last))
                    (let ((mid (truncate (+ first last) 2)))
                      (declare (type fixnum mid))
                      (when (< (abs (- (aref uu mid) um)) *ptol*)
                        (return))
                      (if (< um (aref uu mid))
                          (setf last (1- mid))
                          (setf first (1+ mid)))))
                  (when (> first last)
                    (do ((i tu (1- i))) ((< i (1+ last)))
                      (declare (type fixnum i))
                      (setf (aref uu (1+ i)) (aref uu i)))
                    (setf (aref uu (1+ last)) um)
                    (incf tu))
                  (setf (aref s-ul (+ ts 1)) ul
                        (aref s-ur (+ ts 1)) um
                        (aref s-vb (+ ts 1)) vb
                        (aref s-vt (+ ts 1)) vt
                        (aref s-ul (+ ts 2)) um
                        (aref s-ur (+ ts 2)) ur
                        (aref s-vb (+ ts 2)) vb
                        (aref s-vt (+ ts 2)) vt)))
            (incf ts 2)))))
    ;; the output
    (when (or (= rfl +points+) (= rfl +both+))
      (let ((q (make-point-net (1+ tu) (1+ tv)))
            (s (make-point)))
        (declare (type point-net q) (dynamic-extent s))
        (dotimes (i (1+ tu))
          (dotimes (j (1+ tv))
            (srf-eval-pt sur (aref uu i) (aref vv j) +left+ +left+ s)
            (%pt-net s q i j)))
        (setf pts q)))
    (when (or (= rfl +parameters+) (= rfl +both+))
      (setf uo (make-real-vector (1+ tu))
            vo (make-real-vector (1+ tv)))
      (replace uo uu :end2 (1+ tu))
      (replace vo vv :end2 (1+ tv)))
    (values pts uo vo tu tv)))
