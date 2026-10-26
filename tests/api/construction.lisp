;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Tests of the API in src/api/construction.lisp: each function against the C routines
;;;; it wraps (or the port routines, where the C has no single counterpart).
;;;;
;;;; The inputs come from the generators of the port's own tests (tests/nl/),
;;;; which keep them inside the C's domain.  Inputs go in now as Batten arrays,
;;;; now as plain lists, and every test checks that the API left its arguments
;;;; as they were (the C side answers T there).

(in-package :batten-tests)

;;; ------------------------------------------------------------------ helpers

(defun %apc-rows (a)
  "The 2D array A (a point-array) as a list of lists of its rows."
  (let (rows)
    (dotimes (i (array-dimension a 0) (nreverse rows))
      (let (row)
        (dotimes (k (array-dimension a 1))
          (push (aref a i k) row))
        (push (nreverse row) rows)))))

(defun %apc-net-rows (a)
  "The 3D array A (a point-net) as a list of rows of lists of points."
  (let (rows)
    (dotimes (i (array-dimension a 0) (nreverse rows))
      (let (row)
        (dotimes (j (array-dimension a 1))
          (push (list (aref a i j 0) (aref a i j 1) (aref a i j 2)) row))
        (push (nreverse row) rows)))))

(defun %apc-maybe-list (x)
  "X, a point, array or vector, or now and then the same as lists."
  (if (zerop (rnd 2))
      x
      (case (array-rank x)
        (1 (coerce x 'list))
        (2 (%apc-rows x))
        (t (%apc-net-rows x)))))

(defun %apc-copy (x)
  "A deep copy of X: an array, a curve, or a list of them."
  (cond ((consp x) (mapcar #'%apc-copy x))
        ((typep x 'curve) (batten::crv-copy x))
        ((arrayp x) (%crv-fit-copy-array x))
        (t x)))

(defun %apc-unchanged (before after)
  "T when AFTER is still bit for bit BEFORE."
  (values (same before after)))

(defun %apc-par (key)
  (ecase key (:chord-length +chordlength+) (:centripetal +centripetal+) (:uniform +uniform+)))

(defun %apc-random-par ()
  (elt '(:chord-length :centripetal :uniform) (rnd 3)))

(defmacro %apc-c-surface ((sur sg) call)
  "Run CALL, a C routine filling the NULL surface SUR from stack SG; (values
error-p (list surface))."
  `(with-stacks (,sg)
     (with-foreign-out ((,sur :surface))
       (let ((err (c-flag-error-p ,call)))
         (values err (unless err (list (from-foreign ,sur :surface))))))))

(defun %apc-c-knots (u n p)
  "The C's averaged knot vector (N_FitCrvCalcKnotVector) of the parameters U[0..N]."
  (with-foreign ((uf u :real-vector) (kf (make-real-vector (+ n p 2)) :knots))
    (c-fit-crv-calc-knot-vector uf n p kf)
    (from-foreign kf :knots)))

;;; ---------------------------------------------------------------- from parts

;;; No C counterpart: the result is the port's %make-surface of the control
;;; points NLib's convention makes, (w*x w*y w*z w).
(define-oracle-test api-make-surface (:file "api-construction" :trials 1000)
  (let* ((rational (zerop (rnd 2)))
         (s (random-surface :rational nil))
         (p (surface-degree-u s)) (q (surface-degree-v s))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (n+1 (array-dimension (surface-cpoints s) 0))
         (m+1 (array-dimension (surface-cpoints s) 1))
         (net (make-point-net n+1 m+1))
         (w (make-real-matrix n+1 m+1))
         (expected (make-cpoint-net n+1 m+1))
         (bad (zerop (rnd 15))))
    (dotimes (i n+1)
      (dotimes (j m+1)
        (let ((wij (+ 0.2d0 (rnd 4.8d0))))
          (setf (aref w i j) wij)
          (dotimes (k 3)
            (setf (aref net i j k) (aref (surface-cpoints s) i j k)
                  (aref expected i j k) (if rational (* wij (aref net i j k)) (aref net i j k))))
          (setf (aref expected i j 3) (if rational wij +no-w+)))))
    (let ((weights (when rational
                     (let (rows)
                       (dotimes (i n+1 (nreverse rows))
                         (let (row) (dotimes (j m+1) (push (aref w i j) row)) (push (nreverse row) rows))))))
          (ku-in (if bad (subseq ku 1) (%apc-maybe-list ku)))
          (net-in (%apc-maybe-list net)))
      (agree (list :p p :q q :ku ku :kv kv :net net :weights weights :bad bad)
             (let ((before (%apc-copy net-in)))
               (list (batten::make-surface p q ku-in kv net-in :weights weights)
                     (%apc-unchanged before net-in)))
             (if bad
                 (values t nil)
                 (values nil (list (%make-surface p q (copy-seq ku) (copy-seq kv) expected) t)))))))

;;; ------------------------------------------------------- lines, arcs, conics

(define-oracle-test api-line-segment (:file "api-construction" :trials 1000)
  (let* ((p (random-point))
         (q (case (rnd 10) (0 (copy-seq p)) (1 (random-point 1d-3)) (t (random-point))))
         (v (make-point (- (aref q 0) (aref p 0)) (- (aref q 1) (aref p 1)) (- (aref q 2) (aref p 2)))))
    (agree (list :p p :q q)
           (let ((pin (%apc-maybe-list p)) (qin (%apc-maybe-list q)))
             (list (batten::line-segment pin qin)))
           ;; NLib's line from P along Q - P, its end control point then made Q
           (with-stacks (s)
             (with-foreign ((fp p :point) (fv v :point))
               (with-foreign-out ((cur :curve))
                 (let ((rc (c-crv-line-from-pt-and-vector fp fv cur s)))
                   (values (c-flag-error-p rc)
                           (unless (c-flag-error-p rc)
                             (let ((c (from-foreign cur :curve)))
                               (dotimes (k 3) (setf (aref (curve-cpoints c) 1 k) (aref q k)))
                               (list c)))))))))))

(defun %apc-circle-degree ()
  (elt '(2 4 5) (rnd 3)))

(define-oracle-test api-circular-arc (:file "api-construction" :trials 1000)
  (destructuring-bind (c x y r as ae at) (%cc-arc-input)
    (let ((degree (%apc-circle-degree))
          (xin (%apc-maybe-list x)) (yin (%apc-maybe-list y)))
      (agree (list :c c :x x :y y :r r :as as :ae ae :at at :degree degree)
             (let ((before (%apc-copy (list xin yin))))
               (list (batten::circular-arc c r :x-axis xin :y-axis yin :start-angle as :end-angle ae
                                               :degree degree :angle-tolerance at)
                     (%apc-unchanged before (list xin yin))))
             (with-stacks (s)
               (with-foreign ((fc c :point) (fx x :point) (fy y :point))
                 (with-foreign-out ((cur :curve))
                   (let ((rc (c-create-circ-arc-tol fc fx fy r as ae at degree cur s)))
                     (values (c-flag-error-p rc)
                             (unless (c-flag-error-p rc)
                               (list (from-foreign cur :curve) t)))))))))))

(define-oracle-test api-elliptical-arc (:file "api-construction" :trials 1000)
  (destructuring-bind (c x y r1 as ae at) (%cc-arc-input)
    ;; no zero radius (see the create-elliptical-arc-ang-type-tol test)
    (let* ((r1 (if (and (zerop (rnd 10)) (/= r1 0d0)) r1 (rnd-range 0.5d0 80d0)))
           (r2 (case (rnd 10)
                 (0 r1)
                 (1 (* r1 (rnd-range 1d-4 1d-2)))
                 (2 (* r1 (rnd-range 1d0 3d0)))
                 (3 (- -0.1d0 (rnd 10d0)))
                 (t (* r1 (rnd-range 0.05d0 1d0)))))
           (on (if (zerop (rnd 2)) :circle :ellipse))
           (degree (%apc-circle-degree))
           (xin (%apc-maybe-list x)) (yin (%apc-maybe-list y)))
      (agree (list :c c :x x :y y :r1 r1 :r2 r2 :as as :ae ae :at at :on on :degree degree)
             (let ((before (%apc-copy (list xin yin))))
               (list (batten::elliptical-arc c r1 r2 :x-axis xin :y-axis yin :start-angle as
                                                     :end-angle ae :angles-on on :degree degree
                                                     :angle-tolerance at)
                     (%apc-unchanged before (list xin yin))))
             (with-stacks (s)
               (with-foreign ((fc c :point) (fx x :point) (fy y :point))
                 (with-foreign-out ((cur :curve))
                   (let ((rc (c-create-elliptical-arc-ang-type-tol fc fx fy r1 r2 as ae at
                                                                   (if (eq on :circle) 0 1) degree cur s)))
                     (values (c-flag-error-p rc)
                             (unless (c-flag-error-p rc)
                               (list (from-foreign cur :curve) t)))))))))))

(define-oracle-test api-conic-arc (:file "api-construction" :trials 1000)
  (destructuring-bind (p0 t0 p2 t2 p) (%cc-conic-input)
    (let ((args (mapcar #'%apc-maybe-list (list p0 t0 p2 t2 p))))
      (agree (list :p0 p0 :t0 t0 :p2 p2 :t2 t2 :p p)
             (let ((before (%apc-copy args)))
               (list (apply #'batten::conic-arc args) (%apc-unchanged before args)))
             (with-stacks (s)
               (with-foreign ((a p0 :point) (b t0 :point) (c p2 :point) (d t2 :point) (e p :point))
                 (with-foreign-out ((cur :curve))
                   (let ((rc (c-create-conic-arc a b c d e cur s)))
                     (values (c-flag-error-p rc)
                             (unless (c-flag-error-p rc)
                               (list (from-foreign cur :curve) t)))))))))))

;;; ------------------------------------------------------------- curve fitting

;;; Plain interpolation is N_FitCrvInterp; with end conditions or parameters,
;;; N_FitCrvKnotsAndTangents at the parameters N_FitCalcCrvParamValues gives.
(define-oracle-test api-interpolate-curve (:file "api-construction" :trials 1000)
  (let* ((k (+ 2 (rnd 12)))
         (mode (rnd 4))                 ; 0 plain, 1 parameters, 2 tangents, 3 derivatives
         (tan (if (>= mode 2) (1+ (rnd 3)) 0))
         (n (case tan (0 k) ((1 2) (1+ k)) (t (+ k 2))))
         (p (if (= mode 0)
                (if (zerop (rnd 20)) (+ k 1) (1+ (rnd (min k 8))))
                (if (zerop (rnd 25)) 1 (+ 2 (rnd (max 1 (min 5 (1- n))))))))
         (par (%apc-random-par))
         (pts (%cf2-planar (random-points (1+ k))))
         (params (when (or (= mode 1) (and (>= mode 2) (zerop (rnd 3)))) (%cf2-params (1+ k))))
         (ts (when (member tan '(1 3)) (random-point 50d0)))
         (te (when (member tan '(2 3)) (random-point 50d0)))
         (pin (%apc-maybe-list pts)))
    (agree (list :k k :mode mode :p p :par par :params params :ts ts :te te :pts pts)
           (let ((before (%apc-copy (list pin ts te))))
             (list (apply #'batten::interpolate-curve pin :degree p :parameterization par
                          (append (when params (list :parameters params))
                                  (case mode
                                    (2 (list :start-tangent ts :end-tangent te))
                                    (3 (list :start-derivative ts :end-derivative te)))))
                   (%apc-unchanged before (list pin ts te))))
           (with-foreign ((pp pts :point-array)
                          (uf (or params (make-real-vector (1+ k))) :real-vector)
                          (tsf (or ts (make-point)) :point)
                          (tef (or te (make-point)) :point))
             (block c
               (when (= mode 0)
                 (return-from c
                   (with-stacks (sg)
                     (with-foreign-out ((cur :curve))
                       (let ((err (c-flag-error-p (c-fit-crv-interp pp k p (%apc-par par) cur sg))))
                         (values err (unless err (list (from-foreign cur :curve) t))))))))
               (unless params
                 (when (c-flag-error-p (c-fit-calc-crv-param-values pp k +epoint+ (%apc-par par) uf))
                   (return-from c (values t nil))))
               (with-stacks (sc)
                 (with-stacks (sk)
                   (with-scalars ((kc :pointer))
                     (with-foreign-out ((cur :curve))
                       (let ((err (c-flag-error-p
                                   (c-fit-crv-knots-and-tangents
                                    pp k uf (if ts tsf (cffi:null-pointer)) (if te tef (cffi:null-pointer))
                                    (if (= mode 3) +derivative+ +tangent+) kc 0d0 p cur sc sk))))
                         (values err (unless err (list (from-foreign cur :curve) t)))))))))))))

;;; To a tolerance, N_FitCrvApprox; with a control point count, the
;;; least-squares fit: N_FitCalcCrvParamValues, N_FitCalcKnotVectorCrvApprox,
;;; N_FitCrvApproxKnots.
(define-oracle-test api-approximate-curve (:file "api-construction" :trials 600)
  (let* ((m (+ 3 (rnd 25)))
         (lsq (zerop (rnd 2)))
         (closed (zerop (rnd 3)))
         (pts (case (rnd 4)
                (0 (let ((p (%cf2-planar (random-points (1+ m)))))
                     (when closed (dotimes (j 3) (setf (aref p m j) (aref p 0 j))))
                     p))
                (t (%cf3-smooth-points (1+ m) :closed closed))))
         (pin (%apc-maybe-list pts)))
    (if lsq
        (let* ((p (1+ (rnd (min m 5))))
               (n (if (zerop (rnd 15)) (1+ m) (+ p (rnd (1+ (- m p))))))
               (par (%apc-random-par)))
          (agree (list :m m :n n :p p :par par :pts pts)
                 (let ((before (%apc-copy pin)))
                   (list (batten::approximate-curve pin :degree p :control-points (1+ n) :parameterization par)
                         (%apc-unchanged before pin)))
                 (with-foreign ((pp pts :point-array)
                                (uf (make-real-vector (1+ m)) :real-vector)
                                (kf (make-real-vector (+ n p 2)) :knots))
                   (cond ((c-flag-error-p (c-fit-calc-crv-param-values pp m +epoint+ (%apc-par par) uf))
                          (values t nil))
                         ((c-flag-error-p (c-fit-calc-knot-vector-crv-approx uf m n p kf))
                          (values t nil))
                         (t (with-stacks (sg)
                              (with-foreign-out ((cur :curve))
                                (let ((err (c-flag-error-p (c-fit-crv-approx-knots pp m +epoint+ uf kf n p cur sg))))
                                  (values err (unless err (list (from-foreign cur :curve) t)))))))))))
        (let* ((close-p (plusp (rnd 4)))
               ;; closing closed data needs a start degree of 2 or more
               (ps (if (and closed close-p) (+ 2 (rnd 2)) (1+ (rnd 3))))
               (pr (case (rnd 15) (0 (max 1 (1- ps))) (t (min m (+ ps (rnd 3))))))
               (ext (batten::pts-1d-get-max-extent pts m))
               (e (* ext (case (rnd 6) (0 1d-5) (1 1d-3) (2 1d-2) (3 0.05d0) (4 0.3d0) (t 1d30))))
               (knots (if (zerop (rnd 2)) :single :multiple)))
          (agree (list :m m :close close-p :ps ps :pr pr :e e :knots knots :pts pts)
                 (let ((before (%apc-copy pin)))
                   (list (batten::approximate-curve pin :degree pr :tolerance e :start-degree ps
                                                        :closed close-p :knots knots)
                         (%apc-unchanged before pin)))
                 (with-foreign ((pp pts :point-array))
                   (with-stacks (sg)
                     (with-foreign-out ((cur :curve))
                       (let ((err (c-flag-error-p
                                   (c-fit-crv-approx pp m (if close-p batten::+yes+ batten::+no+) ps pr e
                                                     (if (eq knots :single) +single+ +multiple+) cur sg))))
                         (values err (unless err (list (from-foreign cur :curve) t))))))))))))

(define-oracle-test api-fit-arc (:file "api-construction" :trials 1000)
  (let* ((n (if (zerop (rnd 30)) 0 (1+ (rnd 15))))
         (pts (%cf6-arc-points (1+ n)))
         (tol (case (rnd 4) (0 1d-3) (1 1d-6) (2 0.1d0) (t (expt 10d0 (- (rnd-range 0d0 6d0))))))
         (measure (if (zerop (rnd 2)) :maximum :average))
         (arc (zerop (rnd 3)))
         (degree (%apc-circle-degree))
         (pin (%apc-maybe-list pts)))
    (agree (list :n n :tol tol :measure measure :arc arc :degree degree :pts pts)
           (let ((before (%apc-copy pin)))
             (multiple-value-bind (cur kind era erm)
                 (batten::fit-arc pin :tolerance tol :error-measure measure :arc arc :degree degree)
               (list cur kind era erm (%apc-unchanged before pin))))
           (with-foreign ((pp pts :point-array))
             (with-stacks (sg)
               (with-scalars ((cut :short) (era :double) (erm :double))
                 (with-foreign-out ((cur :curve))
                   (let ((err (c-flag-error-p
                               (c-fit-arc-to-pts pp n tol (if (eq measure :maximum) +maximum+ +average+)
                                                 (if arc batten::+yes+ batten::+no+) degree cur cut era erm sg))))
                     (values err
                             (unless err
                               (list (from-foreign cur :curve)
                                     (if (= (cffi:mem-ref cut :short) +nline+) :line :circle)
                                     (cffi:mem-ref era :double) (cffi:mem-ref erm :double) t)))))))))))

(define-oracle-test api-fit-line (:file "api-construction" :trials 1000)
  (let* ((count (1+ (rnd 20)))
         (pts (random-points count (case (rnd 3) (0 1d0) (1 100d0) (t 1d4)))))
    (case (rnd 5)
      (0                                ; collinear
       (let ((a (random-point)) (d (random-point 1d0)))
         (dotimes (i count)
           (let ((s (rnd-range -10d0 10d0)))
             (dotimes (k 3) (setf (aref pts i k) (+ (aref a k) (* s (aref d k)))))))))
      (1                                ; nearly collinear
       (let ((a (random-point)) (d (random-point 1d0)))
         (dotimes (i count)
           (let ((s (rnd-range -10d0 10d0)))
             (dotimes (k 3) (setf (aref pts i k) (+ (aref a k) (* s (aref d k)) (rnd-range -1d-3 1d-3))))))))
      (2                               ; all the same point
       (dotimes (i count) (dotimes (k 3) (setf (aref pts i k) (aref pts 0 k))))))
    (let ((pin (%apc-maybe-list pts)))
      (agree (list :pts pts)
             (let ((before (%apc-copy pin)))
               (multiple-value-bind (c v era erm) (batten::fit-line pin)
                 (list c v era erm (%apc-unchanged before pin))))
             (with-foreign ((pp pts :point-array))
               (with-foreign-out ((c :point) (v :point))
                 (with-scalars ((era :double) (erm :double))
                   (let ((err (c-flag-error-p (c-line-fit-pts pp (1- count) c v era erm))))
                     (values err
                             (unless err
                               (list (from-foreign c :point) (from-foreign v :point)
                                     (cffi:mem-ref era :double) (cffi:mem-ref erm :double) t)))))))))))

;;; ----------------------------------------------------------- surface fitting

;;; Plain interpolation is N_FitSrfToPts; with parameters or knots given,
;;; N_FitSrfToPtsKnots on the parameters N_FitCalcSrfParamValues gives (where
;;; not given) and the knots N_FitCrvCalcKnotVector averages (where not given).
(define-oracle-test api-interpolate-surface (:file "api-construction" :trials 1000)
  (let* ((n (1+ (rnd 8))) (m (1+ (rnd 8)))
         (p (if (zerop (rnd 15)) (+ n 1) (1+ (rnd (min n 4)))))
         (q (if (zerop (rnd 15)) (+ m 1) (1+ (rnd (min m 4)))))
         (par (%apc-random-par))
         (pts (%random-grid (1+ n) (1+ m)))
         (mode (rnd 5))       ; 0 plain, 1 both parameters, 2 u parameters, 3 u knots, 4 all
         (pu (when (member mode '(1 2 4)) (%increasing-params n)))
         (pv (when (member mode '(1 4)) (%increasing-params m)))
         (ku (when (and (member mode '(3 4)) (<= p n))
               (if pu
                   (batten::fit-crv-calc-knot-vector pu n p (make-real-vector (+ n p 2)))
                   (random-knots p (1+ n)))))
         (kv (when (and (= mode 4) (<= q m))
               (batten::fit-crv-calc-knot-vector pv m q (make-real-vector (+ m q 2)))))
         (net-in (%apc-maybe-list pts)))
    (agree (list :n n :m m :p p :q q :par par :pu pu :pv pv :ku ku :kv kv :pts pts)
           (let ((before (%apc-copy net-in)))
             (list (batten::interpolate-surface net-in :degree-u p :degree-v q :parameterization par
                                                       :parameters-u pu :parameters-v pv
                                                       :knots-u ku :knots-v kv)
                   (%apc-unchanged before net-in)))
           (with-foreign ((ap pts :point-net)
                          (up (make-real-vector (1+ n)) :real-vector)
                          (vp (make-real-vector (1+ m)) :real-vector))
             (block c
               (unless (or pu pv ku kv)
                 (return-from c
                   (multiple-value-bind (err r) (%apc-c-surface (sur sg) (c-fit-srf-to-pts ap n m p q (%apc-par par) sur sg))
                     (values err (when r (list (first r) t))))))
               (unless (and pu pv)
                 (when (c-flag-error-p (c-fit-calc-srf-param-values ap n m +epoint+ (%apc-par par) up vp))
                   (return-from c (values t nil))))
               (let* ((u (or pu (from-foreign up :real-vector (1+ n))))
                      (v (or pv (from-foreign vp :real-vector (1+ m))))
                      (knu (or ku (%apc-c-knots u n p)))
                      (knv (or kv (%apc-c-knots v m q))))
                 (with-foreign ((uf u :real-vector) (vf v :real-vector) (kuf knu :knots) (kvf knv :knots))
                   (multiple-value-bind (err r)
                       (%apc-c-surface (sur sg) (c-fit-srf-to-pts-knots ap n m +epoint+ uf vf kuf kvf p q sur sg))
                     (values err (when r (list (first r) t)))))))))))

;;; To a tolerance, N_FitSrfApproxTol; with control point counts,
;;; N_FitSrfLstSqApprox.
(define-oracle-test api-approximate-surface (:file "api-construction" :trials 600)
  (if (zerop (rnd 2))
      (multiple-value-bind (r s n m p q) (%random-fit-sizes)
        (let* ((pts (%random-grid (1+ r) (1+ s)))
               (par (%apc-random-par))
               (which (rnd 4))          ; 0 u count only, 1 v count only, else both
               (net-in (%apc-maybe-list pts)))
          (when (= which 0) (setf m s))
          (when (= which 1) (setf n r))
          (when (and (= which 0) (> q m)) (setf q m))
          (when (and (= which 1) (> p n)) (setf p n))
          (agree (list :r r :s s :n n :m m :p p :q q :par par :which which :pts pts)
                 (let ((before (%apc-copy net-in)))
                   (list (batten::approximate-surface net-in :degree-u p :degree-v q :parameterization par
                                                             :control-points-u (unless (= which 1) (1+ n))
                                                             :control-points-v (unless (= which 0) (1+ m)))
                         (%apc-unchanged before net-in)))
                 (with-foreign ((pp pts :point-net))
                   (multiple-value-bind (err res)
                       (%apc-c-surface (sur sg) (c-fit-srf-lst-sq-approx pp r s n m p q (%apc-par par) sur sg))
                     (values err (when res (list (first res) t))))))))
      (let* ((mu (+ 1 (rnd 10))) (mv (+ 1 (rnd 10)))
             (ps (1+ (rnd (min mu 2)))) (qs (1+ (rnd (min mv 2))))
             (pr (min mu (+ ps (rnd 3)))) (qr (min mv (+ qs (rnd 3))))
             (pts (%random-grid (1+ mu) (1+ mv)))
             (e (case (rnd 6) (0 1d-6) (1 1d-2) (2 0.5d0) (3 (rnd 5d0)) (4 50d0) (t 1d4)))
             (knots (if (zerop (rnd 2)) :single :multiple))
             (from-removal (zerop (rnd 2)))
             (net-in (%apc-maybe-list pts)))
        (when (zerop (rnd 30)) (setf pr (1+ mu)))
        (agree (list :mu mu :mv mv :ps ps :pr pr :qs qs :qr qr :e e :knots knots :from-removal from-removal :pts pts)
               (let ((before (%apc-copy net-in)))
                 (list (batten::approximate-surface net-in :degree-u pr :degree-v qr :tolerance e
                                                           :start-degree-u ps :start-degree-v qs
                                                           :knots knots :knots-from-removal from-removal)
                       (%apc-unchanged before net-in)))
               (with-foreign ((pp pts :point-net))
                 (multiple-value-bind (err res)
                     (%apc-c-surface (sur sg) (c-fit-srf-approx-tol pp mu mv ps pr qs qr e
                                                                    (if (eq knots :single) +single+ +multiple+)
                                                                    (if from-removal batten::+yes+ batten::+no+)
                                                                    sur sg))
                   (values err (when res (list (first res) t)))))))))

;;; ---------------------------------------------------- classic constructions

(defun %apc-across ()
  (if (zerop (rnd 2)) :u :v))

(defun %apc-dir (across)
  (if (eq across :u) +udir+ +vdir+))

(define-oracle-test api-ruled-surface (:file "api-construction" :trials 1000)
  (let* ((two-d (zerop (rnd 5)))
         (rat (rnd 4))
         (cp (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (oddp rat)))
         (cq (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (>= rat 2)))
         (across (%apc-across)))
    (agree (list :p cp :q cq :across across)
           (let ((before (%apc-copy (list cp cq))))
             (list (batten::ruled-surface cp cq :across across) (%apc-unchanged before (list cp cq))))
           (with-foreign ((fp cp :curve) (fq cq :curve))
             (multiple-value-bind (err res)
                 (%apc-c-surface (sur sg) (c-create-ruled-srf-from-boundary-crvs fp fq (%apc-dir across) sur sg))
               (values err (when res (list (first res) t))))))))

(define-oracle-test api-extruded-surface (:file "api-construction" :trials 1000)
  (let* ((c (random-curve))
         (w (case (rnd 20)
              (0 (make-point 0d0 0d0 0d0))
              (1 (make-point 0d0 0d0 (rnd-range 0.1d0 3d0)))
              (t (random-point 5d0))))
         (d (case (rnd 10) (0 0d0) (1 (- (rnd 10d0))) (t (rnd 50d0))))
         (across (%apc-across))
         (win (%apc-maybe-list w)))
    (agree (list :curve c :w w :d d :across across)
           (let ((before (%apc-copy (list c win))))
             (list (batten::extruded-surface c win d :across across) (%apc-unchanged before (list c win))))
           (with-foreign ((cur c :curve) (wp w :point))
             (multiple-value-bind (err res)
                 (%apc-c-surface (sur sg) (c-create-srf-extrude-crv cur wp d (%apc-dir across) sur sg))
               (values err (when res (list (first res) t))))))))

(define-oracle-test api-revolved-surface (:file "api-construction" :trials 1000)
  (multiple-value-bind (c s tv al ctp) (%sc-revolve-inputs)
    ;; the API takes the circle kinds by degree only
    (let ((degree (if (member ctp '(2 4 5)) ctp 2))
          (s-in (%apc-maybe-list s)) (t-in (%apc-maybe-list tv)))
      (agree (list :curve c :s s :t tv :al al :degree degree)
             (let ((before (%apc-copy (list c s-in t-in))))
               (list (batten::revolved-surface c s-in t-in :angle al :degree degree)
                     (%apc-unchanged before (list c s-in t-in))))
             (with-foreign ((cur c :curve) (sp s :point) (tp tv :point))
               (multiple-value-bind (err res)
                   (%apc-c-surface (sur sg) (c-create-revolved-srf cur sp tp al degree sur sg))
                 (values err (when res (list (first res) t)))))))))

;;; Without parameters N_CreateSkinSrf; with them N_CreateSkinSrfParams, on
;;; the given knots or those N_FitCrvCalcKnotVector averages.
(define-oracle-test api-skinned-surface (:file "api-construction" :trials 600)
  (multiple-value-bind (curves k cmp cld deg dir closed) (%srf-adv-skin-inputs)
    (let* ((across (cond ((= dir +udir+) :u) ((= dir +vdir+) :v) (t (%apc-across))))
           (dir (%apc-dir across))
           (cs (coerce curves 'list))
           (with-params (zerop (rnd 2)))
           (default-degree (and (= deg (min 3 k)) (zerop (rnd 2)))))
      (if with-params
          (multiple-value-bind (lo hi) (if (zerop (rnd 3)) (values 0d0 1d0) (batten-random::random-domain))
            (let* ((tt (%srf-adv-ascending lo hi k))
                   (knt (when (or closed (zerop (rnd 4)))
                          (batten::with-c-float-environment
                            (if closed
                                (batten::fit-calc-knot-vector-end-derivs tt k deg (make-real-vector (+ k deg 4)))
                                (random-knots deg (1+ k) :lo lo :hi hi))))))
              (agree (list :cur curves :k k :cmp cmp :cld cld :deg deg :across across :t tt :knt knt)
                     (let ((before (%apc-copy cs)))
                       (list (apply #'batten::skinned-surface cs :across across
                                    :make-compatible (= cmp batten::+yes+) :closed (= cld batten::+yes+)
                                    :parameters tt :knots knt
                                    (unless default-degree (list :degree deg)))
                             (%apc-unchanged before cs)))
                     (let ((knt (or knt (%apc-c-knots tt k deg))))
                       (with-foreign ((ft tt :real-vector) (fk knt :knots))
                         (with-stacks (sc)
                           (with-stacks (ss)
                             (let ((fc (%srf-adv-curves-to-foreign curves)))
                               (with-foreign-out ((sur :surface))
                                 (let ((failed (c-flag-error-p
                                                (c-create-skin-srf-params fc k cmp cld deg dir ft fk sur sc ss))))
                                   (values failed (unless failed (list (from-foreign sur :surface) t)))))))))))))
          (let ((knot-tol (case (rnd 6) (0 1d-3) (1 1d-6) (2 0.05d0) (t 0d0))))
            (agree (list :cur curves :k k :cmp cmp :cld cld :deg deg :across across :knot-tol knot-tol)
                   (let ((before (%apc-copy cs)))
                     (list (apply #'batten::skinned-surface cs :across across
                                  :make-compatible (= cmp batten::+yes+) :closed (= cld batten::+yes+)
                                  :knot-tolerance knot-tol
                                  (unless default-degree (list :degree deg)))
                           (%apc-unchanged before cs)))
                   (with-foreign ((dummy (make-real-vector 1) :real-vector))
                     (declare (ignore dummy))
                     (with-stacks (sc)
                       (with-stacks (ss)
                         (let ((fc (%srf-adv-curves-to-foreign curves)))
                           (with-foreign-out ((sur :surface))
                             (let ((failed (c-flag-error-p
                                            (c-create-skin-srf fc k cmp cld deg dir knot-tol sur sc ss))))
                               (values failed (unless failed (list (from-foreign sur :surface) t)))))))))))))))

(define-oracle-test api-gordon-surface (:file "api-construction" :trials 600)
  (let* ((k (1+ (rnd 5))) (l (1+ (rnd 5)))
         (defaults (zerop (rnd 3)))
         (deg (lambda (top) (if (zerop (rnd 40)) (1+ top) (1+ (rnd (min top 3))))))
         (dlu (if defaults (min 3 l) (funcall deg l)))
         (dlv (if defaults (min 3 k) (funcall deg k)))
         (dtu (if defaults dlu (funcall deg l)))
         (dtv (if defaults dlv (funcall deg k)))
         (preserve (zerop (rnd 2))))
    (multiple-value-bind (cur-u v-prms cur-v u-prms) (%srf-adv-gordon-network k l)
      (let ((cu (coerce cur-u 'list)) (cv (coerce cur-v 'list))
            (vs (%apc-maybe-list v-prms)) (us (%apc-maybe-list u-prms)))
        (agree (list :cur-u cur-u :v-prms v-prms :cur-v cur-v :u-prms u-prms
                     :degs (list dlu dlv dtu dtv) :defaults defaults :preserve preserve)
               (let ((before (%apc-copy (list cu vs cv us))))
                 (list (apply #'batten::gordon-surface cu vs cv us :preserve-parameters preserve
                              (unless defaults
                                (list :degree-u dlu :degree-v dlv :tensor-degree-u dtu :tensor-degree-v dtv)))
                       (%apc-unchanged before (list cu vs cv us))))
               (with-foreign ((fv v-prms :real-vector) (fu u-prms :real-vector))
                 (with-stacks (sc)
                   (with-stacks (ss)
                     (let ((au (%srf-adv-curves-to-foreign cur-u))
                           (av (%srf-adv-curves-to-foreign cur-v)))
                       (with-foreign-out ((sur :surface))
                         (let ((failed (c-flag-error-p
                                        (c-create-gordon-srf au fv k av fu l dlu dlv dtu dtv
                                                             (if preserve batten::+yes+ batten::+no+)
                                                             sur sc ss))))
                           (values failed (unless failed (list (from-foreign sur :surface) t))))))))))))))

(define-oracle-test api-coons-surface (:file "api-construction" :trials 600)
  ;; the generator's own edits (an elevation after a knot was moved) can fail:
  ;; draw again
  (destructuring-bind (cl dl cr dr cb db ct dt)
      (do ((data nil)) (data data)
        (setf data (handler-case (batten::with-c-float-environment (%srf-adv-coons-data))
                     (nurbs-error () nil))))
    (let ((check (zerop (rnd 2)))
          (boundary (list cl cr cb ct))
          (derivatives (list dl dr db dt)))
      (agree (list :boundary boundary :derivatives derivatives :check check)
             (let ((before (%apc-copy (list boundary derivatives))))
               (list (batten::coons-surface boundary derivatives :check-corners check)
                     (%apc-unchanged before (list boundary derivatives))))
             (with-foreign ((dummy (make-real-vector 1) :real-vector))
               (declare (ignore dummy))
               (let ((fc (mapcar (lambda (c) (to-foreign c :curve)) (list cl dl cr dr cb db ct dt))))
                 (with-stacks (sc)
                   (with-stacks (ss)
                     (with-foreign-out ((sur :surface))
                       (let ((failed (c-flag-error-p
                                      (apply #'c-create-coons-boundary-crvs
                                             (append fc (list (if check batten::+yes+ batten::+no+) sur sc ss))))))
                         (values failed (unless failed (list (from-foreign sur :surface) t)))))))))))))

(define-oracle-test api-sphere (:file "api-construction" :trials 1000)
  ;; as the create-sphere test: angles mostly valid, now and then not; quartic
  ;; revolutions only below 241 degrees (the C's weight array holds 9)
  (let* ((c (if (zerop (rnd 5)) (make-point 0d0 0d0 0d0) (random-point 50d0)))
         (r (case (rnd 15) (0 0d0) (1 (- (rnd 10d0))) (t (rnd-range 0.01d0 100d0))))
         (as (case (rnd 6) ((0 1) 0d0) (2 (float (* 30 (rnd 6)) 1d0)) (t (rnd 180d0))))
         (ae (case (rnd 6) ((0 1) 180d0) (2 (min 180d0 (+ as (rnd 1d-3)))) (t (rnd-range as 180d0))))
         (as (case (rnd 50) (0 ae) (1 (- (rnd 5d0))) (t as)))
         (ae (case (rnd 50) (0 (+ 180d0 (rnd 5d0))) (1 (- as (rnd 5d0))) (t ae)))
         (al (case (rnd 10)
               (0 90d0) (1 180d0) (2 270d0) (3 360d0)
               (4 (rnd-range 0.1d0 5d0))
               (5 (rnd-range 360d0 500d0))
               (6 (- (rnd 100d0)))
               (t (rnd-range 0.1d0 360d0))))
         (dv (%apc-circle-degree))
         (du (case (rnd 5)
               (0 (if (and (> al 0d0) (< al 240d0)) 4 2))
               (t (if (zerop (rnd 2)) 2 5))))
         (defaults (and (= dv 2) (= du 2) (= as 0d0) (= ae 180d0) (= al 360d0)))
         (cin (%apc-maybe-list c)))
    (agree (list :c c :r r :as as :ae ae :al al :du du :dv dv)
           (let ((before (%apc-copy cin)))
             (list (if defaults
                       (batten::sphere cin r)
                       (batten::sphere cin r :start-angle as :end-angle ae :sweep-angle al
                                             :degree-u du :degree-v dv))
                   (%apc-unchanged before cin)))
           (with-foreign ((fc c :point))
             (multiple-value-bind (err res)
                 (%apc-c-surface (sur sg) (c-create-sphere fc r as ae al dv du sur sg))
               (values err (when res (list (first res) t))))))))
