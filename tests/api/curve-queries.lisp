;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Tests of the API in src/api/curve-queries.lisp: each function against the C routines
;;;; it wraps (or the port routines, where the C has no single counterpart).
;;;;
;;;; Every Lisp result list ends with whether the API left its curve argument
;;;; as it was (T), against a T on the C side: the API is pure, the port's
;;;; routines need not be.  The curve generators are the port tests' own
;;;; (tests/nl/crv-basic.lisp, crv-nurb.lisp, crv-geom.lisp).

(in-package :batten-tests)

(defun %api-untouched (c snapshot)
  "T when curve C is still identical to SNAPSHOT, a copy taken before the call."
  (if (same c snapshot) t nil))

(defun %api-side () (if (zerop (rnd 2)) :left :right))

(defun %api-flag (side) (if (eq side :left) +left+ +right+))

(defun %api-real (x)
  "X, or now and then the same value as a rational (the API takes any real)."
  (if (zerop (rnd 4)) (rational x) x))

(defun %api-point (p)
  "Point P, or now and then the same point as a list or a simple-vector (the
API takes any sequence of three reals)."
  (case (rnd 4)
    (0 (coerce p 'list))
    (1 (map 'simple-vector #'identity p))
    (t p)))

(defun %api-ordered-parameters (knots)
  "Two parameters in KNOTS' domain, (values a b), mostly with a <= b."
  (let ((a (random-parameter knots)) (b (random-parameter knots)))
    (if (and (> a b) (/= 0 (rnd 4)))
        (values b a)
        (values a b))))

;;; ------------------------------------------------------- the shared generics

(define-oracle-test api-domain (:file "api-curve-queries" :trials 500)
  (let* ((c (random-curve)) (snap (batten::crv-copy c)))
    (agree c
           (append (multiple-value-list (domain c)) (list (%api-untouched c snap)))
           (with-foreign ((cp c :curve))
             (with-scalars ((ul :double) (ur :double))
               (c-crv-get-param-bounds cp ul ur)
               (values nil (list (cffi:mem-ref ul :double) (cffi:mem-ref ur :double) t)))))))

(define-oracle-test api-bounding-box (:file "api-curve-queries" :trials 1000)
  (let* ((c (random-curve)) (snap (batten::crv-copy c)))
    (agree c
           (append (multiple-value-list (bounding-box c)) (list (%api-untouched c snap)))
           (with-foreign ((cp c :curve))
             (with-foreign-out ((box :box))
               (if (c-flag-error-p (c-crv-get-b-box cp box))
                   (values t nil)
                   (let ((b (from-foreign box :box)))
                     (values nil (list (make-point (aref b 0) (aref b 2) (aref b 4))
                                       (make-point (aref b 1) (aref b 3) (aref b 5))
                                       t)))))))))

;;; The copy also shares no array with the original.
(define-oracle-test api-copy (:file "api-curve-queries" :trials 500)
  (let* ((c (random-curve)) (snap (batten::crv-copy c)))
    (agree c
           (let ((d (copy c)))
             (list d
                   (and (not (eq (curve-knots d) (curve-knots c)))
                        (not (eq (curve-cpoints d) (curve-cpoints c))))
                   (%api-untouched c snap)))
           (with-stacks (s)
             (with-foreign ((cp c :curve))
               (with-foreign-out ((cq :curve))
                 (if (c-flag-error-p (c-crv-copy cp cq s))
                     (values t nil)
                     (values nil (list (from-foreign cq :curve) t t)))))))))

;;; ------------------------------------------------- tangent, curvature, frame

(define-oracle-test api-curve-tangent (:file "api-curve-queries" :trials 1000)
  (let* ((c (if (zerop (rnd 6)) (%curve-with-flat-start) (random-curve)))
         (snap (batten::crv-copy c))
         (u (if (zerop (rnd 6)) (aref (curve-knots c) 0) (random-parameter (curve-knots c))))
         (u-arg (%api-real u))
         (side (%api-side)))
    (agree (list :curve c :u u :side side)
           (append (multiple-value-list (batten::curve-tangent c u-arg :side side))
                   (list (%api-untouched c snap)))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((p :point) (tv :point))
               (if (c-flag-error-p (c-crv-eval-tangent cur u (%api-flag side) p tv))
                   (values t nil)
                   (values nil (list (from-foreign tv :point) (from-foreign p :point) t))))))))

(define-oracle-test api-curve-curvature (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 8)
              (0 (%curve-with-flat-start))
              (1 (random-curve :max-degree 1))
              (2 (random-curve :max-degree 2))
              (3 (random-curve :two-d t))
              (t (random-curve))))
         (snap (batten::crv-copy c))
         (k (curve-knots c))
         (u (case (rnd 8)
              (0 (aref k 0))
              (1 (aref k (1- (length k))))
              (t (random-parameter k))))
         (u-arg (%api-real u))
         (side (%api-side)))
    (agree (list :curve c :u u :side side)
           (list (batten::curve-curvature c u-arg :side side) (%api-untouched c snap))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((p :point) (tv :point) (cc :point) (nv :point))
               (with-scalars ((kap :double))
                 (if (c-flag-error-p (c-crv-eval-curvature cur u (%api-flag side) kap p tv +no+ cc nv))
                     (values t nil)
                     (values nil (list (cffi:mem-ref kap :double) t)))))))))

(define-oracle-test api-curve-frenet-frame (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 8)
              (0 (%curve-with-flat-start))
              (1 (random-curve :max-degree 1))
              (2 (random-curve :max-degree 2))
              (t (random-curve))))
         (snap (batten::crv-copy c))
         (u (if (zerop (rnd 6)) (aref (curve-knots c) 0) (random-parameter (curve-knots c))))
         (u-arg (%api-real u))
         (side (%api-side)))
    (agree (list :curve c :u u :side side)
           (append (multiple-value-list (batten::curve-frenet-frame c u-arg :side side))
                   (list (%api-untouched c snap)))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((p :point) (tv :point) (nv :point) (bv :point))
               (if (c-flag-error-p (c-crv-eval-frenet-frame cur u (%api-flag side) p tv nv bv))
                   (values t nil)
                   (values nil (list (from-foreign p :point) (from-foreign tv :point)
                                     (from-foreign nv :point) (from-foreign bv :point) t))))))))

;;; ------------------------------------------------------------------- length

;;; A tolerance of 0 is refused by the API (the C would subdivide for ever):
;;; there the C side is not called and counts as failing.
(define-oracle-test api-curve-length (:file "api-curve-queries" :trials 1000)
  (let* ((c (%crv-geom-arc-curve))
         (snap (batten::crv-copy c))
         (kn (curve-knots c))
         (lo (aref kn 0)) (hi (aref kn (1- (length kn)))))
    (multiple-value-bind (a b) (%api-ordered-parameters kn)
      (let* ((from (if (zerop (rnd 4)) nil a))
             (to (if (zerop (rnd 4)) nil b))
             (relative (zerop (rnd 2)))
             (tol (cond ((zerop (rnd 25)) 0d0)
                        (relative (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 2d0 7d0)))) (t 1d-5)))
                        (t (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 0d0 6d0)))) (t 0.1d0)))))
             (default-tol (zerop (rnd 8)))
             (tol (if default-tol 1d-6 tol))
             (relative (if default-tol t relative))
             (from-arg (and from (%api-real from)))
             (to-arg (and to (%api-real to))))
        (agree (list :curve c :from from :to to :tol tol :relative relative)
               (list (if default-tol
                         (batten::curve-length c :from from-arg :to to-arg)
                         (batten::curve-length c :from from-arg :to to-arg :tolerance tol :relative relative))
                     (%api-untouched c snap))
               (if (<= tol 0d0)
                   (values t nil)
                   (with-foreign ((cp c :curve))
                     (with-scalars ((len :double))
                       (if (c-flag-error-p (c-crv-arc-length cp (or from lo) (or to hi) tol
                                                             (if relative +relative+ +absolute+) len))
                           (values t nil)
                           (values nil (list (cffi:mem-ref len :double) t)))))))))))

;;; No curve with an interior knot of multiplicity degree + 1 (the C then
;;; writes Bw[-1]); fractions outside [0, 1] fail on both sides.
(define-oracle-test api-curve-point-at-length (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 12)
              ((0 1) (%crv-geom-line-curve))
              ((2 3) (%crv-geom-circle-curve))
              (t (random-curve :max-degree 6 :max-extra 10))))
         (snap (batten::crv-copy c))
         (per (case (rnd 12)
                (0 0d0) (1 1d0) (2 0.5d0)
                (3 (if (zerop (rnd 2)) -0.1d0 1.1d0))
                (t (rnd 1d0))))
         (per-arg (%api-real per))
         (relative (/= 0 (rnd 3)))
         (atl (if relative
                  (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 2d0 6d0)))) (t 1d-4))
                  (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 1d0 5d0)))) (t 0.05d0))))
         (len (case (rnd 6)
                (0 (handler-case (batten::crv-arc-length c (aref (curve-knots c) 0)
                                                         (aref (curve-knots c) (curve-m c))
                                                         1d-6 +relative+)
                     (nurbs-error () nil)))
                (t nil))))
    (agree (list :curve c :per per :atl atl :relative relative :len len)
           (append (multiple-value-list
                    (batten::curve-point-at-length c per-arg :tolerance atl :relative relative :length len))
                   (list (%api-untouched c snap)))
           (with-foreign ((cp c :curve))
             (with-foreign-out ((pt :point))
               (with-scalars ((l :double) (up :double))
                 (setf (cffi:mem-ref l :double) (or len -1d0))
                 (if (c-flag-error-p (c-crv-percentage-along-pt cp per atl batten::*mtol*
                                                                (if relative +relative+ +absolute+)
                                                                l up pt))
                     (values t nil)
                     (values nil (list (from-foreign pt :point) (cffi:mem-ref up :double)
                                       (cffi:mem-ref l :double) t)))))))))

;;; ----------------------------------------------------------- closest points

(defun %api-c-polyline-nearest-parameter (pts pars p)
  "What batten::%polyline-nearest-parameter computes, through NLib's
N_ProjectPtLineParam and N_DistPtPt: the parameter of the point of polyline
PTS (a point-array, PARS its vertices' parameters) nearest point P, or nil
when N_ProjectPtLineParam fails."
  (with-foreign ((pp p :point))
    (with-foreign-out ((a :point) (b :point) (q :point))
      (with-scalars ((alf :double) (d :double))
        (flet ((put (ptr i)
                 (dotimes (k 3) (setf (cffi:mem-aref ptr :double k) (aref pts i k))))
               (copy (from to)
                 (dotimes (k 3) (setf (cffi:mem-aref to :double k) (cffi:mem-aref from :double k)))))
          (put a 0)
          (c-dist-pt-pt pp a d)
          (let ((best-u (aref pars 0)) (best-d (cffi:mem-ref d :double)))
            (dotimes (j (1- (array-dimension pts 0)) best-u)
              (put a j)
              (put b (1+ j))
              (when (c-flag-error-p (c-project-pt-line-param pp a b q alf))
                (return nil))
              (let ((al (cffi:mem-ref alf :double)))
                (cond ((< al 0d0) (setf al 0d0) (copy a q))
                      ((> al 1d0) (setf al 1d0) (copy b q)))
                (c-dist-pt-pt pp q d)
                (when (< (cffi:mem-ref d :double) best-d)
                  (setf best-d (cffi:mem-ref d :double)
                        best-u (+ (aref pars j) (* al (- (aref pars (1+ j)) (aref pars j))))))))))))))

;;; The C side composes the same NLib routines: the box diagonal sizes a
;;; polyline, its point nearest P (N_ProjectPtLineParam on each leg) starts
;;; N_CrvClosestPt, and N_DistPtPt gives the distance.
(define-oracle-test api-curve-closest-point (:file "api-curve-queries" :trials 1000)
  (let* ((c (random-curve :max-degree 5 :max-extra 8))
         (c (if (zerop (rnd 3)) (%crv-geom-close c) c))
         (snap (batten::crv-copy c))
         (kn (curve-knots c))
         (us (random-parameter kn))
         (on (batten::crv-eval c us +left+ (make-point)))
         (p (case (rnd 5)
              (0 on)
              (1 (%crv-geom-offset on 1d0))
              (2 (%crv-geom-offset on 10d0))
              (3 (%crv-geom-offset on 100d0))
              (t (random-point))))
         (p-arg (%api-point p))
         (start (when (zerop (rnd 4)) (random-parameter kn)))
         (tol (case (rnd 5) (0 nil) (1 0d0) (2 1d-7) (3 1d-4) (t 1d-2)))
         (ctol (case (rnd 3) (0 nil) (1 1d-6) (t 1d-9)))
         (top (or tol batten::*mtol*))
         (toc (or ctol top)))
    (agree (list :curve c :p p :start start :tol tol :ctol ctol)
           (append (multiple-value-list
                    (apply #'batten::curve-closest-point c p-arg
                           (append (when tol (list :tolerance tol))
                                   (when ctol (list :cosine-tolerance ctol))
                                   (when start (list :start start)))))
                   (list (%api-untouched c snap)))
           (with-stacks (s)
             (with-foreign ((cp c :curve) (pp p :point))
               (with-foreign-out ((q :point))
                 (with-scalars ((d :double) (pts :pointer) (pars :pointer) (n :long)
                                (u :double) (dist :double))
                   (block c-side
                     (let ((u0 start))
                       (unless u0
                         (when (c-flag-error-p (c-crv-get-b-box-max-diag-dist cp d))
                           (return-from c-side (values t nil)))
                         (let* ((a (* 1d-3 (cffi:mem-ref d :double)))
                                (ptol (if (> a batten::*mtol*) a batten::*mtol*)))
                           (when (c-flag-error-p (c-approx-crv-with-polyline cp ptol +absolute+ +both+
                                                                             pts pars n s))
                             (return-from c-side (values t nil)))
                           (let ((nn (cffi:mem-ref n :long)))
                             (setf u0 (%api-c-polyline-nearest-parameter
                                       (from-foreign (cffi:mem-ref pts :pointer) :point-array (1+ nn))
                                       (from-foreign (cffi:mem-ref pars :pointer) :real-vector (1+ nn))
                                       p))
                             (unless u0
                               (return-from c-side (values t nil))))))
                       (when (c-flag-error-p (c-crv-closest-pt cp pp u0 top toc u q))
                         (return-from c-side (values t nil)))
                       (c-dist-pt-pt pp q dist)
                       (values nil (list (cffi:mem-ref u :double) (from-foreign q :point)
                                         (cffi:mem-ref dist :double) t)))))))))))

(define-oracle-test api-curve-project-points (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 10)
              (0 (%crv-geom-line-curve))
              (1 (%crv-geom-circle-curve))
              (t (random-curve :max-degree 5 :max-extra 8))))
         (snap (batten::crv-copy c))
         (k (rnd 8))
         (pts (%crv-geom-project-points c k))
         (pts-copy (copy-seq (make-array (* 3 (1+ k)) :element-type 'double-float :displaced-to pts)))
         (pts-arg (case (rnd 3)
                    (0 pts)
                    (1 (let (l) (dotimes (i (1+ k) (nreverse l))
                                  (push (list (aref pts i 0) (aref pts i 1) (aref pts i 2)) l))))
                    (t (let ((v (make-array (1+ k))))
                         (dotimes (i (1+ k) v)
                           (setf (svref v i) (make-point (aref pts i 0) (aref pts i 1) (aref pts i 2))))))))
         (default-tol (zerop (rnd 6)))
         (tol (cond (default-tol 1d-3)
                    ((zerop (rnd 25)) 0d0)
                    (t (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 1.5d0 4d0)))) (t 1d-2))))))
    (agree (list :curve c :pts pts :k k :tol tol)
           (append (multiple-value-list
                    (if default-tol
                        (batten::curve-project-points c pts-arg)
                        (batten::curve-project-points c pts-arg :tolerance tol)))
                   (list (and (%api-untouched c snap)
                              (equalp pts-copy (make-array (* 3 (1+ k)) :element-type 'double-float
                                                                        :displaced-to pts)))))
           (if (<= tol 0d0)
               (values t nil)
               (with-foreign ((cp c :curve) (pp pts :point-array))
                 (with-stacks (s)
                   (with-scalars ((tt :pointer) (u :pointer))
                     (if (c-flag-error-p (c-crv-project-pts cp pp k tol tt u s))
                         (values t nil)
                         (values nil
                                 (list (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ k))
                                       (from-foreign (cffi:mem-ref tt :pointer) :point-array (1+ k))
                                       t))))))))))

;;; A degenerate curve gets an absolute tolerance (its length is 0).
(define-oracle-test api-curve-polyline (:file "api-curve-queries" :trials 1000)
  (let* ((degen (zerop (rnd 12)))
         (c (if degen (%crv-geom-degenerate (random-curve :max-degree 5 :max-extra 6)) (%crv-geom-arc-curve)))
         (snap (batten::crv-copy c))
         (relative (and (not degen) (zerop (rnd 3))))
         (tol (cond ((zerop (rnd 25)) 0d0)
                    (relative (expt 10d0 (- (rnd-range 1.5d0 4d0))))
                    (t (case (rnd 4) (0 1d0) (t (expt 10d0 (- (rnd-range -0.5d0 1.7d0))))))))
         (tol-arg (%api-real tol)))
    (agree (list :curve c :tol tol :relative relative)
           (append (multiple-value-list
                    (if relative
                        (batten::curve-polyline c tol-arg :relative t)
                        (batten::curve-polyline c tol-arg)))
                   (list (%api-untouched c snap)))
           (if (<= tol 0d0)
               (values t nil)
               (with-foreign ((cp c :curve))
                 (with-stacks (s)
                   (with-scalars ((pts :pointer) (pars :pointer) (n :long))
                     (if (c-flag-error-p (c-approx-crv-with-polyline cp tol (if relative +relative+ +absolute+)
                                                                     +both+ pts pars n s))
                         (values t nil)
                         (let ((nn (cffi:mem-ref n :long)))
                           (values nil
                                   (list (from-foreign (cffi:mem-ref pts :pointer) :point-array (1+ nn))
                                         (from-foreign (cffi:mem-ref pars :pointer) :real-vector (1+ nn))
                                         t)))))))))))

;;; --------------------------------------------------------------- predicates

(define-oracle-test api-curve-closed-p (:file "api-curve-queries" :trials 1000)
  (let* ((c (if (zerop (rnd 2)) (random-curve) (%closed-curve))) (snap (batten::crv-copy c)))
    (agree c
           (list (if (batten::curve-closed-p c) 1 0) (%api-untouched c snap))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-closed cp) t))))))

(define-oracle-test api-curve-planar-p (:file "api-curve-queries" :trials 1000)
  (let* ((tol (case (rnd 5) (0 nil) (1 1d-3) (2 1d-6) (3 0.5d0) (t 0d0)))
         (gen-tol (if (or (null tol) (zerop tol)) 1d-6 tol))
         (c (case (rnd 6)
              (0 (random-curve))
              (1 (%line-curve gen-tol))
              (2 (%degenerate-curve))
              (t (%planar-curve gen-tol))))
         (snap (batten::crv-copy c)))
    (agree (list :curve c :tol tol)
           (multiple-value-bind (planar pp nn)
               (if tol (batten::curve-planar-p c :tolerance tol) (batten::curve-planar-p c))
             (list (if planar 1 0) pp nn (%api-untouched c snap)))
           (with-stacks (s)
             (with-foreign ((cp c :curve))
               (with-foreign-out ((pp :point) (nn :point))
                 (with-scalars ((flt :short))
                   (if (c-flag-error-p (c-crv-is-planar cp (or tol batten::*mtol*) pp nn flt s))
                       (values t nil)
                       (values nil (if (= (cffi:mem-ref flt :short) +true+)
                                       (list 1 (from-foreign pp :point) (from-foreign nn :point) t)
                                       (list 0 nil nil t)))))))))))

(define-oracle-test api-curve-line-p (:file "api-curve-queries" :trials 1000)
  (let* ((tol (case (rnd 5) (0 nil) (1 1d-3) (2 1d-6) (3 0.5d0) (t 0d0)))
         (c (if (zerop (rnd 5))
                (random-curve)
                (%line-curve (if (or (null tol) (zerop tol)) 1d-6 tol))))
         (snap (batten::crv-copy c)))
    (agree (list :curve c :tol tol)
           (list (if (if tol (batten::curve-line-p c :tolerance tol) (batten::curve-line-p c)) 1 0)
                 (%api-untouched c snap))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-line cp (or tol batten::*mtol*)) t))))))

(define-oracle-test api-curve-degenerate-p (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 4)
              (0 (random-curve))
              (1 (%degenerate-curve))
              (2 (%degenerate-curve :spread 4d-8))
              (t (%degenerate-curve :spread 1d-7))))
         (snap (batten::crv-copy c)))
    (agree c
           (list (if (batten::curve-degenerate-p c) 1 0) (%api-untouched c snap))
           (with-foreign ((cp c :curve))
             (values nil (list (c-crv-is-degen cp) t))))))

;;; N_CrvsAreCoincident makes its two curves compatible in place; the API
;;; works on copies and leaves both arguments alone.
(define-oracle-test api-curves-coincident-p (:file "api-curve-queries" :trials 1000)
  (let* ((c (random-curve :max-degree 6 :max-extra 8))
         (d (%coincident-partner c))
         (snap-c (batten::crv-copy c)) (snap-d (batten::crv-copy d))
         (tol (case (rnd 6) (0 nil) (1 0d0) (2 1d-9) (3 1d-6) (4 1d-2) (t 10d0))))
    (agree (list :curve-1 c :curve-2 d :tol tol)
           (list (if (if tol
                         (batten::curves-coincident-p c d :tolerance tol)
                         (batten::curves-coincident-p c d))
                     1 0)
                 (%api-untouched c snap-c) (%api-untouched d snap-d))
           (with-foreign ((cp c :curve) (cq d :curve))
             (with-stacks (s)
               (values nil (list (c-crvs-are-coincident cp cq (or tol batten::*mtol*) s) t t)))))))

;;; --------------------------------------------------------------------- type

(define-oracle-test api-curve-type (:file "api-curve-queries" :trials 1000)
  (let* ((c (case (rnd 6)
              ((0 1 2) (%cn-circle-curve))
              (3 (%cn-line-curve))
              (4 (%curve-with-flat-start))
              (t (random-curve))))
         (snap (batten::crv-copy c))
         (tol (case (rnd 5) (0 nil) (1 1d-8) (2 1d-5) (3 1d-3) (t 0.1d0))))
    (agree (list :curve c :tol tol)
           (append (multiple-value-list
                    (if tol (batten::curve-type c :tolerance tol) (batten::curve-type c)))
                   (list (%api-untouched c snap)))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((aa :point) (bb :point) (xx :point) (yy :point))
               (with-scalars ((r :double) (al :double) (typ :short))
                 (if (c-flag-error-p (c-crv-get-type cur (or tol batten::*mtol*) aa bb xx yy r al typ))
                     (values t nil)
                     (let ((typ (cffi:mem-ref typ :short)))
                       (values nil
                               (append (cond ((= typ +nline+)
                                              (list :line (from-foreign aa :point) (from-foreign bb :point)))
                                             ((= typ +ncircle+)
                                              (list :circle (from-foreign aa :point) (from-foreign bb :point)
                                                    (cffi:mem-ref r :double) (cffi:mem-ref al :double)
                                                    (from-foreign xx :point) (from-foreign yy :point)))
                                             (t (list :freeform)))
                                       (list t)))))))))))
