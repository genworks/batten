;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-adv.lisp against NLib.

(in-package :batten-tests)

(defun %collapse-first-row (s)
  "Surface S with its first u-row of control points collapsed onto one point,
so that S_v and the normal vanish along u = u_min (a pole)."
  (let ((net (surface-cpoints s)))
    (dotimes (j (array-dimension net 1) s)
      (dotimes (c 4)
        (setf (aref net 0 j c) (aref net 0 0 c))))))

(define-oracle-test srf-eval-pt-normal-deriv (:file "srf-adv" :trials 8000)
  (let* ((s (random-surface :two-d (zerop (rnd 5))))
         (pole (zerop (rnd 8)))
         (s (if pole (%collapse-first-row s) s))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (u (if (and pole (zerop (rnd 2))) (aref ku 0) (random-parameter ku)))
         (v (random-parameter kv))
         (ufl (random-side)) (vfl (random-side)))
    (agree (list :surface s :u u :v v :ufl ufl :vfl vfl)
           (list (batten::srf-eval-pt-normal-deriv s u v ufl vfl (make-point-net 3 3)))
           (with-foreign ((sp s :surface))
             (with-foreign-out ((nd :point-net 3 3))
               (values (c-flag-error-p (c-srf-eval-pt-normal-deriv sp u v ufl vfl nd))
                       (list (from-foreign nd :point-net 3 3))))))))

(define-oracle-test srf-offset-get-max-deriv (:file "srf-adv" :trials 3000)
  (let* ((s (random-surface :two-d (zerop (rnd 5)) :max-degree 4 :max-extra 4))
         (s (if (zerop (rnd 6)) (%collapse-first-row s) s))
         (d (case (rnd 6) (0 0d0) (t (- (rnd 20d0) 10d0)))))
    (agree (list :surface s :d d)
           (multiple-value-list (batten::srf-offset-get-max-deriv s d))
           (with-foreign ((sp s :surface))
             (with-scalars ((muu :double) (muv :double) (mvv :double))
               (values (c-flag-error-p (c-srf-offset-get-max-deriv sp d muu muv mvv))
                       (list (cffi:mem-ref muu :double) (cffi:mem-ref muv :double)
                             (cffi:mem-ref mvv :double))))))))

(defun %srf-adv-grid (knots count)
  "COUNT+1 ascending parameters in the domain of KNOTS: the first now and then
the start, the last now and then the end, the ones between strictly inside
(NLib's span walk runs off the knot vector for an inner parameter at the end),
and some of them knots.  Now and then the first or last leaves the domain."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (inside (remove-if (lambda (x) (or (= x lo) (= x hi))) (coerce knots 'list)))
         (xs nil))
    (dotimes (i (1+ count))
      (push (if (and inside (zerop (rnd 4)))
                (nth (rnd (length inside)) inside)
                (rnd-range lo hi))
            xs))
    (let ((v (coerce (sort xs #'<) 'real-vector)))
      (when (zerop (rnd 4)) (setf (aref v 0) lo))
      (when (zerop (rnd 4)) (setf (aref v count) hi))
      (case (rnd 30)
        (0 (setf (aref v 0) (- lo 0.25d0)))
        (1 (setf (aref v count) (+ hi 0.25d0))))
      v)))

(define-oracle-test srf-get-grid-pts-normals (:file "srf-adv" :trials 5000)
  (let* ((s (random-surface :two-d (zerop (rnd 5))))
         (s (if (zerop (rnd 6)) (%collapse-first-row s) s))
         (mu (if (zerop (rnd 40)) -1 (rnd 6)))
         (mv (if (zerop (rnd 40)) -1 (rnd 6)))
         (u (%srf-adv-grid (surface-knots-u s) (max mu 0)))
         (v (%srf-adv-grid (surface-knots-v s) (max mv 0)))
         (ufl (random-side)) (vfl (random-side))
         (rows (1+ (max mu 0))) (cols (1+ (max mv 0))))
    (agree (list :surface s :u u :v v :mu mu :mv mv :ufl ufl :vfl vfl)
           (let ((p (make-point-net rows cols)) (n (make-point-net rows cols)))
             (batten::srf-get-grid-pts-normals s u v mu mv ufl vfl p n)
             (list p n))
           (with-foreign ((sp s :surface) (cu u :real-vector) (cv v :real-vector))
             (with-foreign-out ((p :point-net rows cols) (n :point-net rows cols))
               (if (c-flag-error-p (c-srf-get-grid-pts-normals sp cu cv mu mv ufl vfl p n))
                   (values t nil)
                   (values nil (list (from-foreign p :point-net rows cols)
                                     (from-foreign n :point-net rows cols)))))))))

(define-oracle-test srf-eval-pt-grid (:file "srf-adv" :trials 5000)
  (let* ((s (random-surface :two-d (zerop (rnd 5))))
         (nn (if (zerop (rnd 40)) -1 (rnd 9)))
         (mm (if (zerop (rnd 40)) -1 (rnd 9)))
         (u (%srf-adv-grid (surface-knots-u s) (max nn 0)))
         (v (%srf-adv-grid (surface-knots-v s) (max mm 0)))
         (ufl (random-side)) (vfl (random-side))
         (rows (1+ (max nn 0))) (cols (1+ (max mm 0))))
    (agree (list :surface s :u u :v v :nn nn :mm mm :ufl ufl :vfl vfl)
           (list (batten::srf-eval-pt-grid s u v nn mm ufl vfl (make-point-net rows cols)))
           (with-foreign ((sp s :surface) (cu u :real-vector) (cv v :real-vector))
             (with-foreign-out ((pts :point-net rows cols))
               (if (c-flag-error-p (c-srf-eval-pt-grid sp cu cv nn mm ufl vfl pts))
                   (values t nil)
                   (values nil (list (from-foreign pts :point-net rows cols)))))))))

(define-oracle-test create-b-vectors (:file "srf-adv" :trials 4000)
  (let* ((c (random-curve))
         (pw (curve-cpoints c))
         (closed (zerop (rnd 3)))
         (k (rnd 9))
         (knots (curve-knots c))
         (lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (tt (make-array (1+ k) :element-type 'double-float :initial-element 0d0))
         (b0 (if (zerop (rnd 40)) (make-point 0d0 0d0 0d0) (random-point))))
    (when closed
      (dotimes (j 4)
        (setf (aref pw (1- (array-dimension pw 0)) j) (aref pw 0 j))))
    (dotimes (i (1+ k))
      (setf (aref tt i) (case (rnd 10) (0 lo) (1 hi) (t (rnd-range lo hi)))))
    (setf tt (sort tt #'<))
    (agree (list :curve c :b0 b0 :t tt :k k)
           (list (batten::create-b-vectors c b0 tt k (make-point-array (1+ k))))
           (with-foreign ((cp c :curve) (cb0 b0 :point) (ct tt :real-vector))
             (with-foreign-out ((b :point-array (1+ k)))
               (values (c-flag-error-p (c-create-b-vectors cp cb0 ct k b))
                       (list (from-foreign b :point-array (1+ k)))))))))

;;; ------------------------------------------------------------- L7-g2

(defun %srf-adv-copy-surface (s)
  "A deep copy of surface S."
  (%make-surface (surface-degree-u s) (surface-degree-v s)
                 (copy-seq (surface-knots-u s)) (copy-seq (surface-knots-v s))
                 (let* ((net (surface-cpoints s))
                        (c (make-array (array-dimensions net) :element-type 'double-float)))
                   (dotimes (i (array-total-size net) c)
                     (setf (row-major-aref c i) (row-major-aref net i))))))

(defun %srf-adv-uv (knots)
  "A parameter for N_SrfModifyBoundaryCrv in the domain of KNOTS: usually
inside, now and then a knot, an end, a hair (relative 1e-7 or 1e-5) off a
knot, or outside the domain."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (kn (aref knots (rnd (length knots)))))
    (case (rnd 12)
      (0 lo)
      (1 hi)
      ((2 3) kn)
      (4 (min hi (max lo (+ kn (* (- hi lo) (if (zerop (rnd 2)) 5d-7 -5d-7))))))
      (5 (min hi (max lo (+ kn (* (- hi lo) (if (zerop (rnd 2)) 1d-5 -1d-5))))))
      (6 (if (zerop (rnd 2)) (- lo 0.1d0) (+ hi 0.1d0)))
      (t (rnd-range lo hi)))))

(define-oracle-test srf-modify-boundary-crv (:file "srf-adv" :trials 1500)
  (let* ((two-d (zerop (rnd 8)))
         (s (random-surface :max-degree 4 :max-extra 4 :two-d two-d))
         (c (random-curve :max-degree 4 :max-extra 5
                          :two-d (if (zerop (rnd 10)) (not two-d) two-d)))
         (dirflg (if (zerop (rnd 2)) +udir+ +vdir+))
         (end (if (zerop (rnd 2)) +start+ +end+))
         (uv (%srf-adv-uv (if (= dirflg +udir+) (surface-knots-u s) (surface-knots-v s))))
         (cont (case (rnd 3) (0 +c-1+) (1 +c-2+) (t +cmax+)))
         (in-place (zerop (rnd 5))))
    (agree (list :sur-p s :newbdy c :dirflg dirflg :end end :uv uv :cont cont :in-place in-place)
           (list (batten::srf-modify-boundary-crv (%srf-adv-copy-surface s) (batten::crv-copy c)
                                                  dirflg end uv cont))
           (with-foreign ((sp s :surface) (cp c :curve))
             (with-stacks (s1)
               (with-stacks (s2)
                 (if in-place
                     (let ((failed (c-flag-error-p
                                    (c-srf-modify-boundary-crv sp cp dirflg end uv cont sp s1 s2))))
                       (values failed (unless failed (list (from-foreign sp :surface)))))
                     (with-foreign-out ((sq :surface))
                       (let ((failed (c-flag-error-p
                                      (c-srf-modify-boundary-crv sp cp dirflg end uv cont sq s1 s2))))
                         (values failed (unless failed (list (from-foreign sq :surface)))))))))))))

(defun %srf-adv-ascending (lo hi count)
  "COUNT+1 ascending parameters from LO to HI, the ends exactly LO and HI."
  (let ((xs (sort (let (l) (dotimes (i (max 0 (1- count)) l)
                             (push (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0)))) l)))
                  #'<)))
    (coerce (append (list lo) xs (list hi)) 'real-vector)))

(defun %srf-adv-params (lo hi base wild)
  "Ascending parameters from LO to HI, one per entry of BASE (ascending
fractions from 0 to 1): the fractions jittered by a fifth of the smallest gap,
or, when WILD, drawn afresh."
  (let* ((count (1- (length base)))
         (ps (if wild (%srf-adv-ascending 0d0 1d0 count) (copy-seq base)))
         (gap 1d0))
    (dotimes (i count) (setf gap (min gap (- (aref base (1+ i)) (aref base i)))))
    (unless wild
      (do ((i 1 (1+ i))) ((>= i count))
        (incf (aref ps i) (* gap (- (rnd 0.4d0) 0.2d0)))))
    (dotimes (i (1+ count) ps)
      (setf (aref ps i) (cond ((= i 0) lo) ((= i count) hi) (t (+ lo (* (- hi lo) (aref ps i)))))))))

(defun %srf-adv-network (k l)
  "A curve network for N_ReparmCrvsIsectPt, as (values cur-u cur-v uu vv):
K+1 u-curves and L+1 v-curves, and the (L+1) x (K+1) matrices of their
intersection parameters, each curve's running from its start to its end,
usually near common fractions of the curves' domains."
  (let ((cur-u (make-array (1+ k))) (cur-v (make-array (1+ l)))
        (uu (make-real-matrix (1+ l) (1+ k))) (vv (make-real-matrix (1+ l) (1+ k)))
        (rat (case (rnd 3) (0 nil) (1 t) (t :random)))
        (two-d (case (rnd 4) (0 t) (t nil)))
        (wild (zerop (rnd 6)))
        (ubase (%srf-adv-ascending 0d0 1d0 l))
        (vbase (%srf-adv-ascending 0d0 1d0 k)))
    (dotimes (j (1+ k))
      (let* ((c (random-curve :max-degree 3 :max-extra 5 :rational rat :two-d two-d))
             (kn (curve-knots c))
             (ps (%srf-adv-params (aref kn 0) (aref kn (1- (length kn))) ubase wild)))
        (setf (svref cur-u j) c)
        (dotimes (i (1+ l)) (setf (aref uu i j) (aref ps i)))))
    (dotimes (i (1+ l))
      (let* ((c (random-curve :max-degree 3 :max-extra 5 :rational rat :two-d two-d))
             (kn (curve-knots c))
             (ps (%srf-adv-params (aref kn 0) (aref kn (1- (length kn))) vbase wild)))
        (setf (svref cur-v i) c)
        (dotimes (j (1+ k)) (setf (aref vv i j) (aref ps j)))))
    (values cur-u cur-v uu vv)))

(defun %srf-adv-curves-to-foreign (curves)
  "A foreign array of pointers to foreign copies of CURVES (inside WITH-FOREIGN)."
  (let ((arr (batten-oracle::falloc :pointer (length curves))))
    (dotimes (i (length curves) arr)
      (setf (cffi:mem-aref arr :pointer i) (to-foreign (svref curves i) :curve)))))

(defun %srf-adv-curves-from-foreign (arr count)
  "COUNT curves read back from the foreign pointer array ARR."
  (let (out)
    (dotimes (i count (nreverse out))
      (push (from-foreign (cffi:mem-aref arr :pointer i) :curve) out))))

(define-oracle-test reparm-crvs-isect-pt (:file "srf-adv" :trials 800)
  (let* ((rep (case (rnd 20) (0 3) (t (if (zerop (rnd 2)) +linear+ +quadratic+))))
         ;; quadratic needs two intervals; now and then it gets one
         (kmin (if (and (= rep +quadratic+) (plusp (rnd 10))) 2 1))
         (k (+ kmin (rnd (- 5 kmin)))) (l (+ kmin (rnd (- 5 kmin)))))
    (multiple-value-bind (cur-u cur-v uu vv) (%srf-adv-network k l)
      (agree (list :cur-u cur-u :cur-v cur-v :k k :l l :uu uu :vv vv :rep rep)
             (let ((cu (map 'simple-vector #'batten::crv-copy cur-u))
                   (cv (map 'simple-vector #'batten::crv-copy cur-v)))
               (multiple-value-bind (u v) (batten::reparm-crvs-isect-pt cu cv k l uu vv rep)
                 (list u v (coerce cu 'list) (coerce cv 'list))))
             (with-foreign ((fuu uu :real-matrix) (fvv vv :real-matrix))
               (with-stacks (sc)
                 (with-stacks (sp)
                   (with-scalars ((u :pointer) (v :pointer))
                     (let* ((au (%srf-adv-curves-to-foreign cur-u))
                            (av (%srf-adv-curves-to-foreign cur-v))
                            (failed (c-flag-error-p
                                     (c-reparm-crvs-isect-pt au av k l fuu fvv rep u v sc sp))))
                       (values failed
                               (unless failed
                                 (list (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ l))
                                       (from-foreign (cffi:mem-ref v :pointer) :real-vector (1+ k))
                                       (%srf-adv-curves-from-foreign au (1+ k))
                                       (%srf-adv-curves-from-foreign av (1+ l))))))))))))))

(defun %srf-adv-iso-curve (net dir idx p knots two-d)
  "The boundary curve of control net NET at row (DIR :u, fixed first index)
or column (:v) IDX, degree P on KNOTS."
  (let* ((len (if (eq dir :u) (array-dimension net 1) (array-dimension net 0)))
         (pw (make-cpoint-array len)))
    (dotimes (i len)
      (dotimes (c 3)
        (setf (aref pw i c) (if (eq dir :u) (aref net idx i c) (aref net i idx c))))
      (when two-d (setf (aref pw i 2) +no-z+))
      (setf (aref pw i 3) +no-w+))
    (%make-curve p (copy-seq knots) pw)))

(defun %srf-adv-cross-curve (net dir i0 i1 fac p knots two-d)
  "The cross-boundary derivative curve FAC * (row I1 - row I0) of net NET
(DIR :u: rows of the first index; :v: columns), degree P on KNOTS."
  (let* ((len (if (eq dir :u) (array-dimension net 1) (array-dimension net 0)))
         (pw (make-cpoint-array len)))
    (dotimes (i len)
      (dotimes (c 3)
        (setf (aref pw i c)
              (if (eq dir :u)
                  (* fac (- (aref net i1 i c) (aref net i0 i c)))
                  (* fac (- (aref net i i1 c) (aref net i i0 c))))))
      (when two-d (setf (aref pw i 2) +no-z+))
      (setf (aref pw i 3) +no-w+))
    (%make-curve p (copy-seq knots) pw)))

(defun %srf-adv-coons-data ()
  "Boundary curves and cross-boundary derivatives for N_CreateCoonsBoundaryCrvs,
as a list (curL derL curR derR curB derB curT derT), taken from a random
non-rational surface (so they satisfy the Coons conditions), then now and
then refined, elevated, nudged, moved to another domain, made rational, or
replaced by unrelated curves on the same domain."
  (let* ((two-d (zerop (rnd 6)))
         (s (random-surface :rational nil :max-degree 4 :max-extra 4 :two-d two-d))
         (net (surface-cpoints s))
         (p (surface-degree-u s)) (q (surface-degree-v s))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (n (batten::surface-n s)) (m (batten::surface-m s))
         (r (1- (length ku))) (ss (1- (length kv)))
         (cur-l (%srf-adv-iso-curve net :u 0 q kv two-d))
         (cur-r (%srf-adv-iso-curve net :u n q kv two-d))
         (der-l (%srf-adv-cross-curve net :u 0 1 (/ p (- (aref ku (1+ p)) (aref ku 0))) q kv two-d))
         (der-r (%srf-adv-cross-curve net :u (1- n) n (/ p (- (aref ku r) (aref ku (- r p 1)))) q kv two-d))
         (cur-b (%srf-adv-iso-curve net :v 0 p ku two-d))
         (cur-t (%srf-adv-iso-curve net :v m p ku two-d))
         (der-b (%srf-adv-cross-curve net :v 0 1 (/ q (- (aref kv (1+ q)) (aref kv 0))) p ku two-d))
         (der-t (%srf-adv-cross-curve net :v (1- m) m (/ q (- (aref kv ss) (aref kv (- ss q 1)))) p ku two-d))
         (curves (vector cur-l der-l cur-r der-r cur-b der-b cur-t der-t)))
    ;; unrelated curves on the same domains
    (when (zerop (rnd 6))
      (dotimes (i 8)
        (when (zerop (rnd 2))
          (let* ((kn (curve-knots (svref curves i)))
                 (c (random-curve :max-degree 4 :max-extra 4 :rational nil :two-d two-d))
                 (ck (curve-knots c))
                 (lo (aref kn 0)) (hi (aref kn (1- (length kn))))
                 (clo (aref ck 0)) (chi (aref ck (1- (length ck)))))
            (dotimes (j (length ck))
              (setf (aref ck j) (cond ((= (aref ck j) clo) lo) ((= (aref ck j) chi) hi)
                                      (t (+ lo (* (- hi lo) (/ (- (aref ck j) clo) (- chi clo))))))))
            (setf (svref curves i) c)))))
    (dotimes (pass (rnd 3))
      (let* ((i (rnd 8)) (c (svref curves i)))
        (setf (svref curves i)
              (case (rnd 12)
                ((0 1 2) (let* ((kn (curve-knots c)) (lo (aref kn 0)) (hi (aref kn (1- (length kn))))
                                (x (make-real-vector (1+ (rnd 2)))))
                           (dotimes (j (length x)) (setf (aref x j) (rnd-range (+ lo (* 0.1d0 (- hi lo))) (- hi (* 0.1d0 (- hi lo))))))
                           (batten::crv-refine c (sort x #'<))))
                ((3 4) (batten::crv-elevate-degree c 1))
                ((5 6) (let ((pw (curve-cpoints c)))
                         (incf (aref pw (rnd (array-dimension pw 0)) 0) (case (rnd 3) (0 1d-9) (1 1d-6) (t 1d-2)))
                         c))
                (7 (let ((kn (curve-knots c)))
                     (incf (aref kn (1- (length kn))) 0.5d0)
                     c))
                (8 (batten::crv-non-rat-to-rat c))
                (t c)))))
    (coerce curves 'list)))

(define-oracle-test create-coons-boundary-crvs (:file "srf-adv" :trials 1500)
  (let* ((curves (%srf-adv-coons-data))
         (cpf (if (zerop (rnd 2)) +yes+ +no+)))
    (agree (list :curves curves :cpf cpf)
           (let* ((cs (mapcar #'batten::crv-copy curves))
                  (sur (apply #'batten::create-coons-boundary-crvs (append cs (list cpf)))))
             (list sur cs))
           (with-foreign ((dummy (make-real-vector 1) :real-vector))
             (declare (ignore dummy))
             (with-stacks (sc)
               (with-stacks (ss)
                 (let ((fc (mapcar (lambda (c) (to-foreign c :curve)) curves)))
                   (with-foreign-out ((sur :surface))
                     (let ((failed (c-flag-error-p
                                    (apply #'c-create-coons-boundary-crvs
                                           (append fc (list cpf sur sc ss))))))
                       (values failed
                               (unless failed
                                 (list (from-foreign sur :surface)
                                       (mapcar (lambda (c) (from-foreign c :curve)) fc)))))))))))))

;;; ------------------------------------------------------------- L8-g0

(defun %srf-adv-skin-sections (k)
  "K+1 cross-sectional curves for skinning, as (values curves closed-p): a
random base curve moved along a random direction with its control points
jittered, or (now and then) unrelated random curves; mixed degrees,
rationality and dimension now and then; half the time the last curve a copy
of the first (a closed set)."
  (let* ((rat (case (rnd 3) (0 nil) (1 t) (t :random)))
         (two-d (case (rnd 6) (0 t) (1 :random) (t nil)))
         (closed (zerop (rnd 2)))
         (curves (make-array (1+ k))))
    (if (zerop (rnd 5))
        (dotimes (i (1+ k))
          (setf (svref curves i) (random-curve :max-degree 3 :max-extra 4 :rational rat :two-d two-d)))
        (let* ((base (random-curve :max-degree 3 :max-extra 4
                                   :rational (if (eq rat :random) (zerop (rnd 2)) rat)
                                   :two-d (if (eq two-d :random) (zerop (rnd 2)) two-d)))
               (dir (random-point 30d0))
               (jit (case (rnd 3) (0 0d0) (1 1d0) (t 10d0))))
          (dotimes (i (1+ k))
            (let* ((c (batten::crv-copy base))
                   (pw (curve-cpoints c)))
              (dotimes (j (array-dimension pw 0))
                (let ((w (if (= (aref pw j 3) +no-w+) 1d0 (aref pw j 3))))
                  (dotimes (d 3)
                    (unless (= (aref pw j d) +no-z+)
                      (incf (aref pw j d)
                            (* w (+ (* i (aref dir d)) (- (rnd (* 2 jit)) jit))))))))
              (setf (svref curves i) c)))
          ;; now and then a section of another degree
          (when (zerop (rnd 4))
            (setf (svref curves (rnd (1+ k)))
                  (random-curve :max-degree 4 :max-extra 3
                                :rational (/= (aref (curve-cpoints base) 0 3) +no-w+) :two-d (= (aref (curve-cpoints base) 0 2) +no-z+))))))
    (when closed
      (setf (svref curves k) (batten::crv-copy (svref curves 0))))
    (values curves closed)))

(defun %srf-adv-skin-inputs ()
  "Random arguments for N_CreateSkinSrf / N_CreateSkinSrfParams, as (values
curves k cmp cld deg dir closed-p): K from 1 to 5, the degree at most K
(the C's interpolation needs it); CMP +no+ only for curves made compatible
beforehand; DIR now and then neither direction."
  (let ((k (1+ (rnd 5))))
    (multiple-value-bind (curves closed) (%srf-adv-skin-sections k)
      (let ((cmp (if (zerop (rnd 4)) +no+ +yes+))
            (cld (if (zerop (rnd 4)) +no+ +yes+))
            ;; a closed set needs degree 2 or more (N_FitCrvKnotsAndDerivs)
            (deg (if (and closed (> k 1) (plusp (rnd 8)))
                     (+ 2 (rnd (1- (min k 4))))
                     (1+ (rnd (min k 4)))))
            (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+)))))
        (when (= cmp +no+)
          (handler-case (batten::with-c-float-environment (batten::crvs-make-compatible curves k))
            (nurbs-error () (setf cmp +yes+))))
        (values curves k cmp cld deg dir (and closed (= cld +yes+)))))))

(define-oracle-test create-skin-srf (:file "srf-adv" :trials 1500)
  (multiple-value-bind (curves k cmp cld deg dir) (%srf-adv-skin-inputs)
    (let ((knot-tol (case (rnd 6) (0 1d-3) (1 1d-6) (2 0.05d0) (t 0d0))))
      (agree (list :cur curves :k k :cmp cmp :cld cld :deg deg :dir dir :knot-tol knot-tol)
             (let ((cs (map 'simple-vector #'batten::crv-copy curves)))
               (list (batten::create-skin-srf cs k cmp cld deg dir knot-tol) (coerce cs 'list)))
             (with-foreign ((dummy (make-real-vector 1) :real-vector))
               (declare (ignore dummy))
               (with-stacks (sc)
                 (with-stacks (ss)
                   (let ((fc (%srf-adv-curves-to-foreign curves)))
                     (with-foreign-out ((sur :surface))
                       (let ((failed (c-flag-error-p
                                      (c-create-skin-srf fc k cmp cld deg dir knot-tol sur sc ss))))
                         (values failed
                                 (unless failed
                                   (list (from-foreign sur :surface)
                                         (%srf-adv-curves-from-foreign fc (1+ k)))))))))))))))

(define-oracle-test create-skin-srf-params (:file "srf-adv" :trials 1500)
  (multiple-value-bind (curves k cmp cld deg dir closed) (%srf-adv-skin-inputs)
    ;; parameters ascending over [0,1] or another domain; the knots averaged
    ;; from them (end-derivative knots for a closed set), now and then random
    ;; clamped knots or the wrong count
    (multiple-value-bind (lo hi) (if (zerop (rnd 3)) (values 0d0 1d0) (batten-random::random-domain))
      (let* ((tt (%srf-adv-ascending lo hi k))
             (count (case (rnd 15) (0 (+ k deg 3)) (1 (+ k deg 1)) (t (if closed (+ k deg 4) (+ k deg 2)))))
             (knt (batten::with-c-float-environment
                    (cond ((zerop (rnd 6))
                           (random-knots deg (- count deg 1) :lo lo :hi hi))
                          ((= count (+ k deg 4))
                           (batten::fit-calc-knot-vector-end-derivs tt k deg (make-real-vector count)))
                          (t (let ((kv (make-real-vector count)))
                               (if (= count (+ k deg 2))
                                   (batten::fit-crv-calc-knot-vector tt k deg kv)
                                   (random-knots deg (- count deg 1) :lo lo :hi hi))))))))
        (agree (list :cur curves :k k :cmp cmp :cld cld :deg deg :dir dir :t tt :knt knt)
               (let ((cs (map 'simple-vector #'batten::crv-copy curves)))
                 (list (batten::create-skin-srf-params cs k cmp cld deg dir tt knt) (coerce cs 'list)))
               (with-foreign ((ft tt :real-vector) (fk knt :knots))
                 (with-stacks (sc)
                   (with-stacks (ss)
                     (let ((fc (%srf-adv-curves-to-foreign curves)))
                       (with-foreign-out ((sur :surface))
                         (let ((failed (c-flag-error-p
                                        (c-create-skin-srf-params fc k cmp cld deg dir ft fk sur sc ss))))
                           (values failed
                                   (unless failed
                                     (list (from-foreign sur :surface)
                                           (%srf-adv-curves-from-foreign fc (1+ k))))))))))))))))

;;; ------------------------------------------------------------- L9-g0

(defun %srf-adv-spine-params (knots k)
  "K+1 ascending spine parameters in the domain of KNOTS: usually from end to
end, now and then over an interior stretch, some of them knots now and then."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (a (if (zerop (rnd 4)) (rnd-range lo (+ lo (* 0.3d0 (- hi lo)))) lo))
         (b (if (zerop (rnd 4)) (rnd-range (- hi (* 0.3d0 (- hi lo))) hi) hi))
         (ps (%srf-adv-ascending a b k)))
    (when (zerop (rnd 4))
      (let ((inside (remove-if (lambda (x) (or (<= x a) (>= x b))) (coerce knots 'list))))
        (when (and inside (> k 1))
          (let ((i (1+ (rnd (1- k)))) (x (nth (rnd (length inside)) inside)))
            (when (and (< (aref ps (1- i)) x) (< x (aref ps (1+ i))))
              (setf (aref ps i) x))))))
    ps))

(define-oracle-test create-skin-spine (:file "srf-adv" :trials 1200)
  (let* ((k (1+ (rnd 5)))
         (deriv (case (rnd 40) (0 7) (t (case (rnd 3) (0 +noder+) (1 +endder+) (t +allder+)))))
         ;; the degree in range, now and then out of it
         (degree (cond ((= deriv +noder+) (if (zerop (rnd 30)) (1+ k) (1+ (rnd (min k 3)))))
                       ((= deriv +endder+) (case (rnd 40) (0 1) (1 (+ k 3)) (t (+ 2 (rnd (min (1+ k) 3))))))
                       ((= deriv +allder+) (case (rnd 40) (0 1) (1 4) (t (+ 2 (rnd 2)))))
                       (t 2)))
         (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+))))
         (spine (random-curve :max-degree 4 :max-extra 5 :two-d (zerop (rnd 10))))
         (params (%srf-adv-spine-params (curve-knots spine) k))
         (align (if (zerop (rnd 3)) +no+ +yes+))
         (z0 (if (zerop (rnd 40)) (make-point 0d0 0d0 0d0) (random-point)))
         (subdiv (- (rnd 7) 1))
         (cmp (if (zerop (rnd 4)) +no+ +yes+)))
    (multiple-value-bind (curves closed) (%srf-adv-skin-sections k)
      (declare (ignore closed))
      (when (= cmp +no+)
        (handler-case (batten::with-c-float-environment (batten::crvs-make-compatible curves k))
          (nurbs-error () (setf cmp +yes+))))
      (agree (list :profiles curves :k k :cmp cmp :spine spine :params params :align align :z0 z0
                   :subdiv subdiv :deriv deriv :degree degree :dir dir)
             (let ((cs (map 'simple-vector #'batten::crv-copy curves)))
               (list (batten::create-skin-spine cs k cmp spine params align (copy-seq z0) subdiv
                                                deriv degree dir)
                     (coerce cs 'list)))
             (with-foreign ((sp spine :curve) (fp params :real-vector) (fz z0 :point))
               (with-stacks (sc)
                 (with-stacks (ss)
                   (let ((fc (%srf-adv-curves-to-foreign curves)))
                     (with-foreign-out ((sur :surface))
                       (let ((failed (c-flag-error-p
                                      (c-create-skin-spine fc k cmp sp fp align fz subdiv deriv degree
                                                           dir sur sc ss))))
                         (values failed
                                 (unless failed
                                   (list (from-foreign sur :surface)
                                         (%srf-adv-curves-from-foreign fc (1+ k)))))))))))))))

(defun %srf-adv-jitter (c amount)
  "Curve C with its control points moved by up to AMOUNT (z left alone in 2D)."
  (let ((pw (curve-cpoints c)))
    (dotimes (i (array-dimension pw 0) c)
      (dotimes (d 3)
        (unless (= (aref pw i d) +no-z+)
          (incf (aref pw i d) (- (rnd (* 2 amount)) amount)))))))

(defun %srf-adv-gordon-network (k l)
  "A curve network for N_CreateGordonSrf, as (values cur-u v-prms cur-v
u-prms): the iso-curves of a random surface, K+1 u-curves at the ascending
V-PRMS and L+1 v-curves at the ascending U-PRMS, now and then jittered off
the surface (still compatible), the parameters now and then jittered off the
curves' intersections; the surface rational (an input error) now and then."
  (let* ((s (random-surface :max-degree 3 :max-extra 4 :rational (zerop (rnd 15))
                            :two-d (zerop (rnd 10))))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (u-prms (%srf-adv-ascending (aref ku 0) (aref ku (1- (length ku))) l))
         (v-prms (%srf-adv-ascending (aref kv 0) (aref kv (1- (length kv))) k))
         (cur-u (make-array (1+ k))) (cur-v (make-array (1+ l)))
         (jit (case (rnd 4) (0 1d0) (1 10d0) (t 0d0))))
    (batten::with-c-float-environment
      (dotimes (i (1+ k))
        (setf (svref cur-u i) (%srf-adv-jitter (batten::srf-extract-iso-crv s (aref v-prms i) +udir+) jit)))
      (dotimes (j (1+ l))
        (setf (svref cur-v j) (%srf-adv-jitter (batten::srf-extract-iso-crv s (aref u-prms j) +vdir+) jit))))
    (when (zerop (rnd 5))
      (let ((up (copy-seq u-prms)) (vp (copy-seq v-prms)))
        (do ((j 1 (1+ j))) ((>= j l)) (incf (aref up j) (* 0.2d0 (- (rnd 1d0) 0.5d0) (- (aref up (1+ j)) (aref up j)))))
        (do ((i 1 (1+ i))) ((>= i k)) (incf (aref vp i) (* 0.2d0 (- (rnd 1d0) 0.5d0) (- (aref vp (1+ i)) (aref vp i)))))
        (setf u-prms up v-prms vp)))
    (values cur-u v-prms cur-v u-prms)))

(define-oracle-test create-gordon-srf (:file "srf-adv" :trials 1500)
  (let* ((k (1+ (rnd 5))) (l (1+ (rnd 5)))
         ;; degrees in range (at most the highest curve index), now and then not
         (deg (lambda (top) (if (zerop (rnd 40)) (1+ top) (1+ (rnd (min top 3))))))
         (dlu (funcall deg l)) (dlv (funcall deg k))
         (dtu (funcall deg l)) (dtv (funcall deg k))
         (preserve (case (rnd 40) (0 7) (t (if (zerop (rnd 2)) +yes+ +no+)))))
    (multiple-value-bind (cur-u v-prms cur-v u-prms) (%srf-adv-gordon-network k l)
      ;; now and then a u-curve set that is not compatible
      (when (zerop (rnd 30))
        (setf (svref cur-u (rnd (1+ k)))
              (batten::with-c-float-environment (batten::crv-elevate-degree (svref cur-u 0) 1))))
      (agree (list :cur-u cur-u :v-prms v-prms :k k :cur-v cur-v :u-prms u-prms :l l
                   :degs (list dlu dlv dtu dtv) :preserve preserve)
             (let ((cu (map 'simple-vector #'batten::crv-copy cur-u))
                   (cv (map 'simple-vector #'batten::crv-copy cur-v)))
               (list (batten::create-gordon-srf cu v-prms k cv u-prms l dlu dlv dtu dtv preserve)
                     (coerce cu 'list) (coerce cv 'list)))
             (with-foreign ((fv v-prms :real-vector) (fu u-prms :real-vector))
               (with-stacks (sc)
                 (with-stacks (ss)
                   (let ((au (%srf-adv-curves-to-foreign cur-u))
                         (av (%srf-adv-curves-to-foreign cur-v)))
                     (with-foreign-out ((sur :surface))
                       (let ((failed (c-flag-error-p
                                      (c-create-gordon-srf au fv k av fu l dlu dlv dtu dtv preserve
                                                           sur sc ss))))
                         (values failed
                                 (unless failed
                                   (list (from-foreign sur :surface)
                                         (%srf-adv-curves-from-foreign au (1+ k))
                                         (%srf-adv-curves-from-foreign av (1+ l)))))))))))))))

;;; ------------------------------------------------------------- L11-g0

(defun %srf-ext-boundary-offset (s dflg end)
  "A curve for N_SrfExtendToCrv: SURFACE S's boundary row across which DFLG
and END extend, pushed outward by a random multiple of the step from the
next row in (Euclidean, keeping the boundary's weights), now and then
jittered; its knots and degree are the surface's along the boundary."
  (let* ((vd (/= dflg +udir+))
         (net (surface-cpoints s))
         (n (1- (array-dimension net 0))) (m (1- (array-dimension net 1)))
         (cnt (if vd (1+ n) (1+ m)))
         (f (rnd-range 0.1d0 2d0))
         (jit (case (rnd 4) (0 0d0) (1 1d-4) (t (rnd-range 0d0 0.5d0))))
         (cp (make-cpoint-array cnt)))
    (dotimes (k cnt)
      (multiple-value-bind (bi bj ii ij)
          (if vd
              (if (= end +start+) (values k 0 k 1) (values k m k (1- m)))
              (if (= end +start+) (values 0 k 1 k) (values n k (1- n) k)))
        (let* ((wb (aref net bi bj 3)) (wi (aref net ii ij 3))
               (sb (if (= wb +no-w+) 1d0 wb)) (si (if (= wi +no-w+) 1d0 wi)))
          (dotimes (c 3)
            (let* ((b (/ (aref net bi bj c) sb)) (i (/ (aref net ii ij c) si))
                   (x (+ b (* f (- b i)) (* jit (- (rnd 2d0) 1d0)))))
              (setf (aref cp k c) (if (= wb +no-w+) x (* x wb)))))
          (setf (aref cp k 3) wb))))
    (batten::%make-curve (if vd (surface-degree-u s) (surface-degree-v s))
                         (copy-seq (if vd (surface-knots-u s) (surface-knots-v s)))
                         cp)))

(define-oracle-test srf-extend-to-crv (:file "srf-adv" :trials 1500)
  (let* ((s (random-surface :max-degree 4 :max-extra 4))
         (dflg (case (rnd 25) (0 3) (t (if (zerop (rnd 2)) +udir+ +vdir+))))
         (end (if (zerop (rnd 2)) +start+ +end+))
         (cont (if (zerop (rnd 2)) +g-1+ +cmax+))
         (c (case (rnd 6)
              (0 (random-curve :max-degree 4 :max-extra 5 :two-d nil))
              (t (%srf-ext-boundary-offset s dflg end))))
         (in-place (zerop (rnd 5))))
    (when (zerop (rnd 4))                 ; the other rationality
      (setf c (if (batten::is-crv-rat c)
                  (let ((cp (curve-cpoints c)))
                    (dotimes (k (array-dimension cp 0))
                      (let ((w (aref cp k 3)))
                        (dotimes (i 3) (setf (aref cp k i) (/ (aref cp k i) w)))
                        (setf (aref cp k 3) +no-w+)))
                    c)
                  (let ((cp (curve-cpoints c)))
                    (dotimes (k (array-dimension cp 0))
                      (let ((w (rnd-range 0.5d0 2d0)))
                        (dotimes (i 3) (setf (aref cp k i) (* (aref cp k i) w)))
                        (setf (aref cp k 3) w)))
                    c))))
    (agree (list :sur-p s :newbdy c :dflg dflg :end end :cont cont :in-place in-place)
           (list (batten::srf-extend-to-crv s c dflg end cont))
           (with-foreign ((sp s :surface) (cp c :curve))
             (with-stacks (s1)
               (with-stacks (s2)
                 (if in-place
                     (let ((failed (c-flag-error-p (c-srf-extend-to-crv sp cp dflg end cont sp s1 s2))))
                       (values failed (unless failed (list (from-foreign sp :surface)))))
                     (with-foreign-out ((sq :surface))
                       (let ((failed (c-flag-error-p (c-srf-extend-to-crv sp cp dflg end cont sq s1 s2))))
                         (values failed (unless failed (list (from-foreign sq :surface)))))))))))))
