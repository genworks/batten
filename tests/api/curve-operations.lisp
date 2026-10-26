;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Tests of the API in src/api/curve-operations.lisp: each function against the C routines
;;;; it wraps (or the port routines, where the C has no single counterpart).

(in-package :batten-tests)

;;; The inputs come from the generators of the port's own tests
;;; (tests/nl/crv-tool.lisp, crv-adv.lisp, crv-offset.lisp, crv-approx.lisp,
;;; crv-nurb.lisp), which keep them inside the C's domain.

(defun %api-co-deep-copy (x)
  "A copy of X sharing nothing with it: curves, arrays (simple-vectors of
curves too) and lists copied all the way down."
  (typecase x
    (curve (batten::crv-copy x))
    (cons (mapcar #'%api-co-deep-copy x))
    (string x)
    (array (let ((a (make-array (array-dimensions x) :element-type (array-element-type x))))
             (dotimes (i (array-total-size x) a)
               (setf (row-major-aref a i) (%api-co-deep-copy (row-major-aref x i))))))
    (t x)))

(defmacro %api-co-agree (case inputs lisp-form c-form)
  "AGREE, with one more result on each side: on the Lisp side :unchanged when
every object in INPUTS (a list of forms) is bit for bit what it was before
LISP-FORM ran, else :changed; on the C side :unchanged."
  (let ((ins (gensym "INPUTS")) (snapshot (gensym "SNAPSHOT"))
        (failed (gensym "FAILED")) (results (gensym "RESULTS")))
    `(let* ((,ins (list ,@inputs))
            (,snapshot (%api-co-deep-copy ,ins)))
       (agree ,case
              (append ,lisp-form
                      (list (if (values (same ,ins ,snapshot)) :unchanged :changed)))
              (multiple-value-bind (,failed ,results) ,c-form
                (values ,failed (unless ,failed (append ,results (list :unchanged)))))))))

(defun %api-co-reals (x)
  "The real-vector X as one of the sequences the API takes: itself, a list,
or a list in reverse order."
  (case (rnd 3)
    (0 x)
    (1 (coerce x 'list))
    (t (reverse (coerce x 'list)))))

(defun %api-co-point (p)
  "The point P as itself or as a list of its coordinates."
  (if (zerop (rnd 2)) p (coerce p 'list)))

;;; ------------------------------------------------------------ knots

(define-oracle-test api-insert-knots-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (if (zerop (rnd 2)) (random-curve) (%crv-tool-curve)))
         (x (%crv-tool-refinement c))
         (knots (%api-co-reals x)))
    (%api-co-agree (list :curve c :knots knots) (c)
                   (list (insert-knots c knots))
                   (if (zerop (length x))
                       ;; no knots: a copy (the C's N_CrvRefine refuses an empty vector)
                       (values nil (list (batten::crv-copy c)))
                       (with-foreign ((cp c :curve) (kx x :knots))
                         (with-foreign-out ((cq :curve))
                           (with-stacks (sp)
                             (with-stacks (sq)
                               (if (c-flag-error-p (c-crv-refine cp kx cq sp sq))
                                   (values t nil)
                                   (values nil (list (from-foreign cq :curve))))))))))))

(define-oracle-test api-remove-knots (:file "api-curve-operations" :trials 600)
  (multiple-value-bind (c tol) (%crv-tool-removals-input)
    (%api-co-agree (list :curve c :tol tol) (c)
                   (list (batten::remove-knots c :tolerance tol))
                   (with-foreign ((cp c :curve))
                     (with-foreign-out ((cq :curve))
                       (with-stacks (sq)
                         (if (c-flag-error-p (c-crv-remove-knots cp tol cq sq))
                             (values t nil)
                             (values nil (list (from-foreign cq :curve))))))))))

;;; ------------------------------------------------------------ degree

(define-oracle-test api-elevate-degree-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (if (zerop (rnd 2))
                (random-curve :max-extra 8)
                (%crv-tool-curve :full (+ 1 (rnd 5)))))
         (p (curve-degree c))
         (by (case (rnd 20)
               (0 -1)
               (1 (- 18 p))
               (2 0)
               (t (1+ (rnd 4))))))
    (when (zerop (rnd 25))
      ;; a knot vector that is not valid: start multiplicity p
      (let ((u (curve-knots c)))
        (setf (aref u p) (* 0.5d0 (+ (aref u p) (aref u (1+ p)))))))
    (%api-co-agree (list :curve c :by by) (c)
                   (list (elevate-degree c by))
                   (with-foreign ((cp c :curve))
                     (with-foreign-out ((cq :curve))
                       (with-stacks (sp)
                         (with-stacks (sq)
                           (if (c-flag-error-p (c-crv-elevate-degree cp by cq sp sq))
                               (values t nil)
                               (values nil (list (from-foreign cq :curve)))))))))))

;;; ------------------------------------------------------------ pieces

(defun %api-co-extract (c ul ur)
  "N_CrvExtractCrvSeg of curve C over [UL, UR]: (values failed curve)."
  (with-foreign ((cp c :curve))
    (with-foreign-out ((cq :curve))
      (with-stacks (sp)
        (with-stacks (sq)
          (if (c-flag-error-p (c-crv-extract-crv-seg cp ul ur cq sp sq))
              (values t nil)
              (values nil (from-foreign cq :curve))))))))

(define-oracle-test api-subcurve (:file "api-curve-operations" :trials 1000)
  (let* ((c (if (zerop (rnd 2)) (random-curve) (%crv-tool-curve)))
         (knots (curve-knots c))
         (runs (%crv-tool-knot-runs knots (curve-degree c)))
         (a (%crv-tool-segment-end knots runs))
         (b (%crv-tool-segment-end knots runs))
         ;; mostly ordered; equal or reversed now and then, for the error path
         (u0 (if (zerop (rnd 25)) (max a b) (min a b)))
         (u1 (if (zerop (rnd 25)) (min a b) (max a b))))
    (%api-co-agree (list :curve c :u0 u0 :u1 u1) (c)
                   (list (batten::subcurve c u0 u1))
                   (multiple-value-bind (failed q) (%api-co-extract c u0 u1)
                     (values failed (list q))))))

(define-oracle-test api-split-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (if (zerop (rnd 2)) (random-curve) (%crv-tool-curve)))
         (knots (curve-knots c))
         (u0 (aref knots 0))
         (u1 (aref knots (1- (length knots))))
         ;; an end now and then, for the error path
         (at (%crv-tool-segment-end knots (%crv-tool-knot-runs knots (curve-degree c)))))
    (%api-co-agree (list :curve c :at at) (c)
                   (multiple-value-list (split c at))
                   (multiple-value-bind (failed-0 q0) (%api-co-extract c u0 at)
                     (if failed-0
                         (values t nil)
                         (multiple-value-bind (failed-1 q1) (%api-co-extract c at u1)
                           (values failed-1 (list q0 q1))))))))

(define-oracle-test api-bezier-segments (:file "api-curve-operations" :trials 1000)
  (let ((c (if (zerop (rnd 2)) (%crv-tool-curve :full (+ 2 (rnd 4))) (random-curve))))
    (%api-co-agree (list :curve c) (c)
                   (list (batten::bezier-segments c))
                   (with-foreign ((cp c :curve))
                     (with-stacks (sq)
                       (with-scalars ((k :long) (cq :pointer))
                         (if (c-flag-error-p (c-crv-decompose-bez cp cq k sq))
                             (values t nil)
                             (let ((kk (cffi:mem-ref k :long))
                                   (arr (cffi:mem-ref cq :pointer))
                                   (curves nil))
                               (dotimes (i (1+ kk))
                                 (push (from-foreign (cffi:mem-aref arr :pointer i) :curve) curves))
                               (values nil (list (nreverse curves)))))))))))

;;; ------------------------------------------------------------ direction and domain

(define-oracle-test api-reverse-curve (:file "api-curve-operations" :trials 1000)
  (let ((c (random-curve :max-extra (if (zerop (rnd 3)) 2 12))))
    (%api-co-agree (list :curve c) (c)
                   (list (batten::reverse-curve c))
                   (with-stacks (s)
                     (with-foreign ((cur c :curve))
                       (with-foreign-out ((q :curve))
                         (let ((rc (c-crv-reverse cur q s)))
                           (values (c-flag-error-p rc)
                                   (unless (c-flag-error-p rc) (list (from-foreign q :curve)))))))))))

(define-oracle-test api-reparameterize-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (random-curve))
         (k (curve-knots c))
         (ul (if (zerop (rnd 4)) (aref k 0) (- (rnd 20d0) 10d0)))
         (ur (if (zerop (rnd 4)) (aref k (1- (length k))) (+ ul 0.1d0 (rnd 10d0))))
         ;; integers now and then: any real is a parameter
         (ul* (if (zerop (rnd 6)) (round ul) ul))
         (ur* (if (zerop (rnd 6)) (+ (round ul) 1 (rnd 5)) ur))
         (ul (float ul* 1d0))
         (ur (float ur* 1d0)))
    (%api-co-agree (list :curve c :u0 ul* :u1 ur*) (c)
                   (list (reparameterize c ul* ur*))
                   (with-foreign ((cur c :curve))
                     (c-crv-reparam cur ul ur)
                     (values nil (list (from-foreign cur :curve)))))))

;;; ------------------------------------------------------------ extension and offset

(define-oracle-test api-extend-curve (:file "api-curve-operations" :trials 300)
  (let* ((c (%crv-adv-extend-curve))
         (relative (zerop (rnd 2)))
         (dist (if relative
                   (case (rnd 30)
                     (0 30d0)                   ; past the G2 recursion limit
                     (1 (rnd-range 1d0 3d0))
                     (t (rnd-range 0.01d0 1.2d0)))
                   (case (rnd 10)
                     (0 (rnd-range 300d0 900d0))
                     (t (rnd-range 0.5d0 150d0)))))
         (end (if (zerop (rnd 2)) :start :end))
         (continuity (case (rnd 3) (0 :g1) (1 :g2) (t :extrapolate))))
    (%api-co-agree (list :curve c :dist dist :relative relative :end end :continuity continuity) (c)
                   (list (batten::extend-curve c dist :end end :continuity continuity :relative relative))
                   (with-foreign ((cp c :curve))
                     (with-stacks (sp)
                       (with-stacks (sq)
                         (with-foreign-out ((q :curve))
                           (let ((failed (c-flag-error-p
                                          (c-crv-extend-by-dist cp dist (if relative +relative+ +absolute+)
                                                                (if (eq end :start) +start+ +end+)
                                                                (ecase continuity
                                                                  (:g1 +g-1+) (:g2 +g-2+) (:extrapolate +cmax+))
                                                                q sp sq))))
                             (values failed (unless failed (list (from-foreign q :curve))))))))))))

(define-oracle-test api-offset-curve (:file "api-curve-operations" :trials 400)
  (multiple-value-bind (c nrm)
      (case (rnd 10)
        (0 (values (%crv-geom-line-curve) (random-point 1d0)))
        (1 (let ((c (%crv-geom-circle-curve))) (values c (random-point 1d0))))
        (t (%co6-planar-curve)))
    (let* ((b (case (rnd 20)
                (0 (make-point 0d0 0d0 0d0))
                (t (%co6-scaled nrm (if (zerop (rnd 2)) (rnd-range 0.5d0 3d0) (- (rnd-range 0.5d0 3d0)))))))
           (normal (%api-co-point b))
           (d (case (rnd 6) (0 0d0) (t (* (if (zerop (rnd 2)) 1d0 -1d0) (rnd-range 0.05d0 2d0)))))
           (deg (if (zerop (rnd 12)) 1 (+ 2 (rnd 4))))
           (keep (zerop (rnd 2)))
           (par (case (rnd 4) (0 :uniform) (1 :chord-length) (2 :centripetal) (t :inherited)))
           (tol (case (rnd 10) (0 0d0) ((1 2 3) 0.05d0) ((4 5) 0.2d0) (6 1d-2) (t 0.1d0))))
      (%api-co-agree (list :curve c :normal normal :d d :deg deg :keep keep :par par :tol tol)
                     (c normal)
                     (list (batten::offset-curve c d :normal normal :tolerance tol :degree deg
                                                     :parameterization par :keep-end-tangents keep))
                     (with-foreign ((cp c :curve) (bp b :point))
                       (with-stacks (sg)
                         (with-foreign-out ((cur :curve))
                           (let ((err (c-flag-error-p
                                       (c-crv-offset-pt-sampling
                                        cp bp d deg (if keep batten::+yes+ batten::+no+)
                                        (ecase par
                                          (:uniform +uniform+) (:chord-length batten::+chordlength+)
                                          (:centripetal batten::+centripetal+) (:inherited +inherited+))
                                        tol cur sg))))
                             (values err (unless err (list (from-foreign cur :curve))))))))))))

;;; ------------------------------------------------------------ compatibility

(define-oracle-test api-make-compatible (:file "api-curve-operations" :trials 600)
  (multiple-value-bind (curves k) (%crv-tool-curve-set)
    (let ((tol (case (rnd 9)
                 ((0 1 2) nil) (3 0d0) (4 1d-9) (5 1d-6) (6 1d-4) (7 1d-3) (t 1d-2)))
          (arg (if (zerop (rnd 2)) curves (coerce curves 'list))))
      (%api-co-agree (list :curves curves :k k :tol tol) (curves)
                     (list (batten::make-compatible arg :knot-tolerance tol))
                     (with-foreign ()
                       (with-stacks (s)
                         (let ((arr (batten-oracle::falloc :pointer (1+ k))))
                           (dotimes (i (1+ k))
                             (setf (cffi:mem-aref arr :pointer i)
                                   (batten-oracle::to-foreign (svref curves i) :curve)))
                           (let ((failed (c-flag-error-p
                                          (if tol
                                              (c-crvs-make-compatible-knot-tol arr k tol s)
                                              (c-crvs-make-compatible arr k s)))))
                             (values failed
                                     (unless failed
                                       (list (let (out)
                                               (dotimes (i (1+ k) (nreverse out))
                                                 (push (from-foreign (cffi:mem-aref arr :pointer i) :curve)
                                                       out))))))))))))))

;;; ------------------------------------------------------------ moving

(define-oracle-test api-translate-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (random-curve))
         (tv (random-point (if (zerop (rnd 4)) 1d-3 100d0)))
         (vector (%api-co-point tv)))
    (%api-co-agree (list :curve c :vector vector) (c vector)
                   (list (translate c vector))
                   (with-foreign ((cur c :curve) (tp tv :point))
                     (c-crv-translate cur tp)
                     (values nil (list (from-foreign cur :curve)))))))

(define-oracle-test api-transform-curve (:file "api-curve-operations" :trials 1000)
  (let* ((c (random-curve))
         (rma (%cn-random-transform))
         (rm (batten::rmatrix-rm rma))
         (matrix (if (zerop (rnd 2))
                     rm
                     (let (rows)
                       (dotimes (i 4 (nreverse rows))
                         (push (list (aref rm i 0) (aref rm i 1) (aref rm i 2) (aref rm i 3)) rows))))))
    (%api-co-agree (list :curve c :matrix matrix) (c matrix)
                   (list (transform c matrix))
                   (with-foreign ((cur c :curve) (rp rma :rmatrix))
                     (c-crv-transform cur rp)
                     (values nil (list (from-foreign cur :curve)))))))

;;; ------------------------------------------------------------ 2D, 3D, rational, polynomial

(define-oracle-test api-curve-to-2d (:file "api-curve-operations" :trials 500)
  (let ((c (random-curve)))
    (%api-co-agree (list :curve c) (c)
                   (list (batten::curve-to-2d c))
                   (with-foreign ((cp c :curve))
                     (c-crv-3d-to-2d cp)
                     (values nil (list (from-foreign cp :curve)))))))

(define-oracle-test api-curve-to-3d (:file "api-curve-operations" :trials 500)
  (let ((c (random-curve :two-d (zerop (rnd 2)))))
    (%api-co-agree (list :curve c) (c)
                   (list (batten::curve-to-3d c))
                   (with-foreign ((cp c :curve))
                     ;; a 3D curve comes back as it is; N_Crv2dTo3d would set z to 0
                     (unless (batten::crv-is-3d c)
                       (c-crv-2d-to-3d cp))
                     (values nil (list (from-foreign cp :curve)))))))

(define-oracle-test api-curve-to-rational (:file "api-curve-operations" :trials 500)
  (let ((c (random-curve)))
    (%api-co-agree (list :curve c) (c)
                   (list (batten::curve-to-rational c))
                   (with-foreign ((cp c :curve))
                     ;; a rational curve comes back as it is; N_CrvNonRatToRat
                     ;; would set its weights to 1
                     (unless (batten::is-crv-rat c)
                       (c-crv-non-rat-to-rat cp))
                     (values nil (list (from-foreign cp :curve)))))))

(define-oracle-test api-curve-to-polynomial (:file "api-curve-operations" :trials 500)
  ;; TOL positive (a negative one makes pow() a NaN the C casts to an index),
  ;; or now and then 0, which the division check rejects
  (let* ((c (%crv-approx-curve))
         (tol (case (rnd 50)
                (0 0d0)
                (t (case (rnd 6) (0 1d-3) (1 1d-2) (2 0.1d0) (3 1d0) (4 5d0) (t (rnd-range 1d-2 1d0))))))
         (deg (case (rnd 10) (0 1) (1 2) (2 2) (3 3) (4 3) (5 4) (6 5) (t (1+ (rnd 3)))))
         (par (case (rnd 4) (0 :uniform) (1 :chord-length) (2 :centripetal) (t :inherited)))
         (keep (if (zerop (rnd 3)) :default (zerop (rnd 2))))
         (tn (if (if (eq keep :default) (/= deg 1) keep) +yes+ +no+)))
    (%api-co-agree (list :curve c :tol tol :deg deg :par par :keep keep) (c)
                   (list (if (eq keep :default)
                             (batten::curve-to-polynomial c :tolerance tol :degree deg :parameterization par)
                             (batten::curve-to-polynomial c :tolerance tol :degree deg :parameterization par
                                                            :keep-end-tangents keep)))
                   (if (batten::is-crv-rat c)
                       (with-foreign ((cp c :curve))
                         (with-stacks (sg)
                           (with-foreign-out ((q :curve))
                             (let ((failed (c-flag-error-p
                                            (c-approx-nurbs-with-non-rat-crv
                                             cp tol deg
                                             (ecase par
                                               (:uniform +uniform+) (:chord-length batten::+chordlength+)
                                               (:centripetal batten::+centripetal+) (:inherited +inherited+))
                                             tn q sg))))
                               (values failed (unless failed (list (from-foreign q :curve))))))))
                       ;; a polynomial curve comes back as it is
                       (values nil (list (batten::crv-copy c)))))))
