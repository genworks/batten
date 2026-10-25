;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-adv.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test crv-derivs-at-knot (:file "crv-adv" :trials 3000)
  (let* ((c (random-curve))
         (knots (curve-knots c))
         ;; now and then a parameter just outside the domain, which both reject
         (u (if (zerop (rnd 25))
                (if (zerop (rnd 2)) (- (aref knots 0) 0.5d0) (+ (aref knots (1- (length knots))) 0.5d0))
                (random-parameter knots)))
         (flg (random-side))
         (der (rnd 5)))
    (agree (list :curve c :u u :flg flg :der der)
           (list (batten::crv-derivs-at-knot c u flg der (make-cpoint-array (1+ der))))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((dw :cpoint-array (1+ der)))
               (values (c-flag-error-p (c-crv-derivs-at-knot cur u flg der dw))
                       (list (from-foreign dw :cpoint-array (1+ der)))))))))

(define-oracle-test crv-get-max-2nd-deriv (:file "crv-adv" :trials 2000)
  (let ((c (random-curve)))
    (agree (list :curve c)
           (list (batten::crv-get-max-2nd-deriv c))
           (with-foreign ((cur c :curve))
             (with-scalars ((muu :double))
               (values (c-flag-error-p (c-crv-get-max-2nd-deriv cur muu))
                       (list (cffi:mem-ref muu :double))))))))

(define-oracle-test crv-non-rat-get-first-deriv (:file "crv-adv" :trials 3000)
  (let* ((c (random-curve :rational (zerop (rnd 10))))
         (p (curve-degree c))
         ;; d in [0, p] mostly; -1 and p+1 now and then for the input check
         (d (case (rnd 20) (0 -1) (1 (1+ p)) (t (rnd (1+ p))))))
    (agree (list :curve c :d d)
           (list (batten::crv-non-rat-get-first-deriv c d))
           (with-foreign ((cur c :curve))
             (with-stacks (sg)
               (with-foreign-out ((der :curve))
                 (let ((failed (c-flag-error-p (c-crv-non-rat-get-first-deriv cur d der sg))))
                   (values failed (unless failed (list (from-foreign der :curve)))))))))))

;;; ------------------------------------------------------------- L2-g7

(defun %crv-adv-planar-curve ()
  "A random planar curve and the normal to its plane, as (values cur nv): 2D,
or 3D in the z = 0 plane; now and then a straight line."
  (let* ((two-d (zerop (rnd 2)))
         (c (random-curve :two-d two-d :max-degree (if (zerop (rnd 15)) 1 5) :max-extra 8)))
    (unless two-d
      (let ((pw (curve-cpoints c)))
        (dotimes (i (array-dimension pw 0)) (setf (aref pw i 2) 0d0))))
    ;; now and then a degenerate curve: coincident control points (zero first
    ;; derivative), which crv-planar-offset-get-deriv rejects
    (when (zerop (rnd 25))
      (let ((pw (curve-cpoints c)))
        (dotimes (i (array-dimension pw 0))
          (dotimes (k 4) (setf (aref pw i k) (aref pw 0 k))))))
    (values c (case (rnd 3)
                (0 (make-point 0d0 0d0 1d0))
                (1 (make-point 0d0 0d0 -1d0))
                (t (random-point 1d0))))))

(define-oracle-test crv-offset-get-max-2nd-deriv (:file "crv-adv" :trials 1500)
  (multiple-value-bind (c nv) (%crv-adv-planar-curve)
    (let ((d (case (rnd 10) (0 0d0) (t (- (rnd 40d0) 20d0)))))
      (agree (list :curve c :d d :n nv)
             (list (batten::crv-offset-get-max-2nd-deriv c d nv))
             (with-foreign ((cur c :curve) (np nv :point))
               (with-scalars ((muu :double))
                 (let ((failed (c-flag-error-p (c-crv-offset-get-max-2nd-deriv cur d np muu))))
                   (values failed (unless failed (list (cffi:mem-ref muu :double)))))))))))

(defun %crv-adv-removable-curve ()
  "A random curve with knots inserted (so some are removable exactly), its
control points now and then perturbed a little (removable within some
tolerance)."
  (let* ((c (random-curve :max-degree 5 :max-extra 6))
         (p (curve-degree c))
         (knots (curve-knots c))
         (lo (aref knots 0))
         (hi (aref knots (1- (length knots)))))
    (when (plusp (rnd 3))
      (let ((xs (sort (let (l)
                        (dotimes (i (1+ (rnd 6)) l)
                          (let ((x (+ lo (* (- hi lo) (+ 0.03d0 (rnd 0.94d0))))))
                            (push x l)
                            (when (and (> p 1) (zerop (rnd 4))) (push x l)))))
                      #'<)))
        (setf c (batten::crv-refine c (coerce xs 'batten::real-vector)))))
    (when (zerop (rnd 3))
      (let ((eps (expt 10d0 (- (rnd 10)))) (pw (curve-cpoints c)))
        (dotimes (i (array-dimension pw 0))
          (dotimes (k 2) (incf (aref pw i k) (* eps (- (rnd 2d0) 1d0)))))))
    c))

(define-oracle-test crv-remove-all-knots (:file "crv-adv" :trials 3000)
  (let ((c (%crv-adv-removable-curve))
        (tol (case (rnd 8) (0 0d0) (1 1d-9) (2 1d-6) (3 1d-3) (4 1d-1) (5 1d0) (6 30d0) (t 1d6))))
    (agree (list :curve c :tol tol)
           (list (batten::crv-remove-all-knots c tol))
           (with-foreign ((cur c :curve))
             (with-stacks (sq)
               (with-foreign-out ((q :curve))
                 (let ((failed (c-flag-error-p (c-crv-remove-all-knots cur tol q sq))))
                   (values failed (unless failed (list (from-foreign q :curve)))))))))))

;;; ------------------------------------------------------------- L3-g6

(define-oracle-test crv-remove-knots-deriv-constraints (:file "crv-adv" :trials 3000)
  (let* ((c (%crv-adv-removable-curve))
         (tol (case (rnd 8) (0 0d0) (1 1d-9) (2 1d-6) (3 1d-3) (4 1d-1) (5 1d0) (6 30d0) (t 1d6)))
         (whr (case (rnd 4) (0 +no+) (1 +start+) (2 +end+) (t +both+)))
         (der (rnd (+ (curve-degree c) 2))))
    (agree (list :curve c :tol tol :whr whr :der der)
           (list (batten::crv-remove-knots-deriv-constraints c tol whr der))
           (with-foreign ((cur c :curve))
             (with-stacks (sq)
               (with-foreign-out ((q :curve))
                 (let ((failed (c-flag-error-p
                                (c-crv-remove-knots-deriv-constraints cur tol whr der q sq))))
                   (values failed (unless failed (list (from-foreign q :curve)))))))))))

(define-oracle-test crv-planar-offset-get-deriv (:file "crv-adv" :trials 8000)
  (let* ((two-d (zerop (rnd 2)))
         ;; now and then a straight line (zero curvature) or a degenerate curve
         ;; whose control points coincide (zero first derivative)
         (c (case (rnd 20)
              (0 (random-curve :two-d two-d :max-degree 1))
              (1 (let* ((c (random-curve :two-d two-d :rational nil))
                        (pw (curve-cpoints c)))
                   (dotimes (i (array-dimension pw 0) c)
                     (dotimes (k 3) (setf (aref pw i k) (aref pw 0 k))))))
              (t (random-curve :two-d two-d))))
         (u (random-parameter (curve-knots c)))
         (d (case (rnd 10) (0 0d0) (t (- (rnd 40d0) 20d0))))
         (nv (case (rnd 3)
               (0 (make-point 0d0 0d0 1d0))
               (1 (make-point 0d0 0d0 -1d0))
               (t (random-point 1d0))))
         (ufl (random-side)))
    (agree (list :curve c :u u :d d :n nv :ufl ufl)
           (list (batten::crv-planar-offset-get-deriv c u d nv ufl (make-point-array 3)))
           (with-foreign ((cur c :curve) (np nv :point))
             (with-foreign-out ((od :point-array 3))
               (let ((failed (c-flag-error-p (c-crv-planar-offset-get-deriv cur u d np ufl od))))
                 (values failed (unless failed (list (from-foreign od :point-array 3))))))))))

;;; ------------------------------------------------------------- L4-g0

(defun %crv-adv-curve-on (lo hi &rest keys)
  "A random curve (RANDOM-CURVE with KEYS) on the domain [LO, HI]."
  (let ((c (apply #'random-curve keys)))
    (setf (curve-knots c)
          (random-knots (curve-degree c) (array-dimension (curve-cpoints c) 0) :lo lo :hi hi))
    c))

(defun %crv-adv-curve-pair ()
  "Two random curves for N_CrvSumDiffCrv, as (values p q): on one domain
(now and then two), rational or not each, 2D or 3D together; now and then Q
with P's degree and knots, or P's knots refined."
  (let* ((lo (if (zerop (rnd 3)) 0d0 (rnd-range -10d0 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0)))
         (two-d (zerop (rnd 5)))
         (rp (case (rnd 3) (0 t) (1 nil) (t :random)))
         (rq (case (rnd 3) (0 t) (1 nil) (t :random)))
         (p (%crv-adv-curve-on lo hi :two-d two-d :rational rp :max-degree 6 :max-extra 8))
         (q (case (rnd 8)
              (0 ;; same degree and knots
               (let ((c (%crv-adv-curve-on lo hi :two-d two-d :rational rq)))
                 (%make-curve (curve-degree p) (copy-seq (curve-knots p))
                              (let ((pw (make-cpoint-array (array-dimension (curve-cpoints p) 0)))
                                    (cw (curve-cpoints c)))
                                (dotimes (i (array-dimension pw 0) pw)
                                  (dotimes (k 4)
                                    (setf (aref pw i k) (aref cw (mod i (array-dimension cw 0)) k))))))))
              (1 ;; P refined
               (let ((xs (sort (let (l) (dotimes (i (1+ (rnd 4)) l)
                                          (push (+ lo (* (- hi lo) (+ 0.05d0 (rnd 0.9d0)))) l)))
                               #'<)))
                 (batten::crv-refine (batten::crv-copy p) (coerce xs 'batten::real-vector))))
              (2 ;; another domain: the products and the merge rescale it
               (%crv-adv-curve-on (+ lo (rnd 1d0)) (+ hi (rnd 2d0)) :two-d two-d :rational rq
                                  :max-degree 6 :max-extra 8))
              (t (%crv-adv-curve-on lo hi :two-d two-d :rational rq :max-degree 6 :max-extra 8)))))
    (if (zerop (rnd 2)) (values p q) (values q p))))

(define-oracle-test crv-sum-diff-crv (:file "crv-adv" :trials 3000)
  (multiple-value-bind (cp cq) (%crv-adv-curve-pair)
    (let ((opr (case (rnd 25) (0 0) (1 3) (t (if (zerop (rnd 2)) +plus+ +minus+)))))
      (agree (list :p cp :q cq :opr opr)
             (let ((pc (batten::crv-copy cp)) (qc (batten::crv-copy cq)))
               ;; the inputs are returned to their rationality
               (list (batten::crv-sum-diff-crv pc qc opr) pc qc))
             (with-foreign ((p cp :curve) (q cq :curve))
               (with-stacks (sg)
                 (with-foreign-out ((r :curve))
                   (let ((failed (c-flag-error-p (c-crv-sum-diff-crv p q opr r sg))))
                     (values failed (unless failed
                                      (list (from-foreign r :curve)
                                            (from-foreign p :curve) (from-foreign q :curve))))))))))))

(defun %crv-adv-reparam-func (a b)
  "A random increasing B-spline function from A to B (f' > 0 for strictly
increasing control values): degree 1..4 on a random domain, now and then a
single linear segment, or control values with flat runs."
  (let* ((lo (if (zerop (rnd 3)) 0d0 (rnd-range -10d0 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0)))
         (q (if (zerop (rnd 8)) 1 (1+ (rnd 4))))
         (n+1 (if (and (= q 1) (zerop (rnd 2))) 2 (+ q 1 (rnd 7))))
         (vals (make-real-vector n+1))
         (inner (sort (let (l) (dotimes (i (- n+1 2) l) (push (+ a (* (- b a) (+ 0.01d0 (rnd 0.98d0)))) l)))
                      #'<)))
    (setf (aref vals 0) a
          (aref vals (1- n+1)) b)
    (let ((i 1))
      (dolist (v inner) (setf (aref vals i) v) (incf i)))
    ;; now and then two equal control values in the middle
    (when (and (> n+1 3) (zerop (rnd 10)))
      (setf (aref vals 2) (aref vals 1)))
    (%make-cfun q (random-knots q n+1 :lo lo :hi hi) vals)))

(define-oracle-test srf-reparam-func (:file "crv-adv" :trials 3000)
  (let* ((c (random-curve :max-degree (if (zerop (rnd 10)) 7 4) :max-extra 8))
         (knots (curve-knots c))
         (a (aref knots 0))
         (b (aref knots (1- (length knots))))
         (f (%crv-adv-reparam-func a b)))
    ;; now and then an end value off the curve's end (the mapping rejects it)
    (when (zerop (rnd 25))
      (incf (aref (cfun-vals f) 0) (* (- b a) (if (zerop (rnd 2)) 1d-2 1d-6))))
    (agree (list :curve c :cfn f)
           (list (batten::srf-reparam-func c f))
           (with-foreign ((cp c :curve) (fp f :cfun))
             (with-stacks (sp)
               (with-stacks (sq)
                 (with-foreign-out ((q :curve))
                   (let ((failed (c-flag-error-p (c-srf-reparam-func cp fp q sp sq))))
                     (values failed (unless failed (list (from-foreign q :curve))))))))))))
