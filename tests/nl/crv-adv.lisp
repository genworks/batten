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
