;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-common.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------ level 2, group 5

(defun %cc-arc-input ()
  "Random arc data (list c x y r as ae at): axes orthogonal, nearly so, or
just not parallel, now and then unnormalized or zero; angles anywhere, with
sweeps near the 90/120/180/240/270/360 degree breaks, tiny sweeps and full
circles."
  (let* ((c (random-point))
         (x (random-point (if (zerop (rnd 3)) 1d0 50d0)))
         (y (let ((v (random-point 1d0)))
              ;; Gram-Schmidt against X, fully or partly
              (let* ((xx (+ (* (aref x 0) (aref x 0)) (* (aref x 1) (aref x 1)) (* (aref x 2) (aref x 2))))
                     (dot (+ (* (aref x 0) (aref v 0)) (* (aref x 1) (aref v 1)) (* (aref x 2) (aref v 2))))
                     (f (case (rnd 3) (0 0.7d0) (t 1d0)))
                     (s (* f (/ dot xx))))
                (dotimes (k 3) (setf (aref v k) (- (aref v k) (* s (aref x k)))))
                (when (zerop (rnd 2))
                  (dotimes (k 3) (setf (aref v k) (* 30d0 (aref v k)))))
                v)))
         (r (case (rnd 8) (0 0d0) (1 (- (rnd 10d0))) (t (rnd 100d0))))
         (as (case (rnd 4) (0 0d0) (1 (float (* 90 (- (rnd 9) 4)) 1d0)) (t (- (rnd 1440d0) 720d0))))
         (sweep (case (rnd 12)
                  (0 0d0)
                  (1 360d0)
                  (2 (+ 1d-6 (rnd 1d-4)))
                  (3 (+ 360d0 (- (rnd 0.05d0) 0.025d0)))
                  (4 (+ (float (* 30 (1+ (rnd 12))) 1d0) (- (rnd 2d0) 1d0)))
                  (5 (- (rnd 1440d0) 720d0))
                  (t (rnd 360d0))))
         (ae (+ as sweep))
         (at (case (rnd 4) (0 0.01d0) (1 0d0) (2 1d0) (t (rnd 5d0)))))
    (case (rnd 40)
      (0 (setf x (make-point)))
      (1 (setf y (make-point))))
    (list c x y r as ae at)))

(defmacro %cc-arc-test (lisp-fn c-fn)
  `(destructuring-bind (c x y r as ae at) (%cc-arc-input)
     (agree (list :c c :x x :y y :r r :as as :ae ae :at at)
            (let ((xc (copy-seq x)) (yc (copy-seq y)))
              (list (,lisp-fn c xc yc r as ae at)
                    ;; the axes are the caller's: they come back untouched
                    (equalp xc x) (equalp yc y)))
            (with-stacks (s)
              (with-foreign ((fc c :point) (fx x :point) (fy y :point))
                (with-foreign-out ((cur :curve))
                  (let ((rc (,c-fn fc fx fy r as ae at cur s)))
                    (values (c-flag-error-p rc)
                            (unless (c-flag-error-p rc)
                              (list (from-foreign cur :curve) t t))))))))))

(define-oracle-test create-quadratic-arc-tol (:file "crv-common" :trials 4000)
  (%cc-arc-test batten::create-quadratic-arc-tol c-create-quadratic-arc-tol))

(define-oracle-test create-quartic-arc-tol (:file "crv-common" :trials 4000)
  (%cc-arc-test batten::create-quartic-arc-tol c-create-quartic-arc-tol))
