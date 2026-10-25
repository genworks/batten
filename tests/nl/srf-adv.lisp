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
