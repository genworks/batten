;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-nurb.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test crv-eval (:file "crv-nurb" :trials 4000)
  (let* ((c (random-curve)) (u (random-parameter (curve-knots c))) (flg (random-side)))
    (agree (list :curve c :u u :flg flg)
           (list (batten::crv-eval c u flg (make-point)))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((pt :point))
               (values (c-flag-error-p (c-crv-eval cur u flg pt))
                       (list (from-foreign pt :point))))))))

(define-oracle-test crv-derivs (:file "crv-nurb" :trials 4000)
  (let* ((c (random-curve)) (u (random-parameter (curve-knots c))) (flg (random-side))
         (der (rnd (1+ +maxder+))))
    (agree (list :curve c :u u :flg flg :der der)
           (list (batten::crv-derivs c u flg der (make-point-array (1+ der))))
           (with-foreign ((cur c :curve))
             (with-foreign-out ((cd :point-array (1+ der)))
               (values (c-flag-error-p (c-crv-derivs cur u flg der cd))
                       (list (from-foreign cd :point-array (1+ der)))))))))
