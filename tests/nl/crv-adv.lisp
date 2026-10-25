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
