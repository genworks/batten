;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-basic.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test is-crv-rat (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (list (if (batten::is-crv-rat c) 1 0))
           (with-foreign ((cp c :curve))
             (values nil (list (c-is-crv-rat cp)))))))

(define-oracle-test crv-get-param-bounds (:file "crv-basic" :trials 1000)
  (let ((c (random-curve)))
    (agree c
           (multiple-value-list (batten::crv-get-param-bounds c))
           (with-foreign ((cp c :curve))
             (with-scalars ((ul :double) (ur :double))
               (c-crv-get-param-bounds cp ul ur)
               (values nil (list (cffi:mem-ref ul :double) (cffi:mem-ref ur :double))))))))
