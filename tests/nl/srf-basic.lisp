;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-basic.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test is-srf-rat (:file "srf-basic" :trials 1000)
  (let ((s (random-surface)))
    (agree s
           (list (if (batten::is-srf-rat s) 1 0))
           (with-foreign ((sp s :surface))
             (values nil (list (c-is-srf-rat sp)))))))

(define-oracle-test srf-get-parameter-bounds (:file "srf-basic" :trials 1000)
  (let ((s (random-surface)))
    (agree s
           (multiple-value-list (batten::srf-get-parameter-bounds s))
           (with-foreign ((sp s :surface))
             (with-scalars ((ul :double) (ur :double) (vb :double) (vt :double))
               (c-srf-get-parameter-bounds sp ul ur vb vt)
               (values nil (list (cffi:mem-ref ul :double) (cffi:mem-ref ur :double)
                                 (cffi:mem-ref vb :double) (cffi:mem-ref vt :double))))))))
