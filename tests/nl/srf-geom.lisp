;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-geom.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test srf-get-average-len (:file "srf-geom" :trials 3000)
  (let ((s (random-surface :two-d (zerop (rnd 5)))))
    (when (zerop (rnd 6))
      ;; a degenerate row: every control point of row 0 the same (a pole)
      (let ((net (surface-cpoints s)))
        (dotimes (j (array-dimension net 1))
          (dotimes (k 4) (setf (aref net 0 j k) (aref net 0 0 k))))))
    (agree (list :surface s)
           (multiple-value-list (batten::srf-get-average-len s))
           (with-foreign ((sp s :surface))
             (with-scalars ((a :double) (lu :double) (lv :double))
               (values (c-flag-error-p (c-srf-get-average-len sp a lu lv))
                       (list (cffi:mem-ref a :double) (cffi:mem-ref lu :double)
                             (cffi:mem-ref lv :double))))))))
