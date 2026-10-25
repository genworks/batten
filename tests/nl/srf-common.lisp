;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-common.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test create-srf-corner-pts (:file "srf-common" :trials 1000)
  (let* ((pts (random-points 4))
         (corners (let (cs)
                    (dotimes (i 4 (nreverse cs))
                      (push (make-point (aref pts i 0) (aref pts i 1) (aref pts i 2)) cs)))))
    (destructuring-bind (p00 p10 p01 p11) corners
      (agree (list :corners corners)
             (list (batten::create-srf-corner-pts p00 p10 p01 p11))
             (with-foreign ((a p00 :point) (b p10 :point) (c p01 :point) (d p11 :point))
               (with-stacks (sg)
                 (with-foreign-out ((sur :surface))
                   (let ((failed (c-flag-error-p (c-create-srf-corner-pts a b c d sur sg))))
                     (values failed (unless failed (list (from-foreign sur :surface))))))))))))
