;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/spiral.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L4-g0

;;; Only valid inputs: on every input error the C jumps to EXIT, which runs
;;; N_EndNurbs on its local stacks before N_InitNurbs has initialized them,
;;; freeing garbage pointers (the test process aborts in free()).  So the
;;; height and the number of turns stay at or above mtol, the radii at or
;;; above -mtol, and one radius at or above mtol.
(define-oracle-test crv-approx-spiral (:file "spiral" :trials 1500)
  (let* ((height (case (rnd 20) (0 1d-7) (1 (rnd-range 1d-7 1d-3)) (t (rnd-range 0.01d0 50d0))))
         (rad (lambda ()
                (case (rnd 12)
                  (0 0d0) (1 5d-8) (2 -5d-8) (3 -1d-7)
                  (t (rnd-range 0.01d0 20d0)))))
         (rad0 (funcall rad))
         (rad1 (if (zerop (rnd 3)) rad0 (funcall rad)))
         (num-turns (case (rnd 20) (0 1d-7) (1 (rnd-range 1d-7 1d-3)) (2 (rnd-range 1d-3 0.2d0))
                      (3 (float (1+ (rnd 4)) 1d0))
                      ;; many turns: the span count reaches *ccplim* - 6,
                      ;; or starts above it
                      (4 (rnd-range 40d0 130d0))
                      (t (rnd-range 0.1d0 6d0))))
         (rad1 (if (< (+ (if (< rad0 1d-7) 0d0 rad0) (if (< rad1 1d-7) 0d0 rad1)) 1d-7)
                   (rnd-range 0.01d0 20d0)
                   rad1))
         (ror-l (case (rnd 5) (0 0) (1 (- (rnd 5) 2)) (t 1)))
         (num-points (case (rnd 6) ((0 1 2) 0) (3 (- (rnd 4) 2)) (t (+ 2 (rnd 60)))))
         (tol (case (rnd 8) (0 0d0) (1 1d-7) (2 (rnd 1d0)) (t (expt 10d0 (- (rnd 6d0)))))))
    (agree (list :height height :rad0 rad0 :rad1 rad1 :num-turns num-turns :ror-l ror-l
                 :num-points num-points :tol tol)
           (multiple-value-list
            (batten::crv-approx-spiral height rad0 rad1 num-turns ror-l num-points tol))
           (with-stacks (s)
             (with-scalars ((np :long) (tl :double))
               (setf (cffi:mem-ref np :long) num-points
                     (cffi:mem-ref tl :double) tol)
               (with-foreign-out ((cur :curve))
                 (let ((failed (c-flag-error-p
                                (c-crv-approx-spiral height rad0 rad1 num-turns ror-l np tl cur s))))
                   (values failed
                           (unless failed
                             (list (cffi:mem-ref np :long) (cffi:mem-ref tl :double)
                                   (from-foreign cur :curve)))))))))))
