;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-geom.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------ level 2, group 6

(defun %crv-geom-close (c)
  "Make curve C closed by copying its first control point onto its last."
  (let ((pw (curve-cpoints c)) (n (curve-n c)))
    (dotimes (k 4) (setf (aref pw n k) (aref pw 0 k)))
    c))

(defun %crv-geom-offset (p scale)
  "A fresh point P plus a random vector of size up to SCALE."
  (let ((d (random-point scale)))
    (make-point (+ (aref p 0) (aref d 0)) (+ (aref p 1) (aref d 1)) (+ (aref p 2) (aref d 2)))))

;;; The C returns the best u and Q found even when it reports an error, so the
;;; test compares them in every case, with the error flag beside them.
(define-oracle-test crv-closest-pt (:file "crv-geom" :trials 3000)
  (let* ((c (random-curve :max-degree 5 :max-extra 8))
         (c (if (zerop (rnd 3)) (%crv-geom-close c) c))
         (kn (curve-knots c))
         (lo (aref kn 0)) (hi (aref kn (1- (length kn))))
         (us (random-parameter kn))
         (on (batten::crv-eval c us +left+ (make-point)))
         (p (case (rnd 5)
              (0 on)
              (1 (%crv-geom-offset on 1d0))
              (2 (%crv-geom-offset on 10d0))
              (3 (%crv-geom-offset on 100d0))
              (t (random-point))))
         (u0 (case (rnd 12)
               (0 (+ hi (rnd-range 0d0 1d0)))  ; outside: the derivatives fail
               ((1 2 3 4 5)
                (min hi (max lo (+ us (* 0.1d0 (- (rnd 2d0) 1d0) (- hi lo))))))
               (t (random-parameter kn))))
         (top (case (rnd 4) (0 0d0) (1 1d-7) (2 1d-4) (t 1d-2)))
         (toc (case (rnd 4) (0 0d0) (1 1d-6) (2 1d-3) (t 1d-9))))
    (agree (list :curve c :p p :u0 u0 :top top :toc toc)
           (let* ((q (make-point)) (ucell (make-real-vector 1)) (err nil)
                  (u (handler-case (batten::crv-closest-pt c p u0 top toc q ucell)
                       (nurbs-error () (setf err t) (aref ucell 0)))))
             (list err u (aref ucell 0) q))
           (with-foreign ((cp c :curve) (pp p :point))
             (with-foreign-out ((q :point))
               (with-scalars ((u :double))
                 (let ((err (c-flag-error-p (c-crv-closest-pt cp pp u0 top toc u q))))
                   (values nil
                           (list err (cffi:mem-ref u :double) (cffi:mem-ref u :double)
                                 (from-foreign q :point))))))))))
