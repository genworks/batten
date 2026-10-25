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

;;; ------------------------------------------------------------ level 1, group 4

(defun %close-surface (s dir)
  "Make surface S closed in DIR by copying its first row (+udir+) or column
(+vdir+) of control points onto the last."
  (let* ((net (surface-cpoints s)) (n (surface-n s)) (m (surface-m s)))
    (if (= dir +udir+)
        (dotimes (j (1+ m)) (dotimes (k 4) (setf (aref net n j k) (aref net 0 j k))))
        (dotimes (i (1+ n)) (dotimes (k 4) (setf (aref net i m k) (aref net i 0 k)))))
    s))

(define-oracle-test get-closest-pt-on-srf (:file "srf-geom" :trials 3000)
  (let* ((s (random-surface :max-degree 4 :max-extra 5 :two-d (zerop (rnd 6))))
         (s (progn (when (zerop (rnd 4)) (%close-surface s +udir+))
                   (when (zerop (rnd 4)) (%close-surface s +vdir+))
                   s))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (us (random-parameter ku)) (vs (random-parameter kv))
         (on (batten::srf-eval-pt s us vs +left+ +left+ (make-point)))
         (p (case (rnd 4)
              (0 on)
              (1 (%pt+ on 1d0 (random-point 5d0)))
              (2 (%pt+ on 1d0 (random-point 50d0)))
              (t (random-point))))
         (near (zerop (rnd 2)))
         (u0 (if near
                 (min (aref ku (1- (length ku)))
                      (max (aref ku 0) (+ us (* 0.05d0 (- (rnd 2d0) 1d0) (- (aref ku (1- (length ku))) (aref ku 0))))))
                 (random-parameter ku)))
         (v0 (if near
                 (min (aref kv (1- (length kv)))
                      (max (aref kv 0) (+ vs (* 0.05d0 (- (rnd 2d0) 1d0) (- (aref kv (1- (length kv))) (aref kv 0))))))
                 (random-parameter kv)))
         (top (case (rnd 3) (0 1d-7) (1 1d-4) (t 1d-2)))
         (toc (case (rnd 3) (0 1d-6) (1 1d-3) (t 1d-9))))
    (agree (list :surface s :p p :u0 u0 :v0 v0 :top top :toc toc)
           (let ((q (make-point)) (d (make-point-net 3 3)))
             (multiple-value-bind (u v) (batten::get-closest-pt-on-srf s p u0 v0 top toc q d)
               (list u v q d)))
           (with-foreign ((sp s :surface) (pp p :point))
             (with-foreign-out ((q :point) (d :point-net 3 3))
               (with-scalars ((u :double) (v :double))
                 (values (c-flag-error-p (c-get-closest-pt-on-srf sp pp u0 v0 top toc u v q d))
                         (list (cffi:mem-ref u :double) (cffi:mem-ref v :double)
                               (from-foreign q :point) (from-foreign d :point-net 3 3)))))))))

;;; ------------------------------------------------------------- L4-g1

(define-oracle-test approx-srf-with-quad-srf (:file "srf-geom" :trials 400)
  (let* ((s (random-surface :max-degree 3 :max-extra 4 :two-d (zerop (rnd 8))))
         (s (progn (when (zerop (rnd 8)) (%close-surface s +udir+))
                   (when (zerop (rnd 8)) (%close-surface s +vdir+))
                   s))
         (tfl (if (zerop (rnd 2)) +relative+ +absolute+))
         (tol (if (= tfl +relative+)
                  (case (rnd 4) (0 0.2d0) (1 0.05d0) (2 0.02d0) (t 0.01d0))
                  (case (rnd 4) (0 50d0) (1 20d0) (2 5d0) (t 2d0))))
         (rfl (case (rnd 7) (0 0) ((1 2) +points+) ((3 4) +parameters+) (t +both+))))
    (agree (list :surface s :tol tol :tfl tfl :rfl rfl)
           (multiple-value-list (batten::approx-srf-with-quad-srf s tol tfl rfl))
           (with-foreign ((sp s :surface))
             (with-stacks (sg)
               (with-scalars ((p :pointer) (u :pointer) (v :pointer) (n :long) (m :long))
                 (let ((failed (c-flag-error-p (c-approx-srf-with-quad-srf sp tol tfl rfl p u v n m sg))))
                   (values failed
                           (unless failed
                             (let ((n (cffi:mem-ref n :long)) (m (cffi:mem-ref m :long))
                                   (pp (cffi:mem-ref p :pointer))
                                   (up (cffi:mem-ref u :pointer))
                                   (vp (cffi:mem-ref v :pointer)))
                               (list (if (cffi:null-pointer-p pp) nil
                                         (from-foreign pp :point-net (1+ n) (1+ m)))
                                     (if (cffi:null-pointer-p up) nil
                                         (from-foreign up :real-vector (1+ n)))
                                     (if (cffi:null-pointer-p vp) nil
                                         (from-foreign vp :real-vector (1+ m)))
                                     n m)))))))))))
