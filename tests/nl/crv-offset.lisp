;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/crv-offset.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L6-g0

(defun %co6-frame ()
  "A random orthonormal frame (values x y z) of points."
  (let* ((a (random-point 1d0))
         (a (if (< (sqrt (reduce #'+ (map 'list #'* a a))) 1d-3) (make-point 1d0 0d0 0d0) a))
         (la (sqrt (reduce #'+ (map 'list #'* a a))))
         (x (map 'point (lambda (c) (/ c la)) a))
         (h (if (< (abs (aref x 0)) 0.9d0) (make-point 1d0 0d0 0d0) (make-point 0d0 1d0 0d0)))
         (z (make-point (- (* (aref x 1) (aref h 2)) (* (aref x 2) (aref h 1)))
                        (- (* (aref x 2) (aref h 0)) (* (aref x 0) (aref h 2)))
                        (- (* (aref x 0) (aref h 1)) (* (aref x 1) (aref h 0)))))
         (lz (sqrt (reduce #'+ (map 'list #'* z z))))
         (z (map 'point (lambda (c) (/ c lz)) z))
         (y (make-point (- (* (aref z 1) (aref x 2)) (* (aref z 2) (aref x 1)))
                        (- (* (aref z 2) (aref x 0)) (* (aref z 0) (aref x 2)))
                        (- (* (aref z 0) (aref x 1)) (* (aref z 1) (aref x 0))))))
    (values x y z)))

(defun %co6-planar-curve ()
  "A planar curve of modest size and span count and the normal of its plane
(values curve normal): 2D, in the xy plane, or in a random plane."
  (let* ((kind (rnd 3))
         (c (random-curve :max-degree 4 :max-extra 4 :two-d (= kind 0)))
         (pw (curve-cpoints c))
         (n (curve-n c))
         (scale (case (rnd 3) (0 0.1d0) (1 0.05d0) (t 0.02d0))))
    (multiple-value-bind (fx fy fz) (if (= kind 2) (%co6-frame)
                                        (values (make-point 1d0 0d0 0d0) (make-point 0d0 1d0 0d0)
                                                (make-point 0d0 0d0 1d0)))
      (dotimes (i (1+ n))
        (let* ((w (if (= (aref pw i 3) +no-w+) 1d0 (aref pw i 3)))
               (a (* scale (/ (aref pw i 0) w)))
               (b (* scale (/ (aref pw i 1) w))))
          (if (= kind 0)
              (setf (aref pw i 0) (* w a) (aref pw i 1) (* w b))
              (dotimes (k 3)
                (setf (aref pw i k) (* w (+ (* a (aref fx k)) (* b (aref fy k)))))))))
      (values c fz))))

(defun %co6-scaled (v s)
  (map 'point (lambda (c) (* s c)) v))

(define-oracle-test crv-offset-pt-sampling (:file "crv-offset" :trials 2000)
  (multiple-value-bind (c nrm)
      (case (rnd 10)
        (0 (values (%crv-geom-line-curve) (random-point 1d0)))
        (1 (let ((c (%crv-geom-circle-curve))) (values c (random-point 1d0))))
        (t (%co6-planar-curve)))
    (let* ((b (case (rnd 20)
                (0 (make-point 0d0 0d0 0d0))
                (t (%co6-scaled nrm (if (zerop (rnd 2)) (rnd-range 0.5d0 3d0) (- (rnd-range 0.5d0 3d0)))))))
           (d (case (rnd 6) (0 0d0) (t (* (if (zerop (rnd 2)) 1d0 -1d0) (rnd-range 0.05d0 2d0)))))
           (deg (if (zerop (rnd 12)) 1 (+ 2 (rnd 4))))
           (tan (if (zerop (rnd 2)) batten::+yes+ batten::+no+))
           (par (case (rnd 20) (0 0) ((1 2) +uniform+) ((3 4) batten::+chordlength+)
                  ((5 6) batten::+centripetal+) (t +inherited+)))
           (tol (case (rnd 10) (0 0d0) ((1 2 3) 0.05d0) ((4 5) 0.2d0) (6 1d-2) (t 0.1d0))))
      (agree (list :curve c :b b :d d :deg deg :tan tan :par par :tol tol)
             (list (batten::crv-offset-pt-sampling c b d deg tan par tol) b)
             (with-foreign ((cp c :curve) (bp b :point))
               (with-stacks (sg)
                 (with-foreign-out ((cur :curve))
                   (let ((err (c-flag-error-p (c-crv-offset-pt-sampling cp bp d deg tan par tol cur sg))))
                     (values err
                             (unless err
                               (list (from-foreign cur :curve) (from-foreign bp :point))))))))))))
