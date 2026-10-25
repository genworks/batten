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

;;; ------------------------------------------------------------- L2-g3

(define-oracle-test create-srf-extrude-crv (:file "srf-common" :trials 3000)
  (let* ((c (random-curve))
         (w (case (rnd 20)
              (0 (make-point 0d0 0d0 0d0))
              (1 (make-point 0d0 0d0 (rnd-range 0.1d0 3d0)))
              (t (random-point 5d0))))
         (d (case (rnd 10) (0 0d0) (1 (- (rnd 10d0))) (t (rnd 50d0))))
         (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+)))))
    (agree (list :curve c :w w :d d :dir dir)
           (let ((cc (%make-curve (curve-degree c) (copy-seq (curve-knots c))
                                  (let* ((a (curve-cpoints c))
                                         (b (make-array (array-dimensions a) :element-type 'double-float)))
                                    (dotimes (i (array-total-size a) b)
                                      (setf (row-major-aref b i) (row-major-aref a i))))))
                 (ww (copy-seq w)))
             (let ((s (batten::create-srf-extrude-crv cc ww d dir)))
               ;; W is the C's by-value argument: it must come back unchanged
               (list s cc ww)))
           (with-foreign ((cur c :curve) (wp w :point))
             (with-stacks (sg)
               (with-foreign-out ((sur :surface))
                 (let ((failed (c-flag-error-p (c-create-srf-extrude-crv cur wp d dir sur sg))))
                   (values failed (unless failed
                                    (list (from-foreign sur :surface) (from-foreign cur :curve)
                                          (from-foreign wp :point)))))))))))

;;; ------------------------------------------------------------- L4-g1

(defun %sc-copy-curve (c)
  "A fresh copy of curve C."
  (%make-curve (curve-degree c) (copy-seq (curve-knots c))
               (let* ((a (curve-cpoints c))
                      (b (make-array (array-dimensions a) :element-type 'double-float)))
                 (dotimes (i (array-total-size a) b)
                   (setf (row-major-aref b i) (row-major-aref a i))))))

(define-oracle-test create-ruled-srf-from-boundary-crvs (:file "srf-common" :trials 2000)
  (let* ((two-d (zerop (rnd 5)))
         (rat (rnd 4))                    ; 0 neither, 1 P, 2 Q, 3 both rational
         (cp (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (oddp rat)))
         (cq (case (rnd 4)
               ;; the same knots as P
               (0 (let ((c (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (>= rat 2))))
                    (if (and (= (curve-degree c) (curve-degree cp)))
                        c
                        (let ((d (%sc-copy-curve cp)))
                          (dotimes (i (1+ (curve-n d)) d)
                            (let ((w (if (>= rat 2) (+ 0.2d0 (rnd 4.8d0)) 1d0)))
                              (dotimes (k 3)
                                (unless (= (aref (curve-cpoints d) i k) +no-z+)
                                  (setf (aref (curve-cpoints d) i k) (* w (- (rnd 200d0) 100d0)))))
                              (setf (aref (curve-cpoints d) i 3) (if (>= rat 2) w +no-w+))))))))
               (t (random-curve :max-degree 4 :max-extra 6 :two-d two-d :rational (>= rat 2)))))
         (dir (case (rnd 30) (0 7) (t (if (zerop (rnd 2)) +udir+ +vdir+)))))
    (agree (list :p cp :q cq :dir dir)
           (let ((lp (%sc-copy-curve cp)) (lq (%sc-copy-curve cq)))
             (list (batten::create-ruled-srf-from-boundary-crvs lp lq dir) lp lq))
           (with-foreign ((fp cp :curve) (fq cq :curve))
             (with-stacks (sg)
               (with-foreign-out ((sur :surface))
                 (let ((failed (c-flag-error-p (c-create-ruled-srf-from-boundary-crvs fp fq dir sur sg))))
                   (values failed (unless failed
                                    (list (from-foreign sur :surface) (from-foreign fp :curve)
                                          (from-foreign fq :curve)))))))))))

;;; ------------------------------------------------------------- L7-g1

(defun %sc-revolve-inputs ()
  "A random profile curve, axis and angle for N_CreateRevolvedSrf, as (values
cur s tv al ctp): the curve rational or not, 2D or not, now and then with
control points on the axis (all of them, rarely); the axis direction now and
then a coordinate axis, rarely zero; the angle mostly in (0, 360], now and
then a right angle multiple or outside; quadratic or quintic circles, quartic
ones only below 241 degrees (the C's weight array holds 9), rarely an unknown
kind."
  (let* ((c (random-curve :max-degree 5 :max-extra 6))
         (pw (curve-cpoints c))
         (n+1 (array-dimension pw 0))
         (s (if (zerop (rnd 4)) (make-point 0d0 0d0 0d0) (random-point 50d0)))
         (tv (case (rnd 40)
               (0 (make-point 0d0 0d0 0d0))
               ((1 2 3 4 5 6) (make-point 0d0 0d0 (if (zerop (rnd 2)) 1d0 -3d0)))
               ((7 8 9) (make-point 0d0 1d0 0d0))
               (t (random-point 5d0))))
         (al (case (rnd 12)
               (0 90d0) (1 180d0) (2 270d0) (3 360d0)
               (4 (rnd-range 0.1d0 5d0))
               (5 (rnd-range 360d0 500d0))
               (6 (- (rnd 100d0)))
               (t (rnd-range 0.1d0 360d0))))
         (ctp (case (rnd 40)
                (0 1)
                ((1 2 3 4) (if (and (> al 0d0) (< al 240d0)) +quartic+ +quadratic+))
                (t (if (zerop (rnd 2)) +quadratic+ +quintic+)))))
    ;; control points on the axis (in a 2D curve's z = 0 plane only when the
    ;; axis lies in it)
    (let ((on (case (rnd 12) (0 :all) ((1 2 3) :some) (t nil))))
      (when on
        (let ((two-d (= (aref pw 0 2) +no-z+)))
          (dotimes (i n+1)
            (when (or (eq on :all) (zerop (rnd 3)))
              (let* ((k (rnd-range -20d0 20d0))
                     (x (+ (aref s 0) (* k (aref tv 0))))
                     (y (+ (aref s 1) (* k (aref tv 1))))
                     (z (+ (aref s 2) (* k (aref tv 2))))
                     (w (aref pw i 3))
                     (ww (if (= w +no-w+) 1d0 w)))
                (unless (and two-d (/= z 0d0))
                  (setf (aref pw i 0) (* ww x)
                        (aref pw i 1) (* ww y))
                  (unless two-d (setf (aref pw i 2) (* ww z))))))))))
    (values c s tv al ctp)))

(define-oracle-test create-revolved-srf (:file "srf-common" :trials 1000)
  (multiple-value-bind (c s tv al ctp) (%sc-revolve-inputs)
    (agree (list :curve c :s s :t tv :al al :ctp ctp)
           (let ((cc (%sc-copy-curve c)) (ss (copy-seq s)) (tt (copy-seq tv)))
             ;; CUR is made rational and 3D in place; S and T come back unchanged
             (list (batten::create-revolved-srf cc ss tt al ctp) cc ss tt))
           (with-foreign ((cur c :curve) (sp s :point) (tp tv :point))
             (with-stacks (sg)
               (with-foreign-out ((sur :surface))
                 (let ((failed (c-flag-error-p (c-create-revolved-srf cur sp tp al ctp sur sg))))
                   (values failed (unless failed
                                    (list (from-foreign sur :surface) (from-foreign cur :curve)
                                          (from-foreign sp :point) (from-foreign tp :point)))))))))))
