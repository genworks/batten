;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-approx.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------ subsurface bounded by 4 curves

(defun %approx-tame-surface ()
  "A non-rational surface of degrees 2 or 3 (1 now and then) over a perturbed
lattice of control points, a wavy height on it."
  (let* ((p (if (zerop (rnd 15)) 1 (+ 2 (rnd 2)))) (q (if (zerop (rnd 15)) 1 (+ 2 (rnd 2))))
         ;; many control points now and then (the C thins out its samples
         ;; for the parameter functions beyond 70)
         (many (rnd 20))
         (n+1 (+ p 1 (if (= many 0) (+ 34 (rnd 8)) (rnd (if (< many 4) 20 6)))))
         (m+1 (+ q 1 (if (= many 1) (+ 34 (rnd 8)) (rnd (if (< many 4) 20 6)))))
         (net (make-cpoint-net n+1 m+1))
         (dx (+ 1d0 (rnd 10d0))) (dy (+ 1d0 (rnd 10d0)))
         (noise (* 0.2d0 (min dx dy) (rnd 1d0)))
         (amp (rnd 5d0)) (fu (rnd 1.5d0)) (fv (rnd 1.5d0)))
    (dotimes (i n+1)
      (dotimes (j m+1)
        (setf (aref net i j 0) (+ (* i dx) (- (rnd (* 2 noise)) noise))
              (aref net i j 1) (+ (* j dy) (- (rnd (* 2 noise)) noise))
              (aref net i j 2) (* amp (sin (+ (* fu i) (* fv j))))
              (aref net i j 3) +no-w+)))
    (multiple-value-bind (ulo uhi) (batten-random::random-domain)
      (multiple-value-bind (vlo vhi) (batten-random::random-domain)
        (batten::%make-surface p q (random-knots p n+1 :lo ulo :hi uhi)
                               (random-knots q m+1 :lo vlo :hi vhi) net)))))

(defun %uv-line-curve (a0 a1 b0 b1)
  "A non-rational curve in the (u,v) plane on the segment (A0,A1)-(B0,B1),
of degree 1 to 3 on a random domain, parametrized nearly in proportion to
the segment: its control points at the fractions of its Greville abscissae
along it (linear precision), now and then moved a little; 2D (no z) now and
then."
  (let* ((deg (1+ (rnd 3)))
         (n+1 (+ deg 1 (rnd 4)))
         (cp (make-cpoint-array n+1))
         (fr (make-real-vector n+1))
         (two-d (zerop (rnd 3)))
         (uneven (zerop (rnd 3))))
    (multiple-value-bind (lo hi) (if (zerop (rnd 2)) (values 0d0 1d0) (batten-random::random-domain))
      (let ((knots (random-knots deg n+1 :lo lo :hi hi)))
        (dotimes (i n+1)
          (let ((g 0d0))
            (dotimes (j deg) (incf g (aref knots (+ i j 1))))
            (setf (aref fr i) (/ (- (/ g deg) lo) (- hi lo)))))
        (when uneven
          (do ((i 1 (1+ i))) ((>= i (1- n+1)))
            (let ((gap (min (- (aref fr i) (aref fr (1- i))) (- (aref fr (1+ i)) (aref fr i)))))
              (incf (aref fr i) (* gap (- (rnd 0.6d0) 0.3d0))))))
        (setf (aref fr 0) 0d0 (aref fr (1- n+1)) 1d0)
        (dotimes (i n+1)
          (let ((f (aref fr i)))
            (setf (aref cp i 0) (+ a0 (* f (- b0 a0)))
                  (aref cp i 1) (+ a1 (* f (- b1 a1)))
                  (aref cp i 2) (if two-d +no-z+ 0d0)
                  (aref cp i 3) +no-w+)))
        (batten::%make-curve deg knots cp)))))

(defun %approx-sub-srf-case ()
  "A surface, a rectangle [U0,U1] x [V0,V1] in its domain, and the four 3D
boundary curves of that subsurface with their (u,v) curves: (values sur
cur-s-3d cur-t-3d cur-s-2d cur-t-2d)."
  (let* ((sur (%approx-tame-surface))
         (ku (surface-knots-u sur)) (kv (surface-knots-v sur))
         (ua (aref ku 0)) (ub (aref ku (1- (length ku))))
         (va (aref kv 0)) (vb (aref kv (1- (length kv))))
         (whole (zerop (rnd 3)))
         (collapse (and whole (rnd 6)))
         (u0 (if whole ua (+ ua (* (- ub ua) (rnd 0.4d0)))))
         (u1 (if whole ub (- ub (* (- ub ua) (rnd 0.4d0)))))
         (v0 (if whole va (+ va (* (- vb va) (rnd 0.4d0)))))
         (v1 (if whole vb (- vb (* (- vb va) (rnd 0.4d0)))))
         (sub (progn
                ;; now and then one edge of the whole surface collapsed to a
                ;; point: a boundary curve of length 0
                (let ((net (surface-cpoints sur)))
                  (case collapse
                    (0 (dotimes (i (array-dimension net 0))
                         (dotimes (c 3) (setf (aref net i 0 c) (aref net 0 0 c)))))
                    (1 (dotimes (j (array-dimension net 1))
                         (dotimes (c 3) (setf (aref net 0 j c) (aref net 0 0 c)))))))
                (if whole sur (batten::srf-extract-patch sur u0 u1 v0 v1))))
         (p (surface-degree-u sub)) (q (surface-degree-v sub))
         (net (surface-cpoints sub))
         (n0 (surface-n sub)) (m0 (surface-m sub))
         (su (surface-knots-u sub)) (sv (surface-knots-v sub))
         (cs3 (list (%net-row-curve p su net 0 t) (%net-row-curve p su net m0 t)))
         (ct3 (list (%net-row-curve q sv net 0 nil) (%net-row-curve q sv net n0 nil))))
    (flet ((vary (curves)
             ;; the same curves, not compatible: one degree raised, or knots
             ;; inserted into one, or weights of 1 on both
             (let ((k (rnd 2)))
               (case (rnd 8)
                 (0 (setf (nth k curves) (batten::crv-elevate-degree (nth k curves) 1)))
                 (1 (let* ((c (nth k curves)) (kn (curve-knots c))
                           (lo (aref kn 0)) (hi (aref kn (1- (length kn))))
                           (x (sort (let (l) (dotimes (i (1+ (rnd 3)) l)
                                               (push (+ lo (* (- hi lo) (+ 0.05d0 (rnd 0.9d0)))) l)))
                                    #'<)))
                      (setf (nth k curves) (batten::crv-refine c (coerce x 'real-vector)))))
                 (2 (dolist (c curves)
                      (let ((cp (curve-cpoints c)))
                        (dotimes (i (array-dimension cp 0)) (setf (aref cp i 3) 1d0))))))
               curves)))
      (values sur (vary cs3) (vary ct3)
              (list (%uv-line-curve u0 v0 u1 v0) (%uv-line-curve u0 v1 u1 v1))
              (list (%uv-line-curve u0 v0 u0 v1) (%uv-line-curve u1 v0 u1 v1))))))

(defun %curve-pair (curves)
  "A foreign NL_CURVE *[2] holding the two curves of the list CURVES."
  (let ((cp (batten-oracle::falloc :pointer 2)))
    (setf (cffi:mem-aref cp :pointer 0) (batten-oracle::to-foreign (first curves) :curve)
          (cffi:mem-aref cp :pointer 1) (batten-oracle::to-foreign (second curves) :curve))
    cp))

(define-oracle-test approx-sub-srf-with-srf (:file "srf-approx" :trials 3000)
  (multiple-value-bind (sur cs3 ct3 cs2 ct2) (%approx-sub-srf-case)
    ;; both functions, neither, or (now and then) one alone
    (let* ((fns (case (rnd 20) (0 :u) (1 :v) (t (if (zerop (rnd 2)) :both :none))))
           (ufn (member fns '(:u :both)))
           (vgn (member fns '(:v :both))))
      (agree (list :surface sur :cur-s-3d cs3 :cur-t-3d ct3 :cur-s-2d cs2 :cur-t-2d ct2 :fns fns)
             (multiple-value-list
              (batten::approx-sub-srf-with-srf cs3 ct3 cs2 ct2 sur (and ufn t) (and vgn t)))
             (with-foreign ((sp sur :surface))
               (let ((a (%curve-pair cs3)) (b (%curve-pair ct3))
                     (c (%curve-pair cs2)) (d (%curve-pair ct2)))
                 (with-stacks (sg)
                   (with-foreign-out ((sub :surface) (uf :sfun) (vg :sfun))
                     (if (c-flag-error-p
                          (c-approx-sub-srf-with-srf a b c d sp sub
                                                     (if ufn uf (cffi:null-pointer))
                                                     (if vgn vg (cffi:null-pointer)) sg))
                         (values t nil)
                         (values nil (list (from-foreign sub :surface)
                                           (if ufn (from-foreign uf :sfun) nil)
                                           (if vgn (from-foreign vg :sfun) nil))))))))))))
