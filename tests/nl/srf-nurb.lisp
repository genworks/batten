;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-nurb.lisp against NLib.

(in-package :batten-tests)

(defun some-surface () (random-surface :two-d (zerop (rnd 5))))

(define-oracle-test srf-eval-pt (:file "srf-nurb" :trials 3000)
  (let* ((s (some-surface))
         (u (random-parameter (surface-knots-u s))) (v (random-parameter (surface-knots-v s)))
         (ufl (random-side)) (vfl (random-side)))
    (agree (list :surface s :u u :v v)
           (list (batten::srf-eval-pt s u v ufl vfl (make-point)))
           (with-foreign ((sp s :surface))
             (with-foreign-out ((pt :point))
               (values (c-flag-error-p (c-srf-eval-pt sp u v ufl vfl pt))
                       (list (from-foreign pt :point))))))))

(define-oracle-test srf-derivs (:file "srf-nurb" :trials 3000)
  (let* ((s (some-surface))
         (u (random-parameter (surface-knots-u s))) (v (random-parameter (surface-knots-v s)))
         (ufl (random-side)) (vfl (random-side))
         (udr (rnd 3)) (vdr (if (zerop (rnd 2)) udr (rnd 3)))
         (mfl (if (zerop (rnd 2)) +true+ +false+))
         (upper (and (= mfl +true+) (= udr vdr))))
    (flet ((only-computed (net)
             ;; with the upper-half flag the C leaves the rest unwritten
             (when upper
               (dotimes (k (1+ udr))
                 (dotimes (l (1+ vdr))
                   (when (> (+ k l) udr) (dotimes (c 3) (setf (aref net k l c) 0d0))))))
             net))
      (agree (list :surface s :u u :v v :udr udr :vdr vdr :mfl mfl)
             (list (only-computed (batten::srf-derivs s u v ufl vfl mfl udr vdr
                                                      (make-point-net (1+ udr) (1+ vdr)))))
             (with-foreign ((sp s :surface))
               (with-foreign-out ((sd :point-net (1+ udr) (1+ vdr)))
                 (values (c-flag-error-p (c-srf-derivs sp u v ufl vfl mfl udr vdr sd))
                         (list (only-computed (from-foreign sd :point-net (1+ udr) (1+ vdr)))))))))))

(define-oracle-test srf-basis-derivs (:file "srf-nurb" :trials 2000)
  (let* ((s (some-surface))
         (p (surface-degree-u s)) (q (surface-degree-v s))
         (u (random-parameter (surface-knots-u s))) (v (random-parameter (surface-knots-v s)))
         (ufl (random-side)) (vfl (random-side))
         (udr (rnd 3)) (vdr (rnd 3)))
    (agree (list :surface s :u u :v v :udr udr :vdr vdr)
           (let ((bd (make-array (list (1+ udr) (1+ vdr) (1+ p) (1+ q)) :element-type 'double-float
                                                                     :initial-element 0d0)))
             (multiple-value-bind (usp vsp) (batten::srf-basis-derivs s u v ufl vfl +false+ udr vdr bd)
               (list usp vsp bd)))
           (with-foreign ((sp s :surface))
             (with-foreign-out ((bd :basis-derivatives))
               (with-scalars ((usp :long) (vsp :long))
                 (values (c-flag-error-p (c-srf-basis-derivs sp u v ufl vfl +false+ udr vdr bd usp vsp))
                         (list (cffi:mem-ref usp :long) (cffi:mem-ref vsp :long)
                               (from-foreign bd :basis-derivatives (1+ udr) (1+ vdr) (1+ p) (1+ q))))))))))

(define-oracle-test srf-eval-pt-pt-deriv-normal (:file "srf-nurb" :trials 3000)
  (let* ((s (random-surface))
         (u (random-parameter (surface-knots-u s))) (v (random-parameter (surface-knots-v s)))
         (ufl (random-side)) (vfl (random-side)))
    (agree (list :surface s :u u :v v)
           (let ((p (make-point)) (su (make-point)) (sv (make-point)) (n (make-point)))
             (batten::srf-eval-pt-pt-deriv-normal s u v ufl vfl p su sv n)
             (list p su sv n))
           (with-foreign ((sp s :surface))
             (with-foreign-out ((p :point) (su :point) (sv :point) (n :point))
               (values (c-flag-error-p (c-srf-eval-pt-pt-deriv-normal sp u v ufl vfl p su sv n))
                       (list (from-foreign p :point) (from-foreign su :point)
                             (from-foreign sv :point) (from-foreign n :point))))))))

(define-oracle-test srf-min-max-weight-pos-vectors (:file "srf-nurb" :trials 3000)
  (let ((s (some-surface)))
    (agree s
           (multiple-value-list (batten::srf-min-max-weight-pos-vectors s))
           (with-foreign ((sp s :surface))
             (with-scalars ((wmin :double) (wmax :double) (pmin :double) (pmax :double))
               (c-srf-min-max-weight-pos-vectors sp wmin wmax pmin pmax)
               (values nil (list (cffi:mem-ref wmin :double) (cffi:mem-ref wmax :double)
                                 (cffi:mem-ref pmin :double) (cffi:mem-ref pmax :double))))))))

(define-oracle-test srf-translate (:file "srf-nurb" :trials 3000)
  (let ((s (some-surface)) (tv (random-point)))
    (agree (list s tv)
           (list (batten::srf-translate (batten::srf-copy s) tv))
           (with-foreign ((sp s :surface) (tp tv :point))
             (c-srf-translate sp tp)
             (values nil (list (from-foreign sp :surface)))))))

(define-oracle-test srf-get-b-box (:file "srf-nurb" :trials 3000)
  (let ((s (random-surface :two-d (zerop (rnd 4)))))
    (agree s
           (list (batten::srf-get-b-box s (make-real-vector 6)))
           (with-foreign ((sp s :surface))
             (with-foreign-out ((box :box))
               (values (c-flag-error-p (c-srf-get-b-box sp box))
                       (list (from-foreign box :box))))))))

(define-oracle-test srf-reparam (:file "srf-nurb" :trials 3000)
  (let* ((s (some-surface))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (u1 (aref ku 0)) (u2 (aref ku (1- (length ku))))
         (v1 (aref kv 0)) (v2 (aref kv (1- (length kv)))))
    (flet ((ends (a b)
             ;; keep the old ends now and then, or one of them
             (case (rnd 5)
               (0 (values a b))
               (1 (values a (+ a 0.5d0 (rnd 10d0))))
               (2 (values (- b 0.5d0 (rnd 10d0)) b))
               (t (let ((lo (- (rnd 20d0) 10d0))) (values lo (+ lo 0.5d0 (rnd 10d0))))))))
      (multiple-value-bind (us ue) (ends u1 u2)
        (multiple-value-bind (vs ve) (ends v1 v2)
          (agree (list :surface s :us us :ue ue :vs vs :ve ve)
                 (list (batten::srf-reparam (batten::srf-copy s) us ue vs ve))
                 (with-foreign ((sp s :surface))
                   (c-srf-reparam sp us ue vs ve)
                   (values nil (list (from-foreign sp :surface))))))))))

;;; ------------------------------------------------------------ level 3, group 4

(define-oracle-test srf-max-diag-dist-b-box (:file "srf-nurb" :trials 3000)
  (let ((s (random-surface :two-d (zerop (rnd 4)))))
    (when (zerop (rnd 8))                 ; all control points the same
      (let ((net (surface-cpoints s)))
        (dotimes (i (array-dimension net 0))
          (dotimes (j (array-dimension net 1))
            (dotimes (k 4) (setf (aref net i j k) (aref net 0 0 k)))))))
    (agree s
           (list (batten::srf-max-diag-dist-b-box s))
           (with-foreign ((sp s :surface))
             (with-scalars ((d :double))
               (values (c-flag-error-p (c-srf-max-diag-dist-b-box sp d))
                       (list (cffi:mem-ref d :double))))))))

;;; ------------------------------------------------------------- L8-g0

(defun %st-transpose (s)
  "Surface S with u and v swapped."
  (let* ((net (surface-cpoints s))
         (n+1 (array-dimension net 0)) (m+1 (array-dimension net 1))
         (out (make-cpoint-net m+1 n+1)))
    (dotimes (i n+1)
      (dotimes (j m+1)
        (dotimes (k 4) (setf (aref out j i k) (aref net i j k)))))
    (%make-surface (surface-degree-v s) (surface-degree-u s)
                   (copy-seq (surface-knots-v s)) (copy-seq (surface-knots-u s)) out)))

(defun %st-unit (v)
  (let ((l (sqrt (+ (* (aref v 0) (aref v 0)) (* (aref v 1) (aref v 1)) (* (aref v 2) (aref v 2))))))
    (make-point (/ (aref v 0) l) (/ (aref v 1) l) (/ (aref v 2) l))))

(defun %st-cross (a b)
  (make-point (- (* (aref a 1) (aref b 2)) (* (aref a 2) (aref b 1)))
              (- (* (aref a 2) (aref b 0)) (* (aref a 0) (aref b 2)))
              (- (* (aref a 0) (aref b 1)) (* (aref a 1) (aref b 0)))))

(defun %st-axis ()
  "A random axis (values s tv x): a point, a unit direction (now and then a
coordinate axis) and a unit vector perpendicular to it."
  (let* ((s (if (zerop (rnd 4)) (make-point 0d0 0d0 0d0) (random-point 30d0)))
         (tv (case (rnd 6)
               (0 (make-point 0d0 0d0 1d0))
               (1 (make-point 0d0 1d0 0d0))
               (t (%st-unit (random-point 1d0)))))
         (x (%st-unit (%st-cross tv (%st-unit (random-point 1d0))))))
    (values s tv x)))

(defun %st-at (s tv x rad hgt)
  "The point S + RAD*X + HGT*TV."
  (make-point (+ (aref s 0) (* rad (aref x 0)) (* hgt (aref tv 0)))
              (+ (aref s 1) (* rad (aref x 1)) (* hgt (aref tv 1)))
              (+ (aref s 2) (* rad (aref x 2)) (* hgt (aref tv 2)))))

(defun %st-conic () (elt (list +quadratic+ +quintic+ +quartic+) (rnd 3)))

(defun %st-angle ()
  (case (rnd 5) (0 360d0) (1 90d0) (2 180d0) (t (rnd-range 10d0 239d0))))

(defun %st-plane ()
  "A bilinear patch over four coplanar corners, rational or not, 2D now and
then; now and then one control point moved off the plane or off a boundary."
  (let* ((p (1+ (rnd 3))) (q (1+ (rnd 3)))
         (n+1 (+ p 1 (rnd 3))) (m+1 (+ q 1 (rnd 3)))
         (two-d (zerop (rnd 5)))
         (rational (zerop (rnd 2)))
         (o (random-point 50d0))
         (a (random-point 40d0)) (b (random-point 40d0))
         (c1 (rnd-range 0.5d0 1.5d0)) (c2 (rnd-range 0.5d0 1.5d0))
         (net (make-cpoint-net n+1 m+1)))
    (when two-d (setf (aref o 2) 0d0 (aref a 2) 0d0 (aref b 2) 0d0))
    (dotimes (i n+1)
      (dotimes (j m+1)
        (let* ((u (/ (float i 1d0) (1- n+1)))
               (v (/ (float j 1d0) (1- m+1)))
               ;; corners o, o+a, o+c1*a+c2*b, o+b
               (fa (+ (* u (- 1d0 v)) (* u v c1)))
               (fb (+ (* (- 1d0 u) v) (* u v c2)))
               (w (if rational (rnd-range 0.3d0 3d0) 1d0)))
          (dotimes (k 3)
            (setf (aref net i j k) (* w (+ (aref o k) (* fa (aref a k)) (* fb (aref b k))))))
          (when two-d (setf (aref net i j 2) +no-z+))
          (setf (aref net i j 3) (if rational w +no-w+)))))
    (when (zerop (rnd 4))
      (let ((i (rnd n+1)) (j (rnd m+1)) (k (if two-d (rnd 2) (rnd 3)))
            (e (case (rnd 3) (0 1d-8) (1 1d-5) (t 1d-2))))
        (incf (aref net i j k) e)))
    (%make-surface p q (random-knots p n+1) (random-knots q m+1) net)))

(defun %st-revolved ()
  "A surface of revolution (circle in u): a sphere, a torus, a cylinder or
cone (now and then a hyperboloid), or a random profile curve revolved."
  (multiple-value-bind (s tv x) (%st-axis)
    (case (rnd 6)
      (0 (let* ((as (case (rnd 3) (0 0d0) (t (rnd 90d0))))
                (ae (case (rnd 3) (0 180d0) (t (rnd-range (+ as 10d0) 180d0)))))
           (batten::create-sphere (random-point 30d0) (rnd-range 0.5d0 60d0) as ae (%st-angle)
                                  (%st-conic) (elt (list +quadratic+ +quintic+) (rnd 2)))))
      (1 (let* ((big (rnd-range 5d0 50d0)) (small (rnd-range 0.5d0 (* 0.9d0 big)))
                (hgt (rnd-range -10d0 10d0))
                (cen (%st-at s tv x big hgt))
                (as (rnd 360d0)) (sw (case (rnd 3) (0 360d0) (t (rnd-range 10d0 300d0))))
                (circ (batten::create-circ-arc cen x tv small as (+ as sw) (%st-conic))))
           (batten::create-revolved-srf circ s tv (%st-angle) (elt (list +quadratic+ +quintic+) (rnd 2)))))
      ((2 3) (let* ((r1 (rnd-range 0.5d0 30d0))
                    (r2 (if (zerop (rnd 2)) r1 (rnd-range 0d0 30d0)))
                    (h1 (rnd-range -20d0 20d0)) (h2 (+ h1 (rnd-range 1d0 40d0)))
                    (pa (%st-at s tv x r1 h1))
                    (pb (%st-at s tv x r2 h2)))
               (when (zerop (rnd 6))       ; skew line: a hyperboloid
                 (let ((y (%st-cross tv x)))
                   (dotimes (k 3) (incf (aref pb k) (* 5d0 (aref y k))))))
               (let ((line (batten::crv-line-from-pt-and-vector
                            pa (make-point (- (aref pb 0) (aref pa 0)) (- (aref pb 1) (aref pa 1))
                                           (- (aref pb 2) (aref pa 2))))))
                 (batten::create-revolved-srf line s tv (%st-angle)
                                              (elt (list +quadratic+ +quintic+) (rnd 2))))))
      (t (batten::create-revolved-srf (random-curve :max-degree 3 :max-extra 3 :two-d nil)
                                      s tv (%st-angle) (elt (list +quadratic+ +quintic+) (rnd 2)))))))

(defun %st-ruled ()
  "A ruled surface between two random curves (linear in u), or a random
curve extruded (linear in u)."
  (if (zerop (rnd 2))
      (let ((a (random-curve :max-degree 3 :max-extra 3))
            (b (random-curve :max-degree 3 :max-extra 3)))
        (values (batten::create-ruled-srf-from-boundary-crvs a b +udir+)))
      (let ((c (random-curve :max-degree 3 :max-extra 3)))
        (batten::create-srf-extrude-crv c (random-point 1d0) (rnd-range 0.5d0 40d0) +udir+))))

(defun %st-degenerate ()
  "A random surface whose u-iso-curves are points (each column of the net
one point), now and then only some of them.  All of them only with three
columns or more: N_SrfType then averages over no circles (a NaN axis), and
with two columns the profile is a line, which sends the NaN axis into
N_ProjectPtLine; Batten's vector-normalize fails on a NaN vector where the C
carries the NaN on."
  (let* ((s (random-surface :max-degree 3 :max-extra 3))
         (net (surface-cpoints s))
         (some (or (< (array-dimension net 1) 3) (zerop (rnd 3))))
         (keep (rnd (array-dimension net 1))))  ; a column left alone when SOME
    (dotimes (j (array-dimension net 1))
      (unless (and some (or (= j keep) (zerop (rnd 2))))
        (dotimes (i (array-dimension net 0))
          (dotimes (k 4) (setf (aref net i j k) (aref net 0 j k))))))
    s))

(defun %st-surface ()
  "A random input for N_SrfType: planes, surfaces of revolution, ruled and
extruded surfaces, degenerate and free-form ones, half the time with u and v
swapped."
  (let ((s (handler-case
               (batten::with-c-float-environment
                 (case (rnd 12)
                   ((0 1) (%st-plane))
                   ((2 3 4 5) (%st-revolved))
                   ((6 7 8) (%st-ruled))
                   (9 (%st-degenerate))
                   (t (random-surface :max-degree 4 :max-extra 4 :two-d (zerop (rnd 5))))))
             (nurbs-error () (random-surface :max-degree 3 :max-extra 3)))))
    (if (zerop (rnd 2)) (%st-transpose s) s)))

(define-oracle-test srf-type (:file "srf-nurb" :trials 1500)
  (let ((s (%st-surface))
        (tol (case (rnd 8) (0 1d-3) (1 1d-9) (2 1d-12) (3 1d-4) (t 1d-6))))
    (agree (list :surface s :tol tol)
           (let ((aa (make-point)) (bb (make-point)) (cc (make-point))
                 (dd (make-point)) (xx (make-point)) (yy (make-point)))
             (multiple-value-bind (ca cb r1 r2 ac al h stp dir)
                 (batten::srf-type s tol aa bb cc dd xx yy)
               (list ca cb aa bb cc dd xx yy r1 r2 ac al h stp dir)))
           (with-foreign ((sp s :surface))
             (with-stacks (sg)
               (with-foreign-out ((faa :point) (fbb :point) (fcc :point)
                                  (fdd :point) (fxx :point) (fyy :point))
                 (with-scalars ((ca :pointer) (cb :pointer) (r1 :double) (r2 :double)
                                (ac :double) (al :double) (h :double) (stp :short) (dir :short))
                   (let ((failed (c-flag-error-p
                                  (c-srf-type sp tol ca cb faa fbb fcc fdd fxx fyy
                                              r1 r2 ac al h stp dir sg))))
                     (flet ((crv (cell)
                              (let ((pc (cffi:mem-ref cell :pointer)))
                                (if (cffi:null-pointer-p pc) nil (from-foreign pc :curve)))))
                       (values failed
                               (unless failed
                                 (list (crv ca) (crv cb)
                                       (from-foreign faa :point) (from-foreign fbb :point)
                                       (from-foreign fcc :point) (from-foreign fdd :point)
                                       (from-foreign fxx :point) (from-foreign fyy :point)
                                       (cffi:mem-ref r1 :double) (cffi:mem-ref r2 :double)
                                       (cffi:mem-ref ac :double) (cffi:mem-ref al :double)
                                       (cffi:mem-ref h :double)
                                       (cffi:mem-ref stp :short) (cffi:mem-ref dir :short)))))))))))))
