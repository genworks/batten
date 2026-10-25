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

;;; ------------------------------------------------------------- L5-g0

(defun %crv-geom-line-curve ()
  "A curve whose control points lie on a line in increasing order (a line to
N_CrvGetType), rational or not, 2D now and then."
  (let* ((c (random-curve :max-degree 5 :max-extra 6))
         (pw (curve-cpoints c))
         (n (curve-n c))
         (a (random-point)) (b (random-point))
         (ts (sort (let (l) (dotimes (i (1+ n) l) (push (rnd 1d0) l))) #'<)))
    (setf (first ts) 0d0 (car (last ts)) 1d0)
    (let ((i 0))
      (dolist (tt ts)
        (let ((w (if (= (aref pw i 3) +no-w+) 1d0 (aref pw i 3))))
          (dotimes (k 3)
            (unless (and (= k 2) (= (aref pw i 2) +no-z+))
              (setf (aref pw i k) (* w (+ (aref a k) (* tt (- (aref b k) (aref a k)))))))))
        (incf i)))
    c))

(defun %crv-geom-circle-curve ()
  "A circle or circular arc made by create-circ-arc-tol, of any of the three
kinds."
  (let* ((x (random-point 1d0))
         (v (random-point 1d0))
         (xx (reduce #'+ (map 'list #'* x x)))
         (s (/ (reduce #'+ (map 'list #'* x v)) xx))
         (y (make-point (- (aref v 0) (* s (aref x 0))) (- (aref v 1) (* s (aref x 1)))
                        (- (aref v 2) (* s (aref x 2)))))
         (as (rnd-range -180d0 180d0))
         (ae (+ as (case (rnd 4) (0 360d0) (t (rnd-range 10d0 350d0))))))
    (batten::create-circ-arc-tol (random-point) x y (rnd-range 0.5d0 50d0) as ae 0.01d0
                                 (elt (list +quadratic+ +quartic+ +quintic+) (rnd 3)))))

(defun %crv-geom-full-mult-curve ()
  "A curve of degree P with one interior knot of multiplicity P+1 (2P+2
control points)."
  (let* ((p (1+ (rnd 4)))
         (n+1 (+ p p 2))
         (rational (zerop (rnd 2)))
         (pw (make-cpoint-array n+1))
         (u (make-real-vector (+ n+1 p 1))))
    (dotimes (i n+1)
      (let ((w (if rational (rnd-range 0.2d0 5d0) 1d0)))
        (dotimes (k 3) (setf (aref pw i k) (* w (- (rnd 200d0) 100d0))))
        (setf (aref pw i 3) (if rational w +no-w+))))
    (dotimes (i (1+ p))
      (setf (aref u i) 0d0
            (aref u (+ i p 1)) 0.4d0
            (aref u (+ i p p 2)) 1d0))
    (%make-curve p u pw)))

(defun %crv-geom-arc-curve ()
  "A curve for the arc length and polyline tests: random, a line, a circle
or a curve with a knot of multiplicity degree + 1."
  (case (rnd 12)
    ((0 1) (%crv-geom-line-curve))
    ((2 3) (%crv-geom-circle-curve))
    (4 (%crv-geom-full-mult-curve))
    (t (random-curve :max-degree 6 :max-extra 10))))

(define-oracle-test crv-arc-length (:file "crv-geom" :trials 2000)
  (let* ((c (%crv-geom-arc-curve))
         (kn (curve-knots c))
         (lo (aref kn 0)) (hi (aref kn (1- (length kn))))
         (ul (case (rnd 4) (0 lo) (t (random-parameter kn))))
         (ur (case (rnd 10)
               ((0 1 2) hi)
               (3 ul)
               (4 (+ ul 1d-13))
               (t (random-parameter kn))))
         ;; mostly ul < ur; now and then reversed (the extraction fails)
         (swap (and (> ul ur) (/= 0 (rnd 4))))
         (ul1 (if swap ur ul))
         (ur (if swap ul ur))
         (ul ul1)
         (flg (case (rnd 12) (0 0) (1 7) ((2 3 4 5 6) +absolute+) (t +relative+)))
         (tol (if (= flg +relative+)
                  (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 2d0 7d0)))) (t 1d-5))
                  (case (rnd 3) (0 1d-3) (1 (expt 10d0 (- (rnd-range 0d0 6d0)))) (t 0.1d0)))))
    (agree (list :curve c :ul ul :ur ur :tol tol :flg flg)
           (list (batten::crv-arc-length c ul ur tol flg))
           (with-foreign ((cp c :curve))
             (with-scalars ((len :double))
               (let ((rc (c-crv-arc-length cp ul ur tol flg len)))
                 (values (c-flag-error-p rc)
                         (list (cffi:mem-ref len :double)))))))))

(defun %crv-geom-degenerate (c)
  "Curve C with every control point (or, half the time, a run of degree + 1
of them) made one point."
  (let* ((pw (curve-cpoints c)) (n (curve-n c)) (p (curve-degree c))
         (all (or (zerop (rnd 2)) (< n (1+ p))))
         (i0 (if all 0 (rnd (- n p -1))))
         (i1 (if all n (+ i0 p))))
    (do ((i (1+ i0) (1+ i))) ((> i i1) c)
      (dotimes (k 4) (setf (aref pw i k) (aref pw i0 k))))))

;;; TOL is kept positive (a zero tolerance subdivides for ever in the C), and
;;; a degenerate curve gets an absolute one.
(define-oracle-test approx-crv-with-polyline (:file "crv-geom" :trials 1500)
  (let* ((degen (zerop (rnd 12)))
         (c (if degen (%crv-geom-degenerate (random-curve :max-degree 5 :max-extra 6)) (%crv-geom-arc-curve)))
         (tflg (case (rnd 8) (0 0) ((1 2 3) (if degen +absolute+ +relative+)) (t +absolute+)))
         (tol (if (= tflg +relative+)
                  (expt 10d0 (- (rnd-range 1.5d0 4d0)))
                  (case (rnd 4) (0 1d0) (t (expt 10d0 (- (rnd-range -0.5d0 1.7d0)))))))
         (rflg (case (rnd 8) (0 0) ((1 2) +points+) ((3 4) +parameters+) (t +both+))))
    (agree (list :curve c :tol tol :tflg tflg :rflg rflg)
           (multiple-value-list (batten::approx-crv-with-polyline c tol tflg rflg))
           (with-foreign ((cp c :curve))
             (with-stacks (s)
               (with-scalars ((pts :pointer) (pars :pointer) (n :long))
                 (let ((rc (c-approx-crv-with-polyline cp tol tflg rflg pts pars n s)))
                   (if (c-flag-error-p rc)
                       (values t nil)
                       (let ((nn (cffi:mem-ref n :long)))
                         (values nil
                                 (list (unless (= rflg +parameters+)
                                         (from-foreign (cffi:mem-ref pts :pointer) :point-array (1+ nn)))
                                       (unless (= rflg +points+)
                                         (from-foreign (cffi:mem-ref pars :pointer) :real-vector (1+ nn)))
                                       nn)))))))))))
