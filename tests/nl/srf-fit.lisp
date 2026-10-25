;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/srf-fit.lisp against NLib.

(in-package :batten-tests)

(defun %copy-real-matrix (a)
  (let ((b (make-real-matrix (array-dimension a 0) (array-dimension a 1))))
    (dotimes (i (array-total-size a) b)
      (setf (row-major-aref b i) (row-major-aref a i)))))

(defun %removable-knot (knots deg count)
  "A random knot index r of an interior knot with knots[r] /= knots[r+1], and
its multiplicity s, for a knot vector of degree DEG with control points
0..COUNT; NIL when there is no interior knot."
  (let ((rs nil))
    (do ((r (1+ deg) (1+ r))) ((> r count))
      (when (/= (aref knots r) (aref knots (1+ r))) (push r rs)))
    (when rs
      (let* ((r (nth (rnd (length rs)) rs))
             (s (do ((j r (1- j)) (c 0 (1+ c))) ((/= (aref knots j) (aref knots r)) c))))
        (values r s)))))

(define-oracle-test fit-srf-removal-boundary (:file "srf-fit" :trials 8000)
  (let ((dir (if (zerop (rnd 50)) 3 (if (zerop (rnd 2)) batten::+udir+ batten::+vdir+)))
        sur r s)
    ;; a surface with an interior knot in the direction of removal
    (do () (r)
      (setf sur (random-surface :two-d (zerop (rnd 5))))
      (multiple-value-setq (r s)
        (if (= dir batten::+vdir+)
            (%removable-knot (surface-knots-v sur) (surface-degree-v sur) (surface-m sur))
            (%removable-knot (surface-knots-u sur) (surface-degree-u sur) (surface-n sur)))))
    (let* ((n (surface-n sur)) (m (surface-m sur))
           (other (if (= dir batten::+vdir+) n m))
           (f (rnd (1+ other)))
           (l (if (zerop (rnd 10)) (max 0 (1- f)) (+ f (rnd (- (1+ other) f)))))
           (br (random-real-matrix (1+ n) (1+ m) 1d0))
           (mr (case (rnd 3) (0 batten::+bigd+) (1 -1d0) (t (rnd 100d0)))))
      (agree (list :surface sur :r r :s s :f f :l l :dir dir :mr mr)
             (let ((b (%copy-real-matrix br)))
               (list (batten::fit-srf-removal-boundary sur r s f l dir b mr) b))
             (with-foreign ((sp sur :surface) (bp br :real-matrix))
               (with-scalars ((mrp :double))
                 (setf (cffi:mem-ref mrp :double) mr)
                 (values (c-flag-error-p (c-fit-srf-removal-boundary sp r s f l dir bp mrp))
                         (list (cffi:mem-ref mrp :double)
                               (from-foreign bp :real-matrix (1+ n) (1+ m))))))))))

;;; ------------------------------------------------ least-squares surface fits

(defun %copy-net (a)
  (let ((b (make-array (array-dimensions a) :element-type 'double-float)))
    (dotimes (i (array-total-size a) b)
      (setf (row-major-aref b i) (row-major-aref a i)))))

(defun %random-grid (rows cols &key (k 3) (rational nil) (two-d nil))
  "A ROWS x COLS net of K-coordinate points sampled roughly from a grid: a
perturbed lattice with a wavy height most of the time, scattered points now
and then, and now and then repeated rows or columns or one point everywhere
(the degenerate cases).  For K = 4 the points are weighted when RATIONAL and
lack z when TWO-D."
  (let ((a (make-array (list rows cols k) :element-type 'double-float :initial-element 0d0))
        (kind (rnd 12))
        (dx (+ 0.5d0 (rnd 10d0))) (dy (+ 0.5d0 (rnd 10d0)))
        (noise (case (rnd 3) (0 0d0) (1 0.1d0) (t (rnd 2d0))))
        (amp (rnd 20d0)) (fu (rnd 2d0)) (fv (rnd 2d0)))
    (dotimes (i rows)
      (dotimes (j cols)
        (let* ((x (if (= kind 0) (- (rnd 200d0) 100d0) (+ (* i dx) (- (rnd (* 2 noise)) noise))))
               (y (if (= kind 0) (- (rnd 200d0) 100d0) (+ (* j dy) (- (rnd (* 2 noise)) noise))))
               (z (if (= kind 0) (- (rnd 200d0) 100d0) (* amp (sin (+ (* fu i) (* fv j))))))
               (w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
          (setf (aref a i j 0) (* w x) (aref a i j 1) (* w y))
          (when (>= k 3)
            (setf (aref a i j 2) (if two-d +no-z+ (* w z))))
          (when (= k 4)
            (setf (aref a i j 3) (if rational w +no-w+))))))
    (case kind
      ;; a repeated row or column
      (1 (when (> rows 2)
           (let ((i (1+ (rnd (- rows 2)))))
             (dotimes (j cols) (dotimes (c k) (setf (aref a i j c) (aref a (1- i) j c)))))))
      (2 (when (> cols 2)
           (let ((j (1+ (rnd (- cols 2)))))
             (dotimes (i rows) (dotimes (c k) (setf (aref a i j c) (aref a i (1- j) c)))))))
      ;; every point the same
      (3 (when (zerop (rnd 4))
           (dotimes (i rows) (dotimes (j cols) (dotimes (c k) (setf (aref a i j c) (aref a 0 0 c))))))))
    a))

(defun %random-fit-sizes ()
  "Highest data indexes R, S, control point indexes N, M and degrees P, Q for
a least-squares surface fit: valid mostly, now and then out of range."
  (let* ((r (1+ (rnd 12))) (s (1+ (rnd 12)))
         (n (1+ (rnd r))) (m (1+ (rnd s)))
         (p (1+ (rnd (min n 5)))) (q (1+ (rnd (min m 5)))))
    ;; ruled and bilinear cases more often
    (case (rnd 12)
      (0 (setf n 1 p 1))
      (1 (setf m 1 q 1))
      (2 (setf n 1 p 1 m 1 q 1))
      (3 (setf n r p (1+ (rnd (min n 5)))))
      (4 (case (rnd 5)
           (0 (setf n (1+ r)))
           (1 (setf n 0))
           (2 (setf p (1+ n)))
           (3 (setf m (1+ s)))
           (t (setf q (1+ m))))))
    (values r s n m p q)))

(define-oracle-test fit-srf-lst-sq-approx (:file "srf-fit" :trials 6000)
  (multiple-value-bind (r s n m p q) (%random-fit-sizes)
    (let ((pts (%random-grid (1+ r) (1+ s)))
          (par (case (rnd 40) (0 7) (t (1+ (rnd 3))))))
      (agree (list :r r :s s :n n :m m :p p :q q :par par :p-net pts)
             (let ((pl (%copy-net pts)))
               (list (batten::fit-srf-lst-sq-approx pl r s n m p q par) pl))
             (with-foreign ((pp pts :point-net))
               (with-stacks (sg)
                 (with-foreign-out ((sur :surface))
                   (if (c-flag-error-p (c-fit-srf-lst-sq-approx pp r s n m p q par sur sg))
                       (values t nil)
                       (values nil (list (from-foreign sur :surface)
                                         (from-foreign pp :point-net (1+ r) (1+ s))))))))))))

(defun %random-fit-params (count lo hi)
  "COUNT+1 nondecreasing parameters from LO to HI, now and then repeated."
  (let ((u (make-real-vector (1+ count))))
    (dotimes (i (1+ count))
      (setf (aref u i) (+ lo (rnd (- hi lo)))))
    (setf u (sort u #'<))
    (setf (aref u 0) lo (aref u count) hi)
    (when (and (> count 2) (zerop (rnd 6)))
      (let ((i (1+ (rnd (1- count)))))
        (setf (aref u i) (aref u (1- i)))))
    u))

(define-oracle-test fit-srf-lst-sq-knots (:file "srf-fit" :trials 6000)
  (multiple-value-bind (r s n m p q) (%random-fit-sizes)
    (let* ((ptp (case (rnd 40) (0 3) (t (if (zerop (rnd 2)) +epoint+ +hpoint+))))
           (a (if (= ptp +hpoint+)
                  (%random-grid (1+ r) (1+ s) :k 4 :rational (zerop (rnd 2)) :two-d (zerop (rnd 5)))
                  (%random-grid (1+ r) (1+ s))))
           (ulo (if (zerop (rnd 2)) 0d0 (- (rnd 10d0) 5d0))) (uhi (+ ulo 0.5d0 (rnd 5d0)))
           (vlo (if (zerop (rnd 2)) 0d0 (- (rnd 10d0) 5d0))) (vhi (+ vlo 0.5d0 (rnd 5d0)))
           (u (%random-fit-params r ulo uhi))
           (v (%random-fit-params s vlo vhi))
           (knu (if (and (<= 1 p n) (<= n r) (zerop (rnd 2)))
                    (batten::fit-calc-knot-vector-crv-approx u r n p (make-real-vector (+ n p 2)))
                    (random-knots (max p 1) (+ n 1) :lo ulo :hi uhi)))
           (knv (if (and (<= 1 q m) (<= m s) (zerop (rnd 2)))
                    (batten::fit-calc-knot-vector-crv-approx v s m q (make-real-vector (+ m q 2)))
                    (random-knots (max q 1) (+ m 1) :lo vlo :hi vhi))))
      ;; now and then a knot vector of the wrong length
      (when (zerop (rnd 30))
        (setf knu (random-knots (max p 1) (+ n 2) :lo ulo :hi uhi)))
      (agree (list :r r :s s :ptp ptp :u u :v v :knu knu :knv knv :n n :m m :p p :q q :a a)
             (let ((al (%copy-net a)))
               (list (batten::fit-srf-lst-sq-knots al r s ptp u v knu knv n m p q) al))
             (with-foreign ((ap a (if (= ptp +hpoint+) :cpoint-net :point-net))
                            (up u :real-vector) (vp v :real-vector)
                            (ku knu :knots) (kv knv :knots))
               (with-stacks (sg)
                 (with-foreign-out ((sur :surface))
                   (if (c-flag-error-p (c-fit-srf-lst-sq-knots ap r s ptp up vp ku kv n m p q sur sg))
                       (values t nil)
                       (values nil (list (from-foreign sur :surface)
                                         (from-foreign ap (if (= ptp +hpoint+) :cpoint-net :point-net)
                                                       (1+ r) (1+ s))))))))))))
