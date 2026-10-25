;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/bezier.lisp against NLib.

(in-package :batten-tests)

;;; NL_RMATRIX, which marshal.lisp does not carry: highest row and column
;;; indexes, the row pointers, the NL_MATRIXTYPE enum and the bandwidth.
(cffi:defcstruct %bezier-nl-rmatrix (n :long) (m :long) (rm :pointer) (mtp :int) (bw :long))

(defun %bezier-call-matrix-routine (rows cols fn)
  "Call FN with a foreign NL_RMATRIX over a zeroed ROWS by COLS matrix and an
NL_STACKS; (values error-p (list matrix)) as AGREE wants."
  (with-foreign ((rp (make-real-matrix rows cols) :real-matrix))
    (let ((dm (batten-oracle::falloc '(:struct %bezier-nl-rmatrix) 1)))
      (setf (cffi:foreign-slot-value dm '(:struct %bezier-nl-rmatrix) 'n) (1- rows)
            (cffi:foreign-slot-value dm '(:struct %bezier-nl-rmatrix) 'm) (1- cols)
            (cffi:foreign-slot-value dm '(:struct %bezier-nl-rmatrix) 'rm) rp
            (cffi:foreign-slot-value dm '(:struct %bezier-nl-rmatrix) 'mtp) 0
            ;; a bandwidth of 0 would mark the matrix unallocated
            (cffi:foreign-slot-value dm '(:struct %bezier-nl-rmatrix) 'bw) (max 1 (1- cols)))
      (with-stacks (sg)
        (if (c-flag-error-p (funcall fn dm sg))
            (values t nil)
            (values nil (list (from-foreign rp :real-matrix rows cols))))))))

(defun %bezier-matrix-size (n m)
  "Rows and columns of a caller's matrix for highest indexes N and M: exact
mostly, larger now and then, too small now and then."
  (case (rnd 6)
    (0 (values (+ n 1 (rnd 3)) (+ m 1 (rnd 3))))
    (1 (if (zerop (rnd 2))
           (values (max 1 n) (1+ m))
           (values (1+ n) (max 1 m))))
    (t (values (1+ n) (1+ m)))))

(define-oracle-test bez-get-degree-elevation-matrix (:file "bezier" :trials 2000)
  (let* ((p (rnd 12)) (tt (rnd (- 18 p))) (q (+ p tt)))
    (multiple-value-bind (rows cols) (%bezier-matrix-size q p)
      (agree (list :p p :t tt :rows rows :cols cols)
             (list (batten::bez-get-degree-elevation-matrix p tt (make-real-matrix rows cols)))
             (%bezier-call-matrix-routine rows cols
                                          (lambda (dm sg) (c-bez-get-degree-elevation-matrix p tt dm sg)))))))

(define-oracle-test bez-func-multiply-bez-crv-matrix (:file "bezier" :trials 2000)
  (let* ((p (rnd 12)) (q (rnd (- 18 p))) (n (+ p q)))
    (multiple-value-bind (rows cols) (%bezier-matrix-size n p)
      (agree (list :p p :q q :rows rows :cols cols)
             (list (batten::bez-func-multiply-bez-crv-matrix p q (make-real-matrix rows cols)))
             (%bezier-call-matrix-routine rows cols
                                          (lambda (pm sg) (c-bez-func-multiply-bez-crv-matrix p q pm sg)))))))

(defun %bezier-triangle ()
  "Three points: random, or with coincident vertices, or isosceles."
  (let ((p0 (random-point)) (p1 (random-point)) (p2 (random-point)))
    (case (rnd 8)
      (0 (setf p1 (copy-seq p0)))
      (1 (setf p2 (copy-seq p0)))
      (2 (setf p1 (copy-seq p2)))
      (3 ;; isosceles: P1 on the bisector of P0 P2
       (let ((h (rnd-range -50d0 50d0)))
         (setf p1 (make-point (+ (* 0.5d0 (+ (aref p0 0) (aref p2 0))) (* h (- (aref p2 1) (aref p0 1))))
                              (- (* 0.5d0 (+ (aref p0 1) (aref p2 1))) (* h (- (aref p2 0) (aref p0 0))))
                              (* 0.5d0 (+ (aref p0 2) (aref p2 2)))))))
      (t nil))
    (values p0 p1 p2)))

(define-oracle-test bez-set-conic-weight (:file "bezier" :trials 4000)
  (multiple-value-bind (p0 p1 p2) (%bezier-triangle)
    (agree (list p0 p1 p2)
           (list (batten::bez-set-conic-weight p0 p1 p2))
           (with-foreign ((a p0 :point) (b p1 :point) (c p2 :point))
             (with-scalars ((w :double))
               (c-bez-set-conic-weight a b c w)
               (values nil (list (cffi:mem-ref w :double))))))))

(define-oracle-test bez-conic-split (:file "bezier" :trials 4000)
  (let ((p (random-points 3))
        (wp (case (rnd 6)
              (0 0d0)
              (1 (rnd-range -0.9d0 0d0))
              (t (rnd-range 0.05d0 6d0)))))
    (agree (list p wp)
           (let ((q (make-point-array 3)) (r (make-point-array 3)))
             (multiple-value-bind (wq wr) (batten::bez-conic-split p wp q r)
               (list q wq r wr)))
           (with-foreign ((pp p :point-array))
             (with-foreign-out ((q :point-array 3) (r :point-array 3))
               (with-scalars ((wq :double) (wr :double))
                 (c-bez-conic-split pp wp q wq r wr)
                 (values nil (list (from-foreign q :point-array 3) (cffi:mem-ref wq :double)
                                   (from-foreign r :point-array 3) (cffi:mem-ref wr :double)))))))))

(defun %bezier-random-cpoints (n+1)
  "N+1 random control points, rational half the time, 2D a fifth of it."
  (let ((rational (zerop (rnd 2))) (two-d (zerop (rnd 5)))
        (a (make-cpoint-array n+1)))
    (dotimes (i n+1 a)
      (let ((cp (random-cpoint :rational rational :two-d two-d)))
        (dotimes (k 4) (setf (aref a i k) (aref cp k)))))))

(define-oracle-test bez-split (:file "bezier" :trials 4000)
  (let* ((p (rnd 13))
         (pw (%bezier-random-cpoints (1+ p)))
         (u (case (rnd 10)
              (0 0d0)
              (1 1d0)
              (2 (rnd-range -0.5d0 0d0))
              (3 (rnd-range 1.000001d0 1.5d0))
              (t (rnd 1d0)))))
    (agree (list :pw pw :p p :u u)
           (let ((qw (make-cpoint-array (1+ p))) (rw (make-cpoint-array (1+ p))))
             (batten::bez-split pw p u qw rw)
             (list qw rw))
           (with-foreign ((a pw :cpoint-array))
             (with-foreign-out ((qw :cpoint-array (1+ p)) (rw :cpoint-array (1+ p)))
               (values (c-flag-error-p (c-bez-split a p u qw rw))
                       (list (from-foreign qw :cpoint-array (1+ p))
                             (from-foreign rw :cpoint-array (1+ p)))))))))

;;; ------------------------------------------------------------- L1-g5

(defun %bezier-rmatrix-spec (n m kind)
  "How to give a routine its NL_RMATRIX for highest indexes N and M: NIL for
the NULL matrix (the routine computes it), :exact for the matrix KIND
(:elevation or :product) computed beforehand, or a random matrix."
  (case (rnd 4)
    ((0 1) nil)
    (2 (list :exact kind n m))
    (t (list :random (random-real-matrix (1+ n) (1+ m) 2d0)))))

(defun %bezier-rmatrix (spec p tt)
  "A fresh rmatrix from SPEC (see %BEZIER-RMATRIX-SPEC) for degree P and
increment (or second degree) TT."
  (cond ((null spec) (batten::make-rmatrix))
        ((eq (first spec) :exact)
         (destructuring-bind (kind n m) (rest spec)
           (batten::make-rmatrix n m (if (eq kind :elevation)
                                         (batten::bez-get-degree-elevation-matrix p tt)
                                         (batten::bez-func-multiply-bez-crv-matrix p tt))
                                 +mt-full+ p)))
        (t (let ((rm (second spec)))
             (batten::make-rmatrix (1- (array-dimension rm 0)) (1- (array-dimension rm 1))
                                   (%bezier-copy-array rm) +mt-full+ p)))))

(defun %bezier-copy-array (a)
  (let ((b (make-array (array-dimensions a) :element-type 'double-float)))
    (dotimes (i (array-total-size a) b) (setf (row-major-aref b i) (row-major-aref a i)))))

(defun %bezier-range (hi)
  "First and last indexes within 0..HI: the whole range now and then."
  (if (zerop (rnd 3))
      (values 0 hi)
      (let* ((a (rnd (1+ hi))) (b (+ a (rnd (- (1+ hi) a))))) (values a b))))

(defun %bezier-bad-range (hi)
  "Start and end indexes, now and then outside 0 <= s <= e <= HI."
  (case (rnd 12)
    (0 (values -1 hi))
    (1 (values 0 (1+ hi)))
    (2 (values (min hi 1) 0))
    (t (%bezier-range hi))))

(define-oracle-test bez-elevate-degree (:file "bezier" :trials 3000)
  (let* ((p (rnd 12)) (tt (rnd (- 18 p))) (q (+ p tt))
         (pw (%bezier-random-cpoints (1+ p)))
         (spec (%bezier-rmatrix-spec q p :elevation))
         (qw0 (%bezier-random-cpoints (1+ q))))
    (multiple-value-bind (f l) (%bezier-range q)
      (agree (list :pw pw :p p :t tt :spec spec :f f :l l)
             (list (batten::bez-elevate-degree pw p tt (%bezier-rmatrix spec p tt) f l (%bezier-copy-array qw0)))
             (with-foreign ((a pw :cpoint-array) (dm (%bezier-rmatrix spec p tt) :rmatrix)
                            (qw qw0 :cpoint-array))
               (values (c-flag-error-p (c-bez-elevate-degree a p tt dm f l qw))
                       (list (from-foreign qw :cpoint-array (1+ q)))))))))

(define-oracle-test bez-func-degree-elevate (:file "bezier" :trials 3000)
  (let* ((p (rnd 12)) (tt (if (zerop (rnd 15)) (- (rnd 3) 3) (rnd (- 18 p))))
         (q (+ p tt))
         (tt (if (zerop (rnd 20)) (+ tt (- 18 q)) tt)) ; now and then past *dmax*
         (q (+ p tt))
         (f (random-real-vector (1+ p)))
         (spec (when (and (<= 0 tt) (<= q 17)) (%bezier-rmatrix-spec q p :elevation)))
         (g0 (random-real-vector (max 1 (1+ q)))))
    (multiple-value-bind (s e) (%bezier-bad-range (max 0 q))
      (agree (list :f f :p p :t tt :spec spec :s s :e e)
             (list (batten::bez-func-degree-elevate f p tt (%bezier-rmatrix spec p tt) s e (copy-seq g0)))
             (with-foreign ((fp f :real-vector) (dm (%bezier-rmatrix spec p tt) :rmatrix)
                            (g g0 :real-vector))
               (values (c-flag-error-p (c-bez-func-degree-elevate fp p tt dm s e g))
                       (list (from-foreign g :real-vector (length g0)))))))))

(define-oracle-test bez-func-multiply-bez-func (:file "bezier" :trials 3000)
  (let* ((p (rnd 12)) (q (rnd (- 18 p))) (n (+ p q))
         (f (random-real-vector (1+ p))) (g (random-real-vector (1+ q)))
         (spec (%bezier-rmatrix-spec n p :product))
         (fg0 (random-real-vector (1+ n))))
    (multiple-value-bind (s e) (%bezier-bad-range n)
      (agree (list :f f :p p :g g :q q :spec spec :s s :e e)
             (list (batten::bez-func-multiply-bez-func f p g q (%bezier-rmatrix spec p q) s e (copy-seq fg0)))
             (with-foreign ((fp f :real-vector) (gp g :real-vector)
                            (pm (%bezier-rmatrix spec p q) :rmatrix) (fg fg0 :real-vector))
               (values (c-flag-error-p (c-bez-func-multiply-bez-func fp p gp q pm s e fg))
                       (list (from-foreign fg :real-vector (1+ n)))))))))

(define-oracle-test bez-func-multiply-bez-crv (:file "bezier" :trials 3000)
  (let* ((p (rnd 12)) (q (rnd (- 18 p))) (n (+ p q))
         (f (random-real-vector (1+ p) 3d0)) (pw (%bezier-random-cpoints (1+ q)))
         (spec (%bezier-rmatrix-spec n p :product))
         (qw0 (%bezier-random-cpoints (1+ n))))
    (multiple-value-bind (s e) (%bezier-bad-range n)
      (agree (list :f f :p p :pw pw :q q :spec spec :s s :e e)
             (list (batten::bez-func-multiply-bez-crv f p pw q (%bezier-rmatrix spec p q) s e (%bezier-copy-array qw0)))
             (with-foreign ((fp f :real-vector) (a pw :cpoint-array)
                            (pm (%bezier-rmatrix spec p q) :rmatrix) (qw qw0 :cpoint-array))
               (values (c-flag-error-p (c-bez-func-multiply-bez-crv fp p a q pm s e qw))
                       (list (from-foreign qw :cpoint-array (1+ n)))))))))

(defun %bezier-random-cnet (rows cols)
  "A random cpoint-net, rational half the time, 2D a fifth of it."
  (let ((rational (zerop (rnd 2))) (two-d (zerop (rnd 5)))
        (a (make-cpoint-net rows cols)))
    (dotimes (i rows a)
      (dotimes (j cols)
        (let ((cp (random-cpoint :rational rational :two-d two-d)))
          (dotimes (k 4) (setf (aref a i j k) (aref cp k))))))))

(define-oracle-test bez-srf-elevate-degree (:file "bezier" :trials 3000)
  (let* ((r (rnd 12)) (tt (rnd (- 18 r))) (q (+ r tt))
         (dir (case (rnd 20) (0 3) (t (if (zerop (rnd 2)) +udir+ +vdir+))))
         (other (1+ (rnd 6)))
         (roc (rnd other))
         (pw (if (= dir +vdir+) (%bezier-random-cnet other (1+ r)) (%bezier-random-cnet (1+ r) other)))
         (qw0 (if (= dir +vdir+) (%bezier-random-cnet other (1+ q)) (%bezier-random-cnet (1+ q) other)))
         (spec (%bezier-rmatrix-spec q r :elevation)))
    (multiple-value-bind (f l) (%bezier-range q)
      (agree (list :pw pw :r r :t tt :spec spec :dir dir :f f :l l :roc roc)
             (list (batten::bez-srf-elevate-degree pw r tt (%bezier-rmatrix spec r tt) dir f l roc
                                                   (%bezier-copy-array qw0)))
             (with-foreign ((a pw :cpoint-net) (dm (%bezier-rmatrix spec r tt) :rmatrix)
                            (qw qw0 :cpoint-net))
               (values (c-flag-error-p (c-bez-srf-elevate-degree a r tt dm dir f l roc qw))
                       (list (from-foreign qw :cpoint-net (array-dimension qw0 0)
                                           (array-dimension qw0 1)))))))))

;;; ------------------------------------------------------------- L2-g7

(defun %bezier-conic-input ()
  "Inputs for N_BezConicArcFromPtsAndTangents, as (list P0 T0 P2 T2 P): a
true conic arc (tangents meeting, P on the arc), parallel tangents (the
middle point at infinity), random points, or degenerate ones; planar a third
of the time."
  (let* ((planar (zerop (rnd 3)))
         (pt (lambda () (let ((q (random-point))) (when planar (setf (aref q 2) 0d0)) q)))
         (p0 (funcall pt)) (p2 (funcall pt)) (p1 (funcall pt))
         (lin (lambda (a x b y) (make-point (+ (* a (aref x 0)) (* b (aref y 0)))
                                            (+ (* a (aref x 1)) (* b (aref y 1)))
                                            (+ (* a (aref x 2)) (* b (aref y 2)))))))
    (case (rnd 10)
      ((0 1 2 3)                         ; a conic arc
       (let* ((w (rnd-range 0.05d0 5d0))
              (u (case (rnd 8) (0 0d0) (1 1d0) (2 0.5d0) (t (rnd 1d0))))
              (a (* (- 1d0 u) (- 1d0 u))) (b (* 2d0 u (- 1d0 u) w)) (c (* u u))
              (d (+ a b c))
              (sc (rnd-range 0.2d0 3d0))
              (q (make-point (/ (+ (* a (aref p0 0)) (* b (aref p1 0)) (* c (aref p2 0))) d)
                             (/ (+ (* a (aref p0 1)) (* b (aref p1 1)) (* c (aref p2 1))) d)
                             (/ (+ (* a (aref p0 2)) (* b (aref p1 2)) (* c (aref p2 2))) d))))
         (list p0 (funcall lin sc p1 (- sc) p0)
               p2 (funcall lin (if (zerop (rnd 2)) 1d0 -1d0) p2 -1d0 p1) q)))
      ((4 5 6)                           ; parallel tangents
       (let* ((tv (funcall pt))
              (s (case (rnd 6) (0 0d0) (1 1d0) (t (rnd 1d0))))
              (k (rnd-range -2d0 2d0)))
         (list p0 tv p2 (funcall lin (if (zerop (rnd 2)) 1d0 -1.5d0) tv 0d0 tv)
               (funcall lin 1d0 (funcall lin (- 1d0 s) p0 s p2) k tv))))
      (7 (list p0 (funcall pt) p2 (funcall pt) (funcall pt)))
      (8 (list p0 (make-point) p2 (funcall pt) (funcall pt)))
      (t (list p0 (funcall lin 1d0 p1 -1d0 p0) p2 (funcall lin 1d0 p2 -1d0 p1)
               (if (zerop (rnd 2)) (copy-seq p1) (copy-seq p0)))))))

(define-oracle-test bez-conic-arc-from-pts-and-tangents (:file "bezier" :trials 6000)
  (destructuring-bind (p0 t0 p2 t2 p) (%bezier-conic-input)
    (agree (list :p0 p0 :t0 t0 :p2 p2 :t2 t2 :p p)
           (let ((p1 (make-point)))
             (let ((w1 (batten::bez-conic-arc-from-pts-and-tangents p0 t0 p2 t2 p p1)))
               (list p1 w1)))
           (with-foreign ((a p0 :point) (b t0 :point) (c p2 :point) (d t2 :point) (e p :point))
             (with-foreign-out ((p1 :point))
               (with-scalars ((w1 :double))
                 (let ((failed (c-flag-error-p (c-bez-conic-arc-from-pts-and-tangents a b c d e p1 w1))))
                   (values failed (unless failed
                                    (list (from-foreign p1 :point) (cffi:mem-ref w1 :double)))))))))))

(define-oracle-test bez-reparam (:file "bezier" :trials 3000)
  (let* ((p (rnd 8)) (q (rnd 8))
         (q (if (and (> (* p q) 17) (plusp (rnd 8))) (floor 17 (max p 1)) q))
         (pq (* p q))
         (pw (%bezier-random-cpoints (1+ p)))
         (f (let ((v (random-real-vector (1+ q) 0.5d0)))
              (dotimes (i (1+ q) (sort v #'<)) (incf (aref v i) 0.5d0))))
         (bin (when (and (<= pq 17) (zerop (rnd 4)))
                (batten::pascal-tri-row (make-array (list (1+ pq) (1+ pq)) :element-type 'fixnum
                                                                            :initial-element 0)
                                        pq)))
         (qw0 (%bezier-random-cpoints (1+ pq))))
    (agree (list :pw pw :p p :f f :q q :bin (and bin t))
           (list (batten::bez-reparam pw p f q bin (%bezier-copy-array qw0)))
           (with-foreign ((a pw :cpoint-array) (fp f :real-vector) (qw qw0 :cpoint-array))
             (with-foreign ((b (or bin (make-array '(1 1) :element-type 'fixnum :initial-element 0))
                               :index-matrix))
               (values (c-flag-error-p (c-bez-reparam a p fp q (if bin b (cffi:null-pointer)) qw))
                       (list (from-foreign qw :cpoint-array (1+ pq)))))))))
