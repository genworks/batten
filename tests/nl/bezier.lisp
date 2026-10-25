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
