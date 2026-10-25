;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/frame-basic.lisp against NLib.

(in-package :batten-tests)

(defun %frame-basic-cpoints (rows)
  "A cpoint-array of ROWS random control points, rational and 2D now and then."
  (let ((pw (make-cpoint-array rows)) (two-d (zerop (rnd 4))))
    (dotimes (i rows pw)
      (let ((c (random-cpoint :two-d two-d)))
        (dotimes (j 4) (setf (aref pw i j) (aref c j)))))))

(defun %frame-basic-net (rows cols k)
  "A ROWS by COLS net of random points with K coordinates (3 or 4)."
  (let ((net (make-array (list rows cols k) :element-type 'double-float)))
    (dotimes (i rows net)
      (dotimes (j cols)
        (dotimes (c k)
          (setf (aref net i j c) (- (rnd 200d0) 100d0)))))))

(defun %frame-basic-highest (rows)
  "A highest index below ROWS: usually ROWS - 1, now and then smaller."
  (if (zerop (rnd 2)) (1- rows) (rnd rows)))

(define-oracle-test c-polygon-from-c-pts (:file "frame-basic" :trials 2000)
  (let* ((rows (1+ (rnd 15))) (pw (%frame-basic-cpoints rows)) (n (%frame-basic-highest rows)))
    (agree (list :pw pw :n n)
           (list (batten::c-polygon-from-c-pts pw n))
           (with-foreign ((a pw :cpoint-array))
             (with-foreign-out ((pol :cpolygon))
               (c-c-polygon-from-c-pts pol a n)
               (values nil (list (from-foreign pol :cpolygon))))))))

(define-oracle-test e-polygon-from-pts (:file "frame-basic" :trials 2000)
  (let* ((rows (1+ (rnd 15))) (p (random-points rows)) (n (%frame-basic-highest rows)))
    (agree (list :n n :p p)
           (list (batten::e-polygon-from-pts n p))
           (with-foreign ((a p :point-array))
             (with-foreign-out ((ppl :epolygon))
               (c-e-polygon-from-pts ppl n a)
               (values nil (list (from-foreign ppl :epolygon))))))))

(define-oracle-test e-polygon-get-pts (:file "frame-basic" :trials 2000)
  (let ((p (random-points (1+ (rnd 15)))))
    (agree (list :p p)
           (multiple-value-list (batten::e-polygon-get-pts p))
           (with-foreign ((ppl p :epolygon))
             (cffi:with-foreign-objects ((nn :long) (pp :pointer))
               (c-e-polygon-get-pts ppl nn pp)
               (let ((n (cffi:mem-ref nn :long)))
                 (values nil (list n (from-foreign (cffi:mem-ref pp :pointer) :point-array (1+ n))))))))))

(define-oracle-test c-net-from-c-pts (:file "frame-basic" :trials 2000)
  (let* ((rows (1+ (rnd 8))) (cols (1+ (rnd 8)))
         (pw (if (zerop (rnd 2))
                 (%frame-basic-net rows cols 4)
                 (let ((s (random-surface :two-d :random)))
                   (setf rows (1+ (batten::surface-n s)) cols (1+ (batten::surface-m s)))
                   (batten::surface-cpoints s))))
         (n (%frame-basic-highest rows)) (m (%frame-basic-highest cols)))
    (agree (list :pw pw :n n :m m)
           (list (batten::c-net-from-c-pts pw n m))
           (with-foreign ((a pw :cpoint-net))
             ;; marshal.lisp has no empty :cnet; allocate the NL_CNET here (it is
             ;; freed with the WITH-FOREIGN allocations) and read its fields back
             (let ((cn (batten-oracle::falloc '(:struct batten-oracle::nl-cnet) 1)))
               (c-c-net-from-c-pts cn a n m)
               (flet ((slot (name) (cffi:foreign-slot-value cn '(:struct batten-oracle::nl-cnet) name)))
                 (values nil (list (from-foreign (slot 'batten-oracle::pw) :cpoint-net
                                                 (1+ (slot 'batten-oracle::n))
                                                 (1+ (slot 'batten-oracle::m)))))))))))

(define-oracle-test e-net-from-pts (:file "frame-basic" :trials 2000)
  (let* ((rows (1+ (rnd 8))) (cols (1+ (rnd 8)))
         (p (%frame-basic-net rows cols 3))
         (n (%frame-basic-highest rows)) (m (%frame-basic-highest cols)))
    (agree (list :n n :m m :p p)
           (list (batten::e-net-from-pts n m p))
           (with-foreign ((a p :point-net))
             (with-foreign-out ((ntl :enet))
               (c-e-net-from-pts ntl n m a)
               (values nil (list (from-foreign ntl :enet))))))))

(define-oracle-test e-net-get-pts (:file "frame-basic" :trials 2000)
  (let ((p (%frame-basic-net (1+ (rnd 8)) (1+ (rnd 8)) 3)))
    (agree (list :p p)
           (multiple-value-list (batten::e-net-get-pts p))
           (with-foreign ((ntl p :enet))
             (cffi:with-foreign-objects ((nn :long) (mm :long) (pp :pointer))
               (c-e-net-get-pts ntl nn mm pp)
               (let ((n (cffi:mem-ref nn :long)) (m (cffi:mem-ref mm :long)))
                 (values nil (list n m (from-foreign (cffi:mem-ref pp :pointer) :point-net (1+ n) (1+ m))))))))))

;;; ------------------------------------------------------------ level 1, group 1

(defun %frame-basic-signed-zeros (a)
  "Now and then set some coordinates of point array or net A to 0 or -0, so
min-max comparisons meet ties."
  (when (zerop (rnd 3))
    (dotimes (i (array-total-size a))
      (when (zerop (rnd 4))
        (setf (row-major-aref a i) (if (zerop (rnd 2)) 0d0 -0d0)))))
  a)

(define-oracle-test polygon-get-b-box (:file "frame-basic" :trials 3000)
  (let ((p (%frame-basic-signed-zeros (random-points (1+ (rnd 15))))))
    (agree (list :p p)
           (list (batten::polygon-get-b-box p (make-real-vector 6)))
           (with-foreign ((ppl p :epolygon))
             (with-foreign-out ((box :box))
               (c-polygon-get-b-box ppl box)
               (values nil (list (from-foreign box :box))))))))

(define-oracle-test e-net-get-b-box (:file "frame-basic" :trials 3000)
  (let ((p (%frame-basic-signed-zeros (%frame-basic-net (1+ (rnd 8)) (1+ (rnd 8)) 3))))
    (agree (list :p p)
           (list (batten::e-net-get-b-box p (make-real-vector 6)))
           (with-foreign ((ntl p :enet))
             (with-foreign-out ((box :box))
               (c-e-net-get-b-box ntl box)
               (values nil (list (from-foreign box :box))))))))

(defun %frame-basic-closed-net (rows cols)
  "A random point net, now and then closed (give or take about *mtol*) in u
or in v."
  (let ((p (%frame-basic-net rows cols 3)))
    (case (rnd 3)
      (0 (dotimes (j cols)
           (dotimes (c 3)
             (setf (aref p (1- rows) j c)
                   (+ (aref p 0 j c) (case (rnd 4) (0 0d0) (1 3d-8) (2 6d-8) (t 2d-7)))))))
      (1 (dotimes (i rows)
           (dotimes (c 3)
             (setf (aref p i (1- cols) c)
                   (+ (aref p i 0 c) (case (rnd 4) (0 0d0) (1 3d-8) (2 6d-8) (t 2d-7))))))))
    p))

(define-oracle-test e-net-is-closed (:file "frame-basic" :trials 3000)
  (let ((p (%frame-basic-closed-net (1+ (rnd 8)) (1+ (rnd 8))))
        (dir (case (rnd 8) (0 0) (1 3) ((2 3 4) +udir+) (t +vdir+))))
    (agree (list :p p :dir dir)
           (list (if (batten::e-net-is-closed p dir) 1 0))
           (with-foreign ((ntl p :enet))
             (values nil (list (c-e-net-is-closed ntl dir)))))))
