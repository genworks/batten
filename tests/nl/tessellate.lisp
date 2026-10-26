;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/tessellate.lisp against NLib.

(in-package :batten-tests)

;;; ------------------------------------------------------------- L14-g0

(defun %tess-edge-count (dt hd m)
  "How many edges N_TessGetTriangleEdges gathers from DT, HD."
  (let ((c 0))
    (dotimes (i (1+ m) c)
      (let ((row (svref dt i)))
        (unless (= (aref row 0) -1)
          (dotimes (r (1+ (aref hd i)))
            (when (< i (aref row r)) (incf c))))))))

(defun %tess-triangle-count (dt hd m)
  "How many triangles N_TessGetTriangles gathers from DT, HD."
  (let ((c 0))
    (dotimes (i (1+ m) c)
      (let ((row (svref dt i)))
        (unless (= (aref row 0) -1)
          (dotimes (r (aref hd i))
            (when (and (< i (aref row r)) (< i (aref row (1+ r)))) (incf c))))))))

(defun %tess-database (kind)
  "A triangulation database (list dt hd m) inside the C's domain for KIND
(:edges or :triangles): a triangulated grid most of the time, random index
lists the rest, kept within the arrays the C gathers into (4m+1 edges; 2m+2
triangles, the last of which overwrites the first's second and third points,
which the lists reach now and then)."
  (case (rnd 12)
    ((0 1 2 3 4 5 6)
     (multiple-value-bind (dt hd m) (random-triangulation :unused (if (zerop (rnd 3)) nil 8))
       (list dt hd m)))
    (7
     ;; exactly at capacity: point 0's list alone gives 4m+1 edges or 2m+2
     ;; triangles; every other point out of the triangulation, or not
     (let* ((m (1+ (rnd 6)))
            (len (if (eq kind :edges) (1+ (* 4 m)) (+ 3 (* 2 m))))
            (dt (make-array (1+ m))) (hd (make-array (1+ m) :element-type 'fixnum)))
       (setf (svref dt 0) (let ((row (make-array len :element-type 'fixnum)))
                            (dotimes (r len row) (setf (aref row r) (1+ (rnd m)))))
             (aref hd 0) (1- len))
       (do ((i 1 (1+ i))) ((> i m))
         (setf (svref dt i) (make-array 1 :element-type 'fixnum :initial-element (if (zerop (rnd 2)) -1 0))
               (aref hd i) 0))
       (list dt hd m)))
    (t
      (do ((db nil)) (db db)
        (let ((m (rnd 8)))
          (multiple-value-bind (dt hd) (random-index-lists m :max-length (+ 1 (rnd 8)))
            (when (if (eq kind :edges)
                      (<= (%tess-edge-count dt hd m) (1+ (* 4 m)))
                      (<= (%tess-triangle-count dt hd m) (+ 2 (* 2 m))))
              (setf db (list dt hd m)))))))))

(define-oracle-test tess-get-triangle-edges (:file "tessellate" :trials 3000)
  (destructuring-bind (dt hd m) (%tess-database :edges)
    (agree (list :dt dt :hd hd :m m)
           (multiple-value-list (batten::tess-get-triangle-edges dt hd m))
           (with-foreign ((dp dt :index-lists) (hp hd :index-vector))
             (with-stacks (sg)
               (with-scalars ((alf :pointer) (bet :pointer) (n :long))
                 (let ((err (c-flag-error-p (c-tess-get-triangle-edges dp hp m alf bet n sg))))
                   (values err
                           (unless err
                             (let ((n (cffi:mem-ref n :long)))
                               (list (from-foreign (cffi:mem-ref alf :pointer) :index-vector (1+ n))
                                     (from-foreign (cffi:mem-ref bet :pointer) :index-vector (1+ n))
                                     n)))))))))))

(define-oracle-test tess-get-triangles (:file "tessellate" :trials 3000)
  (destructuring-bind (dt hd m) (%tess-database :triangles)
    (agree (list :dt dt :hd hd :m m)
           (multiple-value-list (batten::tess-get-triangles dt hd m))
           (with-foreign ((dp dt :index-lists) (hp hd :index-vector))
             (with-stacks (sg)
               (with-scalars ((alf :pointer) (bet :pointer) (gam :pointer) (n :long))
                 (let ((err (c-flag-error-p (c-tess-get-triangles dp hp m alf bet gam n sg))))
                   (values err
                           (unless err
                             (let ((n (cffi:mem-ref n :long)))
                               (list (from-foreign (cffi:mem-ref alf :pointer) :index-vector (1+ n))
                                     (from-foreign (cffi:mem-ref bet :pointer) :index-vector (1+ n))
                                     (from-foreign (cffi:mem-ref gam :pointer) :index-vector (1+ n))
                                     n)))))))))))

(define-oracle-test tess-srf-area (:file "tessellate" :trials 3000)
  (let ((s (random-surface :two-d (zerop (rnd 5)))))
    ;; now and then a column or row collapsed to a point
    (when (zerop (rnd 6))
      (let ((net (surface-cpoints s)))
        (if (zerop (rnd 2))
            (dotimes (j (array-dimension net 1))
              (dotimes (c 4) (setf (aref net 0 j c) (aref net 0 0 c))))
            (dotimes (i (array-dimension net 0))
              (dotimes (c 4) (setf (aref net i 0 c) (aref net 0 0 c)))))))
    (agree (list :surface s)
           (multiple-value-list (batten::tess-srf-area s))
           (with-foreign ((sp s :surface))
             (with-scalars ((a :double) (lu :double) (lv :double))
               (let ((err (c-flag-error-p (c-tess-srf-area sp a lu lv))))
                 (values err (unless err (list (cffi:mem-ref a :double) (cffi:mem-ref lu :double)
                                               (cffi:mem-ref lv :double))))))))))
