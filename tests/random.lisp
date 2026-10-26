;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Random curves, surfaces, functions and points for the differential tests
;;;; and benchmarks, from a seeded generator of our own (splitmix64), so every
;;;; implementation draws the same ones from the same seed.

(defpackage :batten-random
  (:use :common-lisp :batten)
  (:export #:*seed* #:reseed #:rnd #:rnd-range
           #:random-knots #:random-curve #:random-surface #:random-cfun #:random-sfun
           #:random-parameter #:random-point #:random-cpoint #:random-points
           #:random-real-vector #:random-real-matrix
           #:random-triangulation #:random-index-lists
           #:evaluation-digest))

(in-package :batten-random)

(defvar *seed* 20261025)

(defun reseed (&optional (seed 20261025)) (setf *seed* (ldb (byte 64 0) seed)))

(defun next-u64 ()
  (setf *seed* (ldb (byte 64 0) (+ *seed* #x9E3779B97F4A7C15)))
  (let ((z *seed*))
    (setf z (ldb (byte 64 0) (* (logxor z (ash z -30)) #xBF58476D1CE4E5B9))
          z (ldb (byte 64 0) (* (logxor z (ash z -27)) #x94D049BB133111EB)))
    (logxor z (ash z -31))))

(defun rnd (limit)
  "A random integer below an integer LIMIT, or a double in [0, LIMIT)."
  (if (integerp limit)
      (mod (next-u64) limit)
      (* limit (/ (float (ash (next-u64) -11) 1d0) 9007199254740992d0))))

(defun rnd-range (lo hi) "A double in [LO, HI)." (+ lo (rnd (- hi lo))))

(defun random-knots (p n+1 &key (lo 0d0) (hi 1d0))
  "A clamped knot vector for N+1 control points of degree P: P+1 equal knots at
each end, interior knots random, now and then repeated up to multiplicity P."
  (let* ((interior-count (- n+1 p 1))
         (interior (sort (let (ks)
                           (dotimes (i interior-count ks)
                             (push (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0)))) ks)))
                         #'<))
         (v (coerce interior 'vector)))
    (dotimes (i (length v))
      (when (and (> i 0) (zerop (rnd 6)))
        (let ((mult 0))
          (do ((j (1- i) (1- j))) ((or (< j 0) (/= (aref v j) (aref v (1- i)))))
            (incf mult))
          (when (< mult p)
            (setf (aref v i) (aref v (1- i)))))))
    (coerce (append (make-list (1+ p) :initial-element lo)
                    (coerce v 'list)
                    (make-list (1+ p) :initial-element hi))
            'real-vector)))

(defun random-domain ()
  (let* ((lo (if (zerop (rnd 3)) 0d0 (- (rnd 20d0) 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0))))
    (values lo hi)))

(defun random-curve (&key (max-degree 7) (max-extra 12) (rational :random) (two-d :random))
  "A random NURBS curve: degree 1..MAX-DEGREE, rational half the time and 2D a
fifth of the time unless RATIONAL / TWO-D say otherwise (t or nil)."
  (multiple-value-bind (lo hi) (random-domain)
    (let* ((p (1+ (rnd max-degree)))
           (n+1 (+ p 1 (rnd max-extra)))
           (rational (if (eq rational :random) (zerop (rnd 2)) rational))
           (two-d (if (eq two-d :random) (zerop (rnd 5)) two-d))
           (cpoints (make-cpoint-array n+1)))
      (dotimes (i n+1)
        (let ((w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
          (setf (aref cpoints i 0) (* w (- (rnd 200d0) 100d0))
                (aref cpoints i 1) (* w (- (rnd 200d0) 100d0))
                (aref cpoints i 2) (if two-d +no-z+ (* w (- (rnd 200d0) 100d0)))
                (aref cpoints i 3) (if rational w +no-w+))))
      (%make-curve p (random-knots p n+1 :lo lo :hi hi) cpoints))))

(defun random-surface (&key (max-degree 5) (max-extra 6) (rational :random) (two-d nil))
  "A random NURBS surface: degrees 1..MAX-DEGREE in u and v, rational half the
time; 3D unless TWO-D."
  (multiple-value-bind (ulo uhi) (random-domain)
    (multiple-value-bind (vlo vhi) (random-domain)
      (let* ((p (1+ (rnd max-degree))) (q (1+ (rnd max-degree)))
             (n+1 (+ p 1 (rnd max-extra))) (m+1 (+ q 1 (rnd max-extra)))
             (rational (if (eq rational :random) (zerop (rnd 2)) rational))
             (two-d (if (eq two-d :random) (zerop (rnd 5)) two-d))
             (net (make-cpoint-net n+1 m+1)))
        (dotimes (i n+1)
          (dotimes (j m+1)
            (let ((w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
              (setf (aref net i j 0) (* w (- (rnd 200d0) 100d0))
                    (aref net i j 1) (* w (- (rnd 200d0) 100d0))
                    (aref net i j 2) (if two-d +no-z+ (* w (- (rnd 200d0) 100d0)))
                    (aref net i j 3) (if rational w +no-w+)))))
        (%make-surface p q (random-knots p n+1 :lo ulo :hi uhi)
                       (random-knots q m+1 :lo vlo :hi vhi) net)))))

(defun random-cfun (&key (max-degree 7) (max-extra 12))
  "A random scalar B-spline function."
  (multiple-value-bind (lo hi) (random-domain)
    (let* ((p (1+ (rnd max-degree))) (n+1 (+ p 1 (rnd max-extra)))
           (vals (make-real-vector n+1)))
      (dotimes (i n+1) (setf (aref vals i) (- (rnd 200d0) 100d0)))
      (%make-cfun p (random-knots p n+1 :lo lo :hi hi) vals))))

(defun random-sfun (&key (max-degree 5) (max-extra 6))
  "A random scalar tensor-product function."
  (multiple-value-bind (ulo uhi) (random-domain)
    (multiple-value-bind (vlo vhi) (random-domain)
      (let* ((p (1+ (rnd max-degree))) (q (1+ (rnd max-degree)))
             (n+1 (+ p 1 (rnd max-extra))) (m+1 (+ q 1 (rnd max-extra)))
             (vals (make-real-matrix n+1 m+1)))
        (dotimes (i n+1) (dotimes (j m+1) (setf (aref vals i j) (- (rnd 200d0) 100d0))))
        (%make-sfun p q (random-knots p n+1 :lo ulo :hi uhi) (random-knots q m+1 :lo vlo :hi vhi) vals)))))

(defun random-parameter (knots)
  "A parameter in the knot vector's domain: usually interior, now and then
exactly a knot or an end."
  (let ((lo (aref knots 0))
        (hi (aref knots (1- (length knots)))))
    (case (rnd 8)
      (0 lo)
      (1 hi)
      (2 (aref knots (rnd (length knots))))
      (t (+ lo (rnd (- hi lo)))))))

(defun random-point (&optional (scale 100d0))
  (make-point (- (rnd (* 2 scale)) scale) (- (rnd (* 2 scale)) scale) (- (rnd (* 2 scale)) scale)))

(defun random-cpoint (&key (rational :random) (two-d nil) (scale 100d0))
  (let* ((rational (if (eq rational :random) (zerop (rnd 2)) rational))
         (w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
    (make-cpoint (* w (- (rnd (* 2 scale)) scale)) (* w (- (rnd (* 2 scale)) scale))
                 (if two-d +no-z+ (* w (- (rnd (* 2 scale)) scale)))
                 (if rational w +no-w+))))

(defun random-points (n &optional (scale 100d0))
  "A point-array of N random points."
  (let ((a (make-point-array n)))
    (dotimes (i n a)
      (dotimes (k 3) (setf (aref a i k) (- (rnd (* 2 scale)) scale))))))

(defun random-real-vector (n &optional (scale 100d0))
  (let ((a (make-real-vector n)))
    (dotimes (i n a) (setf (aref a i) (- (rnd (* 2 scale)) scale)))))

(defun random-real-matrix (rows cols &optional (scale 100d0))
  (let ((a (make-real-matrix rows cols)))
    (dotimes (i rows a) (dotimes (j cols) (setf (aref a i j) (- (rnd (* 2 scale)) scale))))))

;;; NLib's tessellator keeps a triangulation as NL_INDEX **DT and NL_INDEX
;;; *hd: DT[i] lists the neighbors of point i counterclockwise -- closed by
;;; repeating the first for a point inside, from one boundary neighbor to the
;;; other for a point on the boundary -- so that consecutive neighbors make a
;;; triangle with i; DT[i][0] = -1 marks a point not in the triangulation;
;;; hd[i] is the highest index in DT[i].  Batten: a simple-vector of
;;; index-vectors and an index-vector.

(defun random-triangulation (&key (max-cells 5) (unused 8))
  "A random triangulation database (values dt hd m) over points 0..M: the
grid of points of an a by b cell rectangle (1 <= a, b <= MAX-CELLS), each
cell cut by a random diagonal, its points labeled in a random order; about
one point in UNUSED (an integer, or NIL for none) marked not in the
triangulation (DT[i] = #(-1))."
  (let* ((a (1+ (rnd max-cells))) (b (1+ (rnd max-cells)))
         (count (* (1+ a) (1+ b)))
         (label (let ((v (make-array count)))
                  (dotimes (i count) (setf (svref v i) i))
                  (do ((i (1- count) (1- i))) ((< i 1) v)   ; Fisher-Yates
                    (rotatef (svref v i) (svref v (rnd (1+ i)))))))
         (next (make-array count :initial-element nil))) ; per point: ((p . q) ...)
    (flet ((pt (r c) (svref label (+ (* r (1+ b)) c)))
           (tri (x y z)                     ; counterclockwise
             (push (cons y z) (svref next x))
             (push (cons z x) (svref next y))
             (push (cons x y) (svref next z))))
      (dotimes (r a)
        (dotimes (c b)
          (let ((v00 (pt r c)) (v01 (pt r (1+ c))) (v10 (pt (1+ r) c)) (v11 (pt (1+ r) (1+ c))))
            (if (zerop (rnd 2))
                (progn (tri v00 v01 v11) (tri v00 v11 v10))
                (progn (tri v00 v01 v10) (tri v01 v11 v10)))))))
    (let ((dt (make-array count)) (hd (make-array count :element-type 'fixnum)))
      (dotimes (v count)
        (let* ((pairs (svref next v))
               (start (or (car (find-if (lambda (pq) (not (find (car pq) pairs :key #'cdr))) pairs))
                          (car (first pairs))))
               (chain (list start)))
          (do ((p start (cdr (assoc p pairs))))
              ((null (assoc p pairs)))
            (push (cdr (assoc p pairs)) chain)
            (when (= (cdr (assoc p pairs)) start) (return)))
          (let ((row (if (and unused (zerop (rnd unused)))
                         (list -1)
                         (nreverse chain))))
            (setf (svref dt v) (make-array (length row) :element-type 'fixnum :initial-contents row)
                  (aref hd v) (1- (length row))))))
      (values dt hd (1- count)))))

(defun random-index-lists (m &key (max-length 6) (unused 6))
  "M+1 random index lists in the shape of a triangulation database (values
dt hd): list i of 1..MAX-LENGTH indexes in 0..M (now and then -1), or #(-1)
about one time in UNUSED; HD[i] its highest index."
  (let ((dt (make-array (1+ m))) (hd (make-array (1+ m) :element-type 'fixnum)))
    (dotimes (i (1+ m) (values dt hd))
      (let* ((len (if (zerop (rnd unused)) 1 (1+ (rnd max-length))))
             (row (make-array len :element-type 'fixnum)))
        (dotimes (r len)
          (setf (aref row r) (if (zerop (rnd 15)) -1 (rnd (1+ m)))))
        (when (= len 1) (setf (aref row 0) (if (zerop (rnd 2)) -1 (aref row 0))))
        (setf (svref dt i) row (aref hd i) (1- len))))))
