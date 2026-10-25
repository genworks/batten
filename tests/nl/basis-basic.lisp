;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/basis-basic.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test basis-find-span (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)))
    (agree (list k p u flg)
           (list (batten::basis-find-span k p u flg))
           (with-foreign ((kp k :knots))
             (with-scalars ((spn :long))
               (values (c-flag-error-p (c-basis-find-span kp p u flg spn))
                       (list (cffi:mem-ref spn :long))))))))

(define-oracle-test basis-eval (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)))
    (agree (list k p u flg)
           (let ((n (make-real-vector (1+ p))))
             (list (batten::basis-eval k p u flg n) n))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((n :real-vector (1+ p)))
               (with-scalars ((spn :long))
                 (values (c-flag-error-p (c-basis-eval kp p u flg n spn))
                         (list (cffi:mem-ref spn :long) (from-foreign n :real-vector (1+ p))))))))))

(define-oracle-test basis-derivs (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (random-side)) (der (rnd (1+ +maxder+))))
    (agree (list k p u flg der)
           (let ((nd (make-real-matrix (1+ der) (1+ p))))
             (list (batten::basis-derivs k p u flg der nd) nd))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((nd :real-matrix (1+ der) (1+ p)))
               (with-scalars ((spn :long))
                 (values (c-flag-error-p (c-basis-derivs kp p u flg der nd spn))
                         (list (cffi:mem-ref spn :long) (from-foreign nd :real-matrix (1+ der) (1+ p))))))))))

(define-oracle-test basis-get-span-count (:file "basis-basic" :trials 2000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)))
    (agree (list k p)
           (list (batten::basis-get-span-count k p))
           (with-foreign ((kp k :knots))
             (with-scalars ((nsp :long))
               (c-basis-get-span-count kp p nsp)
               (values nil (list (cffi:mem-ref nsp :long))))))))

(define-oracle-test knot-vector-is-param-out-of-bounds (:file "basis-basic" :trials 2000)
  (let* ((k (curve-knots (random-curve)))
         (span (- (aref k (1- (length k))) (aref k 0)))
         (u (+ (aref k 0) (* span (rnd-range -0.2d0 1.2d0)))))
    (agree (list k u)
           (progn (batten::knot-vector-is-param-out-of-bounds k u 'test) (list t))
           (with-foreign ((kp k :knots))
             (values (c-flag-error-p (c-knot-vector-is-param-out-of-bounds kp u (cffi:null-pointer)))
                     (list t))))))

(define-oracle-test knot-vector-is-end-param (:file "basis-basic" :trials 2000)
  (let* ((k (curve-knots (random-curve)))
         (span (- (aref k (1- (length k))) (aref k 0)))
         (u (case (rnd 4)
              (0 (aref k 0))
              (1 (aref k (1- (length k))))
              (t (+ (aref k 0) (* span (rnd-range -0.2d0 1.2d0)))))))
    (agree (list k u)
           (progn (batten::knot-vector-is-end-param k u 'test) (list t))
           (with-foreign ((kp k :knots))
             (values (c-flag-error-p (c-knot-vector-is-end-param kp u (cffi:null-pointer)))
                     (list t))))))

;;; ------------------------------------------------------------- L0-g2

(defun %spoil-knots (k p)
  "A copy of knot vector K, now and then made invalid for degree P: an end
multiplicity changed, a knot moved out of order, or an interior knot
repeated past P."
  (let* ((k (copy-seq k)) (m (1- (length k))))
    (case (rnd 10)
      (0 (setf (aref k p) (aref k (1+ p))))                       ; start multiplicity off
      (1 (setf (aref k (- m p)) (aref k (- m p 1))))              ; end multiplicity off
      (2 (when (> (- m p) (1+ p))                                 ; out of order
           (let ((i (+ p 1 (rnd (- m p p 1)))))
             (setf (aref k i) (+ (aref k i) (* 0.5d0 (- (aref k m) (aref k 0))))))))
      (3 (when (> (- m p) (1+ p))                                 ; interior multiplicity p+1
           (let ((i (+ p 1 (rnd (- m p p 1)))))
             (do ((j i (1+ j)) (c 0 (1+ c))) ((or (> c p) (>= j (- m p))))
               (setf (aref k j) (aref k i)))))))
    k))

(define-oracle-test knot-vector-is-valid (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (p0 (curve-degree c))
         (k (%spoil-knots (curve-knots c) p0))
         (p (case (rnd 16) (0 (1- p0)) (1 (1+ p0)) (t p0))))
    (agree (list k p)
           (progn (batten::knot-vector-is-valid k p 'test) (list t))
           (with-foreign ((kp k :knots))
             (values (c-flag-error-p (c-knot-vector-is-valid kp p (cffi:null-pointer)))
                     (list t))))))

(define-oracle-test knot-vector-from-real-array (:file "basis-basic" :trials 1000)
  (let* ((len (1+ (rnd 20))) (u (random-real-vector len)) (m (rnd len)))
    (agree (list u m)
           (list (batten::knot-vector-from-real-array u m))
           (with-foreign ((up u :real-vector))
             (with-foreign-out ((kp :knots))
               (c-knot-vector-from-real-array kp up m)
               (values nil (list (from-foreign kp :knots))))))))

(define-oracle-test knots-check (:file "basis-basic" :trials 2000)
  (let* ((fresh (zerop (rnd 4)))
         (k (if fresh nil (random-real-vector (1+ (rnd 20)))))
         (m (if fresh (rnd 20) (rnd (+ (length k) 2)))))
    ;; NLib's allocator leaves a fresh knot array uninitialized: compare only
    ;; its size there
    (agree (list k m)
           (let ((kv (batten::knots-check k m 'test)))
             (list (if fresh (length kv) kv)))
           (with-stacks (s)
             (cffi:with-foreign-object (kpp :pointer)
               (with-foreign ((kp (or k (make-real-vector 1)) :knots))
                 (setf (cffi:mem-ref kpp :pointer) (if fresh (cffi:null-pointer) kp))
                 (let ((err (c-flag-error-p (c-knots-check kpp m (cffi:null-pointer) s))))
                   (values err (unless err
                                 (let ((kv (from-foreign (cffi:mem-ref kpp :pointer) :knots)))
                                   (list (if fresh (length kv) kv))))))))))))

(defun %flag () (if (zerop (rnd 40)) 3 (random-side)))

(define-oracle-test basis-find-span-and-mult (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (u (random-parameter k)) (flg (%flag)))
    (agree (list k p u flg)
           (multiple-value-list (batten::basis-find-span-and-mult k p u flg))
           (with-foreign ((kp k :knots))
             (with-scalars ((spn :long) (mlt :long))
               (values (c-flag-error-p (c-basis-find-span-and-mult kp p u flg spn mlt))
                       (list (cffi:mem-ref spn :long) (cffi:mem-ref mlt :long))))))))

(define-oracle-test basis-find-index-node-array (:file "basis-basic" :trials 3000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (len (- (length k) 1 p)))
    (agree (list k p)
           (list (batten::basis-find-index-node-array k p (make-real-vector len)))
           (with-foreign ((kp k :knots))
             (with-foreign-out ((tt :real-vector len))
               (c-basis-find-index-node-array kp p tt)
               (values nil (list (from-foreign tt :real-vector len))))))))

(defun %index-and-parameter (k p)
  "A basis function index (now and then out of range) and a parameter (now and
then outside the domain) for knot vector K of degree P."
  (let* ((n (- (length k) p 2))
         (i (case (rnd 30) (0 -1) (1 (1+ n)) (t (rnd (1+ n)))))
         (u (if (zerop (rnd 30))
                (+ (aref k (1- (length k))) 0.25d0)
                (random-parameter k))))
    (values i u)))

(define-oracle-test basis-i-eval (:file "basis-basic" :trials 10000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)) (flg (%flag)))
    (multiple-value-bind (i u) (%index-and-parameter k p)
      (agree (list k i p u flg)
             (list (batten::basis-i-eval k i p u flg))
             (with-foreign ((kp k :knots))
               (with-scalars ((n :double))
                 (values (c-flag-error-p (c-basis-i-eval kp i p u flg n))
                         (list (cffi:mem-ref n :double)))))))))

(define-oracle-test basis-i-derivs (:file "basis-basic" :trials 10000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c)) (flg (%flag))
         (der (rnd (+ p 3))))
    (multiple-value-bind (i u) (%index-and-parameter k p)
      (agree (list k i p u flg der)
             (list (batten::basis-i-derivs k i p u flg der (make-real-vector (1+ der))))
             (with-foreign ((kp k :knots))
               (with-foreign-out ((nd :real-vector (1+ der)))
                 (values (c-flag-error-p (c-basis-i-derivs kp i p u flg der nd))
                         (list (from-foreign nd :real-vector (1+ der))))))))))

;;; ------------------------------------------------------------- L1-g5

(define-oracle-test knots-copy (:file "basis-basic" :trials 2000)
  (let* ((kna (curve-knots (random-curve)))
         (ma (1- (length kna)))
         (fresh (zerop (rnd 4)))
         ;; KNB now and then too short
         (knb (if fresh nil (random-real-vector (max 1 (+ ma (rnd 6) -1))))))
    (agree (list kna knb)
           (list (batten::knots-copy kna (and knb (copy-seq knb))))
           (with-stacks (s)
             (cffi:with-foreign-object (kpp :pointer)
               (with-foreign ((ka kna :knots) (kb (or knb (make-real-vector 1)) :knots))
                 (setf (cffi:mem-ref kpp :pointer) (if fresh (cffi:null-pointer) kb))
                 (let ((err (c-flag-error-p (c-knots-copy ka kpp s))))
                   (values err (unless err
                                 (list (from-foreign (cffi:mem-ref kpp :pointer) :knots)))))))))))

(define-oracle-test basis-find-global-max (:file "basis-basic" :trials 4000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (n (- (length k) p 2))
         (i (case (rnd 40) (0 -1) (1 (1+ n)) (t (rnd (1+ n)))))
         (tol (case (rnd 4) (0 1d-6) (1 1d-10) (2 1d-12) (t 1d-8))))
    (agree (list k i p tol)
           (multiple-value-list (batten::basis-find-global-max k i p tol))
           (with-foreign ((kp k :knots))
             (with-scalars ((mx :double) (u :double))
               (values (c-flag-error-p (c-basis-find-global-max kp i p tol mx u))
                       (list (cffi:mem-ref mx :double) (cffi:mem-ref u :double))))))))

(defun %basis-interval (a b)
  (let ((v (make-real-vector 2))) (setf (aref v 0) a (aref v 1) b) v))

(define-oracle-test basis-reparam (:file "basis-basic" :trials 3000)
  (let* ((c (random-curve)) (k (curve-knots c)) (p (curve-degree c))
         (lo (aref k 0)) (hi (aref k (1- (length k))))
         (iv (case (rnd 6)
               (0 (%basis-interval lo hi))
               (1 (%basis-interval lo (+ hi 1d0)))
               (2 (%basis-interval (- lo 0.5d0) hi))
               (t (let ((a (rnd-range -10d0 10d0)))
                    (%basis-interval a (+ a 0.1d0 (rnd 20d0))))))))
    (agree (list k p iv)
           (list (batten::basis-reparam (copy-seq k) p iv))
           (with-foreign ((kp k :knots) (ip iv :interval))
             (c-basis-reparam kp p ip)
             (values nil (list (from-foreign kp :knots)))))))

;;; ------------------------------------------------------------- L2-g3

(defun %p-fold-knots (k p)
  "A copy of knot vector K (degree P) with, now and then, an interior knot
repeated to multiplicity P (when K has room), so that some basis function has
a P-fold inner knot.  Second value: the index of the basis function whose
inner knots are that P-fold knot, or nil."
  (let* ((k (copy-seq k)) (m (1- (length k))) (room (- m p p 1)) (at nil))
    (when (and (>= room p) (zerop (rnd 3)))
      (let ((i (+ p 1 (rnd (- room p -1)))))
        (dotimes (c p)
          (setf (aref k (+ i c)) (aref k i)))
        (setf at (1- i))))
    (values k at)))

(define-oracle-test basis-find-all-span-maxima (:file "basis-basic" :trials 3000)
  (let* ((c (random-curve)) (p (curve-degree c))
         (n (- (length (curve-knots c)) p 2)))
    (multiple-value-bind (k at) (%p-fold-knots (curve-knots c) p)
      (let ((i (if (and at (zerop (rnd 2)))
                   at
                   (case (rnd 40) (0 -1) (1 (1+ n)) (t (rnd (1+ n))))))
            (tol (case (rnd 4) (0 1d-6) (1 1d-10) (2 1d-12) (t 1d-8))))
        (agree (list k i p tol)
               (let ((mn (make-real-vector (1+ p))) (mx (make-real-vector (1+ p))))
                 (let ((u (batten::basis-find-all-span-maxima k i p tol mn mx)))
                   (list u mn mx)))
               (with-foreign ((kp k :knots))
                 (with-foreign-out ((mn :real-vector (1+ p)) (mx :real-vector (1+ p)))
                   (with-scalars ((u :double))
                     (values (c-flag-error-p (c-basis-find-all-span-maxima kp i p tol mn mx u))
                             (list (cffi:mem-ref u :double)
                                   (from-foreign mn :real-vector (1+ p))
                                   (from-foreign mx :real-vector (1+ p))))))))))))

(defun %compatible-knots (p lo hi &optional other)
  "A random clamped degree-P knot vector on [LO, HI], interior knots repeated
now and then (multiplicity at most P); when OTHER (a knot vector on the same
domain) is given, some of its interior knots are reused."
  (let* ((oint (when other
                 (remove-if (lambda (x) (or (<= x lo) (>= x hi))) (coerce other 'list))))
         (int nil))
    (dotimes (c (rnd 10))
      (push (cond ((and oint (zerop (rnd 2))) (nth (rnd (length oint)) oint))
                  ((and int (zerop (rnd 5))) (first int))
                  (t (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0))))))
            int))
    (setf int (sort int #'<))
    (let ((kept nil))
      (dolist (x int)
        (when (< (count x kept :test #'=) p) (push x kept)))
      (coerce (append (make-list (1+ p) :initial-element lo)
                      (nreverse kept)
                      (make-list (1+ p) :initial-element hi))
              'real-vector))))

(define-oracle-test make-knots-compatible (:file "basis-basic" :trials 3000)
  (let* ((p (1+ (rnd 5))) (q (1+ (rnd 5)))
         (lo (if (zerop (rnd 3)) 0d0 (rnd-range -5d0 5d0)))
         (hi (+ lo 0.5d0 (rnd 5d0)))
         (knr (%compatible-knots p lo hi))
         (kns (case (rnd 4)
                ;; another domain: one is rescaled
                (0 (let ((l2 (rnd-range -5d0 5d0))) (%compatible-knots q l2 (+ l2 0.5d0 (rnd 5d0)))))
                ;; ends a hair apart, within *ptol*
                (1 (let ((k (%compatible-knots q lo hi knr))
                         (e (* 0.25d-12 (rnd 1d0))))
                     (dotimes (j (1+ q)) (setf (aref k (- (length k) 1 j)) (+ hi e)))
                     k))
                (t (%compatible-knots q lo hi knr)))))
    (agree (list knr p kns q)
           (let ((r (copy-seq knr)) (s (copy-seq kns)))
             (multiple-value-bind (xr xs) (batten::make-knots-compatible r s p q)
               (list r s xr xs)))
           (let ((cap (+ (length knr) (length kns) 2)))
             (with-foreign ((rp knr :knots) (sp kns :knots)
                            (xrp (make-real-vector cap) :knots) (xsp (make-real-vector cap) :knots))
               (c-make-knots-compatible rp sp p q xrp xsp)
               (values nil (list (from-foreign rp :knots) (from-foreign sp :knots)
                                 (from-foreign xrp :knots) (from-foreign xsp :knots))))))))

;;; ------------------------------------------------------------- L3-g6

(defun %jiggle-knots (kv eps ends)
  "Knot vector KV with each interior knot moved by up to EPS either way and
left unsorted (NLib's [B168] case of slightly out-of-order knots); with ENDS
(a distance), each end run moved as a whole by up to ENDS too."
  (let* ((kv (copy-seq kv)) (m (1- (length kv))) (lo (aref kv 0)) (hi (aref kv m)))
    (when (plusp eps)
      (dotimes (i (1+ m))
        (let ((x (aref kv i)))
          (unless (or (= x lo) (= x hi))
            (setf (aref kv i) (+ x (* eps (- (rnd 2d0) 1d0)))))))
      (when ends
        (let ((dl (* ends (- (rnd 2d0) 1d0))) (dh (* ends (- (rnd 2d0) 1d0))))
          (dotimes (i (1+ m))
            (cond ((= (aref kv i) lo) (setf (aref kv i) (+ lo dl)))
                  ((= (aref kv i) hi) (setf (aref kv i) (+ hi dh))))))))
    kv))

(define-oracle-test get-compatible-knot-vector-to-tol (:file "basis-basic" :trials 3000)
  (let* ((k (rnd 4))
         (lo (if (zerop (rnd 3)) 0d0 (rnd-range -5d0 5d0)))
         (hi (+ lo 0.5d0 (rnd 5d0)))
         (tol (case (rnd 6) (0 1d-12) (1 1d-9) (2 1d-6) (3 1d-4) (4 1d-2) (t 0.3d0)))
         ;; perturbations below, around or above the tolerance
         ;; (interior knots kept off the ends)
         (eps (min (* 0.01d0 (- hi lo))
                   (case (rnd 5) (0 0d0) (1 (* 0.1d0 tol)) (2 (* 0.45d0 tol)) (3 tol) (t (* 3d0 tol)))))
         ;; end runs a hair apart, mostly within the tolerance
         (ends (case (rnd 6) ((0 1) (* 0.45d0 tol)) (2 (* 0.1d0 tol)) (3 (* 2d0 tol)) (t nil)))
         ;; each vector shares some interior knots with the one before it
         (knt (let (l (prev nil))
                (dotimes (r (1+ k) (coerce (nreverse l) 'simple-vector))
                  (setf prev (%compatible-knots (1+ (rnd 5)) lo hi prev))
                  (push (%jiggle-knots prev eps ends) l)))))
    ;; now and then one vector on another domain: the vectors run out unevenly
    (when (and (plusp k) (zerop (rnd 15)))
      (setf (svref knt (rnd (1+ k))) (%compatible-knots (1+ (rnd 5)) lo (+ hi 1d0))))
    (agree (list :knt knt :k k :tol tol)
           (let* ((kc (map 'simple-vector #'copy-seq knt))
                  (kx (batten::get-compatible-knot-vector-to-tol kc k tol)))
             (list (coerce kc 'list) (coerce kx 'list)))
           (with-stacks (sk)
             (cffi:with-foreign-objects ((kp :pointer (1+ k)) (kxp :pointer))
               ;; TO-FOREIGN inside WITH-FOREIGN: the copies are freed with it
               (with-foreign ((dummy (make-real-vector 1) :real-vector))
                 (declare (ignore dummy))
                 (dotimes (r (1+ k))
                   (setf (cffi:mem-aref kp :pointer r) (to-foreign (svref knt r) :knots)))
                 (setf (cffi:mem-ref kxp :pointer) (cffi:null-pointer))
                 (let ((failed (c-flag-error-p
                                (c-get-compatible-knot-vector-to-tol kp k tol kxp sk))))
                   (values failed
                           (unless failed
                             (let ((kna (cffi:mem-ref kxp :pointer)) (ks nil) (xs nil))
                               (dotimes (r (1+ k))
                                 (push (from-foreign (cffi:mem-aref kp :pointer r) :knots) ks)
                                 (push (from-foreign (cffi:mem-aref kna :pointer r) :knots) xs))
                               (list (nreverse ks) (nreverse xs))))))))))))
