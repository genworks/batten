;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/funcs-adv.lisp against NLib.

(in-package :batten-tests)

(defun %fa-copy-cfun (f)
  "A copy of curve function F with its own arrays."
  (%make-cfun (cfun-degree f) (copy-seq (cfun-knots f)) (copy-seq (cfun-vals f))))

(defun %fa-insert-knots (knots p)
  "Knots to insert into KNOTS (degree P): sorted, strictly inside the domain,
some existing interior knots and some repeated, every resulting multiplicity
at most P.  Now and then an end knot or one outside (:par-err), or none at
all (:inp-err)."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (interior (remove-if (lambda (x) (or (= x lo) (= x hi))) (coerce knots 'list)))
         (cands nil))
    (dotimes (i (1+ (rnd 7)))
      (push (case (rnd 6)
              (0 (if interior (nth (rnd (length interior)) interior) (rnd-range lo hi)))
              (1 (if cands (first cands) (rnd-range lo hi)))
              (t (rnd-range lo hi)))
            cands))
    (setf cands (sort (remove-if (lambda (x) (or (<= x lo) (>= x hi))) cands) #'<))
    (let ((accepted nil))
      (dolist (x cands)
        (when (<= (+ (count x knots :test #'=) (count x accepted :test #'=) 1) p)
          (push x accepted)))
      (setf accepted (nreverse accepted))
      (case (rnd 25)
        (0 (setf accepted nil))
        (1 (push lo accepted))
        (2 (setf accepted (append accepted (list hi))))
        (3 (setf accepted (append accepted (list (+ hi 1d0))))))
      (coerce accepted 'real-vector))))

(define-oracle-test crv-func-refine (:file "funcs-adv" :trials 4000)
  (let* ((f (random-cfun))
         (knx (%fa-insert-knots (cfun-knots f) (cfun-degree f)))
         (in-place (zerop (rnd 2))))
    (agree (list :cfun f :knx knx :in-place in-place)
           (list (batten::crv-func-refine f knx))
           (with-foreign ((fp f :cfun) (kx knx :knots))
             (with-stacks (sp)
               (if in-place
                   (if (c-flag-error-p (c-crv-func-refine fp kx fp sp sp))
                       (values t nil)
                       (values nil (list (from-foreign fp :cfun))))
                   (with-stacks (sq)
                     (with-foreign-out ((fq :cfun))
                       (if (c-flag-error-p (c-crv-func-refine fp kx fq sp sq))
                           (values t nil)
                           (values nil (list (from-foreign fq :cfun))))))))))))

(define-oracle-test crv-func-deriv (:file "funcs-adv" :trials 4000)
  (let* ((f (random-cfun))
         (p (cfun-degree f))
         (d (case (rnd 20) (0 -1) (1 (1+ p)) (t (rnd (1+ p))))))
    (agree (list :cfun f :d d)
           (list (batten::crv-func-deriv f d))
           (with-foreign ((fp f :cfun))
             (with-stacks (sg)
               (with-foreign-out ((der :cfun))
                 (if (c-flag-error-p (c-crv-func-deriv fp d der sg))
                     (values t nil)
                     (values nil (list (from-foreign der :cfun))))))))))

(define-oracle-test crv-func-multiply-constant (:file "funcs-adv" :trials 2000)
  (let ((f (random-cfun))
        (alpha (case (rnd 8) (0 0d0) (1 -1d0) (t (rnd-range -10d0 10d0)))))
    (agree (list :cfun f :alpha alpha)
           (list (batten::crv-func-multiply-constant alpha (%fa-copy-cfun f)))
           (with-foreign ((fp f :cfun))
             (c-crv-func-multiply-constant alpha fp)
             (values nil (list (from-foreign fp :cfun)))))))

;;; ------------------------------------------------------------- L2-g3

(defun %fa-monotone-cfun ()
  "A random curve function, two times in three with sorted control values (so
increasing, and invertible)."
  (let ((f (random-cfun)))
    (if (plusp (rnd 3))
        (%make-cfun (cfun-degree f) (cfun-knots f) (sort (copy-seq (cfun-vals f)) #'<))
        f)))

(define-oracle-test crv-func-eval-invert-pt (:file "funcs-adv" :trials 4000)
  (let* ((f (%fa-monotone-cfun))
         (k (cfun-knots f))
         (a (batten::c-func-eval f (aref k 0) +left+))
         (b (batten::c-func-eval f (aref k (1- (length k))) +left+))
         (val (case (rnd 12)
                (0 (+ a (* 1d-13 (rnd 1d0))))
                (1 (- b (* 1d-13 (rnd 1d0))))
                (2 (- (min a b) 1d0))
                (3 (+ (max a b) 1d0))
                (4 (rnd-range (min a b) (max a b)))
                (t (batten::c-func-eval f (random-parameter k) +left+))))
         (u0 (random-parameter k))
         (tol (case (rnd 4) (0 1d-6) (1 1d-10) (2 1d-12) (t 1d-8)))
         (itl (case (rnd 6) (0 0) (1 2) (t (+ 5 (rnd 40))))))
    (agree (list :cfun f :f val :u0 u0 :tol tol :itl itl)
           (list (batten::crv-func-eval-invert-pt f val u0 tol itl))
           (with-foreign ((fp f :cfun))
             (with-scalars ((u :double))
               (values (c-flag-error-p (c-crv-func-eval-invert-pt fp val u0 tol itl u))
                       (list (cffi:mem-ref u :double))))))))

(defun %fa-spoil-knots (k p)
  "A copy of knot vector K, now and then invalid for degree P: an end
multiplicity off or an interior knot out of order."
  (let* ((k (copy-seq k)) (m (1- (length k))))
    (case (rnd 16)
      (0 (setf (aref k p) (aref k (1+ p))))
      (1 (setf (aref k (- m p)) (aref k (- m p 1))))
      (2 (when (> (- m p) (1+ p))
           (let ((i (+ p 1 (rnd (- m p p 1)))))
             (setf (aref k i) (+ (aref k i) (* 0.5d0 (- (aref k m) (aref k 0)))))))))
    k))

(define-oracle-test crv-func-degree-elevate (:file "funcs-adv" :trials 3000)
  (let* ((f0 (random-cfun))
         (p (cfun-degree f0))
         (f (%make-cfun p (%fa-spoil-knots (cfun-knots f0) p) (cfun-vals f0)))
         (tt (case (rnd 25) (0 -1) (1 (- 18 p)) (t (rnd 4))))
         (in-place (zerop (rnd 2))))
    (agree (list :cfun f :t tt :in-place in-place)
           (list (batten::crv-func-degree-elevate f tt))
           (with-foreign ((fp f :cfun))
             (with-stacks (sp)
               (if in-place
                   (if (c-flag-error-p (c-crv-func-degree-elevate fp tt fp sp sp))
                       (values t nil)
                       (values nil (list (from-foreign fp :cfun))))
                   (with-stacks (sq)
                     (with-foreign-out ((fq :cfun))
                       (if (c-flag-error-p (c-crv-func-degree-elevate fp tt fq sp sq))
                           (values t nil)
                           (values nil (list (from-foreign fq :cfun))))))))))))

;;; ------------------------------------------------------------- L3-g1

(defun %fa-domain ()
  "A random parameter domain (values lo hi)."
  (let* ((lo (if (zerop (rnd 3)) 0d0 (- (rnd 20d0) 10d0)))
         (hi (+ lo 0.5d0 (rnd 10d0))))
    (values lo hi)))

(defun %fa-kink-cfun ()
  "A curve function of degree 2..4 with a knot of multiplicity p (a kink)
followed closely by a simple knot; its values are set by the caller."
  (multiple-value-bind (s0 s1) (%fa-domain)
    (let* ((p (+ 2 (rnd 3)))
           (a (+ s0 (* (- s1 s0) (+ 0.2d0 (rnd 0.5d0)))))
           (b (+ a (* (- s1 s0) (+ 0.002d0 (rnd 0.02d0))))))
      (%make-cfun p
                  (coerce (append (make-list (1+ p) :initial-element s0)
                                  (make-list p :initial-element a)
                                  (list b)
                                  (make-list (1+ p) :initial-element s1))
                          'real-vector)
                  (make-real-vector (+ p p 2))))))

(defun %fa-map-knots-case ()
  "A curve function F mapping its domain onto [lo,hi] (its end values at the
ends of [lo,hi], now and then off by less than the routine's tolerance or,
rarely, by more), increasing three times in four, and a knot vector KNU of
degree PU on [lo,hi] whose interior knots are random or F's values at F's
interior knots (now and then exactly, or very nearly).  Returns
(values f knu pu)."
  (multiple-value-bind (lo hi) (%fa-domain)
    (let* ((kink (zerop (rnd 5)))
           (lo (if (and kink (zerop (rnd 2))) 0d0 lo))
           (hi (if (> hi lo) hi (+ lo 0.5d0 (rnd 10d0))))
           (f0 (if kink (%fa-kink-cfun) (random-cfun :max-degree 5 :max-extra 8)))
           (n+1 (length (cfun-vals f0)))
           (vals (make-real-vector n+1))
           (monotone (plusp (rnd 4)))
           (pu (1+ (rnd 4)))
           (base (random-knots pu (+ pu 1 (rnd 8)) :lo lo :hi hi))
           (spans (max 1 (batten::basis-get-span-count base pu)))
           (ptol (/ (- hi lo) (* spans 1000d0))))
      ;; the function's values
      (let ((inner (let (l)
                     (dotimes (i (- n+1 2) l)
                       (push (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0)))) l)))))
        (when (or monotone kink) (setf inner (sort inner #'<)))
        (setf (aref vals 0) lo
              (aref vals (1- n+1)) hi)
        (let ((i 1))
          (dolist (v inner) (setf (aref vals i) v) (incf i)))
        ;; a kink function is flat (at lo) up to its kink
        (when kink
          (dotimes (i (1+ (cfun-degree f0)))
            (setf (aref vals i) lo))))
      (case (rnd 10)
        (0 (incf (aref vals 0) (* (- (rnd 0.8d0) 0.4d0) ptol)))
        (1 (incf (aref vals (1- n+1)) (* (- (rnd 0.8d0) 0.4d0) ptol)))
        (2 (when (zerop (rnd 4))
             (if (zerop (rnd 2))
                 (incf (aref vals 0) (* 3d0 ptol))
                 (decf (aref vals (1- n+1)) (* 3d0 ptol))))))
      (let* ((f (%make-cfun (cfun-degree f0) (cfun-knots f0) vals))
             (fk (cfun-knots f))
             (p (cfun-degree f))
             (images (let (l)
                       (do ((i (1+ p) (1+ i))) ((> i (1- n+1)) l)
                         (push (batten::c-func-eval f (aref fk i) +left+) l))))
             (interior (coerce (subseq base (1+ pu) (- (length base) pu 1)) 'list)))
        ;; replace some interior knots by images of F's knots
        (when (and images (plusp (rnd 2)))
          (dotimes (j (1+ (rnd 3)))
            (let ((x (nth (rnd (length images)) images)))
              (case (rnd 3)
                (0 nil)
                (1 (setf x (+ x (* (- (rnd 0.6d0) 0.3d0) ptol))))
                (t (setf x (+ x (* (- (rnd 0.04d0) 0.02d0) (- hi lo))))))
              (when (and (> x (+ lo (* 0.01d0 (- hi lo)))) (< x (- hi (* 0.01d0 (- hi lo)))))
                (if (and interior (zerop (rnd 2)))
                    (setf (nth (rnd (length interior)) interior) x)
                    (push x interior))))))
        ;; a non-monotone F: now and then a knot just below a local maximum,
        ;; where the inversion's Newton iteration may fail
        (when (and (not monotone) (plusp (rnd 3)))
          (let ((maxima nil) (fa (cfun-knots f)) (prev nil) (cur nil))
            (dotimes (i 400)
              (let ((v (batten::c-func-eval
                        f (+ (aref fa 0) (* (/ i 400d0) (- (aref fa (1- (length fa))) (aref fa 0))))
                        +left+)))
                (when (and prev cur (> cur prev) (>= cur v)) (push cur maxima))
                (setf prev cur cur v)))
            (when maxima
              (let ((x (- (nth (rnd (length maxima)) maxima)
                          (* (case (rnd 3) (0 1d-3) (1 1d-5) (t 1d-7)) (- hi lo)))))
                (when (and (> x (+ lo (* 0.01d0 (- hi lo)))) (< x (- hi (* 0.01d0 (- hi lo)))))
                  (push x interior))))))
        ;; a kink function: a knot just above its flat value
        (when kink
          (push (+ lo (* (case (rnd 3) (0 2d0) (1 5d0) (t 20d0)) ptol)) interior))
        (setf interior (sort interior #'<))
        ;; multiplicities at most pu
        (let ((kept nil))
          (dolist (x interior)
            (when (< (count x kept :test #'=) pu) (push x kept)))
          (setf interior (nreverse kept)))
        (values f
                (coerce (append (make-list (1+ pu) :initial-element lo) interior
                                (make-list (1+ pu) :initial-element hi))
                        'real-vector)
                pu)))))

(define-oracle-test map-knots-between-crv-func-and-knot-vector (:file "funcs-adv" :trials 2000)
  (multiple-value-bind (f knu pu) (%fa-map-knots-case)
    (let ((cap (+ (length (cfun-knots f)) (length knu) 2)))
      (agree (list :cfun f :knu knu :pu pu)
             (multiple-value-list (batten::map-knots-between-crv-func-and-knot-vector f knu pu))
             (with-foreign ((fp f :cfun) (ku knu :knots)
                            (xu (make-real-vector cap) :knots) (xs (make-real-vector cap) :knots))
               (if (c-flag-error-p (c-map-knots-between-crv-func-and-knot-vector fp ku pu xu xs))
                   (values t nil)
                   (values nil (list (from-foreign xu :knots) (from-foreign xs :knots)))))))))

(defun %fa-copy-curve (c)
  "A copy of curve C with its own arrays."
  (let* ((pw (curve-cpoints c))
         (qw (make-cpoint-array (array-dimension pw 0))))
    (dotimes (i (array-dimension pw 0))
      (dotimes (k 4) (setf (aref qw i k) (aref pw i k))))
    (%make-curve (curve-degree c) (copy-seq (curve-knots c)) qw)))

(defun %fa-partner-knots (knots p q)
  "A knot vector of degree Q for a partner of a function with knot vector
KNOTS of degree P: on another random domain, on the same domain, or sharing
some or all of KNOTS's distinct interior knots.  Returns it and the count of
control values it takes."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (distinct (remove-duplicates
                    (coerce (subseq knots (1+ p) (- (length knots) p 1)) 'list) :test #'=))
         (k (case (rnd 4)
              (0 (multiple-value-bind (l h) (%fa-domain)
                   (random-knots q (+ q 1 (rnd 8)) :lo l :hi h)))
              (1 (random-knots q (+ q 1 (rnd 8)) :lo lo :hi hi))
              (t (let ((in (remove-if (lambda (x) (declare (ignore x)) (zerop (rnd 3))) distinct)))
                   (dotimes (i (rnd 3))
                     (push (+ lo (* (- hi lo) (+ 0.02d0 (rnd 0.96d0)))) in))
                   (setf in (sort in #'<))
                   ;; now and then a repeated knot
                   (when (and in (zerop (rnd 3)) (> q 1))
                     (push (first in) in))
                   (coerce (append (make-list (1+ q) :initial-element lo) in
                                   (make-list (1+ q) :initial-element hi))
                           'real-vector))))))
    (values k (- (length k) q 1))))

(defun %fa-degree-0-cfun (segments)
  "A degree-0 curve function with SEGMENTS spans."
  (multiple-value-bind (lo hi) (%fa-domain)
    (let ((knots (make-real-vector (1+ segments)))
          (vals (random-real-vector segments)))
      (dotimes (i (1+ segments))
        (setf (aref knots i) (+ lo (* (- hi lo) (/ (float i 1d0) segments)))))
      (%make-cfun 0 knots vals))))

(defun %fa-bezier-cfun (q)
  "A degree-Q curve function of one segment."
  (multiple-value-bind (lo hi) (%fa-domain)
    (%make-cfun q (random-knots q (1+ q) :lo lo :hi hi) (random-real-vector (1+ q)))))

(defun %fa-product-pair ()
  "Two curve functions to multiply: usually random of degree 1..7 with
partner knots of several kinds, now and then degrees summing past *dmax*, a
degree-0 factor with a single-segment partner, or a degree-0 factor or partner
with more than one segment."
  (case (rnd 16)
    (0 (let ((f (%fa-degree-0-cfun (if (zerop (rnd 4)) 2 1)))
             (g (if (zerop (rnd 4)) (random-cfun) (%fa-bezier-cfun (rnd 6)))))
         (if (zerop (rnd 2)) (values f g) (values g f))))
    (1 (values (random-cfun :max-degree 10) (random-cfun :max-degree 10)))
    (t (let* ((f (random-cfun))
              (q (1+ (rnd 7))))
         (multiple-value-bind (k n+1) (%fa-partner-knots (cfun-knots f) (cfun-degree f) q)
           (let ((g (%make-cfun q k (random-real-vector n+1))))
             (if (zerop (rnd 2)) (values f g) (values g f))))))))

(define-oracle-test crv-func-multiply-crv-func (:file "funcs-adv" :trials 3000)
  (multiple-value-bind (f g) (%fa-product-pair)
    (agree (list :f f :g g)
           (let ((lf (%fa-copy-cfun f)) (lg (%fa-copy-cfun g)))
             (list (batten::crv-func-multiply-crv-func lf lg) lf lg))
           (with-foreign ((fp f :cfun) (gp g :cfun))
             (with-stacks (sh)
               (with-foreign-out ((hp :cfun))
                 (if (c-flag-error-p (c-crv-func-multiply-crv-func fp gp hp sh))
                     (values t nil)
                     (values nil (list (from-foreign hp :cfun) (from-foreign fp :cfun)
                                       (from-foreign gp :cfun))))))))))

(defun %fa-curve-on (q knots n+1)
  "A random curve of degree Q on KNOTS with N+1 control points: rational half
the time, 2D a fifth of the time."
  (let* ((rational (zerop (rnd 2)))
         (two-d (zerop (rnd 5)))
         (pw (make-cpoint-array n+1)))
    (dotimes (i n+1)
      (let ((w (if rational (+ 0.2d0 (rnd 4.8d0)) 1d0)))
        (setf (aref pw i 0) (* w (- (rnd 200d0) 100d0))
              (aref pw i 1) (* w (- (rnd 200d0) 100d0))
              (aref pw i 2) (if two-d +no-z+ (* w (- (rnd 200d0) 100d0)))
              (aref pw i 3) (if rational w +no-w+))))
    (%make-curve q knots pw)))

(defun %fa-function-curve-pair ()
  "A curve function and a curve to multiply, as %FA-PRODUCT-PAIR makes pairs
of functions."
  (case (rnd 16)
    (0 (if (zerop (rnd 2))
           (let ((f (%fa-degree-0-cfun (if (zerop (rnd 4)) 2 1)))
                 (q (rnd 6)))
             (values f (if (zerop (rnd 4))
                           (random-curve)
                           (multiple-value-bind (lo hi) (%fa-domain)
                             (%fa-curve-on q (random-knots q (1+ q) :lo lo :hi hi) (1+ q))))))
           (let* ((f (if (zerop (rnd 4)) (random-cfun) (%fa-bezier-cfun (1+ (rnd 6)))))
                  (segs (if (zerop (rnd 4)) 2 1))
                  (g0 (%fa-degree-0-cfun segs)))
             (values f (%fa-curve-on 0 (cfun-knots g0) segs)))))
    (1 (values (random-cfun :max-degree 10) (random-curve :max-degree 10)))
    (t (let* ((f (random-cfun))
              (q (1+ (rnd 7))))
         (multiple-value-bind (k n+1) (%fa-partner-knots (cfun-knots f) (cfun-degree f) q)
           (values f (%fa-curve-on q k n+1)))))))

(define-oracle-test crv-func-multiply-crv (:file "funcs-adv" :trials 3000)
  (multiple-value-bind (f g) (%fa-function-curve-pair)
    (agree (list :f f :g g)
           (let ((lf (%fa-copy-cfun f)) (lg (%fa-copy-curve g)))
             (list (batten::crv-func-multiply-crv lf lg) lf lg))
           (with-foreign ((fp f :cfun) (gp g :curve))
             (with-stacks (sh)
               (with-foreign-out ((hp :curve))
                 (if (c-flag-error-p (c-crv-func-multiply-crv fp gp hp sh))
                     (values t nil)
                     (values nil (list (from-foreign hp :curve) (from-foreign fp :cfun)
                                       (from-foreign gp :curve))))))))))

(define-oracle-test crv-func-sum-diff-crv-func (:file "funcs-adv" :trials 3000)
  (let* ((p (random-cfun))
         (cmp (if (zerop (rnd 3)) +yes+ +no+))
         (q (if (= cmp +yes+)
                (%make-cfun (cfun-degree p) (copy-seq (cfun-knots p))
                            (random-real-vector (length (cfun-vals p))))
                (let ((qd (1+ (rnd 7))))
                  (multiple-value-bind (k n+1) (%fa-partner-knots (cfun-knots p) (cfun-degree p) qd)
                    (%make-cfun qd k (random-real-vector n+1))))))
         (opr (case (rnd 20) (0 7) (t (if (zerop (rnd 2)) +plus+ +minus+)))))
    (agree (list :p p :q q :opr opr :cmp cmp)
           (let ((lp (%fa-copy-cfun p)) (lq (%fa-copy-cfun q)))
             (list (batten::crv-func-sum-diff-crv-func lp lq opr cmp) lp lq))
           (with-foreign ((pp p :cfun) (qp q :cfun))
             (with-stacks (sg)
               (with-foreign-out ((rp :cfun))
                 (if (c-flag-error-p (c-crv-func-sum-diff-crv-func pp qp opr cmp rp sg))
                     (values t nil)
                     (values nil (list (from-foreign rp :cfun) (from-foreign pp :cfun)
                                       (from-foreign qp :cfun))))))))))
