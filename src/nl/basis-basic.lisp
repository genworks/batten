;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_BasisBasic.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp; scratch arrays on the stack where
;;;; the implementation allows, in place of NL_STACKS.  Each routine keeps its
;;;; C source's arithmetic, operation for operation, so that results agree with
;;;; NLib to the bit.

(in-package :batten)

(defun knot-vector-is-param-out-of-bounds (knots u rname)
  "Signal :par-err in routine RNAME when U lies outside [U[0], U[m]].
NLib: N_KnotVectorIsParamOutOfBounds."
  (declare (type real-vector knots) (type double-float u))
  (when (or (< u (aref knots 0)) (> u (aref knots (1- (length knots)))))
    (nurbs-error :par-err rname)))

(declaim (inline knot-vector-is-end-param))
(defun knot-vector-is-end-param (knots u rname)
  "Signal :par-err in routine RNAME when U is on or outside the ends of the
domain: u <= U[0] or u >= U[m].  NLib: N_KnotVectorIsEndParam."
  (declare (type real-vector knots) (type double-float u))
  (when (or (<= u (aref knots 0)) (>= u (aref knots (1- (length knots)))))
    (nurbs-error :par-err rname)))

(declaim (inline basis-find-span))
(defun basis-find-span (knots p u flg)
  "The knot span containing parameter U: the index i with U[i] < u <= U[i+1]
for FLG +left+, U[i] <= u < U[i+1] for +right+.  A parameter on or beyond an
end of the domain lands in the first or last span, stepping past degenerate
(zero-length) end spans.  NLib: N_BasisFindSpan."
  (declare (type real-vector knots) (type fixnum p flg) (type double-float u))
  (fast
    (let* ((m (1- (length knots)))
           (low p)
           (high (- m p))
           (mid (ash (+ low high) -1))
           (lastmid mid))
      (declare (type fixnum m low high mid lastmid))
      (cond
        ((>= u (aref knots (- m p)))
         (let ((i (- m p 1)))
           (declare (type fixnum i))
           (do ()
               ((not (and (< (aref knots (1+ i)) (+ (aref knots i) 1d-12)) (> i p))) i)
             (decf i))))
        ((<= u (aref knots p))
         (let ((i p))
           (declare (type fixnum i))
           (do ()
               ((not (and (< (aref knots (1+ i)) (+ (aref knots i) 1d-12)) (< i (- m p 1)))) i)
             (incf i))))
        ((= flg +left+)
         (do ()
             ((not (or (< u (aref knots mid)) (>= u (aref knots (1+ mid))))) mid)
           (if (< u (aref knots mid))
               (setf high mid)
               (setf low mid))
           (setf mid (ash (+ low high) -1))
           (when (= mid lastmid) (return mid))
           (setf lastmid mid)))
        ((= flg +right+)
         (do ()
             ((not (or (<= u (aref knots mid)) (> u (aref knots (1+ mid))))) mid)
           (if (> u (aref knots mid))
               (setf low mid)
               (setf high mid))
           (setf mid (ash (+ low high) -1))
           (when (= mid lastmid) (return mid))
           (setf lastmid mid)))
        (t (nurbs-error :cal-err 'basis-find-span))))))

(defun basis-get-span-count (knots p)
  "The number of non-zero-length knot spans of a degree-P knot vector.
NLib: N_BasisGetSpanCount."
  (declare (type real-vector knots) (type fixnum p))
  (let ((m (1- (length knots))) (nsp 0))
    (declare (type fixnum m nsp))
    (do ((i p (1+ i))) ((>= i (- m p)) nsp)
      (declare (type fixnum i))
      (when (/= (aref knots i) (aref knots (1+ i)))
        (incf nsp)))))

(declaim (inline basis-eval))
(defun basis-eval (knots p u flg n)
  "Fill N[0..p] with the non-vanishing B-spline basis functions of degree P at
U, and return the span index i (N[k] is basis function i-p+k).
NLib: N_BasisEval."
  (declare (type real-vector knots n) (type fixnum p flg) (type double-float u))
  (fast
    (let ((m (1- (length knots))))
      (declare (type fixnum m))
      (cond
        ((= u (aref knots p))
         (setf (aref n 0) 1d0)
         (do ((k 1 (1+ k))) ((> k p)) (declare (type fixnum k))
           (setf (aref n k) 0d0))
         p)
        ((= u (aref knots (- m p)))
         (setf (aref n p) 1d0)
         (dotimes (k p) (setf (aref n k) 0d0))
         (- m p 1))
        (t
         (let ((i (basis-find-span knots p u flg))
               (left (make-array (+ p 2) :element-type 'double-float))
               (right (make-array (+ p 2) :element-type 'double-float)))
           (declare (type fixnum i) (dynamic-extent left right))
           (setf (aref n 0) 1d0)
           (do ((k 1 (1+ k))) ((> k p))
             (declare (type fixnum k))
             (setf (aref left k) (- u (aref knots (- (+ i 1) k)))
                   (aref right k) (- (aref knots (+ i k)) u))
             (let ((saved 0d0))
               (declare (type double-float saved))
               (dotimes (l k)
                 (let ((delta (+ (aref right (1+ l)) (aref left (- k l))))
                       (temp (aref n l)))
                   (declare (type double-float delta temp))
                   (when (< delta 1d-24) (nurbs-error :num-err 'basis-eval))
                   (setf (aref n l) (+ saved (* temp (/ (aref right (1+ l)) delta)))
                         saved (* temp (/ (aref left (- k l)) delta)))))
               (setf (aref n k) saved)))
           i))))))

(defun basis-derivs (knots p u flg der nd)
  "Fill ND[k][0..p], k = 0..DER, with the non-vanishing basis functions of
degree P at U and their derivatives (ND[k][j] is the kth derivative of basis
function i-p+j), and return the span index i.  ND is at least (1+ DER) by
(1+ P).  NLib: N_BasisDerivs."
  (declare (type real-vector knots) (type real-matrix nd)
           (type fixnum p flg der) (type double-float u))
  (fast
    (let* ((i (basis-find-span knots p u flg))
           (mder (min p der))
           (p+1 (1+ p))
           (left (make-array p+1 :element-type 'double-float))
           (right (make-array p+1 :element-type 'double-float))
           (ndu (make-array (* p+1 p+1) :element-type 'double-float))
           (a (make-array (* 2 p+1) :element-type 'double-float)))
      (declare (type fixnum i mder p+1) (dynamic-extent left right ndu a))
      (macrolet ((ndu (r c) `(aref ndu (+ (the fixnum (* ,r p+1)) ,c)))
                 (a (r c) `(aref a (+ (the fixnum (* ,r p+1)) ,c))))
        (do ((k (1+ p) (1+ k))) ((> k der))
          (declare (type fixnum k))
          (dotimes (j p+1) (setf (aref nd k j) 0d0)))
        (setf (ndu 0 0) 1d0)
        (do ((j 1 (1+ j))) ((> j p))
          (declare (type fixnum j))
          (setf (aref left j) (- u (aref knots (- (+ i 1) j)))
                (aref right j) (- (aref knots (+ i j)) u))
          (let ((saved 0d0))
            (declare (type double-float saved))
            (dotimes (r j)
              (let ((delta (+ (aref right (1+ r)) (aref left (- j r)))))
                (declare (type double-float delta))
                (setf (ndu j r) delta)
                (when (< delta 1d-24) (nurbs-error :num-err 'basis-derivs))
                (setf (ndu r j) (+ saved (* (ndu r (1- j)) (/ (aref right (1+ r)) delta)))
                      saved (* (ndu r (1- j)) (/ (aref left (- j r)) delta)))))
            (setf (ndu j j) saved)))
        (dotimes (j p+1)
          (setf (aref nd 0 j) (ndu j p)))
        (dotimes (r p+1)
          (let ((s1 0) (s2 1))
            (declare (type fixnum s1 s2))
            (setf (a 0 0) 1d0)
            (do ((k 1 (1+ k))) ((> k mder))
              (declare (type fixnum k))
              (let* ((d 0d0)
                     (rk (- r k))
                     (pk (- p k))
                     (j1 (if (>= rk -1) 1 (- rk)))
                     (j2 (if (<= (1- r) pk) (1- k) (- p r))))
                (declare (type double-float d) (type fixnum rk pk j1 j2))
                (when (>= r k)
                  (setf (a s2 0) (/ (a s1 0) (ndu (1+ pk) rk))
                        d (* (a s2 0) (ndu rk pk))))
                (do ((j j1 (1+ j))) ((> j j2))
                  (declare (type fixnum j))
                  (setf (a s2 j) (/ (- (a s1 j) (a s1 (1- j))) (ndu (1+ pk) (+ rk j))))
                  (incf d (* (a s2 j) (ndu (+ rk j) pk))))
                (when (<= r pk)
                  (setf (a s2 k) (/ (- (a s1 (1- k))) (ndu (1+ pk) r)))
                  (incf d (* (a s2 k) (ndu r pk))))
                (setf (aref nd k r) d)
                (rotatef s1 s2)))))
        (let ((r p))
          (declare (type fixnum r))
          (do ((k 1 (1+ k))) ((> k mder))
            (declare (type fixnum k))
            (dotimes (j p+1)
              (setf (aref nd k j) (* (aref nd k j) (float r 1d0))))
            (setf r (* r (- p k)))))
        i))))

;;; ------------------------------------------------- knot vector construction

(defun knot-vector-is-valid (knots p rname)
  "Check that KNOTS is a proper degree-P knot vector: end multiplicities P+1,
knots non-decreasing, no interior knot of multiplicity above P.  Signal
:knt-err in routine RNAME when it is not; return T when it is.
NLib: N_KnotVectorIsValid."
  (declare (type real-vector knots) (type fixnum p))
  (let ((m (1- (length knots)))
        (i 0))
    (declare (type fixnum m i))
    ;; end knots
    (do () ((not (and (< i m) (= (aref knots i) (aref knots (1+ i))))))
      (incf i))
    (unless (= i p)
      (nurbs-error :knt-err rname))
    (setf i m)
    (do () ((not (and (> i 0) (= (aref knots i) (aref knots (1- i))))))
      (decf i))
    (unless (= (- m i) p)
      (nurbs-error :knt-err rname))
    ;; internal knots
    (let ((is p) (ie (1+ p)))
      (declare (type fixnum is ie))
      (do () ((not (<= ie (- m p))))
        (when (> (aref knots is) (aref knots ie))
          (nurbs-error :knt-err rname))
        (setf i ie)
        (do () ((not (and (< ie (- m p)) (= (aref knots ie) (aref knots (1+ ie))))))
          (incf ie))
        (when (> (1+ (- ie i)) p)
          (nurbs-error :knt-err rname))
        (setf is ie)
        (incf ie)))
    t))

(defun knot-vector-from-real-array (u m)
  "The knot vector whose knots are U[0..M].  The C points a knot vector struct
at U; in Batten a knot vector is the real-vector of its knots, so this returns
U itself when it holds exactly M+1 knots (sharing it as the C does) and a copy
of its first M+1 otherwise.  NLib: N_KnotVectorFromRealArray."
  (declare (type real-vector u) (type fixnum m))
  (if (= (length u) (1+ m))
      u
      (subseq u 0 (1+ m))))

(defun knots-check (knots m rname)
  "A knot vector with highest index M: a fresh zeroed one when KNOTS is NIL,
else KNOTS itself cut down to M+1 knots (KNOTS when it has exactly that many,
a copy of its first M+1 otherwise).  Signal :sto-err in routine RNAME when
KNOTS has fewer than M+1 knots.  NLib: N_KnotsCheck."
  (declare (type (or null real-vector) knots) (type fixnum m))
  (cond ((null knots) (make-real-vector (1+ m)))
        ((< (1- (length knots)) m) (nurbs-error :sto-err rname))
        ((= (length knots) (1+ m)) knots)
        (t (subseq knots 0 (1+ m)))))

(defun knots-copy (kna knb)
  "Copy the knots of knot vector KNA into KNB and return KNB: a fresh knot
vector when KNB is NIL, else KNB itself cut down to the length of KNA (as
KNOTS-CHECK).  Signals :sto-err when KNB is shorter than KNA.
NLib: N_KnotsCopy."
  (declare (type real-vector kna) (type (or null real-vector) knb))
  (let* ((ma (1- (length kna)))
         (ub (knots-check knb ma 'knots-copy)))
    (declare (type fixnum ma) (type real-vector ub))
    (fast (dotimes (i (1+ ma)) (setf (aref ub i) (aref kna i))))
    ub))

;;; ---------------------------------------------------------- spans and nodes

(defun basis-find-span-and-mult (knots p u flg)
  "The knot span containing U (as BASIS-FIND-SPAN) and the multiplicity of U
in KNOTS (0 when U is not a knot; P+1 at the end of the domain on the side
FLG faces).  Returns (values spn mlt).  NLib: N_BasisFindSpanAndMult."
  (declare (type real-vector knots) (type fixnum p flg) (type double-float u))
  (let* ((m (1- (length knots)))
         (k (basis-find-span knots p u flg))
         (s 0))
    (declare (type fixnum m k s))
    (cond
      ((= flg +left+)
       (when (= u (aref knots (- m p)))
         (return-from basis-find-span-and-mult (values k (1+ p))))
       (do ((j k (1- j))) ((not (and (>= j 0) (= u (aref knots j)))))
         (declare (type fixnum j))
         (incf s)))
      ((= flg +right+)
       (when (= u (aref knots p))
         (return-from basis-find-span-and-mult (values k (1+ p))))
       (do ((j (1+ k) (1+ j))) ((not (and (<= j m) (= u (aref knots j)))))
         (declare (type fixnum j))
         (incf s)))
      (t (nurbs-error :cal-err 'basis-find-span-and-mult)))
    (values k s)))

(defun basis-find-index-node-array (knots p tt)
  "Write the nodes (Greville abscissae) of the clamped degree-P knot vector
KNOTS into TT, which holds at least m-p values, and return TT.
NLib: N_BasisFindIndexNodeArray."
  (declare (type real-vector knots tt) (type fixnum p))
  (fast
    (let* ((m (1- (length knots)))
           (n (- m p 1))
           (fact (/ 1d0 (float p 1d0))))
      (declare (type fixnum m n) (type double-float fact))
      (setf (aref tt 0) (aref knots 0))
      (setf (aref tt n) (aref knots m))
      (do ((i 1 (1+ i))) ((> i (1- n)))
        (declare (type fixnum i))
        (let ((sum 0d0))
          (declare (type double-float sum))
          (do ((j 1 (1+ j))) ((> j p))
            (declare (type fixnum j))
            (setf sum (+ sum (aref knots (+ i j)))))
          (setf (aref tt i) (* fact sum))))
      tt)))

;;; ------------------------------------------------ one basis function at a time

(defun basis-i-eval (knots i p u flg)
  "The value at U of the Ith B-spline basis function of degree P on the
clamped knot vector KNOTS.  Signals :par-err for U outside the domain,
:ind-err for I outside 0..n and :cal-err for a bad FLG.  NLib: N_BasisIEval."
  (declare (type real-vector knots) (type fixnum i p flg) (type double-float u))
  (let* ((m (1- (length knots)))
         (n (- m p 1)))
    (declare (type fixnum m n))
    (knot-vector-is-param-out-of-bounds knots u 'basis-i-eval)
    (when (or (< i 0) (> i n))
      (nurbs-error :ind-err 'basis-i-eval))
    ;; special cases for end values
    (when (= u (aref knots p))
      (return-from basis-i-eval (if (= i 0) 1d0 0d0)))
    (when (= u (aref knots (- m p)))
      (return-from basis-i-eval (if (= i n) 1d0 0d0)))
    (fast
      (let ((na (make-array (1+ p) :element-type 'double-float)))
        (declare (dynamic-extent na))
        ;; degree zero B-splines
        (cond
          ((= flg +left+)
           (when (or (< u (aref knots i)) (>= u (aref knots (+ i p 1))))
             (return-from basis-i-eval 0d0))
           (dotimes (j (1+ p))
             (setf (aref na j)
                   (if (and (>= u (aref knots (+ i j))) (< u (aref knots (+ i j 1)))) 1d0 0d0))))
          ((= flg +right+)
           (when (or (<= u (aref knots i)) (> u (aref knots (+ i p 1))))
             (return-from basis-i-eval 0d0))
           (dotimes (j (1+ p))
             (setf (aref na j)
                   (if (and (> u (aref knots (+ i j))) (<= u (aref knots (+ i j 1)))) 1d0 0d0))))
          (t (nurbs-error :cal-err 'basis-i-eval)))
        ;; the triangular table
        (do ((k 1 (1+ k))) ((> k p))
          (declare (type fixnum k))
          (let ((saved (if (= (aref na 0) 0d0)
                           0d0
                           (/ (* (- u (aref knots i)) (aref na 0))
                              (- (aref knots (+ i k)) (aref knots i))))))
            (declare (type double-float saved))
            (dotimes (j (1+ (- p k)))
              (let ((ur (aref knots (+ i j k 1)))
                    (ul (aref knots (+ i j 1))))
                (declare (type double-float ur ul))
                (if (= (aref na (1+ j)) 0d0)
                    (setf (aref na j) saved
                          saved 0d0)
                    (let ((temp (/ (aref na (1+ j)) (- ur ul))))
                      (declare (type double-float temp))
                      (setf (aref na j) (+ saved (* (- ur u) temp))
                            saved (* (- u ul) temp))))))))
        (aref na 0)))))

(defun basis-i-derivs (knots i p u flg der nd)
  "Write the value at U of the Ith B-spline basis function of degree P on the
clamped knot vector KNOTS, and its derivatives through DER, into ND[0..DER];
return ND.  FLG +left+ takes the derivative from the right except at the end
of the domain, +right+ from the left except at its start.  Signals :par-err,
:ind-err or :cal-err as BASIS-I-EVAL.  NLib: N_BasisIDerivs."
  (declare (type real-vector knots nd) (type fixnum i p flg der) (type double-float u))
  (let* ((m (1- (length knots)))
         (n (- m p 1)))
    (declare (type fixnum m n))
    (knot-vector-is-param-out-of-bounds knots u 'basis-i-derivs)
    (when (or (< i 0) (> i n))
      (nurbs-error :ind-err 'basis-i-derivs))
    (fast
      (let* ((p+1 (1+ p))
             (nt (make-array (* p+1 p+1) :element-type 'double-float))
             (ndl (make-array p+1 :element-type 'double-float)))
        (declare (type fixnum p+1) (dynamic-extent nt ndl))
        (macrolet ((nt (r c) `(aref nt (+ (the fixnum (* ,r p+1)) ,c))))
          (flet ((closed-right ()
                   ;; degree zero on (U[j], U[j+1]]
                   (dotimes (j p+1)
                     (setf (nt j 0)
                           (if (and (> u (aref knots (+ i j))) (<= u (aref knots (+ i j 1)))) 1d0 0d0))))
                 (closed-left ()
                   ;; degree zero on [U[j], U[j+1])
                   (dotimes (j p+1)
                     (setf (nt j 0)
                           (if (and (>= u (aref knots (+ i j))) (< u (aref knots (+ i j 1)))) 1d0 0d0))))
                 (zero-out ()
                   (dotimes (j (1+ der)) (setf (aref nd j) 0d0))
                   (return-from basis-i-derivs nd)))
            (declare (inline closed-right closed-left))
            (cond
              ((= flg +left+)
               (if (and (= u (aref knots (- m p))) (>= i (- n p)))
                   (closed-right)
                   (progn
                     (when (or (< u (aref knots i)) (>= u (aref knots (+ i p 1))))
                       (zero-out))
                     (closed-left))))
              ((= flg +right+)
               (if (and (= u (aref knots p)) (<= i p))
                   (closed-left)
                   (progn
                     (when (or (<= u (aref knots i)) (> u (aref knots (+ i p 1))))
                       (zero-out))
                     (closed-right))))
              (t (nurbs-error :cal-err 'basis-i-derivs))))
          ;; the full triangular array
          (do ((k 1 (1+ k))) ((> k p))
            (declare (type fixnum k))
            (let ((saved (if (= (nt 0 (1- k)) 0d0)
                             0d0
                             (/ (* (- u (aref knots i)) (nt 0 (1- k)))
                                (- (aref knots (+ i k)) (aref knots i))))))
              (declare (type double-float saved))
              (dotimes (j (1+ (- p k)))
                (let ((ur (aref knots (+ i j k 1)))
                      (ul (aref knots (+ i j 1))))
                  (declare (type double-float ur ul))
                  (if (= (nt (1+ j) (1- k)) 0d0)
                      (setf (nt j k) saved
                            saved 0d0)
                      (let ((temp (/ (nt (1+ j) (1- k)) (- ur ul))))
                        (declare (type double-float temp))
                        (setf (nt j k) (+ saved (* (- ur u) temp))
                              saved (* (- u ul) temp))))))))
          ;; derivatives
          (setf (aref nd 0) (nt 0 p))
          (let ((mder (nl-min p der)))
            (declare (type fixnum mder))
            (do ((k (1+ p) (1+ k))) ((> k der))
              (declare (type fixnum k))
              (setf (aref nd k) 0d0))
            (do ((k 1 (1+ k))) ((> k mder))
              (declare (type fixnum k))
              (dotimes (j (1+ k))
                (setf (aref ndl j) (nt j (- p k))))
              (do ((l 1 (1+ l))) ((> l k))
                (declare (type fixnum l))
                (let ((saved (if (= (aref ndl 0) 0d0)
                                 0d0
                                 (/ (aref ndl 0) (- (aref knots (+ i (- p k) l)) (aref knots i)))))
                      (fac (+ (- (float p 1d0) (float k 1d0)) (float l 1d0))))
                  (declare (type double-float saved fac))
                  (dotimes (j (1+ (- k l)))
                    (let ((ur (aref knots (+ (- p k) l i j 1)))
                          (ul (aref knots (+ i j 1))))
                      (declare (type double-float ur ul))
                      (if (= (aref ndl (1+ j)) 0d0)
                          (setf (aref ndl j) (* fac saved)
                                saved 0d0)
                          (let ((temp (/ (aref ndl (1+ j)) (- ur ul))))
                            (declare (type double-float temp))
                            (setf (aref ndl j) (* fac (- saved temp))
                                  saved temp)))))))
              (setf (aref nd k) (aref ndl 0))))
          nd)))))

(defun basis-find-global-max (knots i p tol)
  "The global maximum of the Ith B-spline basis function of degree P on
KNOTS, found by bracketing and Newton iteration to tolerance TOL.  Returns
(values max u), U the parameter of the maximum.  Signals :ind-err for I
outside 0..n, :num-err for a bad Newton division and :con-err when the
iteration does not converge.  NLib: N_BasisFindGlobalMax."
  (declare (type real-vector knots) (type fixnum i p) (type double-float tol))
  (let* ((itl 20)
         (nok 10)
         (m (1- (length knots)))
         (conv nil)
         (nd (make-array 3 :element-type 'double-float :initial-element 0d0))
         (u0 0d0))
    (declare (type fixnum itl nok m) (type double-float u0) (dynamic-extent nd))
    (when (or (< i 0) (> i (- m p 1)))
      (nurbs-error :ind-err 'basis-find-global-max))
    ;; linear as special case
    (when (= p 1)
      (return-from basis-find-global-max (values 1d0 (aref knots (1+ i)))))
    ;; p-fold inner multiple knot
    (let ((s 1))
      (declare (type fixnum s))
      (do ((k (1+ i) (1+ k))) ((>= k (+ i p)))
        (declare (type fixnum k))
        (when (= (aref knots k) (aref knots (1+ k)))
          (incf s)))
      (when (= s p)
        (return-from basis-find-global-max (values 1d0 (aref knots (1+ i))))))
    ;; guess parameter
    (let* ((nos (* p nok))
           (uend (aref knots (+ i p 1)))
           (du (/ (- uend (aref knots i)) (float nos 1d0)))
           (ur (aref knots i))
           (ul (aref knots i))
           (it 0))
      (declare (type fixnum nos it) (type double-float uend du ur ul))
      (do () ((not (< ur uend)))
        (setf ul ur
              ur (+ ur du))
        (when (> ur uend) (setf ur uend))
        (basis-i-derivs knots i p ur +left+ 1 nd)
        (when (< (aref nd 1) 0d0) (return)))
      ;; Newton iterations till convergence
      (do () ((not (< it itl)))
        (setf u0 (* 0.5d0 (+ ul ur)))
        (let ((k 0))
          (declare (type fixnum k))
          (do () ((not (< k itl)))
            (basis-i-derivs knots i p u0 +left+ 2 nd)
            (when (and (< (abs (aref nd 1)) tol) (> (aref nd 0) tol))
              (setf conv t)
              (return))
            (when (float-op-is-bad (aref nd 1) (aref nd 2) +division+)
              (nurbs-error :num-err 'basis-find-global-max))
            (setf u0 (- u0 (/ (aref nd 1) (aref nd 2))))
            (when (or (<= u0 (aref knots i)) (>= u0 uend))
              (return))
            (incf k)))
        (when conv (return))
        ;; no convergence: refine [ul,ur] for a better guess
        (setf du (/ (- ur ul) (float nok 1d0))
              ur ul)
        (do ((k 1 (1+ k))) ((> k nok))
          (declare (type fixnum k))
          (setf ul ur
                ur (+ ur du))
          (basis-i-derivs knots i p ur +left+ 1 nd)
          (when (< (aref nd 1) 0d0) (return)))
        (incf it))
      (when (>= it itl)
        (nurbs-error :con-err 'basis-find-global-max))
      (values (aref nd 0) u0))))

(defun basis-reparam (knots p interval)
  "Map the degree-P knot vector KNOTS in place onto INTERVAL [a, b] (a
real-vector of 2): end knots set to a and b, interior knots scaled linearly.
Left alone when its ends already are a and b.  Returns KNOTS.
NLib: N_BasisReparam."
  (declare (type real-vector knots interval) (type fixnum p))
  (let ((m (1- (length knots))))
    (declare (type fixnum m))
    (multiple-value-bind (ul ur) (interval-get-data interval)
      (declare (type double-float ul ur))
      (when (or (/= (aref knots 0) ul) (/= (aref knots m) ur))
        (let* ((u0 (aref knots 0))
               (fac (/ (- ur ul) (- (aref knots m) (aref knots 0)))))
          (declare (type double-float u0 fac))
          (fast
            (do ((i 0 (1+ i))) ((> i p))
              (declare (type fixnum i))
              (setf (aref knots i) ul))
            (do ((i (1+ p) (1+ i))) ((> i (- m p 1)))
              (declare (type fixnum i))
              (setf (aref knots i) (+ ul (* fac (- (aref knots i) u0)))))
            (do ((i (- m p) (1+ i))) ((> i m))
              (declare (type fixnum i))
              (setf (aref knots i) ur))))))
    knots))
