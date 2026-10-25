;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvFit.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.
;;;;
;;;; Knot vector outputs (NL_KNOTVECTOR *knt, "MEMORY FOR knt IS ALLOCATED IN
;;;; THE CALLING ROUTINE") are real-vectors the caller allocates with exactly
;;;; the number of knots the C routine records in knt->m (its N_SetKnotIndex),
;;;; plus one; the routine writes the knots the C writes and leaves the others
;;;; as they were.

(in-package :batten)

;;; ------------------------------------------------- parameters for interpolation

(defun fit-calc-crv-param-values (a n ptp par u)
  "Write into U[0..N] parameter values for interpolating the points of A
(rows 0..N): a point-array when PTP is +epoint+, a cpoint-array when it is
+hpoint+.  PAR is +uniform+, +chordlength+ or +centripetal+.  U is allocated
by the caller; when the total chord length is not above *ptol*, U[1..N-1] are
left as they were.  Signals :cal-err for an unknown PTP or PAR (U[0] and U[N]
already written).  Returns U.
NLib: N_FitCalcCrvParamValues."
  (declare (type (simple-array double-float (* *)) a) (type fixnum n ptp par) (type real-vector u))
  (setf (aref u 0) 0d0
        (aref u n) 1d0)
  (when (= par +uniform+)
    ;; fact = 1.0/n; with n = 0 the C's infinite fact is never used
    (when (> n 0)
      (let ((fact (/ 1d0 (float n 1d0))))
        (declare (type double-float fact))
        (fast (do ((i 1 (1+ i))) ((>= i n))
                (declare (type fixnum i))
                (setf (aref u i) (* (float i 1d0) fact))))))
    (return-from fit-calc-crv-param-values u))
  (let ((dist (make-real-vector (1+ n)))
        (sum 0d0))
    (declare (type real-vector dist) (type double-float sum))
    (cond
      ((= ptp +epoint+)
       (let ((pts a) (p0 (make-point)) (p1 (make-point)))
         (declare (type point-array pts) (type point p0 p1) (dynamic-extent p0 p1))
         (flet ((dist-at (i)
                  (declare (type fixnum i))
                  (set-point p0 (aref pts (1- i) 0) (aref pts (1- i) 1) (aref pts (1- i) 2))
                  (set-point p1 (aref pts i 0) (aref pts i 1) (aref pts i 2))
                  (dist-pt-pt p0 p1)))
           (declare (inline dist-at))
           (cond ((= par +chordlength+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf sum (+ sum (aref dist i)))))
                 ((= par +centripetal+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf (aref dist i) (sqrt (the (double-float 0d0) (aref dist i))))
                    (setf sum (+ sum (aref dist i)))))
                 (t (nurbs-error :cal-err 'fit-calc-crv-param-values))))))
      ((= ptp +hpoint+)
       (let ((pw a) (q0 (make-cpoint)) (q1 (make-cpoint)))
         (declare (type cpoint-array pw) (type cpoint q0 q1) (dynamic-extent q0 q1))
         (flet ((dist-at (i)
                  (declare (type fixnum i))
                  (set-cpoint q0 (aref pw (1- i) 0) (aref pw (1- i) 1) (aref pw (1- i) 2) (aref pw (1- i) 3))
                  (set-cpoint q1 (aref pw i 0) (aref pw i 1) (aref pw i 2) (aref pw i 3))
                  (dist-cpt-cpt-homo q0 q1)))
           (declare (inline dist-at))
           (cond ((= par +chordlength+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf sum (+ sum (aref dist i)))))
                 ((= par +centripetal+)
                  (do ((i 1 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref dist i) (dist-at i))
                    (setf (aref dist i) (sqrt (the (double-float 0d0) (aref dist i))))
                    (setf sum (+ sum (aref dist i)))))
                 (t (nurbs-error :cal-err 'fit-calc-crv-param-values))))))
      (t (nurbs-error :cal-err 'fit-calc-crv-param-values)))
    (when (> sum *ptol*)
      (fast (do ((i 1 (1+ i))) ((>= i n))
              (declare (type fixnum i))
              (setf (aref u i) (+ (aref u (1- i)) (/ (aref dist i) sum))))))
    u))

;;; ------------------------------------------------------------- knot vectors

(defun fit-crv-calc-knot-vector (u n p knt)
  "Write into KNT (N+P+2 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..N]: P+1 end knots
equal to U[0] and U[N], interior knots averages of P consecutive parameters.
Returns KNT.
NLib: N_FitCrvCalcKnotVector."
  (declare (type real-vector u knt) (type fixnum n p))
  (fast
    (dotimes (i (1+ p))
      (setf (aref knt i) (aref u 0)
            (aref knt (+ n i 1)) (aref u n)))
    (do ((i 1 (1+ i))) ((> i (- n p)))
      (declare (type fixnum i))
      (let ((sum 0d0))
        (declare (type double-float sum))
        (do ((j i (1+ j))) ((> j (+ i p -1)))
          (declare (type fixnum j))
          (setf sum (+ sum (aref u j))))
        (setf (aref knt (+ i p)) (/ sum (float p 1d0))))))
  knt)

(defun fit-calc-knot-vector-end-derivs (u k p knt)
  "Write into KNT (K+P+4 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..K] with both end
derivatives given.  Returns KNT.
NLib: N_FitCalcKnotVectorEndDerivs."
  (declare (type real-vector u knt) (type fixnum k p))
  (let ((n (+ k 2)))
    (declare (type fixnum n))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (if (= p 1)
          (progn
            (if (> k 1)
                (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1)))
                      (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k)))))
                (setf (aref knt 2) (+ (* 0.75d0 (aref u 0)) (* 0.25d0 (aref u 1)))
                      (aref knt 3) (+ (* 0.25d0 (aref u 0)) (* 0.75d0 (aref u 1)))))
            (do ((i 1 (1+ i))) ((>= i k))
              (declare (type fixnum i))
              (setf (aref knt (+ i 2)) (aref u i))))
          (do ((i 0 (1+ i))) ((> i (+ (- k p) 1)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p 1)) (/ sum (float p 1d0)))))))
    knt))

(defun fit-calc-knot-vector-deriv (u k p whr knt)
  "Write into KNT (K+P+3 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P at parameters U[0..K] with the end
derivative given at WHR, +start+ or +end+.  Signals :inp-err when P > K+1 or
WHR is neither.  Returns KNT.
NLib: N_FitCalcKnotVectorDeriv."
  (declare (type real-vector u knt) (type fixnum k p whr))
  (when (or (> p (1+ k)) (and (/= whr +start+) (/= whr +end+)))
    (nurbs-error :inp-err 'fit-calc-knot-vector-deriv))
  (let ((n (1+ k)))
    (declare (type fixnum n))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (when (= p 1)                     ; special case out
        (if (= k 1)
            (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1))))
            (progn
              (when (= whr +start+)
                (setf (aref knt 2) (* 0.5d0 (+ (aref u 0) (aref u 1))))
                (do ((i 1 (1+ i))) ((>= i k))
                  (declare (type fixnum i))
                  (setf (aref knt (+ i 2)) (aref u i))))
              (when (= whr +end+)
                (setf (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k)))))
                (do ((i 1 (1+ i))) ((>= i k))
                  (declare (type fixnum i))
                  (setf (aref knt (+ i 1)) (aref u i)))))))
      (when (and (> p 1) (< p n))
        (when (= whr +start+)
          (do ((i 1 (1+ i))) ((> i (- k p)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p 1)) (/ sum (float p 1d0)))))
          (let ((sum 0d0))
            (declare (type double-float sum))
            (do ((j 0 (1+ j))) ((> j (1- p)))
              (declare (type fixnum j))
              (setf sum (+ sum (aref u j))))
            (setf (aref knt (1+ p)) (/ sum (float p 1d0)))))
        (when (= whr +end+)
          (do ((i 1 (1+ i))) ((> i (- k p)))
            (declare (type fixnum i))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j i (1+ j))) ((> j (+ i p -1)))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref knt (+ i p)) (/ sum (float p 1d0)))))
          (let ((sum 0d0))
            (declare (type double-float sum))
            (do ((j (1+ (- k p)) (1+ j))) ((> j k))
              (declare (type fixnum j))
              (setf sum (+ sum (aref u j))))
            (setf (aref knt n) (/ sum (float p 1d0)))))))
    knt))

(defun fit-calc-knot-vector-crv-approx (u m n p knt)
  "Write into KNT (N+P+2 knots, allocated by the caller) the knot vector for
global approximation of samples at ascending parameters U[0..M] by a curve of
degree P with control points 0..N: interior knots average P consecutive
representatives of the N+1 clusters of parameters.  Signals :inp-err when
N > M.  Returns KNT.
NLib: N_FitCalcKnotVectorCrvApprox."
  (declare (type real-vector u knt) (type fixnum m n p))
  (when (> n m)
    (nurbs-error :inp-err 'fit-calc-knot-vector-crv-approx))
  (fast
    (dotimes (i (1+ p))
      (setf (aref knt i) (aref u 0)
            (aref knt (+ n i 1)) (aref u m))))
  (let ((uk (make-real-vector (1+ n)))
        (d (/ (+ (float m 1d0) 1d0) (+ (float n 1d0) 1d0)))
        (il 0) (ih 0) (l -1d0))
    (declare (type real-vector uk) (type double-float d l) (type fixnum il ih))
    (fast
      (do ((i 0 (1+ i))) ((> i n))
        (declare (type fixnum i))
        (setf l (+ l d))
        (setf ih (the fixnum (truncate (+ l 0.5d0)))) ; ROUND(l)
        (if (= il ih)
            (setf (aref uk i) (aref u il))
            (let ((sum 0d0))
              (declare (type double-float sum))
              (do ((j il (1+ j))) ((> j ih))
                (declare (type fixnum j))
                (setf sum (+ sum (aref u j))))
              (setf (aref uk i) (/ sum (+ (- (float ih 1d0) (float il 1d0)) 1d0)))))
        (setf il (1+ ih)))
      (do ((i 1 (1+ i))) ((> i (- n p)))
        (declare (type fixnum i))
        (let ((sum 0d0))
          (declare (type double-float sum))
          (do ((j i (1+ j))) ((> j (+ i p -1)))
            (declare (type fixnum j))
            (setf sum (+ sum (aref uk j))))
          (setf (aref knt (+ i p)) (/ sum (float p 1d0))))))
    knt))

(defun fit-calc-knots-derivs (u k p knt)
  "Write into KNT (2K+P+3 knots, allocated by the caller) the knot vector for
global curve interpolation of degree P (2 or 3; any other degree gets only its
end knots) at parameters U[0..K] with a derivative given at every point.
Returns KNT.
NLib: N_FitCalcKnotsDerivs."
  (declare (type real-vector u knt) (type fixnum k p))
  (let ((n (1+ (* 2 k)))
        (j (1+ p)))
    (declare (type fixnum n j))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ n i 1)) (aref u k)))
      (case p
        (2
         (setf (aref knt (1+ p)) (* 0.5d0 (+ (aref u 1) (aref u 0))))
         (do ((i 1 (1+ i))) ((>= i k))
           (declare (type fixnum i))
           (setf (aref knt (incf j)) (aref u i))
           (setf (aref knt (incf j)) (* 0.5d0 (+ (aref u i) (aref u (1+ i)))))))
        (3
         (when (> k 1)
           (setf (aref knt (1+ p)) (* 0.5d0 (+ (aref u 1) (aref u 0)))
                 (aref knt n) (* 0.5d0 (+ (aref u k) (aref u (1- k))))))
         (let ((alf (/ 1d0 3d0))
               (bet (/ 2d0 3d0)))
           (declare (type double-float alf bet))
           (do ((i 1 (1+ i))) ((>= i (1- k)))
             (declare (type fixnum i))
             (setf (aref knt (incf j)) (+ (* bet (aref u i)) (* alf (aref u (1+ i)))))
             (setf (aref knt (incf j)) (+ (* alf (aref u i)) (* bet (aref u (1+ i))))))))))
    knt))

(defun fit-calc-knots-high-end-derivs (u n p k l knt)
  "Write into KNT (N+K+L+P+2 knots, allocated by the caller) the knot vector
for global curve interpolation of degree P at parameters U[0..N] with the
first K derivatives given at the start and the first L at the end (K, L < P).
Returns KNT.
NLib: N_FitCalcKnotsHighEndDerivs."
  (declare (type real-vector u knt) (type fixnum n p k l))
  (let* ((m (+ n k l))
         (is (- 1 k))
         (ie (+ (- n p) l))
         (r p))
    (declare (type fixnum m is ie r))
    (fast
      (dotimes (i (1+ p))
        (setf (aref knt i) (aref u 0)
              (aref knt (+ m i 1)) (aref u n)))
      (do ((i is (1+ i))) ((> i ie))
        (declare (type fixnum i))
        (let ((js (nl-max 0 i))
              (je (nl-min n (+ i p -1)))
              (sum 0d0))
          (declare (type fixnum js je) (type double-float sum))
          (do ((j js (1+ j))) ((> j je))
            (declare (type fixnum j))
            (setf sum (+ sum (aref u j))))
          (setf (aref knt (incf r)) (/ sum (+ (- (float je 1d0) (float js 1d0)) 1d0))))))
    knt))

;;; ------------------------------------------------------- surface data rows

(defun fit-crv-c-pts-from-srf-data (a n m ptp k dir cm qw)
  "Write into QW (a cpoint-array the caller allocates) the control points of
the curve interpolating row or column K of the surface data A[0..N][0..M]: a
point-net when PTP is +epoint+, a cpoint-net when it is +hpoint+.  DIR +udir+
takes A[i][K], i = 0..N; +vdir+ takes A[K][j], j = 0..M.  CM is the LU
decomposed interpolation coefficient matrix, an rmatrix.  Signals :cal-err
for another PTP or DIR, and whatever REAL-MATRIX-FOR-BACK signals.  Returns
QW.  NLib: N_FitCrvCPtsFromSrfData."
  (declare (type (simple-array double-float (* * *)) a) (type fixnum n m ptp k dir)
           (type rmatrix cm) (type cpoint-array qw))
  (let* ((cnt (1+ (nl-max n m)))
         (dim (cond ((= ptp +epoint+) 3)
                    ((= ptp +hpoint+) 4)
                    (t (nurbs-error :cal-err 'fit-crv-c-pts-from-srf-data))))
         (r (make-array (list cnt dim) :element-type 'double-float :initial-element 0d0)))
    (declare (type fixnum cnt dim) (type (simple-array double-float (* *)) r))
    (cond
      ((= dir +udir+)
       (fast (do ((i 0 (1+ i))) ((> i n))
               (declare (type fixnum i))
               (dotimes (c dim) (setf (aref r i c) (aref a i k c))))))
      ((= dir +vdir+)
       (fast (do ((j 0 (1+ j))) ((> j m))
               (declare (type fixnum j))
               (dotimes (c dim) (setf (aref r j c) (aref a k j c))))))
      (t (nurbs-error :cal-err 'fit-crv-c-pts-from-srf-data)))
    (real-matrix-for-back cm r ptp qw)))

;;; ----------------------------------------------------------- knot removal

(defconstant +crv-fit--norem+ 1d25 "NOREM in N_FitRemoveKnots: a knot found not removable.")

(defun fit-remove-knots (cur u er mu e)
  "Remove from curve CUR, in place, every knot whose removal keeps the
cumulative error at the parameters U[0..MU] within E, cheapest knot first.
ER[0..MU] holds the error bounds so far and is updated in place.  CUR's knot
vector and control points are replaced by shorter arrays when knots go.
Returns CUR.  NLib: N_FitRemoveKnots."
  (declare (type curve cur) (type real-vector u er) (type fixnum mu) (type double-float e))
  (let* ((n (curve-n cur))
         (pw (curve-cpoints cur))
         (p (curve-degree cur))
         (m (curve-m cur))
         (uk (curve-knots cur))
         (knt uk)
         (rw (make-cpoint-array (1+ (* 2 p))))
         (br (make-real-vector (1+ m)))
         (sr (make-index-vector (1+ m)))
         (te (make-real-vector (1+ mu)))
         (left (make-index-vector (1+ n)))
         (right (make-index-vector (1+ n)))
         (n0 n)
         (lam 0d0) (oml 0d0))
    (declare (type fixnum n p m n0) (type cpoint-array pw rw) (type real-vector uk knt br te)
             (type index-vector sr left right) (type double-float lam oml))
    (unwind-protect
         (progn
           ;; initialize
           (dotimes (i (1+ m))
             (setf (aref br i) +bigd+
                   (aref sr i) 0))
           ;; knot removal errors for each distinct knot
           (let ((r (1+ p)))
             (declare (type fixnum r))
             (do () ((not (<= r n)))
               (let ((i r))
                 (declare (type fixnum i))
                 (do () ((not (and (<= r n) (= (aref uk r) (aref uk (1+ r))))))
                   (incf r))
                 (setf (aref sr r) (1+ (- r i)))
                 (setf (aref br r) (crv-remove-knot-max-err cur r (aref sr r)))
                 (incf r))))
           ;; ranges of parameter indexes
           (let ((j 1))
             (declare (type fixnum j))
             (dotimes (i (1+ n))
               (setf (aref left i) j)
               (do () ((not (and (< j mu) (> (aref u j) (aref uk i)) (<= (aref u j) (aref uk (1+ i))))))
                 (incf j))
               (let ((k j))
                 (declare (type fixnum k))
                 (do () ((not (< k mu)))
                   (when (>= (aref u k) (aref uk (+ i p 1))) (return))
                   (incf k))
                 (setf (aref right i) (1- k)))))
           ;; remove knots while the cumulative error stays within E
           (do () (nil)
             (let ((b (aref br (1+ p)))
                   (s (aref sr (1+ p)))
                   (r (1+ p))
                   (rmf t)
                   (i 0))
               (declare (type double-float b) (type fixnum s r i))
               ;; knot with the smallest error
               (do ((ii (+ p 2) (1+ ii))) ((> ii (- m p 1)))
                 (declare (type fixnum ii))
                 (when (< (aref br ii) b)
                   (setf b (aref br ii)
                         s (aref sr ii)
                         r ii)))
               ;; no more removable knots: finished
               (when (or (= b +bigd+) (= b +crv-fit--norem+))
                 (return))
               ;; error of removal
               (if (oddp (+ p s))
                   (let* ((k (truncate (+ p s 1) 2))
                          (i0 (- r k))
                          (alf (/ (- (aref uk r) (aref uk i0)) (- (aref uk (+ i0 p 1)) (aref uk i0))))
                          (bet (/ (- (aref uk r) (aref uk (1+ i0))) (- (aref uk (+ i0 p 2)) (aref uk (1+ i0)))))
                          (omb (- 1d0 bet)))
                     (declare (type fixnum k i0) (type double-float alf bet omb))
                     (setf i i0
                           lam (/ alf (+ alf bet))
                           oml (- 1d0 lam))
                     (do ((j (aref left i) (1+ j))) ((> j (aref right (1+ i))))
                       (declare (type fixnum j))
                       (let ((nb (basis-i-eval knt i p (aref u j) +left+))
                             (nb1 (basis-i-eval knt (1+ i) p (aref u j) +left+)))
                         (declare (type double-float nb nb1))
                         (setf (aref te j) (+ (aref er j)
                                              (* (abs (- (* (* lam alf) nb) (* (* oml omb) nb1))) b)))
                         (when (> (aref te j) e)
                           (setf rmf nil)
                           (return)))))
                   (let ((k (truncate (+ p s) 2)))
                     (declare (type fixnum k))
                     (setf i (- r k))
                     (do ((j (aref left i) (1+ j))) ((> j (aref right i)))
                       (declare (type fixnum j))
                       (let ((nb (basis-i-eval knt i p (aref u j) +left+)))
                         (declare (type double-float nb))
                         (setf (aref te j) (+ (aref er j) (* nb b)))
                         (when (> (aref te j) e)
                           (setf rmf nil)
                           (return))))))
               (cond
                 (rmf
                  ;; removable: update the error and remove the knot
                  (let ((l (if (oddp (+ p s)) (aref right (1+ i)) (aref right i))))
                    (declare (type fixnum l))
                    (do ((j (aref left i) (1+ j))) ((> j l))
                      (declare (type fixnum j))
                      (setf (aref er j) (aref te j))))
                  (let* ((fout (truncate (- (* 2 r) s p) 2))
                         (first (- r p))
                         (last (- r s))
                         (off (1- first))
                         (i first)
                         (j last)
                         (ii 1)
                         (jj (- last off)))
                    (declare (type fixnum fout first last off i j ii jj))
                    (%copy-c-pt-row pw off rw 0)
                    (%copy-c-pt-row pw (1+ last) rw (- (1+ last) off))
                    ;; new control points for one removal step
                    (fast
                      (do () ((not (> (- j i) 0)))
                        (let* ((alf (/ (- (aref uk (+ i p 1)) (aref uk i)) (- (aref uk r) (aref uk i))))
                               (oma (- 1d0 alf))
                               (bet (/ (- (aref uk (+ j p 1)) (aref uk j)) (- (aref uk (+ j p 1)) (aref uk r))))
                               (omb (- 1d0 bet)))
                          (declare (type double-float alf oma bet omb))
                          (%combine-2-c-pt-rows alf pw i oma rw (1- ii) rw ii)
                          (%combine-2-c-pt-rows bet pw j omb rw (1+ jj) rw jj)
                          (incf i) (decf j) (incf ii) (decf jj))))
                    ;; save new control points
                    (when (oddp (+ p s))
                      (%combine-2-c-pt-rows lam rw (1+ jj) oml rw (1- ii) rw (1+ jj)))
                    (setf i first
                          j last)
                    (fast
                      (do () ((not (> (- j i) 0)))
                        (%copy-c-pt-row rw (- i off) pw i)
                        (%copy-c-pt-row rw (- j off) pw j)
                        (incf i) (decf j)))
                    ;; shift down some parameters
                    (when (> s 1)
                      (setf (aref sr (1- r)) (1- (aref sr r))))
                    (do ((i (1+ r) (1+ i))) ((> i m))
                      (declare (type fixnum i))
                      (setf (aref br (1- i)) (aref br i)
                            (aref sr (1- i)) (aref sr i)))
                    ;; shift down knots and control points
                    (fast
                      (do ((i (1+ r) (1+ i))) ((> i m))
                        (declare (type fixnum i))
                        (setf (aref uk (1- i)) (aref uk i)))
                      (do ((i (1+ fout) (1+ i))) ((> i n))
                        (declare (type fixnum i))
                        (%copy-c-pt-row pw i pw (1- i))))
                    (decf n)
                    (decf m)
                    ;; the knot vector N_BasisIEval sees has highest index m
                    (setf knt (subseq uk 0 (1+ m))))
                  ;; no more internal knots: finished
                  (when (= n p)
                    (return))
                  ;; update error bounds
                  (let ((k (nl-max (- r p) (1+ p)))
                        (l (nl-min n (- (+ r p) s))))
                    (declare (type fixnum k l))
                    (do ((i k (1+ i))) ((> i l))
                      (declare (type fixnum i))
                      (when (and (/= (aref uk i) (aref uk (1+ i))) (/= (aref br i) +crv-fit--norem+))
                        (setf (aref br i) (crv-remove-knot-max-err cur i (aref sr i))))))
                  ;; update index ranges
                  (do ((i (- r p 1) (1+ i))) ((> i (- r s)))
                    (declare (type fixnum i))
                    (let ((k (1+ (aref right i))))
                      (declare (type fixnum k))
                      (do () ((not (<= k mu)))
                        (when (>= (aref u k) (aref uk (+ i p 1))) (return))
                        (incf k))
                      (setf (aref right i) (1- k))))
                  (do ((i (1+ (- r s)) (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (setf (aref left i) (aref left (1+ i))
                          (aref right i) (aref right (1+ i)))))
                 (t
                  ;; not removable
                  (setf (aref br r) +crv-fit--norem+))))))
      ;; CUR keeps the knots and control points left, as N_CrvSetSizeIndices
      ;; records them in C
      (when (/= n n0)
        (setf (curve-knots cur) (subseq uk 0 (1+ m))
              (curve-cpoints cur) (let ((q (make-cpoint-array (1+ n))))
                                    (dotimes (i (1+ n) q)
                                      (%copy-c-pt-row pw i q i)))))))
  cur)
