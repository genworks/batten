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
                    (setf (aref dist i) (c-sqrt (aref dist i)))
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
                    (setf (aref dist i) (c-sqrt (aref dist i)))
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

;;; ------------------------------------------------- global curve interpolation
;;;
;;; The output curves and functions (N_CrvSizeArrays / N_CFuncSizeArrays on a
;;; NULL object, "memory is allocated locally") are returned as fresh objects.
;;; NLib keeps a matrix row in N_AllocReal2dArray's one contiguous block, and
;;; several of these routines write or read a row at a column index that may
;;; step past either end of it (N_FitCrvInterp, N_FitCrvFirstDeriv ... check
;;; no band limits); %CRV-FIT-RM-SET and %CRV-FIT-RM-REF address a real-matrix
;;; the same way, so such an index lands where it lands in C.

(declaim (inline %crv-fit-rm-set %crv-fit-rm-ref))
(defun %crv-fit-rm-set (a i j v)
  "A[i][j] = V, addressing real-matrix A as C addresses its contiguous rows: a
column index past either end of row I reaches into the neighboring row.  The
index is checked against the whole array."
  (declare (type real-matrix a) (type fixnum i j) (type double-float v))
  (setf (row-major-aref a (+ (the fixnum (* i (array-dimension a 1))) j)) v))

(defun %crv-fit-rm-ref (a i j)
  "A[i][j], addressing real-matrix A as C addresses its contiguous rows (see
%CRV-FIT-RM-SET)."
  (declare (type real-matrix a) (type fixnum i j))
  (row-major-aref a (+ (the fixnum (* i (array-dimension a 1))) j)))

(defun fit-func-interp-given-params (f n p u)
  "The curve function of degree P interpolating the values F[0..N] at the
parameters U[0..N], its knots averaged from U (FIT-CRV-CALC-KNOT-VECTOR), as
a fresh cfun.  Signals :inp-err when P > N, and whatever BASIS-EVAL and the
banded solver signal.
NLib: N_FitFuncInterpGivenParams."
  (declare (type real-vector f u) (type fixnum n p))
  (when (> p n)
    (nurbs-error :inp-err 'fit-func-interp-given-params))
  (let ((cfn (%make-cfun p (make-real-vector (+ n p 2)) (make-real-vector (1+ n)))))
    (multiple-value-bind (fu knt) (crv-func-cntrl-val-knot-vector cfn)
      (declare (type real-vector fu knt))
      (fit-crv-calc-knot-vector u n p knt)
      ;; linear splines
      (when (= p 1)
        (dotimes (i (1+ n))
          (setf (aref fu i) (aref f i)))
        (return-from fit-func-interp-given-params cfn))
      ;; coefficient matrix
      (let* ((bw (1- (* 2 p)))
             (sbw (1- p))
             (nb (make-real-vector (1+ p)))
             (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
             (a (rmatrix-rm cm)))
        (declare (type fixnum bw sbw) (type real-vector nb) (type real-matrix a))
        (setf (aref a 0 sbw) 1d0
              (aref a n sbw) 1d0)
        (do ((i 1 (1+ i))) ((>= i n))
          (declare (type fixnum i))
          (let* ((j (basis-eval knt p (aref u i) +left+ nb))
                 (l (- j i 1)))
            (declare (type fixnum j l))
            (dotimes (k (1+ p))
              (%crv-fit-rm-set a i (+ l k) (aref nb k)))))
        (real-matrix-lu-decompose cm)
        (real-matrix-right-for-back cm f fu)))
    cfn))

(defun fit-crv-interp (pts n p par)
  "The curve of degree P interpolating the points PTS[0..N] (a point-array),
parametrized by PAR (+uniform+, +chordlength+ or +centripetal+), its knots
averaged from the parameters, as a fresh non-rational curve.  Signals :inp-err
when P > N, :cal-err for another PAR, and whatever BASIS-EVAL and the banded
solver signal.
NLib: N_FitCrvInterp."
  (declare (type point-array pts) (type fixnum n p par))
  (when (> p n)
    (nurbs-error :inp-err 'fit-crv-interp))
  (let* ((pw (make-cpoint-array (1+ n)))
         (knt (make-real-vector (+ n p 2)))
         (cur (%make-curve p knt pw))
         (u (make-real-vector (1+ n))))
    (declare (type cpoint-array pw) (type real-vector knt u))
    (fit-calc-crv-param-values pts n +epoint+ par u)
    (fit-crv-calc-knot-vector u n p knt)
    ;; the control points are the data points of a linear curve
    (when (= p 1)
      (let ((s (make-point)) (c (make-cpoint)))
        (declare (dynamic-extent s c))
        (dotimes (i (1+ n))
          (pt-to-c-pt (%row-pt pts i s) c)
          (%c-pt-row c pw i)))
      (return-from fit-crv-interp cur))
    ;; coefficient matrix
    (let* ((bw (1- (* 2 p)))
           (sbw (1- p))
           (nb (make-real-vector (1+ p)))
           (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
           (a (rmatrix-rm cm)))
      (declare (type fixnum bw sbw) (type real-vector nb) (type real-matrix a))
      (setf (aref a 0 sbw) 1d0
            (aref a n sbw) 1d0)
      (do ((i 1 (1+ i))) ((>= i n))
        (declare (type fixnum i))
        (let* ((j (basis-eval knt p (aref u i) +left+ nb))
               (l (- j i 1)))
          (declare (type fixnum j l))
          (dotimes (k (1+ p))
            (%crv-fit-rm-set a i (+ l k) (aref nb k)))))
      (real-matrix-lu-decompose cm)
      (real-matrix-for-back cm pts +epoint+ pw))
    cur))

(defun %crv-fit-interp-band (n p)
  "The bandwidth and semi-bandwidth N_FitCrvInterpGivenParams and
N_FitCrvKnotsAndDerivs give their n x n matrix for degree P: 2p+1 and p when
2p <= n, else 2p-1 and p-1."
  (declare (type fixnum n p))
  (if (<= (* 2 p) n)
      (values (1+ (* 2 p)) p)
      (values (1- (* 2 p)) (1- p))))

(defun fit-crv-interp-given-params (a n ptp u knt p)
  "The curve of degree P on the knot vector KNT interpolating A[0..N] at the
parameters U[0..N], as a fresh curve: A is a point-array when PTP is
+epoint+ (a non-rational curve), a cpoint-array when it is +hpoint+.  Signals
:inp-err when KNT does not have N+P+2 knots or a parameter falls outside the
matrix band, :cal-err for another PTP, and whatever BASIS-EVAL and the banded
solver signal.
NLib: N_FitCrvInterpGivenParams."
  (declare (type (simple-array double-float (* *)) a) (type fixnum n ptp p)
           (type real-vector u knt))
  (let ((m (1- (length knt)))
        (mc (+ n p 1)))
    (declare (type fixnum m mc))
    (when (/= m mc)
      (nurbs-error :inp-err 'fit-crv-interp-given-params))
    (let* ((pw (make-cpoint-array (1+ n)))
           (cur (%make-curve p (copy-seq knt) pw)))
      (declare (type cpoint-array pw))
      ;; the control points are the data points of a linear curve
      (when (= p 1)
        (cond ((= ptp +epoint+)
               (let ((s (make-point)) (c (make-cpoint)))
                 (declare (dynamic-extent s c))
                 (dotimes (i (1+ n))
                   (pt-to-c-pt (%row-pt a i s) c)
                   (%c-pt-row c pw i))))
              ((= ptp +hpoint+)
               (dotimes (i (1+ n))
                 (%copy-c-pt-row a i pw i)))
              (t (nurbs-error :cal-err 'fit-crv-interp-given-params)))
        (return-from fit-crv-interp-given-params cur))
      ;; coefficient matrix
      (multiple-value-bind (bw sbw) (%crv-fit-interp-band n p)
        (declare (type fixnum bw sbw))
        (let* ((sp (- sbw p))
               (nb (make-real-vector (1+ p)))
               (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
               (b (rmatrix-rm cm)))
          (declare (type fixnum sp) (type real-vector nb) (type real-matrix b))
          (setf (aref b 0 sbw) 1d0
                (aref b n sbw) 1d0)
          (do ((i 1 (1+ i))) ((>= i n))
            (declare (type fixnum i))
            (let* ((j (basis-eval knt p (aref u i) +left+ nb))
                   (l (+ (- sp i) j)))
              (declare (type fixnum j l))
              (when (or (< l 0) (>= (+ l p) bw))
                (nurbs-error :inp-err 'fit-crv-interp-given-params))
              (dotimes (k (1+ p))
                (setf (aref b i (+ l k)) (aref nb k)))))
          (real-matrix-lu-decompose cm)
          (real-matrix-for-back cm a ptp pw)))
      cur)))

;;; ------------------------------------------- interpolation with end derivatives

(defun fit-crv-knots-and-derivs (a kk ptp u knt p as ae)
  "The curve of degree P on the knot vector KNT interpolating A[0..KK] at the
parameters U[0..KK] with start and end derivatives AS and AE, as a fresh
curve of KK+3 control points: A is a point-array and AS, AE points when PTP
is +epoint+ (a non-rational curve), a cpoint-array and cpoints when it is
+hpoint+.  Signals :inp-err when P < 2, P > KK+2, KNT does not have KK+P+4
knots or a parameter falls outside the matrix band, :cal-err for another PTP,
and whatever BASIS-EVAL and the banded solver signal.
NLib: N_FitCrvKnotsAndDerivs."
  (declare (type (simple-array double-float (* *)) a) (type fixnum kk ptp p)
           (type real-vector u knt) (type (simple-array double-float (*)) as ae))
  (let* ((m (1- (length knt)))
         (n (+ kk 2))
         (mc (+ n p 1)))
    (declare (type fixnum m n mc))
    (when (or (< p 2) (> p (+ kk 2)))
      (nurbs-error :inp-err 'fit-crv-knots-and-derivs))
    (when (/= m mc)
      (nurbs-error :inp-err 'fit-crv-knots-and-derivs))
    (let* ((pw (make-cpoint-array (1+ n)))
           (cur (%make-curve p (copy-seq knt) pw)))
      (declare (type cpoint-array pw))
      ;; coefficient matrix
      (multiple-value-bind (bw sbw) (%crv-fit-interp-band n p)
        (declare (type fixnum bw sbw))
        (let* ((sp (- sbw p))
               (nb (make-real-vector (1+ p)))
               (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
               (b (rmatrix-rm cm)))
          (declare (type fixnum sp) (type real-vector nb) (type real-matrix b))
          (setf (aref b 0 sbw) 1d0
                (aref b 1 (1- sbw)) -1d0
                (aref b 1 sbw) 1d0
                (aref b (1- n) sbw) -1d0
                (aref b (1- n) (1+ sbw)) 1d0
                (aref b n sbw) 1d0)
          (do ((i 2 (1+ i))) ((> i kk))
            (declare (type fixnum i))
            (let* ((j (basis-eval knt p (aref u (1- i)) +left+ nb))
                   (l (+ (- sp i) j)))
              (declare (type fixnum j l))
              (when (or (< l 0) (>= (+ l p) bw))
                (nurbs-error :inp-err 'fit-crv-knots-and-derivs))
              (dotimes (k (1+ p))
                (setf (aref b i (+ l k)) (aref nb k)))))
          (real-matrix-lu-decompose cm)
          ;; solve
          (let ((pf (float p 1d0)))
            (declare (type double-float pf))
            (cond
              ((= ptp +epoint+)
               (let ((q (make-point-array (1+ n))) (s (make-point)))
                 (declare (type point-array q) (dynamic-extent s))
                 (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
                   (%copy-pt-row a 0 q 0)
                   (%pt-row (scale-pt fact as s) q 1))
                 (do ((i 2 (1+ i))) ((> i kk))
                   (declare (type fixnum i))
                   (%copy-pt-row a (1- i) q i))
                 (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
                   (%pt-row (scale-pt fact ae s) q (1- n))
                   (%copy-pt-row a kk q n))
                 (real-matrix-for-back cm q +epoint+ pw)))
              ((= ptp +hpoint+)
               (let ((qw (make-cpoint-array (1+ n))) (s (make-cpoint)))
                 (declare (type cpoint-array qw) (dynamic-extent s))
                 (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
                   (%copy-c-pt-row a 0 qw 0)
                   (%c-pt-row (scale-c-pt fact as s) qw 1))
                 (do ((i 2 (1+ i))) ((> i kk))
                   (declare (type fixnum i))
                   (%copy-c-pt-row a (1- i) qw i))
                 (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
                   (%c-pt-row (scale-c-pt fact ae s) qw (1- n))
                   (%copy-c-pt-row a kk qw n))
                 (real-matrix-for-back cm qw +hpoint+ pw)))
              (t (nurbs-error :cal-err 'fit-crv-knots-and-derivs))))))
      cur)))

(defun fit-crv-knots-and-deriv (a kk ptp u knt p tv whr)
  "The curve of degree P on the knot vector KNT interpolating A[0..KK] at the
parameters U[0..KK] with the derivative TV at the start (WHR +start+) or the
end (+end+), as a fresh curve of KK+2 control points: A is a point-array and
TV a point when PTP is +epoint+ (a non-rational curve), a cpoint-array and a
cpoint when it is +hpoint+.  Signals :inp-err when P < 2, P > KK+1, KNT does
not have KK+P+3 knots, PTP or WHR is another value or a parameter falls
outside the matrix band, and whatever BASIS-EVAL and the banded solver signal.
NLib: N_FitCrvKnotsAndDeriv."
  (declare (type (simple-array double-float (* *)) a) (type fixnum kk ptp p whr)
           (type real-vector u knt) (type (simple-array double-float (*)) tv))
  (let* ((m (1- (length knt)))
         (n (1+ kk))
         (mc (+ n p 1)))
    (declare (type fixnum m n mc))
    (when (or (< p 2) (> p n))
      (nurbs-error :inp-err 'fit-crv-knots-and-deriv))
    (when (/= m mc)
      (nurbs-error :inp-err 'fit-crv-knots-and-deriv))
    (unless (or (= ptp +epoint+) (= ptp +hpoint+))
      (nurbs-error :inp-err 'fit-crv-knots-and-deriv))
    (unless (or (= whr +start+) (= whr +end+))
      (nurbs-error :inp-err 'fit-crv-knots-and-deriv))
    (let* ((pw (make-cpoint-array (1+ n)))
           (cur (%make-curve p (copy-seq knt) pw))
           (bw (1- (* 2 p)))
           (sbw (1- p))
           (nb (make-real-vector (1+ p)))
           (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
           (b (rmatrix-rm cm)))
      (declare (type cpoint-array pw) (type fixnum bw sbw) (type real-vector nb)
               (type real-matrix b))
      (flet ((row (i ui)
               (declare (type fixnum i) (type double-float ui))
               (let* ((j (basis-eval knt p ui +left+ nb))
                      (l (if (= p 2) 0 (- j i 1))))
                 (declare (type fixnum j l))
                 (when (or (< l 0) (>= (+ l p) bw))
                   (nurbs-error :inp-err 'fit-crv-knots-and-deriv))
                 (dotimes (k (1+ p))
                   (setf (aref b i (+ l k)) (aref nb k))))))
        (if (= whr +start+)
            (progn
              (setf (aref b 0 sbw) 1d0
                    (aref b 1 (1- sbw)) -1d0
                    (aref b 1 sbw) 1d0
                    (aref b n sbw) 1d0)
              (do ((i 2 (1+ i))) ((> i kk))
                (declare (type fixnum i))
                (row i (aref u (1- i)))))
            (progn
              (setf (aref b 0 sbw) 1d0
                    (aref b (1- n) sbw) -1d0
                    (aref b (1- n) (1+ sbw)) 1d0
                    (aref b n sbw) 1d0)
              (do ((i 1 (1+ i))) ((> i (1- kk)))
                (declare (type fixnum i))
                (row i (aref u i))))))
      (real-matrix-lu-decompose cm)
      ;; solve
      (let ((pf (float p 1d0)))
        (declare (type double-float pf))
        (if (= ptp +epoint+)
            (let ((q (make-point-array (1+ n))) (s (make-point)))
              (declare (type point-array q) (dynamic-extent s))
              (if (= whr +start+)
                  (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
                    (%copy-pt-row a 0 q 0)
                    (%pt-row (scale-pt fact tv s) q 1)
                    (do ((i 2 (1+ i))) ((> i n))
                      (declare (type fixnum i))
                      (%copy-pt-row a (1- i) q i)))
                  (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
                    (%pt-row (scale-pt fact tv s) q (1- n))
                    (%copy-pt-row a kk q n)
                    (do ((i 0 (1+ i))) ((> i (- n 2)))
                      (declare (type fixnum i))
                      (%copy-pt-row a i q i))))
              (real-matrix-for-back cm q +epoint+ pw))
            (let ((qw (make-cpoint-array (1+ n))) (s (make-cpoint)))
              (declare (type cpoint-array qw) (dynamic-extent s))
              (if (= whr +start+)
                  (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
                    (%copy-c-pt-row a 0 qw 0)
                    (%c-pt-row (scale-c-pt fact tv s) qw 1)
                    (do ((i 2 (1+ i))) ((> i n))
                      (declare (type fixnum i))
                      (%copy-c-pt-row a (1- i) qw i)))
                  (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
                    (%c-pt-row (scale-c-pt fact tv s) qw (1- n))
                    (%copy-c-pt-row a kk qw n)
                    (do ((i 0 (1+ i))) ((> i (- n 2)))
                      (declare (type fixnum i))
                      (%copy-c-pt-row a i qw i))))
              (real-matrix-for-back cm qw +hpoint+ pw))))
      cur)))

;;; ------------------------------------------------- least-squares approximation

(defun fit-crv-approx-knots (a m ptp u knt n p)
  "The least-squares curve of degree P with control points 0..N on the knot
vector KNT approximating A[0..M] at the parameters U[0..M], its end points
interpolating A[0] and A[M], as a fresh curve: A is a point-array when PTP is
+epoint+ (a non-rational curve), a cpoint-array otherwise.  KNT holds at
least N+P+2 knots, the curve's.  When N < 2 the curve's control points are
A[0] and A[1] (a cpoint-array unless PTP is +epoint+).  Signals :inp-err when
N < P, KNT does not have N+P+2 knots, M < N or a span has no data, :cal-err
for another PTP, and whatever BASIS-EVAL, BASIS-I-EVAL and the banded solver
signal.
NLib: N_FitCrvApproxKnots."
  (declare (type (simple-array double-float (* *)) a) (type fixnum m ptp n p)
           (type real-vector u knt))
  (when (< n p)
    (nurbs-error :inp-err 'fit-crv-approx-knots))
  (let* ((mk (1- (length knt)))
         (mc (+ n p 1))
         (pw (make-cpoint-array (1+ n)))
         (cur (%make-curve p (subseq knt 0 (1+ mc)) pw)))
    (declare (type fixnum mk mc) (type cpoint-array pw))
    ;; just a line between two points
    (when (< n 2)
      (if (= ptp +epoint+)
          (let ((s (make-point)) (c (make-cpoint)))
            (declare (dynamic-extent s c))
            (%c-pt-row (pt-to-c-pt (%row-pt a 0 s) c) pw 0)
            (%c-pt-row (pt-to-c-pt (%row-pt a 1 s) c) pw 1))
          (progn
            (%copy-c-pt-row a 0 pw 0)
            (%copy-c-pt-row a 1 pw 1)))
      (return-from fit-crv-approx-knots cur))
    (when (/= mk mc)
      (nurbs-error :inp-err 'fit-crv-approx-knots))
    (when (or (< n 2) (> p n) (< m n))
      (nurbs-error :inp-err 'fit-crv-approx-knots))
    ;; coefficient matrix
    (let* ((bw (1+ (* 2 p)))
           (sbw p)
           (rj p)
           (sj (1- p))
           (ej -2)
           (b (make-real-matrix (1- m) (1+ p)))
           (nb (make-real-vector (1+ p)))
           (index (make-index-vector (1- m)))
           (start (make-index-vector (1- n)))
           (end (make-index-vector (1- n)))
           (cm (set-real-matrix (make-rmatrix) (- n 2) (- n 2) +mt-banded+ bw))
           (ntn (rmatrix-rm cm)))
      (declare (type fixnum bw sbw rj sj ej) (type real-matrix b ntn) (type real-vector nb)
               (type index-vector index start end))
      (do ((i 0 (1+ i))) ((> i (nl-min (1- p) (- n 2))))
        (declare (type fixnum i))
        (setf (aref start i) 0))
      (setf (aref end 0) -2)
      (do ((i 1 (1+ i))) ((> i (1- m)))
        (declare (type fixnum i))
        (let* ((j (basis-eval knt p (aref u i) +left+ nb))
               (l (if (= j p) 1 0))
               (hk (if (or (= j p) (= j n)) (1- p) p)))
          (declare (type fixnum j l hk))
          (do ((k 0 (1+ k))) ((> k hk))
            (declare (type fixnum k))
            (setf (aref b (1- i) k) (aref nb (+ l k))))
          (setf (aref index (1- i)) (nl-max 0 (- j p 1)))
          (when (> j rj)
            (do ((k 1 (1+ k))) ((> k (- j rj)))
              (declare (type fixnum k))
              (incf sj)
              (incf ej)
              (when (<= sj (- n 2))
                (setf (aref start sj) (1- i)))
              (when (>= ej 0)
                (setf (aref end ej) (- i 2))))
            (setf rj j))))
      (when (or (< sj (- n 2)) (= (aref end 0) -1))
        (nurbs-error :inp-err 'fit-crv-approx-knots))
      (do ((i (nl-max 0 (1+ ej)) (1+ i))) ((> i (- n 2)))
        (declare (type fixnum i))
        (setf (aref end i) (- m 2)))
      (do ((i 0 (1+ i))) ((> i (- n 2)))
        (declare (type fixnum i))
        (let ((lj (nl-max 0 (- i p)))
              (hj (nl-min (- n 2) (+ i p))))
          (declare (type fixnum lj hj))
          (do ((j lj (1+ j))) ((> j hj))
            (declare (type fixnum j))
            (let ((lk (nl-max (aref start i) (aref start j)))
                  (hk (nl-min (aref end i) (aref end j))))
              (declare (type fixnum lk hk))
              (do ((k lk (1+ k))) ((> k hk))
                (declare (type fixnum k))
                (setf (aref ntn i (+ (- j i) sbw))
                      (+ (aref ntn i (+ (- j i) sbw))
                         (* (%crv-fit-rm-ref b k (- i (aref index k)))
                            (%crv-fit-rm-ref b k (- j (aref index k)))))))))))
      (real-matrix-lu-decompose cm)
      ;; set up and solve the system of equations
      (cond
        ((= ptp +epoint+)
         (let ((s (make-point)) (c (make-cpoint)) (q0 (make-point)) (qm (make-point)) (t1 (make-point))
               (rk (make-point-array (1- m)))
               (r (make-point-array (1- n)))
               (sol (make-cpoint-array (1- n))))
           (declare (dynamic-extent s c q0 qm t1) (type point-array rk r) (type cpoint-array sol))
           (%c-pt-row (pt-to-c-pt (%row-pt a 0 s) c) pw 0)
           (%c-pt-row (pt-to-c-pt (%row-pt a m s) c) pw n)
           (%row-pt a 0 q0)
           (%row-pt a m qm)
           (do ((i 1 (1+ i))) ((> i (1- m)))
             (declare (type fixnum i))
             (let ((n0 (basis-i-eval knt 0 p (aref u i) +left+))
                   (np (basis-i-eval knt n p (aref u i) +left+)))
               (declare (type double-float n0 np))
               (%pt-row (translate-sum-2-pts (%row-pt a i s) (- n0) q0 (- np) qm t1) rk (1- i))))
           (do ((i 1 (1+ i))) ((> i (1- n)))
             (declare (type fixnum i))
             (set-point s 0d0 0d0 0d0)
             (do ((k (aref start (1- i)) (1+ k))) ((> k (aref end (1- i))))
               (declare (type fixnum k))
               (vector-blend-pt (%crv-fit-rm-ref b k (- i (aref index k) 1)) (%row-pt rk k t1) s))
             (%pt-row s r (1- i)))
           (real-matrix-for-back cm r +epoint+ sol)
           (do ((i 0 (1+ i))) ((> i (- n 2)))
             (declare (type fixnum i))
             (%copy-c-pt-row sol i pw (1+ i)))))
        ((= ptp +hpoint+)
         (let ((s (make-cpoint)) (q0 (make-cpoint)) (qm (make-cpoint)) (t1 (make-cpoint))
               (rkw (make-cpoint-array (1- m)))
               (rw (make-cpoint-array (1- n)))
               (sol (make-cpoint-array (1- n))))
           (declare (dynamic-extent s q0 qm t1) (type cpoint-array rkw rw sol))
           (%copy-c-pt-row a 0 pw 0)
           (%copy-c-pt-row a m pw n)
           (%row-c-pt a 0 q0)
           (%row-c-pt a m qm)
           (do ((i 1 (1+ i))) ((> i (1- m)))
             (declare (type fixnum i))
             (let ((n0 (basis-i-eval knt 0 p (aref u i) +left+))
                   (np (basis-i-eval knt n p (aref u i) +left+)))
               (declare (type double-float n0 np))
               (%c-pt-row (translate-sum-2-c-pts (%row-c-pt a i s) (- n0) q0 (- np) qm t1) rkw (1- i))))
           (do ((i 1 (1+ i))) ((> i (1- n)))
             (declare (type fixnum i))
             (set-cpoint s 0d0 0d0 0d0 0d0)
             (do ((k (aref start (1- i)) (1+ k))) ((> k (aref end (1- i))))
               (declare (type fixnum k))
               (vector-blend-c-pt (%crv-fit-rm-ref b k (- i (aref index k) 1)) (%row-c-pt rkw k t1) s))
             (%c-pt-row s rw (1- i)))
           (real-matrix-for-back cm rw +hpoint+ sol)
           (do ((i 0 (1+ i))) ((> i (- n 2)))
             (declare (type fixnum i))
             (%copy-c-pt-row sol i pw (1+ i)))))
        (t (nurbs-error :cal-err 'fit-crv-approx-knots))))
    cur))

;;; ------------------------------------------ interpolation with derivatives given

(defun %crv-fit-deriv-matrix (knt p u kk check rname)
  "The LU-decomposed banded matrix of N_FitCrvFirstDeriv and
N_FitCrvFirstDerivAndKnots: rows 2i and 2i+1 (0 < i < KK) the basis functions
of degree P on KNT at U[i] and their first derivatives.  With CHECK, a band
index out of range signals :inp-err in routine RNAME (N_FitCrvFirstDerivAndKnots
checks; N_FitCrvFirstDeriv does not)."
  (declare (type real-vector knt u) (type fixnum p kk))
  (let* ((n (1+ (* 2 kk)))
         (bw (1- (* 2 p)))
         (sbw (1- p))
         (nd (make-real-matrix 2 (1+ p)))
         (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
         (a (rmatrix-rm cm))
         (hp (if (= p 3) p (1- p))))
    (declare (type fixnum n bw sbw hp) (type real-matrix nd a))
    (%crv-fit-rm-set a 0 sbw 1d0)
    (%crv-fit-rm-set a 1 (1- sbw) -1d0)
    (%crv-fit-rm-set a 1 sbw 1d0)
    (%crv-fit-rm-set a (1- n) sbw -1d0)
    (%crv-fit-rm-set a (1- n) (1+ sbw) 1d0)
    (%crv-fit-rm-set a n sbw 1d0)
    (do ((i 1 (1+ i))) ((>= i kk))
      (declare (type fixnum i))
      (let* ((j (basis-derivs knt p (aref u i) +left+ 1 nd))
             (ii (* 2 i))
             (l (- j ii 1)))
        (declare (type fixnum j ii l))
        (when (and check (or (< l 0) (>= (+ l hp) bw)))
          (nurbs-error :inp-err rname))
        (do ((k 0 (1+ k))) ((> k hp))
          (declare (type fixnum k))
          (%crv-fit-rm-set a ii (+ l k) (aref nd 0 k))
          (%crv-fit-rm-set a (1+ ii) (+ l k -1) (aref nd 1 k)))))
    (real-matrix-lu-decompose cm)))

(defun fit-crv-first-deriv (pts d kk p vec par)
  "The curve of degree P (2 or 3) interpolating the points PTS[0..KK] and the
first derivatives D[0..KK] (point-arrays), parametrized by PAR, as a fresh
non-rational curve of 2KK+2 control points.  With VEC +tangent+ the vectors
of D are directions: each is normalized and scaled by the length of the data
polygon, IN PLACE in D, as the C does; any other VEC takes D as it is.
Signals :deg-err for another degree, :cal-err for another PAR, :geo-err for a
zero tangent, and whatever BASIS-DERIVS and the banded solver signal.
NLib: N_FitCrvFirstDeriv."
  (declare (type point-array pts d) (type fixnum kk p vec par))
  (when (or (< p 2) (> p 3))
    (nurbs-error :deg-err 'fit-crv-first-deriv))
  (let* ((n (1+ (* 2 kk)))
         (pw (make-cpoint-array (1+ n)))
         (knt (make-real-vector (+ n p 2)))
         (cur (%make-curve p knt pw))
         (u (make-real-vector (1+ kk))))
    (declare (type fixnum n) (type cpoint-array pw) (type real-vector knt u))
    (fit-calc-crv-param-values pts kk +epoint+ par u)
    (fit-calc-knots-derivs u kk p knt)
    (let ((cm (%crv-fit-deriv-matrix knt p u kk nil 'fit-crv-first-deriv))
          (s (make-point)))
      (declare (dynamic-extent s))
      ;; magnitudes, if needed
      (when (= vec +tangent+)
        (let ((alf (dist-polygon pts kk)))
          (declare (type double-float alf))
          (do ((i 0 (1+ i))) ((> i kk))
            (declare (type fixnum i))
            (%row-pt d i s)
            (vector-normalize-ref s)
            (%pt-row (vector-scale s alf s) d i))))
      ;; solve
      (let ((q (make-point-array (1+ n)))
            (pf (float p 1d0)))
        (declare (type point-array q) (type double-float pf))
        (let ((fact (/ (aref knt (1+ p)) pf)))
          (%copy-pt-row pts 0 q 0)
          (%pt-row (scale-pt fact (%row-pt d 0 s) s) q 1))
        (do ((i 1 (1+ i))) ((>= i kk))
          (declare (type fixnum i))
          (%copy-pt-row pts i q (* 2 i))
          (%copy-pt-row d i q (1+ (* 2 i))))
        (let ((fact (/ (- 1d0 (aref knt n)) pf)))
          (%pt-row (scale-pt fact (%row-pt d kk s) s) q (1- n))
          (%copy-pt-row pts kk q n))
        (real-matrix-for-back cm q +epoint+ pw)))
    cur))

(defun fit-crv-first-deriv-and-knots (a b kk ptp u knt p)
  "The curve of degree P (2 or 3) on the knot vector KNT interpolating the
points A[0..KK] and first derivatives B[0..KK] at the parameters U[0..KK], as
a fresh curve of 2KK+2 control points: A and B are point-arrays when PTP is
+epoint+ (a non-rational curve), cpoint-arrays when it is +hpoint+.  Signals
:deg-err for another degree, :inp-err when KNT does not have 2KK+P+3 knots or
a parameter falls outside the matrix band, :cal-err for another PTP, and
whatever BASIS-DERIVS and the banded solver signal.
NLib: N_FitCrvFirstDerivAndKnots."
  (declare (type (simple-array double-float (* *)) a b) (type fixnum kk ptp p)
           (type real-vector u knt))
  (let* ((m (1- (length knt)))
         (n (1+ (* 2 kk)))
         (mc (+ n p 1)))
    (declare (type fixnum m n mc))
    (when (or (< p 2) (> p 3))
      (nurbs-error :deg-err 'fit-crv-first-deriv-and-knots))
    (when (/= m mc)
      (nurbs-error :inp-err 'fit-crv-first-deriv-and-knots))
    (let* ((pw (make-cpoint-array (1+ n)))
           (cur (%make-curve p (copy-seq knt) pw))
           (cm (%crv-fit-deriv-matrix knt p u kk t 'fit-crv-first-deriv-and-knots))
           (pf (float p 1d0)))
      (declare (type cpoint-array pw) (type double-float pf))
      (cond
        ((= ptp +epoint+)
         (let ((q (make-point-array (1+ n))) (s (make-point)))
           (declare (type point-array q) (dynamic-extent s))
           (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
             (%copy-pt-row a 0 q 0)
             (%pt-row (scale-pt fact (%row-pt b 0 s) s) q 1))
           (do ((i 1 (1+ i))) ((>= i kk))
             (declare (type fixnum i))
             (%copy-pt-row a i q (* 2 i))
             (%copy-pt-row b i q (1+ (* 2 i))))
           (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
             (%pt-row (scale-pt fact (%row-pt b kk s) s) q (1- n))
             (%copy-pt-row a kk q n))
           (real-matrix-for-back cm q +epoint+ pw)))
        ((= ptp +hpoint+)
         (let ((qw (make-cpoint-array (1+ n))) (s (make-cpoint)))
           (declare (type cpoint-array qw) (dynamic-extent s))
           (let ((fact (/ (- (aref knt (1+ p)) (aref knt 1)) pf)))
             (%copy-c-pt-row a 0 qw 0)
             (%c-pt-row (scale-c-pt fact (%row-c-pt b 0 s) s) qw 1))
           (do ((i 1 (1+ i))) ((>= i kk))
             (declare (type fixnum i))
             (%copy-c-pt-row a i qw (* 2 i))
             (%copy-c-pt-row b i qw (1+ (* 2 i))))
           (let ((fact (/ (- (aref knt (+ n p)) (aref knt n)) pf)))
             (%c-pt-row (scale-c-pt fact (%row-c-pt b kk s) s) qw (1- n))
             (%copy-c-pt-row a kk qw n))
           (real-matrix-for-back cm qw +hpoint+ pw)))
        (t (nurbs-error :cal-err 'fit-crv-first-deriv-and-knots)))
      cur)))

(defun fit-crv-high-derivs (pts nn p ds kk de ll par)
  "The curve of degree P interpolating the points PTS[0..NN], the start
derivatives DS[1..KK] and the end derivatives DE[1..LL] (point-arrays; row 0
of DS and DE is not read), parametrized by PAR, as a fresh non-rational curve
of NN+KK+LL+1 control points.  Signals :inp-err when P < 1, P > NN+KK+LL,
KK or LL >= P or either is negative, :cal-err for another PAR, and whatever
BASIS-DERIVS, BASIS-EVAL and the banded solver signal.
NLib: N_FitCrvHighDerivs."
  (declare (type point-array pts ds de) (type fixnum nn p kk ll par))
  (let ((n (+ nn kk ll)))
    (declare (type fixnum n))
    (when (or (< p 1) (> p n))
      (nurbs-error :inp-err 'fit-crv-high-derivs))
    (when (or (>= kk p) (>= ll p))
      (nurbs-error :inp-err 'fit-crv-high-derivs))
    (when (or (< kk 0) (< ll 0))
      (nurbs-error :inp-err 'fit-crv-high-derivs))
    (let* ((pw (make-cpoint-array (1+ n)))
           (knt (make-real-vector (+ n p 2)))
           (cur (%make-curve p knt pw))
           (u (make-real-vector (1+ nn))))
      (declare (type cpoint-array pw) (type real-vector knt u))
      (fit-calc-crv-param-values pts nn +epoint+ par u)
      (fit-calc-knots-high-end-derivs u nn p kk ll knt)
      ;; coefficient matrix
      (let* ((bw (1- (* 2 p)))
             (sbw (1- p))
             (nb (make-real-vector (1+ p)))
             (nd (make-real-matrix (1+ (nl-max kk ll)) (1+ p)))
             (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
             (a (rmatrix-rm cm))
             (row 1))
        (declare (type fixnum bw sbw row) (type real-vector nb) (type real-matrix nd a))
        (%crv-fit-rm-set a 0 sbw 1d0)
        (basis-derivs knt p (aref u 0) +left+ kk nd)
        (do ((i 1 (1+ i))) ((> i kk))
          (declare (type fixnum i))
          (do ((j 0 (1+ j))) ((> j i))
            (declare (type fixnum j))
            (%crv-fit-rm-set a row (+ (- sbw i) j) (aref nd i j)))
          (incf row))
        (do ((i 1 (1+ i))) ((>= i nn))
          (declare (type fixnum i))
          (let* ((j (basis-eval knt p (aref u i) +left+ nb))
                 (l (- j row 1)))
            (declare (type fixnum j l))
            (dotimes (k (1+ p))
              (%crv-fit-rm-set a row (+ l k) (aref nb k)))
            (incf row)))
        (basis-derivs knt p (aref u nn) +left+ ll nd)
        (do ((i ll (1- i))) ((< i 1))
          (declare (type fixnum i))
          (do ((j 0 (1+ j))) ((> j i))
            (declare (type fixnum j))
            (%crv-fit-rm-set a row (+ sbw j) (aref nd i (+ (- p i) j))))
          (incf row))
        (%crv-fit-rm-set a row sbw 1d0)
        (real-matrix-lu-decompose cm)
        ;; solve
        (let ((q (make-point-array (1+ n))))
          (declare (type point-array q))
          (%copy-pt-row pts 0 q 0)
          (setf row 1)
          (do ((i 1 (1+ i))) ((> i kk))
            (declare (type fixnum i))
            (%copy-pt-row ds i q row)
            (incf row))
          (do ((i 1 (1+ i))) ((>= i nn))
            (declare (type fixnum i))
            (%copy-pt-row pts i q row)
            (incf row))
          (do ((i ll (1- i))) ((< i 1))
            (declare (type fixnum i))
            (%copy-pt-row de i q row)
            (incf row))
          (%copy-pt-row pts nn q row)
          (real-matrix-for-back cm q +epoint+ pw)))
      cur)))

;;; ---------------------------------------------------------- piecewise conics

(defconstant +crv-fit--conics-alfa+ 0.66d0 "ALFA in N_FitCrvConics.")

(defun %crv-fit-points (count)
  "A simple-vector of COUNT fresh points: an NL_POINT array whose elements the
caller hands to point routines directly."
  (declare (type fixnum count))
  (let ((v (make-array count)))
    (dotimes (i count v)
      (setf (svref v i) (make-point)))))

(defun fit-crv-conics (pp kk par tan csp)
  "The G1 piecewise conic (rational quadratic) curve interpolating the points
PP[0..KK] (a point-array), as a fresh curve: tangents by Bessel's (TAN
+bessel+, the parametrization PAR) or Akima's method (+akima+), corners kept
(CSP +cusp+) or smoothed (+nocusp+), each segment one conic arc, two arcs or
a line; interior points collinear with their neighbors are dropped.  The knots
come from the parametrization PAR of the segment end points.  Signals :inp-err
when KK < 2, :cal-err for another TAN or PAR, :geo-err for coincident
consecutive points or a zero tangent, :num-err when a division fails, and
whatever ISECT-LINE-LINE signals.
NLib: N_FitCrvConics."
  (declare (type point-array pp) (type fixnum kk par tan csp))
  (when (< kk 2)
    (nurbs-error :inp-err 'fit-crv-conics))
  (let* ((k kk)
         (pts (%crv-fit-points (1+ kk)))
         (cw (* 0.5d0 (sqrt 2d0)))
         (n (* 4 k))
         (pw (make-cpoint-array (1+ n)))
         (uk (make-real-vector (+ n 4)))
         (u (make-real-vector (1+ (* 2 k))))
         (q (%crv-fit-points (1+ (* 2 k))))
         (s (%crv-fit-points (+ k 4)))
         (su (%crv-fit-points (+ k 4)))
         (tv (%crv-fit-points (1+ k)))
         (closed nil)
         (l 0) (m 0)
         (alf 0d0) (oma 0d0) (a 0d0) (b 0d0) (w 0d0) (gam 0d0) (gam1 0d0)
         (r (make-point)) (mp (make-point)) (ri (make-point)) (ri1 (make-point)) (av (make-point))
         (c (make-cpoint))
         (li1 (%make-lineseg (make-point) (make-point) 0))
         (li (%make-lineseg (make-point) (make-point) 0)))
    (declare (type fixnum k n l m) (type simple-vector pts q s su tv) (type double-float cw)
             (type cpoint-array pw) (type real-vector uk u)
             (type double-float alf oma a b w gam gam1) (type point r mp ri ri1 av) (type cpoint c))
    (macrolet ((pt (i) `(the point (svref pts ,i)))
               (sv (i) `(the point (svref s ,i)))
               (suv (i) `(the point (svref su ,i)))
               (tg (i) `(the point (svref tv ,i)))
               (qv (i) `(the point (svref q ,i))))
      ;; copy of the point array, since it may be overwritten
      (dotimes (i (1+ k))
        (%row-pt pp i (pt i)))
      ;; chords
      (setf closed (< (dist-pt-pt (pt 0) (pt k)) *mtol*))
      (do ((i 1 (1+ i))) ((> i k))
        (declare (type fixnum i))
        (vector-diff (pt i) (pt (1- i)) (sv (1+ i))))
      (if closed
          (progn
            (vector-copy (sv (1+ k)) (sv 1))
            (vector-copy (sv k) (sv 0))
            (vector-copy (sv 2) (sv (+ k 2)))
            (vector-copy (sv 3) (sv (+ k 3))))
          (progn
            (vector-combine 2d0 (sv 2) -1d0 (sv 3) (sv 1))
            (vector-combine 2d0 (sv 1) -1d0 (sv 2) (sv 0))
            (vector-combine 2d0 (sv (1+ k)) -1d0 (sv k) (sv (+ k 2)))
            (vector-combine 2d0 (sv (+ k 2)) -1d0 (sv (1+ k)) (sv (+ k 3)))))
      (do ((i 2 (1+ i))) ((> i (1+ k)))
        (declare (type fixnum i))
        (vector-normalize (sv i) (suv i)))
      ;; unit tangents
      (cond
        ((= tan +bessel+)
         (let ((pa (make-point-array (1+ k))) (dui 0d0) (dui1 0d0))
           (declare (type point-array pa) (type double-float dui dui1))
           (dotimes (i (1+ k))
             (%pt-row (pt i) pa i))
           (fit-calc-crv-param-values pa k +epoint+ par u)
           (do ((i 1 (1+ i))) ((>= i k))
             (declare (type fixnum i))
             (setf dui (- (aref u i) (aref u (1- i)))
                   dui1 (- (aref u (1+ i)) (aref u i)))
             (when (float-op-is-bad dui (+ dui dui1) +division+)
               (nurbs-error :num-err 'fit-crv-conics))
             (setf alf (/ dui (+ dui dui1))
                   oma (- 1d0 alf))
             (when (float-op-is-bad oma dui +division+)
               (nurbs-error :num-err 'fit-crv-conics))
             (when (float-op-is-bad alf dui1 +division+)
               (nurbs-error :num-err 'fit-crv-conics))
             (setf oma (/ oma dui)
                   alf (/ alf dui1))
             (vector-combine oma (sv (1+ i)) alf (sv (+ i 2)) (tg i)))
           (setf dui (- (aref u 1) (aref u 0))
                 dui1 (- (aref u k) (aref u (1- k))))
           (if closed
               (progn
                 (when (float-op-is-bad 1d0 (* 2d0 dui) +division+)
                   (nurbs-error :num-err 'fit-crv-conics))
                 (when (float-op-is-bad 1d0 (* 2d0 dui1) +division+)
                   (nurbs-error :num-err 'fit-crv-conics))
                 (setf oma (/ 1d0 (* 2d0 dui))
                       alf (/ 1d0 (* 2d0 dui1)))
                 (vector-combine oma (sv 2) alf (sv (1+ k)) (tg 0))
                 (vector-copy (tg 0) (tg k)))
               (progn
                 (when (float-op-is-bad 2d0 dui +division+)
                   (nurbs-error :num-err 'fit-crv-conics))
                 (when (float-op-is-bad 2d0 dui1 +division+)
                   (nurbs-error :num-err 'fit-crv-conics))
                 (setf oma (/ 2d0 dui)
                       alf (/ 2d0 dui1))
                 (vector-combine oma (sv 2) -1d0 (tg 1) (tg 0))
                 (vector-combine alf (sv (1+ k)) -1d0 (tg (1- k)) (tg k))))))
        ((= tan +akima+)
         (do ((i 0 (1+ i))) ((> i k))
           (declare (type fixnum i))
           (vector-cross (sv i) (sv (1+ i)) av)
           (setf a (vector-magnitude av))
           (vector-cross (sv (+ i 2)) (sv (+ i 3)) av)
           (setf b (vector-magnitude av))
           (setf alf (if (< (+ a b) *ltol*)
                         (if (= csp +cusp+) 1d0 0.5d0)
                         (/ a (+ a b))))
           (setf oma (- 1d0 alf))
           (vector-combine oma (sv (1+ i)) alf (sv (+ i 2)) (tg i))))
        (t (nurbs-error :cal-err 'fit-crv-conics)))
      (do ((i 0 (1+ i))) ((> i k))
        (declare (type fixnum i))
        (vector-normalize-ref (tg i)))
      ;; interpolate segment by segment
      (setf n 0 l 0)
      (%c-pt-row (weight (pt 0) 1d0 c) pw 0)
      (copy-pt (pt 0) (qv 0))
      (let ((i 1))
        (declare (type fixnum i))
        (do () ((> i k))
          (let ((type 0))
            (declare (type fixnum type))
            ;; establish the case
            (create-line-start-dir-vector li1 (pt (1- i)) (tg (1- i)) +unbounded+)
            (create-line-start-dir-vector li (pt i) (tg i) +unbounded+)
            (multiple-value-bind (its ta tb) (isect-line-line li1 li r)
              (declare (type fixnum its) (type double-float ta tb))
              (setf a ta b tb)
              (if (= its +true+)
                  ;; the lines intersect
                  (setf type (cond ((and (> a 0d0) (< b 0d0)) 1)
                                   ((or (< (abs a) *ptol*) (< (abs b) *ptol*))
                                    (if (= csp +cusp+) 2 3))
                                   (t 3)))
                  ;; they do not
                  (progn
                    (setf a (vector-dot (tg (1- i)) (tg i)))
                    (vector-cross (tg (1- i)) (suv (1+ i)) av)
                    (setf b (vector-magnitude av))
                    (if (> b *ltol*)
                        (setf type (if (> a 0d0) 3 4))
                        ;; collinear points
                        (if (= i k)
                            (setf type 6)
                            (progn
                              (vector-cross (tg i) (tg (1+ i)) av)
                              (setf b (vector-magnitude av))
                              (if (> b *ltol*)
                                  (setf type 6)
                                  (progn
                                    (vector-cross (suv (1+ i)) (suv (+ i 2)) av)
                                    (setf b (vector-magnitude av))
                                    (setf type (if (< b *ltol*) 5 6))))))))))
            ;; the segment
            (case type
              (1                        ; single arc
               (setf w (bez-set-conic-weight (pt (1- i)) r (pt i)))
               (%c-pt-row (weight r w c) pw (incf n))
               (%c-pt-row (weight (pt i) 1d0 c) pw (incf n))
               (copy-pt (pt i) (qv (incf l))))
              (2                        ; corner
               (combine-2-pts 0.5d0 (pt (1- i)) 0.5d0 (pt i) r)
               (%c-pt-row (weight r 1d0 c) pw (incf n))
               (%c-pt-row (weight (pt i) 1d0 c) pw (incf n))
               (copy-pt (pt i) (qv (incf l))))
              (3                        ; split segment
               (let* ((d (dist-pt-pt (pt (1- i)) (pt i)))
                      (cosi1 (vector-dot (tg (1- i)) (suv (1+ i))))
                      (cosi (vector-dot (tg i) (suv (1+ i))))
                      (fa +crv-fit--conics-alfa+))
                 (declare (type double-float d cosi1 cosi fa))
                 (setf gam1 (/ (* 0.25d0 d) (+ (+ 1d0 (* fa (abs cosi))) (* (- 1d0 fa) (abs cosi1))))
                       gam (/ (* 0.25d0 d) (+ (+ 1d0 (* fa (abs cosi1))) (* (- 1d0 fa) (abs cosi))))
                       alf (/ gam1 (+ gam gam1))
                       oma (- 1d0 alf)))
               (vector-pt-along-vector (pt (1- i)) gam1 (tg (1- i)) ri1)
               (vector-pt-along-vector (pt i) (- gam) (tg i) ri)
               (combine-2-pts oma ri1 alf ri mp)
               (setf w (bez-set-conic-weight (pt (1- i)) ri1 mp))
               (%c-pt-row (weight ri1 w c) pw (incf n))
               (%c-pt-row (weight mp 1d0 c) pw (incf n))
               (copy-pt mp (qv (incf l)))
               (setf w (bez-set-conic-weight mp ri (pt i)))
               (%c-pt-row (weight ri w c) pw (incf n))
               (%c-pt-row (weight (pt i) 1d0 c) pw (incf n))
               (copy-pt (pt i) (qv (incf l))))
              (4                        ; antiparallel tangents
               (combine-2-pts 0.5d0 (pt (1- i)) 0.5d0 (pt i) mp)
               (setf gam (dist-pt-pt (pt (1- i)) mp))
               (vector-pt-along-vector (pt (1- i)) gam (tg (1- i)) ri1)
               (vector-pt-along-vector mp gam (tg (1- i)) mp)
               (vector-pt-along-vector (pt i) gam (tg (1- i)) ri)
               (%c-pt-row (weight ri1 cw c) pw (incf n))
               (%c-pt-row (weight mp 1d0 c) pw (incf n))
               (copy-pt mp (qv (incf l)))
               (%c-pt-row (weight ri cw c) pw (incf n))
               (%c-pt-row (weight (pt i) 1d0 c) pw (incf n))
               (copy-pt (pt i) (qv (incf l))))
              (5                        ; discard the point: shift down the data
               (do ((j i (1+ j))) ((>= j k))
                 (declare (type fixnum j))
                 (copy-pt (pt (1+ j)) (pt j))
                 (vector-copy (tg (1+ j)) (tg j)))
               (do ((j i (1+ j))) ((> j (1+ k)))
                 (declare (type fixnum j))
                 (vector-copy (suv (+ j 2)) (suv (1+ j))))
               (decf k)
               (decf i))
              (6                        ; collinear data: a line segment
               (combine-2-pts 0.5d0 (pt (1- i)) 0.5d0 (pt i) r)
               (%c-pt-row (weight r 1d0 c) pw (incf n))
               (%c-pt-row (weight (pt i) 1d0 c) pw (incf n))
               (copy-pt (pt i) (qv (incf l))))
              (t (nurbs-error :num-err 'fit-crv-conics))))
          (incf i)))
      ;; knot vector
      (setf m 2)
      (dotimes (i 3)
        (setf (aref uk i) 0d0))
      (let ((qa (make-point-array (1+ l))))
        (declare (type point-array qa))
        (dotimes (i (1+ l))
          (%pt-row (qv i) qa i))
        (fit-calc-crv-param-values qa l +epoint+ par u))
      (do ((i 1 (1+ i))) ((>= i l))
        (declare (type fixnum i))
        (setf (aref uk (incf m)) (aref u i))
        (setf (aref uk (incf m)) (aref u i)))
      (dotimes (i 3)
        (setf (aref uk (incf m)) 1d0)))
    ;; the curve as N_CrvSetSizeIndices and N_CrvCompress leave it
    (let ((cpts (make-cpoint-array (1+ n))))
      (dotimes (i (1+ n))
        (%copy-c-pt-row pw i cpts i))
      (%make-curve 2 (subseq uk 0 (1+ m)) cpts))))

;;; --------------------------------------- weighted, constrained least squares

(defun fit-crv-weighted-lst-sq-knots (e wp r f wd ii s ptp u knt n p)
  "The weighted and constrained least-squares curve of degree P with control
points 0..N on the knot vector KNT, as a fresh curve.  E[0..R] are the data
points at the parameters U[0..R], with weights WP[0..R] (WP[i] > 0
unconstrained, else interpolated exactly); F[0..S] are derivatives at the
points E[II[j]] (II ascending), with weights WD[0..S] likewise.  E and F are
point-arrays when PTP is +epoint+ (a non-rational curve), cpoint-arrays when
it is +hpoint+.  Constraints are met through Lagrange multipliers.  Signals
:inp-err when there are N or more constraints, too few unconstrained data
(mc+n >= mu+1), N < 1, P > N, KNT does not have N+P+2 knots or a span lacks
unconstrained data, :cal-err for another PTP, and whatever BASIS-DERIVS and
the matrix routines signal.
NLib: N_FitCrvWeightedLstSqKnots."
  (declare (type (simple-array double-float (* *)) e f) (type real-vector wp wd u knt)
           (type index-vector ii) (type fixnum r s ptp n p))
  (let ((rj p) (sj p) (ej -1)
        (ru -1) (su -1) (rc -1) (sc -1))
    (declare (type fixnum rj sj ej ru su rc sc))
    ;; constrained and unconstrained counts
    (do ((i 0 (1+ i))) ((> i r))
      (declare (type fixnum i))
      (if (> (aref wp i) 0d0) (incf ru) (incf rc)))
    (do ((j 0 (1+ j))) ((> j s))
      (declare (type fixnum j))
      (if (> (aref wd j) 0d0) (incf su) (incf sc)))
    (let* ((mknt (1- (length knt)))
           (mu (+ ru su 1))
           (mc (+ rc sc 1))
           (mcur (+ n p 1)))
      (declare (type fixnum mknt mu mc mcur))
      (when (or (>= mc n) (>= (+ mc n) (1+ mu)))
        (nurbs-error :inp-err 'fit-crv-weighted-lst-sq-knots))
      (when (or (< n 1) (> p n))
        (nurbs-error :inp-err 'fit-crv-weighted-lst-sq-knots))
      (when (/= mknt mcur)
        (nurbs-error :inp-err 'fit-crv-weighted-lst-sq-knots))
      (let* ((pw (make-cpoint-array (1+ n)))
             (cur (%make-curve p (copy-seq knt) pw))
             (index (make-index-vector (1+ mu)))
             (start (make-index-vector (1+ n)))
             (end (make-index-vector (1+ n)))
             (w (make-real-vector (1+ mu)))
             (nd (make-real-matrix 2 (1+ p)))
             (nn (make-real-matrix (1+ mu) (1+ p)))
             (wn (make-real-matrix (1+ mu) (1+ p)))
             (ntwn (set-real-matrix (make-rmatrix) n n +mt-full+ n))
             (mm (if (>= mc 0) (set-real-matrix (make-rmatrix) mc n +mt-full+ n) (make-rmatrix)))
             (am (rmatrix-rm ntwn))
             (cmat (rmatrix-rm mm))
             (mu2 0) (mc2 0) (jd 0))
        (declare (type cpoint-array pw) (type index-vector index start end) (type real-vector w)
                 (type real-matrix nd nn wn am) (type fixnum mu2 mc2 jd))
        ;; N, M and W
        (do ((i 0 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (setf (aref start i) 0
                (aref end i) mu))
        (do ((i 0 (1+ i))) ((> i r))
          (declare (type fixnum i))
          (let ((j (basis-derivs knt p (aref u i) +left+ 1 nd))
                (der (and (<= jd s) (= i (aref ii jd)))))
            (declare (type fixnum j))
            (if (> (aref wp i) 0d0)
                ;; unconstrained point data
                (progn
                  (setf (aref w mu2) (aref wp i))
                  (dotimes (k (1+ p))
                    (setf (aref nn mu2 k) (aref nd 0 k)))
                  (setf (aref index mu2) (- j p))
                  (when (> j rj)
                    (do ((k 1 (1+ k))) ((> k (- j rj)))
                      (declare (type fixnum k))
                      (incf sj)
                      (incf ej)
                      (setf (aref start sj) mu2
                            (aref end ej) (1- mu2)))
                    (setf rj j))
                  (incf mu2))
                ;; constrained point data
                (let ((cmat cmat))
                  (declare (type real-matrix cmat))
                  (dotimes (k (1+ p))
                    (setf (aref cmat mc2 (+ (- j p) k)) (aref nd 0 k)))
                  (incf mc2)))
            (when der
              (if (> (aref wd jd) 0d0)
                  ;; unconstrained derivative data
                  (progn
                    (setf (aref w mu2) (aref wd jd))
                    (dotimes (k (1+ p))
                      (setf (aref nn mu2 k) (aref nd 1 k)))
                    (setf (aref index mu2) (- j p))
                    (when (> j rj)
                      (do ((k 1 (1+ k))) ((> k (- j rj)))
                        (declare (type fixnum k))
                        (incf sj)
                        (incf ej)
                        (setf (aref start sj) mu2
                              (aref end ej) (1- mu2)))
                      (setf rj j))
                    (incf mu2))
                  ;; constrained derivative data
                  (let ((cmat cmat))
                    (declare (type real-matrix cmat))
                    (dotimes (k (1+ p))
                      (setf (aref cmat mc2 (+ (- j p) k)) (aref nd 1 k)))
                    (incf mc2)))
              (incf jd))))
        (when (or (< sj n) (= (aref end 0) -1))
          (nurbs-error :inp-err 'fit-crv-weighted-lst-sq-knots))
        ;; WN, NTWN and the LU decomposition of NTWN
        (do ((j 0 (1+ j))) ((> j n))
          (declare (type fixnum j))
          (do ((k (aref start j) (1+ k))) ((> k (aref end j)))
            (declare (type fixnum k))
            (%crv-fit-rm-set wn k (- j (aref index k))
                             (* (aref w k) (%crv-fit-rm-ref nn k (- j (aref index k)))))))
        (do ((i 0 (1+ i))) ((> i n))
          (declare (type fixnum i))
          (do ((j 0 (1+ j))) ((> j n))
            (declare (type fixnum j))
            (let ((lk (nl-max (aref start i) (aref start j)))
                  (hk (nl-min (aref end i) (aref end j))))
              (declare (type fixnum lk hk))
              (do ((k lk (1+ k))) ((> k hk))
                (declare (type fixnum k))
                (setf (aref am i j)
                      (+ (aref am i j)
                         (* (%crv-fit-rm-ref nn k (- i (aref index k)))
                            (%crv-fit-rm-ref wn k (- j (aref index k))))))))))
        (real-matrix-lu-decompose ntwn)
        ;; the right-hand sides, and the unconstrained solution
        (let ((epoint (= ptp +epoint+))
              (sa nil) (ntws nil) (intws nil) (v nil) (aa nil) (ta nil))
          (macrolet ((build (make-array* scale copy-row blend-row zero)
                       `(progn
                          (setf sa (,make-array* (1+ mu))
                                ntws (,make-array* (1+ n))
                                intws (,make-array* (1+ n))
                                v (,make-array* (1+ (nl-max n mc))))
                          (when (>= mc 0)
                            (setf aa (,make-array* (1+ mc))
                                  ta (,make-array* (1+ mc))))
                          (setf mu2 0 mc2 0 jd 0)
                          (do ((i 0 (1+ i))) ((> i r))
                            (declare (type fixnum i))
                            (let ((der (and (<= jd s) (= i (aref ii jd)))))
                              (if (> (aref wp i) 0d0)
                                  (progn (,scale (aref w mu2) e i sa mu2) (incf mu2))
                                  (progn (,copy-row e i ta mc2) (incf mc2)))
                              (when der
                                (if (> (aref wd jd) 0d0)
                                    (progn (,scale (aref w mu2) f jd sa mu2) (incf mu2))
                                    (progn (,copy-row f jd ta mc2) (incf mc2)))
                                (incf jd))))
                          (do ((i 0 (1+ i))) ((> i n))
                            (declare (type fixnum i))
                            (,zero ntws i)
                            (do ((k (aref start i) (1+ k))) ((> k (aref end i)))
                              (declare (type fixnum k))
                              (,blend-row (%crv-fit-rm-ref nn k (- i (aref index k))) sa k ntws i))))))
            (flet ((scale-pt-row (alpha src i dst j)
                     (declare (type double-float alpha) (type point-array src dst) (type fixnum i j))
                     (let ((x (make-point)))
                       (declare (dynamic-extent x))
                       (%pt-row (scale-pt alpha (%row-pt src i x) x) dst j)))
                   (scale-c-pt-row (alpha src i dst j)
                     (declare (type double-float alpha) (type cpoint-array src dst) (type fixnum i j))
                     (let ((x (make-cpoint)))
                       (declare (dynamic-extent x))
                       (%c-pt-row (scale-c-pt alpha (%row-c-pt src i x) x) dst j)))
                   (blend-pt-row (alpha src i dst j)
                     (declare (type double-float alpha) (type point-array src dst) (type fixnum i j))
                     (let ((x (make-point)) (y (make-point)))
                       (declare (dynamic-extent x y))
                       (%pt-row (vector-blend-pt alpha (%row-pt src i x) (%row-pt dst j y)) dst j)))
                   (zero-pt-row (dst j)
                     (declare (type point-array dst) (type fixnum j))
                     (setf (aref dst j 0) 0d0 (aref dst j 1) 0d0 (aref dst j 2) 0d0))
                   (zero-c-pt-row (dst j)
                     (declare (type cpoint-array dst) (type fixnum j))
                     (setf (aref dst j 0) 0d0 (aref dst j 1) 0d0 (aref dst j 2) 0d0 (aref dst j 3) 0d0)))
              (cond
                (epoint
                 (build make-point-array scale-pt-row %copy-pt-row blend-pt-row zero-pt-row)
                 (when (< mc 0)
                   (real-matrix-for-back ntwn ntws +epoint+ pw)
                   (return-from fit-crv-weighted-lst-sq-knots cur)))
                ((= ptp +hpoint+)
                 (build make-cpoint-array scale-c-pt-row %copy-c-pt-row %vector-blend-c-pt-row zero-c-pt-row)
                 (when (< mc 0)
                   (real-matrix-for-back ntwn ntws +hpoint+ pw)
                   (return-from fit-crv-weighted-lst-sq-knots cur)))
                (t (nurbs-error :cal-err 'fit-crv-weighted-lst-sq-knots)))))
          ;; matrices for the Lagrange multipliers
          (let* ((intwn (real-matrix-inverse ntwn +yes+ (make-rmatrix)))
                 (mt (real-matrix-transpose mm (make-rmatrix)))
                 (mit (real-matrix-multiply intwn mt (make-rmatrix)))
                 (mi (real-matrix-multiply mm mit (make-rmatrix)))
                 (imit (real-matrix-inverse mi +no+ (make-rmatrix))))
            ;; control points
            (if epoint
                (let ((x (make-point)) (y (make-point)) (c (make-cpoint)))
                  (declare (dynamic-extent x y c))
                  (real-matrix-multiply-pt-array intwn ntws intws)
                  (real-matrix-multiply-pt-array mm intws v)
                  (do ((i 0 (1+ i))) ((> i mc))
                    (declare (type fixnum i))
                    (%pt-row (diff-2-pts (%row-pt v i x) (%row-pt ta i y) x) v i))
                  (real-matrix-multiply-pt-array imit v aa)
                  (real-matrix-multiply-pt-array mit aa v)
                  (do ((i 0 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (diff-2-pts (%row-pt intws i x) (%row-pt v i y) x)
                    (%c-pt-row (pt-to-c-pt x c) pw i))
                  (when (< (aref wp 0) 0d0)
                    (%c-pt-row (pt-to-c-pt (%row-pt e 0 x) c) pw 0))
                  (when (< (aref wp r) 0d0)
                    (%c-pt-row (pt-to-c-pt (%row-pt e r x) c) pw n)))
                (let ((x (make-cpoint)) (y (make-cpoint)))
                  (declare (dynamic-extent x y))
                  (real-matrix-multiply-c-pt-array intwn ntws intws)
                  (real-matrix-multiply-c-pt-array mm intws v)
                  (do ((i 0 (1+ i))) ((> i mc))
                    (declare (type fixnum i))
                    (%c-pt-row (diff-2-c-pts (%row-c-pt v i x) (%row-c-pt ta i y) x) v i))
                  (real-matrix-multiply-c-pt-array imit v aa)
                  (real-matrix-multiply-c-pt-array mit aa v)
                  (do ((i 0 (1+ i))) ((> i n))
                    (declare (type fixnum i))
                    (%c-pt-row (diff-2-c-pts (%row-c-pt intws i x) (%row-c-pt v i y) x) pw i))
                  (when (< (aref wp 0) 0d0)
                    (%copy-c-pt-row e 0 pw 0))
                  (when (< (aref wp r) 0d0)
                    (%copy-c-pt-row e r pw n)))))
          cur)))))

;;; ------------------------------------------------ L3: interpolation with end derivatives

(defun %crv-fit-estimate-end-deriv (pts kk whr d)
  "Estimate the end derivative of the data PTS[0..KK] at the start (WHR
+start+) or the end (+end+) as N_FitCrvDerivs and N_FitCubicSplineInterp do
when it is the zero vector: the derivative of the cubic Bezier interpolating
the first (or last) four points, normalized and scaled by the length of the
data polygon, written into D.  The C ignores the errors of N_CrvDerivs and
N_VectorNormalizeRef (D is then left as it was), but not those of
N_FitCrvInterp, which propagate here.  KK must be at least 3."
  (declare (type point-array pts) (type fixnum kk whr) (type point d))
  (let ((cd (make-point-array 2)) (s (make-point)))
    (declare (type point-array cd) (dynamic-extent s))
    (if (= whr +start+)
        (let ((cub (fit-crv-interp pts 3 3 +chordlength+)))
          (handler-case (crv-derivs cub 0d0 +left+ 1 cd) (nurbs-error () nil)))
        (let ((sub (make-point-array 4)))
          (declare (type point-array sub))
          (dotimes (i 4)
            (%copy-pt-row pts (+ (- kk 3) i) sub i))
          (let ((cub (fit-crv-interp sub 3 3 +chordlength+)))
            (handler-case (crv-derivs cub 1d0 +right+ 1 cd) (nurbs-error () nil)))))
    (vector-copy (%row-pt cd 1 s) d)
    (let ((len (dist-polygon pts kk)))
      (declare (type double-float len))
      (handler-case (vector-normalize-ref d) (nurbs-error () nil))
      (vector-scale d len d))
    d))

(declaim (inline %crv-fit-zero-vector-p))
(defun %crv-fit-zero-vector-p (v)
  "True when all three coordinates of V are 0.0."
  (declare (type point v))
  (and (= (aref v 0) 0d0) (= (aref v 1) 0d0) (= (aref v 2) 0d0)))

(defun fit-crv-derivs (pts kk p ds de par &optional knots)
  "The curve of degree P interpolating the points PTS[0..KK] (a point-array)
with the start and end derivatives DS and DE (points, not modified), as a
fresh non-rational curve of KK+3 control points.  A zero DS or DE is
estimated from a cubic Bezier through the first or last four points (KK must
then be at least 3).  PAR is +uniform+, +chordlength+ or +centripetal+, the
parameters must then lie in [0,1] as these give them; or +inherited+, which
takes the parameters from the knot vector of the curve passed in: KNOTS,
standing for the knots of the C's non-NULL output curve, of at least KK+P+4
knots, which the curve then keeps (with no KNOTS the C reads its freshly
allocated knot array; here that is zero-filled).  Signals :inp-err when P < 2
or P > KK+2, :sto-err when KNOTS is too short, :cal-err for another PAR, and
whatever BASIS-EVAL and the banded solver signal.
NLib: N_FitCrvDerivs."
  (declare (type point-array pts) (type fixnum kk p par) (type point ds de))
  (when (or (< p 2) (> p (+ kk 2)))
    (nurbs-error :inp-err 'fit-crv-derivs))
  ;; C passes the derivatives by value
  (let ((ds (copy-seq ds))
        (de (copy-seq de)))
    (declare (type point ds de))
    (when (%crv-fit-zero-vector-p ds)
      (%crv-fit-estimate-end-deriv pts kk +start+ ds))
    (when (%crv-fit-zero-vector-p de)
      (%crv-fit-estimate-end-deriv pts kk +end+ de))
    (let ((n (+ kk 2)))
      (declare (type fixnum n))
      (when (and knots (< (length (the real-vector knots)) (+ n p 2)))
        (nurbs-error :sto-err 'fit-crv-derivs))
      (let* ((pw (make-cpoint-array (1+ n)))
             (knt (if knots (subseq (the real-vector knots) 0 (+ n p 2)) (make-real-vector (+ n p 2))))
             (cur (%make-curve p knt pw))
             (u (make-real-vector (1+ kk))))
        (declare (type cpoint-array pw) (type real-vector knt u))
        ;; parameters and knot vector
        (if (= par +inherited+)
            (dotimes (i (1+ kk))
              (setf (aref u i) (aref knt (+ i p))))
            (progn
              (fit-calc-crv-param-values pts kk +epoint+ par u)
              (fit-calc-knot-vector-end-derivs u kk p knt)))
        ;; coefficient matrix
        (let* ((bw (1- (* 2 p)))
               (sbw (1- p))
               (nb (make-real-vector (1+ p)))
               (cm (set-real-matrix (make-rmatrix) n n +mt-banded+ bw))
               (a (rmatrix-rm cm)))
          (declare (type fixnum bw sbw) (type real-vector nb) (type real-matrix a))
          (%crv-fit-rm-set a 0 sbw 1d0)
          (%crv-fit-rm-set a 1 (1- sbw) -1d0)
          (%crv-fit-rm-set a 1 sbw 1d0)
          (%crv-fit-rm-set a (1- n) sbw -1d0)
          (%crv-fit-rm-set a (1- n) (1+ sbw) 1d0)
          (%crv-fit-rm-set a n sbw 1d0)
          (do ((i 2 (1+ i))) ((> i kk))
            (declare (type fixnum i))
            (let* ((j (basis-eval knt p (aref u (1- i)) +left+ nb))
                   (l (- j i 1)))
              (declare (type fixnum j l))
              (dotimes (k (1+ p))
                (%crv-fit-rm-set a i (+ l k) (aref nb k)))))
          (real-matrix-lu-decompose cm)
          ;; solve
          (let ((q (make-point-array (1+ n)))
                (s (make-point))
                (pf (float p 1d0)))
            (declare (type point-array q) (dynamic-extent s) (type double-float pf))
            (let ((fact (/ (aref knt (1+ p)) pf)))
              (%copy-pt-row pts 0 q 0)
              (%pt-row (scale-pt fact ds s) q 1))
            (do ((i 2 (1+ i))) ((> i kk))
              (declare (type fixnum i))
              (%copy-pt-row pts (1- i) q i))
            (let ((fact (/ (- 1d0 (aref knt n)) pf)))
              (%pt-row (scale-pt fact de s) q (1- n))
              (%copy-pt-row pts kk q n))
            (real-matrix-for-back cm q +epoint+ pw)))
        cur))))

(defun fit-cubic-spline-interp (pts k ds de par)
  "The cubic spline interpolating the points PTS[0..K] (a point-array) with
the start and end derivatives DS and DE (points, not modified), parametrized
by PAR, as a fresh non-rational curve of K+3 control points, solved as a
tridiagonal system.  A zero DS or DE is estimated from a cubic Bezier through
the first or last four points (K must then be at least 3).  Signals :inp-err
when K <= 0, :cal-err for another PAR, and whatever BASIS-EVAL signals.
NLib: N_FitCubicSplineInterp."
  (declare (type point-array pts) (type fixnum k par) (type point ds de))
  (when (<= k 0)
    (nurbs-error :inp-err 'fit-cubic-spline-interp))
  (let ((ds (copy-seq ds))
        (de (copy-seq de)))
    (declare (type point ds de))
    (when (%crv-fit-zero-vector-p ds)
      (%crv-fit-estimate-end-deriv pts k +start+ ds))
    (when (%crv-fit-zero-vector-p de)
      (%crv-fit-estimate-end-deriv pts k +end+ de))
    (let* ((n (+ k 2))
           (pw (make-cpoint-array (1+ n)))
           (uk (make-real-vector (+ n 5)))
           (cur (%make-curve 3 uk pw))
           (u (make-real-vector (1+ k)))
           (abc (make-real-vector 4))
           (a (make-point)) (s (make-point))
           (aw (make-cpoint)) (c1 (make-cpoint)) (c2 (make-cpoint))
           (alf 0d0) (bet 0d0) (gam 0d0) (den 0d0))
      (declare (type fixnum n) (type cpoint-array pw) (type real-vector uk u abc)
               (dynamic-extent a s aw c1 c2) (type double-float alf bet gam den))
      ;; parameters
      (fit-calc-crv-param-values pts k +epoint+ par u)
      ;; knot vector
      (dotimes (i 4)
        (setf (aref uk i) 0d0
              (aref uk (+ n i 1)) 1d0))
      (do ((i 1 (1+ i))) ((>= i k))
        (declare (type fixnum i))
        (setf (aref uk (+ i 3)) (aref u i)))
      ;; first two and last two control points
      (setf alf 1d0
            bet (/ (aref uk 4) 3d0))
      (%c-pt-row (pt-to-c-pt (%row-pt pts 0 s) aw) pw 0)
      (combine-2-pts alf (%row-pt pts 0 s) bet ds a)
      (%c-pt-row (pt-to-c-pt a aw) pw 1)
      (setf bet (/ (- (aref uk n) 1d0) 3d0))
      (%c-pt-row (pt-to-c-pt (%row-pt pts k s) aw) pw n)
      (combine-2-pts alf (%row-pt pts k s) bet de a)
      (%c-pt-row (pt-to-c-pt a aw) pw (1- n))
      (when (= k 1)
        (return-from fit-cubic-spline-interp cur))
      ;; k = 2
      (when (= k 2)
        (basis-eval uk 3 (aref uk 4) +left+ abc)
        (setf den (aref abc 1)
              alf (/ 1d0 den)
              bet (/ (- (aref abc 0)) den))
        (pt-to-c-pt (%row-pt pts 1 s) aw)
        (%c-pt-row (combine-2-c-pts alf aw bet (%row-c-pt pw 1 c1) c2) pw 2)
        (setf alf 1d0
              bet (/ (- (aref abc 2)) den))
        (%combine-2-c-pt-rows alf pw 2 bet pw 3 pw 2)
        (return-from fit-cubic-spline-interp cur))
      ;; tridiagonal system
      (let ((rw (make-cpoint-array (1+ k)))
            (dd (make-real-vector (1+ k))))
        (declare (type cpoint-array rw) (type real-vector dd))
        (do ((i 3 (1+ i))) ((>= i k))
          (declare (type fixnum i))
          (%c-pt-row (pt-to-c-pt (%row-pt pts (1- i) s) aw) rw i))
        (basis-eval uk 3 (aref uk 4) +left+ abc)
        (setf den (aref abc 1)
              alf (/ 1d0 den)
              bet (/ (- (aref abc 0)) den))
        (pt-to-c-pt (%row-pt pts 1 s) aw)
        (%c-pt-row (combine-2-c-pts alf aw bet (%row-c-pt pw 1 c1) c2) pw 2)
        (do ((i 3 (1+ i))) ((>= i k))
          (declare (type fixnum i))
          (setf (aref dd i) (/ (aref abc 2) den))
          (basis-eval uk 3 (aref uk (+ i 2)) +left+ abc)
          (setf den (- (aref abc 1) (* (aref abc 0) (aref dd i)))
                alf (/ 1d0 den)
                bet (/ (- (aref abc 0)) den))
          (%combine-2-c-pt-rows alf rw i bet pw (1- i) pw i))
        (setf (aref dd k) (/ (aref abc 2) den))
        (basis-eval uk 3 (aref uk n) +left+ abc)
        (setf den (- (aref abc 1) (* (aref abc 0) (aref dd k)))
              alf (/ 1d0 den)
              bet (/ (- (aref abc 2)) den)
              gam (/ (- (aref abc 0)) den))
        (pt-to-c-pt (%row-pt pts (1- k) s) aw)
        (%c-pt-row (combine-2-c-pts alf aw bet (%row-c-pt pw (1+ k) c1) c2) pw k)
        (%vector-blend-c-pt-row gam pw (1- k) pw k)
        (do ((i (1- k) (1- i))) ((< i 2))
          (declare (type fixnum i))
          (%vector-blend-c-pt-row (- (aref dd (1+ i))) pw (1+ i) pw i)))
      cur)))

;;; ------------------------------------------------- approximation to a tolerance

(defun fit-crv-approx (pts m cld ps pr e ktp)
  "The curve of degree PR approximating the points PTS[0..M] (a point-array)
within the error E, as a fresh non-rational curve: starting from the
interpolant of degree PS (chord length parameters), knots are removed, the
curve refitted by least squares one degree higher on single (KTP +single+)
or multiple knots (any other KTP), and the parameters corrected by point
projection, up to degree PR; when a step fails the process restarts from the
interpolant one degree higher.  With CLD +yes+ and PTS[0] = PTS[M] (within
*mtol*) the curve is closed with matching end derivatives (PS must then be at
least 2).  Signals :inp-err when PR > M, PR < PS, or the data is closed and
PS = 1, and whatever the interpolation and knot removal routines signal.
NLib: N_FitCrvApprox."
  (declare (type point-array pts) (type fixnum m cld ps pr ktp) (type double-float e))
  (let* ((d (pts-1d-get-max-extent pts m))
         (ltop (* d 1.0d-4))
         (ltoc 1.0d-6))
    (declare (type double-float d ltop ltoc))
    (when (> pr m)
      (nurbs-error :inp-err 'fit-crv-approx))
    (when (< pr ps)
      (nurbs-error :inp-err 'fit-crv-approx))
    (let ((uu (make-real-vector (+ m pr 4)))      ; the knot array of the C's cur
          (u (make-real-vector (1+ m)))
          (us (make-real-vector (1+ m)))
          (er (make-real-vector (1+ m)))
          (closed nil)
          (wp (make-real-vector (1+ m)))
          (wd (make-real-vector 2))
          (ii (make-index-vector 2))
          (dv (make-point-array 2))             ; D[0], D[1]
          (vs (make-point))
          (ve (make-point))
          (r (make-point))
          (pk (make-point))
          (cur nil))
      (declare (type real-vector uu u us er wp wd) (type index-vector ii) (type point-array dv)
               (type point vs ve r pk))
      ;; closed data set?
      (when (1d-pt-set-is-closed pts m +epoint+ *mtol*)
        (when (= cld +yes+)
          (setf closed t)))
      (when closed
        (when (= ps 1)
          (nurbs-error :inp-err 'fit-crv-approx))
        ;; derivative at the end points
        (let ((q (make-point-array (1+ m)))
              (nc (truncate m 2)))
          (declare (type point-array q) (type fixnum nc))
          (do ((k 0 (1+ k))) ((> k m))
            (declare (type fixnum k))
            (let ((l (mod (+ nc k) m)))
              (declare (type fixnum l))
              (when (= l 0)
                (setf l m))
              (%copy-pt-row pts l q k)
              (setf (aref wp k) 1d0
                    (aref er k) 0d0)))
          (fit-calc-crv-param-values q m +epoint+ +chordlength+ u)
          (let* ((kna (fit-crv-calc-knot-vector u m pr (make-real-vector (+ m pr 2))))
                 (cura (fit-crv-interp-given-params q m +epoint+ u kna pr)))
            (fit-remove-knots cura u er m e)
            (setf cura (fit-crv-approx-knots q m +epoint+ u (curve-knots cura) (curve-n cura) pr))
            (let ((i (truncate (1+ m) 2)))
              (declare (type fixnum i))
              (crv-derivs cura (aref u i) +left+ 1 dv))
            (%copy-pt-row dv 1 dv 0)
            (%row-pt dv 0 vs)
            (%row-pt dv 1 ve)))
        ;; prepare for constrained fitting
        (setf (aref wp 0) -1d0
              (aref wp m) -1d0
              (aref wd 0) -1d0
              (aref wd 1) -1d0
              (aref ii 0) 0
              (aref ii 1) m))
      ;; parameters
      (fit-calc-crv-param-values pts m +epoint+ +chordlength+ us)
      ;; for i = ps..pr: interpolate with degree i; for j = i..pr: remove knots,
      ;; least-squares fit one degree higher, adjust parameters and errors
      (do ((i ps (1+ i))) ((> i pr))
        (declare (type fixnum i))
        (dotimes (k (1+ m))
          (setf (aref u k) (aref us k)
                (aref er k) 0d0))
        (if closed
            (progn
              (fit-calc-knot-vector-end-derivs u m i uu)
              (setf cur (fit-crv-knots-and-derivs pts m +epoint+ u (subseq uu 0 (+ m i 4)) i vs ve)))
            (progn
              (fit-crv-calc-knot-vector u m i uu)
              (setf cur (fit-crv-interp-given-params pts m +epoint+ u (subseq uu 0 (+ m i 2)) i))))
        (let ((j i))
          (declare (type fixnum j))
          (do () ((> j pr))
            ;; remove as many knots as possible
            (fit-remove-knots cur u er m e)
            (replace uu (the real-vector (curve-knots cur)))
            (when (= j pr)
              (return))
            (let ((nc (curve-n cur))
                  (mc (curve-m cur))
                  (nh 0) (mh 0))
              (declare (type fixnum nc mc nh mh))
              ;; new knot vector, in place in the curve's knot array
              (if (= ktp +single+)
                  (progn
                    (setf nh (1+ nc)
                          mh (+ mc 2))
                    (when (> nh m)
                      (return))
                    (dotimes (k (+ j 2))
                      (setf (aref uu mh) (aref uu mc))
                      (decf mh))
                    (do ((k nc (1- k))) ((< k (1+ j)))
                      (declare (type fixnum k))
                      (setf (aref uu mh) (aref uu k))
                      (decf mh))
                    (dotimes (k (+ j 2))
                      (setf (aref uu mh) (aref uu 0))
                      (decf mh)))
                  (let ((nsp (basis-get-span-count (curve-knots cur) j)))
                    (declare (type fixnum nsp))
                    (setf nh (+ nc nsp)
                          mh (+ mc nsp 1))
                    (when (> nh m)
                      (return))
                    (do () ((not (> mc 0)))
                      (let ((k mc))
                        (declare (type fixnum k))
                        (do () ((not (and (> mc 0) (= (aref uu mc) (aref uu (1- mc))))))
                          (decf mc))
                        (let ((mlt (1+ (- k mc))))
                          (declare (type fixnum mlt))
                          (dotimes (l (1+ mlt))
                            (setf (aref uu mh) (aref uu mc))
                            (decf mh)))
                        (decf mc)))))
              ;; approximate by least squares
              (let ((knt (subseq uu 0 (+ nh j 3)))
                    (q (1+ j)))
                (declare (type real-vector knt) (type fixnum q))
                (handler-case
                    (setf cur (if closed
                                  (fit-crv-weighted-lst-sq-knots pts wp m dv wd ii 1 +epoint+ u knt nh q)
                                  (fit-crv-approx-knots pts m +epoint+ u knt nh q)))
                  (nurbs-error () (return))))
              ;; adjust parameters and error vector
              (let ((apr t))
                (do ((k 1 (1+ k))) ((>= k m))
                  (declare (type fixnum k))
                  (%row-pt pts k pk)
                  (handler-case
                      (setf (aref u k) (crv-closest-pt cur pk (aref u k) ltop ltoc r))
                    (nurbs-error ()
                      (setf apr nil)
                      (return)))
                  (setf (aref er k) (dist-pt-pt pk r))
                  (when (> (aref er k) e)
                    (setf apr nil)
                    (return)))
                (unless apr
                  (return))))
            (incf j))
          (when (= j pr)
            (return))))
      ;; compact curve
      (crv-compress cur))))

;;; ---------------------------------------- interpolation with knots and tangents

(defun %crv-fit-tan-knot-vector (ul k p tan knv)
  "Write into KNV the knot vector N_FitCrvKnotsAndTangents uses for the
tangency case TAN (+no+, +start+, +end+ or +both+) at the parameters
UL[0..K] for degree P.  The C ignores the error of N_FitCalcKnotVectorDeriv,
which then leaves KNV as it was.  Returns KNV."
  (declare (type real-vector ul knv) (type fixnum k p tan))
  (cond ((= tan +no+) (fit-crv-calc-knot-vector ul k p knv))
        ((= tan +start+)
         (handler-case (fit-calc-knot-vector-deriv ul k p +start+ knv) (nurbs-error () nil)))
        ((= tan +end+)
         (handler-case (fit-calc-knot-vector-deriv ul k p +end+ knv) (nurbs-error () nil)))
        ((= tan +both+) (fit-calc-knot-vector-end-derivs ul k p knv)))
  knv)

(defun fit-crv-knots-and-tangents (pts k u ts te der knt per p)
  "The curve of degree P interpolating the points PTS[0..K] (a point-array)
at the parameters U[0..K] (or, when U is nil, chord length parameters scaled
to KNT's domain), with the start tangent TS and/or end tangent TE when given
(points; nil for none), as a fresh non-rational curve.  With DER +tangent+
TS and TE are directions, normalized and scaled IN PLACE by the polygon
length times the parameter range; with +derivative+ they are used as given.
KNT, a knot vector or nil, offers the knots: for each interior knot the
averaged one is replaced by the KNT knot closest to it within a window set
by PER in [0,1], and averaged knots for which KNT has none are added to it.
Returns (values curve knot-vector): the knot vector is KNT itself when
nothing was added, a fresh merged one when knots were added, and the curve's
own knots (a fresh copy) when KNT is nil.  Signals :inp-err for a KNT with
no distinct knots, P < 2, PER outside [0,1] or U not spanning KNT's domain,
:knt-err for an invalid KNT, :geo-err for a zero tangent, and whatever the
interpolation routines signal.
NLib: N_FitCrvKnotsAndTangents."
  (declare (type point-array pts) (type fixnum k der p) (type double-float per))
  (let ((mt 0) (pt 0) (n 0) (tan +no+) (mo 0) (ms -1)
        (ut (or knt (make-real-vector 1))))
    (declare (type fixnum mt pt n tan mo ms) (type real-vector ut))
    ;; incoming data
    (when knt
      (setf mt (1- (length ut)))
      (setf pt (do ((q 0 (1+ q))) ((or (>= q mt) (/= (aref ut q) (aref ut (1+ q)))) q)
                 (declare (type fixnum q))))
      (when (>= pt mt)
        (nurbs-error :inp-err 'fit-crv-knots-and-tangents))
      (knot-vector-is-valid ut pt 'fit-crv-knots-and-tangents))
    (when (< p 2)
      (nurbs-error :inp-err 'fit-crv-knots-and-tangents))
    (when (or (< per 0d0) (> per 1d0))
      (nurbs-error :inp-err 'fit-crv-knots-and-tangents))
    (when (and u knt)
      (when (or (/= (aref (the real-vector u) 0) (aref ut 0))
                (/= (aref (the real-vector u) k) (aref ut mt)))
        (nurbs-error :inp-err 'fit-crv-knots-and-tangents)))
    ;; tangency conditions
    (cond ((and (null ts) (null te)) (setf tan +no+ n k))
          ((and ts (null te)) (setf tan +start+ n (1+ k)))
          ((and (null ts) te) (setf tan +end+ n (1+ k)))
          (t (setf tan +both+ n (+ k 2))))
    ;; parameters, if not passed in
    (let ((ul (or u
                  (let ((ul (fit-calc-crv-param-values pts k +epoint+ +chordlength+ (make-real-vector (1+ k)))))
                    (declare (type real-vector ul))
                    (when knt
                      (when (or (/= (aref ul 0) (aref ut 0)) (/= (aref ul k) (aref ut mt)))
                        (let ((fac (/ (- (aref ut mt) (aref ut 0)) (- (aref ul k) (aref ul 0)))))
                          (declare (type double-float fac))
                          (do ((i 1 (1+ i))) ((>= i k))
                            (declare (type fixnum i))
                            (setf (aref ul i) (+ (aref ut 0) (* fac (- (aref ul i) (aref ul 0))))))
                          (setf (aref ul 0) (aref ut 0)
                                (aref ul k) (aref ut mt)))))
                    ul)))
          (ui (make-real-vector (+ n p 2)))
          (us nil))
      (declare (type real-vector ul ui))
      ;; knot vector, if none passed in
      (unless knt
        (%crv-fit-tan-knot-vector ul k p tan ui)
        (setf mo (+ n p 1)))
      ;; tangent magnitudes
      (when (and (/= tan +no+) (= der +tangent+))
        (let ((fac (dist-polygon pts k)))
          (declare (type double-float fac))
          (setf fac (* fac (- (aref ul k) (aref ul 0))))
          (when ts
            (vector-normalize-ref ts)
            (vector-scale ts fac ts))
          (when te
            (vector-normalize-ref te)
            (vector-scale te fac te))))
      ;; select knots for interpolation
      (when knt
        (let ((uss (make-real-vector (+ n p 2)))
              (ua (make-real-vector (+ n p 2)))
              (up (make-real-vector (+ n p 1))))
          (declare (type real-vector uss ua up))
          (setf us uss)
          (dotimes (i (1+ p))
            (setf (aref ui i) (aref ul 0)
                  (aref ui (+ n i 1)) (aref ul k)))
          ;; ideal knot vectors for degrees p and p-1
          (%crv-fit-tan-knot-vector ul k p tan ua)
          (%crv-fit-tan-knot-vector ul k (1- p) tan up)
          ;; select knots from knt
          (let ((j (1+ pt)))
            (declare (type fixnum j))
            (setf ms -1)
            (do ((i (1+ p) (1+ i))) ((> i n))
              (declare (type fixnum i))
              (let* ((a (+ (* (- 1d0 per) (aref ua i)) (* per (aref up (1- i)))))
                     (b (+ (* (- 1d0 per) (aref ua i)) (* per (aref up i))))
                     (fac (- b a))
                     (jc -1))
                (declare (type double-float a b fac) (type fixnum jc))
                (do () ((not (<= (aref ut j) a)))
                  (incf j))
                (do () ((not (and (> (aref ut j) a) (< (aref ut j) b))))
                  (let ((d (abs (- (aref ua i) (aref ut j)))))
                    (declare (type double-float d))
                    (when (< d fac)
                      (setf fac d
                            jc j)))
                  (incf j))
                (if (= jc -1)
                    (setf (aref ui i) (aref ua i)
                          (aref uss (incf ms)) (aref ua i))
                    (setf (aref ui i) (aref ut jc))))))))
      ;; interpolate the data points
      (let ((cur (cond ((= tan +no+)
                        (fit-crv-interp-given-params pts k +epoint+ ul ui p))
                       ((= tan +start+)
                        (fit-crv-knots-and-deriv pts k +epoint+ ul ui p ts +start+))
                       ((= tan +end+)
                        (fit-crv-knots-and-deriv pts k +epoint+ ul ui p te +end+))
                       (t
                        (fit-crv-knots-and-derivs pts k +epoint+ ul ui p ts te)))))
        ;; the knot vector handed back
        (values cur
                (cond
                  ((null knt) (subseq ui 0 (1+ mo)))
                  ((> ms -1)
                   (let ((uo (make-real-vector (+ mt ms 2)))
                         (us us)
                         (o -1) (i 0) (j 0))
                     (declare (type real-vector uo us) (type fixnum o i j))
                     (do () ((not (and (<= i mt) (<= j ms))))
                       (if (< (aref ut i) (aref us j))
                           (progn (setf (aref uo (incf o)) (aref ut i)) (incf i))
                           (progn (setf (aref uo (incf o)) (aref us j)) (incf j))))
                     (do () ((not (<= i mt)))
                       (setf (aref uo (incf o)) (aref ut i))
                       (incf i))
                     uo))
                  (t knt)))))))

;;; ------------------------------------------------------ closed interpolation

(defun fit-crv-approx-closed-conditions (pts n p d k par)
  "The curve of degree P interpolating the closed data PTS[0..N] (a
point-array, PTS[0] = PTS[N]) with the derivatives of orders 1..K the same
at both ends, parametrized by PAR, as a fresh non-rational curve
(FIT-CRV-HIGH-DERIVS).  D, a point-array whose rows 1..K are the
derivatives, or nil: they are then taken off a local interpolant of degree P
through the data around the start point.  Signals whatever
FIT-CALC-CRV-PARAM-VALUES, FIT-CRV-INTERP-GIVEN-PARAMS, CRV-DERIVS and
FIT-CRV-HIGH-DERIVS signal.
NLib: N_FitCrvApproxClosedConditions."
  (declare (type point-array pts) (type fixnum n p k par))
  (let ((ds (make-point-array (1+ k)))
        (de (make-point-array (1+ k))))
    (declare (type point-array ds de))
    (if d
        ;; derivatives passed in
        (do ((i 1 (1+ i))) ((> i k))
          (declare (type fixnum i))
          (%copy-pt-row d i ds i)
          (%copy-pt-row d i de i))
        ;; none: compute them
        (let* ((nm (truncate n 2))
               (is (nl-max nm (- n p 1)))
               (ie (nl-min nm (1+ p)))
               (nl (+ (- n is) ie))
               (q (make-point-array (1+ nl)))
               (dl (make-point-array (1+ k)))
               (u (make-real-vector (1+ nl)))
               (v (make-real-vector (1+ n)))
               (knt (make-real-vector (+ nl p 2)))
               (id 0))
          (declare (type fixnum nm is ie nl id) (type point-array q dl) (type real-vector u v knt))
          ;; data points for local fitting
          (fit-calc-crv-param-values pts n +epoint+ par v)
          (%copy-pt-row pts is q 0)
          (setf (aref u 0) 0d0)
          (let ((i (1+ is))
                (j 1))
            (declare (type fixnum i j))
            (do () ((not (<= j nl)))
              (%copy-pt-row pts i q j)
              (setf (aref u j) (+ (aref u (1- j)) (- (aref v i) (aref v (1- i)))))
              (when (= i n)
                (setf i 0
                      id j))
              (incf i)
              (incf j)))
          ;; interpolate local data through the start point
          (fit-crv-calc-knot-vector u nl p knt)
          (let ((curi (fit-crv-interp-given-params q nl +epoint+ u knt p)))
            ;; derivatives off the local interpolant
            (crv-derivs curi (aref u id) +left+ k dl))
          (do ((i 1 (1+ i))) ((> i k))
            (declare (type fixnum i))
            (%copy-pt-row dl i ds i)
            (%copy-pt-row dl i de i))))
    ;; interpolate the data set
    (fit-crv-high-derivs pts n p ds k de k par)))

;;; ----------------------------------------------- interpolation with normals

(defun fit-crv-pts-normals (pts nrm kk p)
  "The curve of degree P (2 or 3) interpolating the points PTS[0..KK] with
tangents perpendicular to the normals NRM[0..KK] (point-arrays), as a fresh
non-rational curve: the derivatives of the chord length interpolant at the
data points (found by projection), projected into the planes normal to NRM,
are the tangent directions of FIT-CRV-FIRST-DERIV (centripetal parameters).
For another degree the C returns without error and without a curve: NIL.
Signals whatever FIT-CRV-INTERP, CRV-CLOSEST-PT, CRV-DERIVS and
FIT-CRV-FIRST-DERIV signal.
NLib: N_FitCrvPtsNormals."
  (declare (type point-array pts nrm) (type fixnum kk p))
  (when (or (< p 2) (> p 3))
    (return-from fit-crv-pts-normals nil))
  (let ((crv (fit-crv-interp pts kk p +chordlength+))
        (tdir (make-point-array (1+ kk)))
        (pntder (make-point-array 2))
        (q (make-point)) (pi* (make-point)) (ni (make-point)) (d1 (make-point)) (mv (make-point)) (tv (make-point))
        (ui 0d0))
    (declare (type point-array tdir pntder) (type point q pi* ni d1 mv tv) (type double-float ui))
    (multiple-value-bind (u0 u1) (crv-get-param-bounds crv)
      (declare (type double-float u0 u1))
      (do ((i 0 (1+ i))) ((> i kk))
        (declare (type fixnum i))
        (cond ((= i 0) (setf ui u0))
              ((= i kk) (setf ui u1))
              (t
               ;; parameter value for this point, from the previous one
               (setf u0 ui)
               (setf ui (crv-closest-pt crv (%row-pt pts i pi*) u0 1.0d-4 1.0d-4 q))))
        ;; point and first derivative
        (crv-derivs crv ui +left+ 1 pntder)
        ;; the derivative projected into the plane normal to nrm[i]
        (%row-pt nrm i ni)
        (vector-cross ni (%row-pt pntder 1 d1) mv)
        (vector-cross mv ni tv)
        (%pt-row tv tdir i)))
    (fit-crv-first-deriv pts tdir kk p +tangent+ +centripetal+)))

;;; ------------------------------------------------------------- L6-g0

(defun fit-arc-to-pts (p n tol etp arc cit)
  "The least-squares circle or circular arc to the 2D points P[0..N] (a
point-array; N at least 1, else :inp-err), or the least-squares line segment
when the points are collinear -- the maximum (ETP +maximum+) or average
(+average+; another ETP signals :cal-err) distance from the line below TOL --
or when the line's average error is below the circle's.  ARC +no+ gives the
full circle, +yes+ the arc spanning the points (leaving out the largest gap
in their angles); CIT is create-circ-arc's kind (+quadratic+, +quartic+ or
+quintic+).  As in line-fit-pts, the z of 2D points is set to +no-z+ in P.
Returns (values cur cut era erm): the fresh curve, its type +nline+ or
+ncircle+, and the average and maximum absolute errors.  Signals what
line-fit-pts, create-circ-from-2d-pts, project-pt-line-param,
vector-directed-angle and create-circ-arc signal.
NLib: N_FitArcToPts."
  (declare (type point-array p) (type fixnum n etp arc cit) (type double-float tol))
  ;; check input
  (when (< n 1)
    (nurbs-error :inp-err 'fit-arc-to-pts))
  (let ((cc (make-point)) (vv (make-point)) (dd (make-point))
        (a (make-point)) (b (make-point)) (pi* (make-point))
        (xx (make-point)) (yy (make-point))
        (pli +no+) (lfi +no+) (cfi +no+)
        (cea 0d0) (cem 0d0) (cx 0d0) (cy 0d0) (rr 0d0))
    (declare (dynamic-extent cc vv dd a b pi* xx yy)
             (type fixnum pli lfi cfi) (type double-float cea cem cx cy rr))
    ;; fit a line first and check whether the points are collinear
    (multiple-value-bind (lea lem) (line-fit-pts p n cc vv)
      (declare (type double-float lea lem))
      (cond ((= etp +maximum+) (when (< lem tol) (setf pli +yes+)))
            ((= etp +average+) (when (< lea tol) (setf pli +yes+)))
            (t (nurbs-error :cal-err 'fit-arc-to-pts)))
      ;; if the points are not collinear, attempt to fit a circle
      (if (= pli +no+)
          (multiple-value-bind (cx1 cy1 rr1 cea1 cem1 cfl) (create-circ-from-2d-pts p n tol etp)
            (declare (type double-float cx1 cy1 rr1 cea1 cem1) (type fixnum cfl))
            (setf cx cx1 cy cy1 rr rr1 cea cea1 cem cem1)
            (if (= cfl +true+)
                (if (< lea cea)
                    (setf lfi +yes+)
                    (setf cfi +yes+))
                (setf lfi +yes+)))
          (setf lfi +yes+))
      ;; the line: the points are collinear, or its average error is smaller
      (when (= lfi +yes+)
        (vector-sum cc vv dd)
        (let ((tmin +bigd+) (tmax -1.0d0) (imin 0) (imax 0))
          (declare (type double-float tmin tmax) (type fixnum imin imax))
          (dotimes (i (1+ n))
            (let ((tt (project-pt-line-param (%row-pt p i pi*) cc dd a)))
              (declare (type double-float tt))
              (when (< tt tmin)
                (setf tmin tt
                      imin i))
              (when (> tt tmax)
                (setf tmax tt
                      imax i))))
          (project-pt-line-param (%row-pt p imin pi*) cc dd a)
          (project-pt-line-param (%row-pt p imax pi*) cc dd b)
          (vector-diff b a vv)
          (return-from fit-arc-to-pts
            (values (crv-line-from-pt-and-vector a vv) +nline+ lea lem))))
      ;; the circle: its average error is smaller
      (when (= cfi +yes+)
        (vector-create 1.0d0 0.0d0 0.0d0 xx)
        (vector-create 0.0d0 1.0d0 0.0d0 yy)
        (vector-create cx cy 0.0d0 cc)
        (let ((as 0d0) (ae 0d0))
          (declare (type double-float as ae))
          (if (= arc +no+)
              (setf as 0.0d0
                    ae 360.0d0)
              (let ((alf (make-real-vector (+ n 2)))
                    (gap -1.0d0)
                    (ig 0))
                (declare (type real-vector alf) (type double-float gap) (type fixnum ig))
                ;; the sweep angles
                (dotimes (i (1+ n))
                  (vector-diff (%row-pt p i pi*) cc vv)
                  (setf (aref alf i) (vector-directed-angle xx vv)))
                ;; sort them
                (shell-sort-real alf n)
                (setf (aref alf (1+ n)) (+ 360.0d0 (aref alf 0)))
                ;; the largest gap
                (dotimes (i (1+ n))
                  (let ((dal (abs (- (aref alf (1+ i)) (aref alf i)))))
                    (declare (type double-float dal))
                    (when (> dal gap)
                      (setf gap dal
                            ig i))))
                (setf as (aref alf (1+ ig))
                      ae (aref alf ig))
                (when (> as 360.0d0) (decf as 360.0d0))
                (when (> ae 360.0d0) (decf ae 360.0d0))))
          (values (create-circ-arc cc xx yy rr as ae cit) +ncircle+ cea cem))))))
