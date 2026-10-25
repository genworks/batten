;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvAdv.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.

(in-package :batten)

;;; Control point rows: %copy-c-pt-row, %combine-2-c-pt-rows and
;;; %vector-blend-c-pt-row in math.lisp.  Bezier decomposition:
;;; crv-decompose-bez in crv-tool.lisp.

;;; ------------------------------------------------------------------ routines

(defun crv-derivs-at-knot (cur u flg der dw)
  "Write the derivatives of curve CUR at U in homogeneous space, orders 0
through DER, into DW, a cpoint-array of at least (1+ DER) rows: DW[k] is the
kth derivative.  FLG +left+ means u is in [u_j, u_j+1) (right derivatives),
+right+ u in (u_j, u_j+1] (left derivatives).  Return DW.
NLib: N_CrvDerivsAtKnot."
  (declare (type curve cur) (type double-float u) (type fixnum flg der) (type cpoint-array dw))
  (let ((knots (curve-knots cur))
        (pw (curve-cpoints cur))
        (p (curve-degree cur)))
    (declare (type real-vector knots) (type cpoint-array pw) (type fixnum p))
    (knot-vector-is-param-out-of-bounds knots u 'crv-derivs-at-knot)
    (let* ((bd (make-real-matrix (1+ der) (1+ p)))
           (spn (basis-derivs knots p u flg der bd)))
      (declare (type real-matrix bd) (type fixnum spn))
      (fast
        (dotimes (k (1+ der))
          (setf (aref dw k 0) 0d0 (aref dw k 1) 0d0 (aref dw k 2) 0d0 (aref dw k 3) 0d0)
          (dotimes (i (1+ p))
            (%vector-blend-c-pt-row (aref bd k i) pw (+ (- spn p) i) dw k))))
      dw)))

(defun crv-get-max-2nd-deriv (cur)
  "An upper bound on the magnitude of the second derivative of curve CUR,
from sampling each Bezier segment at its nodes.  NLib: N_CrvGetMax2ndDeriv."
  (declare (type curve cur))
  (let* ((p (curve-degree cur))
         (bez (crv-decompose-bez cur))
         (k (1- (length bez)))
         (cd (make-point-array 3))
         (c2 (make-point))
         (u (make-real-vector (1+ p)))
         (duu 0d0))
    (declare (type fixnum p k) (type simple-vector bez) (type point-array cd)
             (type point c2) (type real-vector u) (type double-float duu))
    (do ((i 0 (1+ i))) ((> i k))
      (declare (type fixnum i))
      (let* ((al (if (= i 0) 0 1))
             (b (svref bez i))
             (uk (curve-knots b))
             (r (1- (length uk)))
             (ui (/ (- (aref uk r) (aref uk 0)) (float p 1d0))))
        (declare (type fixnum al r) (type real-vector uk) (type double-float ui))
        (setf (aref u 0) (aref uk 0)
              (aref u p) (aref uk r))
        (do ((a 1 (1+ a))) ((> a (1- p)))
          (declare (type fixnum a))
          (setf (aref u a) (+ (aref uk 0) (* (float a 1d0) ui))))
        (do ((a al (1+ a))) ((> a p))
          (declare (type fixnum a))
          (crv-derivs b (aref u a) +left+ 2 cd)
          (setf (aref c2 0) (aref cd 2 0) (aref c2 1) (aref cd 2 1) (aref c2 2) (aref cd 2 2))
          (let ((muu (pt-magnitude c2)))
            (declare (type double-float muu))
            (when (> muu duu)
              (setf duu muu))))))
    duu))

(defun crv-non-rat-get-first-deriv (cur d)
  "The Dth derivative curve of NON-RATIONAL curve CUR, 0 <= D <= degree, as a
fresh curve of degree p-D.  Signals :inp-err for a D out of range or a rational
curve, :der-err for a zero-length knot span in the difference, :num-err when a
coefficient would overflow.  NLib: N_CrvNonRatGetFirstDeriv."
  (declare (type curve cur) (type fixnum d))
  (let* ((n (curve-n cur))
         (pw (curve-cpoints cur))
         (p (curve-degree cur))
         (m (curve-m cur))
         (u (curve-knots cur)))
    (declare (type fixnum n p m) (type cpoint-array pw) (type real-vector u))
    (when (or (< d 0) (> d p))
      (nurbs-error :inp-err 'crv-non-rat-get-first-deriv))
    (when (is-crv-rat cur)
      (nurbs-error :inp-err 'crv-non-rat-get-first-deriv))
    (let* ((nd (- n d))
           (pd (- p d))
           (md (+ nd pd 1))
           (dw (make-cpoint-array (1+ nd)))
           (v (make-real-vector (1+ md)))
           (tw (make-cpoint-array (1+ n))))
      (declare (type fixnum nd pd md) (type cpoint-array dw tw) (type real-vector v))
      (fast
        (dotimes (i (1+ n))
          (%copy-c-pt-row pw i tw i))
        ;; control points of the derivative curve
        (do ((k 1 (1+ k))) ((> k d))
          (declare (type fixnum k))
          (do ((i 0 (1+ i))) ((> i (- n k)))
            (declare (type fixnum i))
            (when (= (aref u (+ i p 1)) (aref u (+ i k)))
              (nurbs-error :der-err 'crv-non-rat-get-first-deriv))
            (when (float-op-is-bad (+ (- (float p 1d0) (float k 1d0)) 1d0)
                                   (- (aref u (+ i p 1)) (aref u (+ i k)))
                                   +division+)
              (nurbs-error :num-err 'crv-non-rat-get-first-deriv))
            (let ((alf (/ (+ (- (float p 1d0) (float k 1d0)) 1d0)
                          (- (aref u (+ i p 1)) (aref u (+ i k))))))
              (declare (type double-float alf))
              (%combine-2-c-pt-rows alf tw (1+ i) (- alf) tw i tw i))))
        (do ((i 0 (1+ i))) ((> i nd))
          (declare (type fixnum i))
          (%copy-c-pt-row tw i dw i))
        ;; the knot vector
        (let ((j -1))
          (declare (type fixnum j))
          (do ((i 0 (1+ i))) ((> i pd))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u 0)))
          (do ((i (1+ p) (1+ i))) ((> i n))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u i)))
          (do ((i 0 (1+ i))) ((> i pd))
            (declare (type fixnum i))
            (setf (aref v (incf j)) (aref u m)))))
      (%make-curve pd v dw))))

(defun crv-planar-offset-get-deriv (cur u d nv ufl od)
  "Write the first and second derivatives at U of the offset by distance D of
planar curve CUR into rows 1 and 2 of OD, a point-array of at least 3 rows
allocated by the caller (row 0 is left alone).  NV is the normal to the
curve's plane: the curve is offset in the direction NV x T.  UFL +left+ means
u is in [u_j, u_j+1) (right derivatives), +right+ u in (u_j, u_j+1] (left
derivatives).  Signals as crv-derivs and crv-get-curvature-deriv do, and
:num-err when the curvature's division would overflow.  Returns OD.
NLib: N_CrvPlanarOffsetGetDeriv."
  (declare (type curve cur) (type double-float u d) (type point nv) (type fixnum ufl)
           (type point-array od))
  (let ((dd (make-point-array 3))
        (d1 (make-point)) (d2 (make-point)) (b (make-point)) (o (make-point))
        (num 0d0) (den 0d0) (kap 0d0) (kpp 0d0) (f1 0d0) (f2 0d0) (sgn 0d0))
    (declare (type point-array dd) (type point d1 d2 b o)
             (type double-float num den kap kpp f1 f2 sgn)
             (dynamic-extent dd d1 d2 b o))
    ;; derivatives, the curvature and its derivative
    (crv-derivs cur u ufl 2 dd)
    (%point-array-row dd 1 d1)
    (%point-array-row dd 2 d2)
    (vector-cross d1 d2 b)
    (setf num (vector-magnitude b))
    (setf den (vector-magnitude d1))
    (setf den (* (* den den) den))
    (when (float-op-is-bad num den +division+)
      (nurbs-error :num-err 'crv-planar-offset-get-deriv))
    (setf kap (/ num den))
    (setf den (vector-dot b nv))
    (setf sgn (if (> den 0d0) -1d0 1d0))
    (setf kpp (crv-get-curvature-deriv cur u ufl))
    ;; offset derivatives
    (setf f1 (+ 1d0 (* (* sgn kap) d)))
    (vector-scale d1 f1 o)
    (setf (aref od 1 0) (aref o 0) (aref od 1 1) (aref o 1) (aref od 1 2) (aref o 2))
    (setf f2 (* (* sgn kpp) d))
    (vector-combine f1 d2 f2 d1 o)
    (setf (aref od 2 0) (aref o 0) (aref od 2 1) (aref o 1) (aref od 2 2) (aref o 2))
    od))

;;; ------------------------------------------------------------- L2-g7

(defun crv-offset-get-max-2nd-deriv (cur d nv)
  "An upper bound on the magnitude of the second derivative of the offset by
distance D of planar curve CUR (NV the normal to its plane), from sampling
each Bezier segment at 10p+1 equally spaced parameters.  Signals as
crv-decompose-bez and crv-planar-offset-get-deriv do.
NLib: N_CrvOffsetGetMax2ndDeriv."
  (declare (type curve cur) (type double-float d) (type point nv))
  (let ((p (curve-degree cur)))
    (declare (type fixnum p))
    (multiple-value-bind (bez k) (crv-decompose-bez cur)
      (declare (type simple-vector bez) (type fixnum k))
      (let* ((duu 0d0)
             (mu (* 10 p))
             (u (make-real-vector (1+ mu)))
             (od (make-point-array 3))
             (c2 (make-point)))
        (declare (type double-float duu) (type fixnum mu) (type real-vector u)
                 (type point-array od) (type point c2))
        (do ((i 0 (1+ i))) ((> i k))
          (declare (type fixnum i))
          (let* ((al (if (= i 0) 0 1))
                 (b (svref bez i))
                 (uk (curve-knots b))
                 (r (1- (length uk)))
                 (ui (/ (- (aref uk r) (aref uk 0)) (float mu 1d0))))
            (declare (type fixnum al r) (type real-vector uk) (type double-float ui))
            (setf (aref u 0) (aref uk 0)
                  (aref u mu) (aref uk r))
            (do ((a 1 (1+ a))) ((> a (1- mu)))
              (declare (type fixnum a))
              (setf (aref u a) (+ (aref uk 0) (* (float a 1d0) ui))))
            (do ((a al (1+ a))) ((> a mu))
              (declare (type fixnum a))
              (crv-planar-offset-get-deriv b (aref u a) d nv +left+ od)
              (setf (aref c2 0) (aref od 2 0) (aref c2 1) (aref od 2 1) (aref c2 2) (aref od 2 2))
              (let ((muu (pt-magnitude c2)))
                (declare (type double-float muu))
                (when (> muu duu)
                  (setf duu muu))))))
        duu))))

(defconstant +crv-adv--norem+ 1d25 "NOREM in NL_CrvAdv.cpp: a knot found not removable.")

(defun crv-remove-all-knots (curp tol)
  "A fresh copy of curve CURP with every knot removed whose removal keeps the
accumulated error bound within TOL (scaled by the minimum weight for a
rational curve), cheapest knot first; knots whose removal would push a weight
outside [*wmin*, *wmax*] stay.  The C writes into a caller's curve, or into
CURP itself; Batten returns the new curve.  NLib: N_CrvRemoveAllKnots."
  (declare (type curve curp) (type double-float tol))
  (let* ((n (curve-n curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (ns n)
         (rat nil)
         (lam 0d0) (oml 0d0))
    (declare (type fixnum n p m ns) (type double-float lam oml))
    ;; removal tolerance for rational curves
    (when (is-crv-rat curp)
      (multiple-value-bind (wmin wmax pmin pmax) (crv-get-min-max-weights-and-pts curp)
        (declare (type double-float wmin pmax) (ignore wmax pmin))
        (setf tol (/ (* tol wmin) (+ 1d0 pmax))
              rat t)))
    (let* ((curq (crv-copy curp))
           (qw (curve-cpoints curq))
           (uq (curve-knots curq))
           (rw (make-cpoint-array (1+ (* 2 p))))
           (br (make-real-vector (1+ m)))
           (sr (make-index-vector (1+ m)))
           (er (make-real-vector (1+ m)))
           (te (make-real-vector (1+ m))))
      (declare (type curve curq) (type cpoint-array qw rw) (type real-vector uq br er te)
               (type index-vector sr))
      (dotimes (i (1+ m))
        (setf (aref br i) +bigd+
              (aref sr i) 0
              (aref er i) 0d0))
      ;; knot removal errors for each distinct knot
      (let ((r (1+ p)))
        (declare (type fixnum r))
        (do () ((not (<= r n)))
          (let ((i r))
            (declare (type fixnum i))
            (do () ((not (and (<= r n) (= (aref uq r) (aref uq (1+ r)))))) (incf r))
            (setf (aref sr r) (1+ (- r i)))
            (setf (aref br r) (crv-remove-knot-max-err curp r (aref sr r)))
            (incf r))))
      ;; try to remove each knot
      (do () (nil)
        (let ((b (aref br (1+ p)))
              (s (aref sr (1+ p)))
              (r (1+ p))
              (k 0) (l 0)
              (rmf t))
          (declare (type double-float b) (type fixnum s r k l))
          ;; knot with the smallest error
          (do ((i (+ p 2) (1+ i))) ((> i (- m p 1)))
            (declare (type fixnum i))
            (when (< (aref br i) b)
              (setf b (aref br i)
                    s (aref sr i)
                    r i)))
          ;; no more removable knots: finished
          (when (or (= b +bigd+) (= b +crv-adv--norem+))
            (return))
          ;; error of removal
          (cond ((oddp (+ p s))
                 (setf k (truncate (+ p s 1) 2)
                       l (+ (- r k) p 1))
                 (let ((alf (/ (- (aref uq r) (aref uq (- r k)))
                               (- (aref uq (+ (- r k) p 1)) (aref uq (- r k)))))
                       (bet (/ (- (aref uq r) (aref uq (+ (- r k) 1)))
                               (- (aref uq (+ (- r k) p 2)) (aref uq (+ (- r k) 1))))))
                   (declare (type double-float alf bet))
                   (setf lam (/ alf (+ alf bet))
                         oml (- 1d0 lam))))
                (t
                 (setf k (truncate (+ p s) 2)
                       l (+ (- r k) p))))
          (do ((i (- r k) (1+ i))) ((> i l))
            (declare (type fixnum i))
            (when (/= (aref uq i) (aref uq (1+ i)))
              (setf (aref te i) (+ (aref er i) b))
              (when (> (aref te i) tol)
                (setf rmf nil)
                (return))))
          (cond
            (rmf
             ;; error test passed: update the error vector and remove the knot
             (do ((i (- r k) (1+ i))) ((> i l))
               (declare (type fixnum i))
               (when (/= (aref uq i) (aref uq (1+ i)))
                 (setf (aref er i) (aref te i))))
             (let* ((fout (truncate (- (* 2 r) s p) 2))
                    (first (- r p))
                    (last (- r s))
                    (off (1- first))
                    (i first)
                    (j last)
                    (ii 1)
                    (jj (- last off))
                    (skip nil))
               (declare (type fixnum fout first last off i j ii jj))
               (%copy-c-pt-row qw off rw 0)
               (%copy-c-pt-row qw (1+ last) rw (- (1+ last) off))
               ;; new control points for one removal step
               (fast
                 (do () ((not (> (- j i) 0)))
                   (let* ((alf (/ (- (aref uq (+ i p 1)) (aref uq i)) (- (aref uq r) (aref uq i))))
                          (oma (- 1d0 alf))
                          (bet (/ (- (aref uq (+ j p 1)) (aref uq j)) (- (aref uq (+ j p 1)) (aref uq r))))
                          (omb (- 1d0 bet)))
                     (declare (type double-float alf oma bet omb))
                     (%combine-2-c-pt-rows alf qw i oma rw (1- ii) rw ii)
                     (%combine-2-c-pt-rows bet qw j omb rw (1+ jj) rw jj)
                     (incf i) (decf j) (incf ii) (decf jj))))
               ;; disallowed weights
               (when rat
                 (let ((i first) (j last) (wmin +bigd+) (wmax +smad+))
                   (declare (type fixnum i j) (type double-float wmin wmax))
                   (do () ((not (> (- j i) 0)))
                     (let ((wi (aref rw (- i off) 3))
                           (wj (aref rw (- j off) 3)))
                       (declare (type double-float wi wj))
                       (when (< wi wmin) (setf wmin wi))
                       (when (< wj wmin) (setf wmin wj))
                       (when (> wi wmax) (setf wmax wi))
                       (when (> wj wmax) (setf wmax wj))
                       (incf i) (decf j)))
                   (when (or (< wmin (the double-float *wmin*)) (> wmax (the double-float *wmax*)))
                     (setf (aref br r) +crv-adv--norem+
                           skip t))))
               (unless skip
                 ;; save new control points
                 (when (oddp (+ p s))
                   (%combine-2-c-pt-rows lam rw (1+ jj) oml rw (1- ii) rw (1+ jj)))
                 (setf i first
                       j last)
                 (fast
                   (do () ((not (> (- j i) 0)))
                     (%copy-c-pt-row rw (- i off) qw i)
                     (%copy-c-pt-row rw (- j off) qw j)
                     (incf i) (decf j)))
                 ;; shift down some parameters
                 (when (= s 1)
                   (setf (aref er (1- r)) (nl-max (aref er (1- r)) (aref er r))))
                 (when (> s 1)
                   (setf (aref sr (1- r)) (1- (aref sr r))))
                 (do ((i (1+ r) (1+ i))) ((> i m))
                   (declare (type fixnum i))
                   (setf (aref br (1- i)) (aref br i)
                         (aref sr (1- i)) (aref sr i)
                         (aref er (1- i)) (aref er i)))
                 ;; shift down knots and control points
                 (fast
                   (do ((i (1+ r) (1+ i))) ((> i m))
                     (declare (type fixnum i))
                     (setf (aref uq (1- i)) (aref uq i)))
                   (do ((i (1+ fout) (1+ i))) ((> i n))
                     (declare (type fixnum i))
                     (%copy-c-pt-row qw i qw (1- i))))
                 (decf n)
                 (decf m)
                 ;; no more internal knots: finished
                 (when (= n p)
                   (return))
                 ;; update error bounds
                 (let ((k (nl-max (- r p) (1+ p)))
                       (l (nl-min n (- (+ r p) s))))
                   (declare (type fixnum k l))
                   (do ((i k (1+ i))) ((> i l))
                     (declare (type fixnum i))
                     (when (and (/= (aref uq i) (aref uq (1+ i))) (/= (aref br i) +crv-adv--norem+))
                       (setf (aref br i) (crv-remove-knot-max-err curq i (aref sr i)))))))))
            (t
             ;; knot not removable
             (setf (aref br r) +crv-adv--norem+)))))
      ;; compact the output curve: it keeps the knots and control points left,
      ;; as N_CrvSetSizeIndices records them in C
      (when (< n ns)
        (setf (curve-knots curq) (subseq uq 0 (1+ m))
              (curve-cpoints curq) (let ((q (make-cpoint-array (1+ n))))
                                     (dotimes (i (1+ n) q)
                                       (%copy-c-pt-row qw i q i))))
        (crv-compress curq))
      curq)))

;;; ------------------------------------------------------------- L3-g6

(defconstant +crv-adv--cto+ 1d-5 "cto in NL_CrvAdv.cpp: relative tolerance for basis maxima.")

(defun crv-remove-knots-deriv-constraints (curp tol whr der)
  "A fresh copy of curve CURP with every knot removed whose removal keeps the
accumulated error bound within TOL (scaled by the minimum weight for a
rational curve), cheapest knot first, the bound weighted by the maxima of the
basis functions involved.  WHR (+no+, +start+, +end+ or +both+) keeps the
knots that carry the 1st to DERth derivatives at the start and/or end;
knots whose removal would push a weight outside [*wmin*, *wmax*] stay.  The
C writes into a caller's curve, or into CURP itself; Batten returns the new
curve.  Signals what the C reports: the error of the last failed
basis-find-all-span-maxima when no knot at all was removed (the C carries on
past that failure with fallback maxima, and its error flag survives only
then).
NLib: N_CrvRemoveKnotsDerivConstraints."
  (declare (type curve curp) (type double-float tol) (type fixnum whr der))
  (let* ((n (curve-n curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (ns n)
         (rat nil)
         (lam 0d0) (oml 0d0)
         (lto 0d0)
         ;; the C's error flag, left set by a failed basis-find-all-span-maxima
         (pending nil))
    (declare (type fixnum n p m ns) (type real-vector up) (type double-float lam oml lto))
    ;; removal tolerance for rational curves
    (when (is-crv-rat curp)
      (multiple-value-bind (wmin wmax pmin pmax) (crv-get-min-max-weights-and-pts curp)
        (declare (type double-float wmin pmax) (ignore wmax pmin))
        (setf tol (/ (* tol wmin) (+ 1d0 pmax))
              rat t)))
    (setf lto (* +crv-adv--cto+ (abs (- (aref up m) (aref up 0)))))
    (let* ((curq (crv-copy curp))
           (qw (curve-cpoints curq))
           (uq (curve-knots curq))
           (rw (make-cpoint-array (1+ (* 2 p))))
           (br (make-real-vector (1+ m)))
           (sr (make-index-vector (1+ m)))
           (er (make-real-vector (1+ m)))
           (te (make-real-vector (1+ m)))
           (minl (make-real-vector (1+ p)))
           (maxl (make-real-vector (1+ p)))
           (minr (make-real-vector (1+ p)))
           (maxr (make-real-vector (1+ p)))
           (mx (make-real-vector (+ p 2))))
      (declare (type curve curq) (type cpoint-array qw rw)
               (type real-vector uq br er te minl maxl minr maxr mx)
               (type index-vector sr))
      (dotimes (i (1+ m))
        (setf (aref br i) +bigd+
              (aref sr i) 0
              (aref er i) 0d0))
      ;; knot removal errors for each distinct knot
      (let ((r (1+ p)))
        (declare (type fixnum r))
        (do () ((not (<= r n)))
          (let ((i r))
            (declare (type fixnum i))
            (do () ((not (and (<= r n) (= (aref uq r) (aref uq (1+ r)))))) (incf r))
            (setf (aref sr r) (1+ (- r i)))
            (setf (aref br r) (crv-remove-knot-max-err curp r (aref sr r)))
            (incf r))))
      (flet ((span-maxima (knt i mins maxs)
               ;; the C's knt is curQ's knot vector, whose highest index
               ;; follows the removals; a failure is recorded, not raised
               (handler-case
                   (progn (basis-find-all-span-maxima knt i p lto mins maxs)
                          (setf pending nil)
                          t)
                 (nurbs-error (c)
                   (setf pending c)
                   nil))))
        ;; try to remove each knot
        (do () (nil)
          (let ((b (aref br (1+ p)))
                (s (aref sr (1+ p)))
                (r (1+ p))
                (k 0) (l 0)
                (rmf t))
            (declare (type double-float b) (type fixnum s r k l))
            ;; knot with the smallest error
            (do ((i (+ p 2) (1+ i))) ((> i (- m p 1)))
              (declare (type fixnum i))
              (when (< (aref br i) b)
                (setf b (aref br i)
                      s (aref sr i)
                      r i)))
            ;; no more removable knots: finished
            (when (or (= b +bigd+) (= b +crv-adv--norem+))
              (return))
            ;; derivative constraints
            (when (or (= whr +start+) (= whr +both+))
              (when (<= r (+ p der))
                (setf rmf nil)))
            (when (or (= whr +end+) (= whr +both+))
              (when (>= r (+ (- n der) 1))
                (setf rmf nil)))
            (cond
              ((not rmf)
               (setf (aref br r) +crv-adv--norem+))
              (t
               ;; error of removal
               (let ((knt (subseq uq 0 (1+ m))))
                 (declare (type real-vector knt))
                 (cond
                   ((oddp (+ p s))
                    ;; maxima of the basis functions over each span
                    (setf k (truncate (+ p s 1) 2)
                          l (+ (- r k) p 1))
                    (let* ((alf (/ (- (aref uq r) (aref uq (- r k)))
                                   (- (aref uq (+ (- r k) p 1)) (aref uq (- r k)))))
                           (bet (/ (- (aref uq r) (aref uq (+ (- r k) 1)))
                                   (- (aref uq (+ (- r k) p 2)) (aref uq (+ (- r k) 1)))))
                           (omb (- 1d0 bet)))
                      (declare (type double-float alf bet omb))
                      (setf lam (/ alf (+ alf bet))
                            oml (- 1d0 lam))
                      (unless (span-maxima knt (- r k) minl maxl)
                        (dotimes (i (1+ p))
                          (setf (aref minl i) 0d0
                                (aref maxl i) 1d0)))
                      (unless (span-maxima knt (+ (- r k) 1) minr maxr)
                        (dotimes (i (1+ p))
                          (setf (aref minr i) 0d0
                                (aref maxr i) 1d0)))
                      (setf (aref mx 0) (* (* lam alf) (aref maxl 0)))
                      (do ((i 1 (1+ i))) ((> i p))
                        (declare (type fixnum i))
                        (setf (aref minl i) (* (aref minl i) (* lam alf))
                              (aref minr (1- i)) (* (aref minr (1- i)) (* oml omb))
                              (aref maxl i) (* (aref maxl i) (* lam alf))
                              (aref maxr (1- i)) (* (aref maxr (1- i)) (* oml omb)))
                        (setf (aref mx i) (nl-max (abs (- (aref maxl i) (aref minr (1- i))))
                                                  (abs (- (aref maxr (1- i)) (aref minl i))))))
                      (setf (aref mx (1+ p)) (* (* oml omb) (aref maxr p)))))
                   (t
                    ;; maximum of the basis function
                    (setf k (truncate (+ p s) 2)
                          l (+ (- r k) p))
                    (unless (span-maxima knt (- r k) minl mx)
                      (dotimes (i (1+ p))
                        (setf (aref mx i) 1d0))))))
               ;; check the error
               (do ((i (- r k) (1+ i))) ((> i l))
                 (declare (type fixnum i))
                 (when (/= (aref uq i) (aref uq (1+ i)))
                   (setf (aref te i) (+ (aref er i) (* (aref mx (+ (- i r) k)) b)))
                   (when (> (aref te i) tol)
                     (setf rmf nil)
                     (return))))
               (cond
                 (rmf
                  ;; error test passed: update the error vector and remove the knot
                  (do ((i (- r k) (1+ i))) ((> i l))
                    (declare (type fixnum i))
                    (when (/= (aref uq i) (aref uq (1+ i)))
                      (setf (aref er i) (aref te i))))
                  (let* ((fout (truncate (- (* 2 r) s p) 2))
                         (first (- r p))
                         (last (- r s))
                         (off (1- first))
                         (i first)
                         (j last)
                         (ii 1)
                         (jj (- last off))
                         (skip nil))
                    (declare (type fixnum fout first last off i j ii jj))
                    (%copy-c-pt-row qw off rw 0)
                    (%copy-c-pt-row qw (1+ last) rw (- (1+ last) off))
                    ;; new control points for one removal step
                    (fast
                      (do () ((not (> (- j i) 0)))
                        (let* ((alf (/ (- (aref uq (+ i p 1)) (aref uq i)) (- (aref uq r) (aref uq i))))
                               (oma (- 1d0 alf))
                               (bet (/ (- (aref uq (+ j p 1)) (aref uq j)) (- (aref uq (+ j p 1)) (aref uq r))))
                               (omb (- 1d0 bet)))
                          (declare (type double-float alf oma bet omb))
                          (%combine-2-c-pt-rows alf qw i oma rw (1- ii) rw ii)
                          (%combine-2-c-pt-rows bet qw j omb rw (1+ jj) rw jj)
                          (incf i) (decf j) (incf ii) (decf jj))))
                    ;; disallowed weights
                    (when rat
                      (let ((i first) (j last) (wmin +bigd+) (wmax +smad+))
                        (declare (type fixnum i j) (type double-float wmin wmax))
                        (do () ((not (> (- j i) 0)))
                          (let ((wi (aref rw (- i off) 3))
                                (wj (aref rw (- j off) 3)))
                            (declare (type double-float wi wj))
                            (when (< wi wmin) (setf wmin wi))
                            (when (< wj wmin) (setf wmin wj))
                            (when (> wi wmax) (setf wmax wi))
                            (when (> wj wmax) (setf wmax wj))
                            (incf i) (decf j)))
                        (when (or (< wmin (the double-float *wmin*)) (> wmax (the double-float *wmax*)))
                          (setf (aref br r) +crv-adv--norem+
                                skip t))))
                    (unless skip
                      ;; save new control points
                      (when (oddp (+ p s))
                        (%combine-2-c-pt-rows lam rw (1+ jj) oml rw (1- ii) rw (1+ jj)))
                      (setf i first
                            j last)
                      (fast
                        (do () ((not (> (- j i) 0)))
                          (%copy-c-pt-row rw (- i off) qw i)
                          (%copy-c-pt-row rw (- j off) qw j)
                          (incf i) (decf j)))
                      ;; shift down some parameters
                      (when (= s 1)
                        (setf (aref er (1- r)) (nl-max (aref er (1- r)) (aref er r))))
                      (when (> s 1)
                        (setf (aref sr (1- r)) (1- (aref sr r))))
                      (do ((i (1+ r) (1+ i))) ((> i m))
                        (declare (type fixnum i))
                        (setf (aref br (1- i)) (aref br i)
                              (aref sr (1- i)) (aref sr i)
                              (aref er (1- i)) (aref er i)))
                      ;; shift down knots and control points
                      (fast
                        (do ((i (1+ r) (1+ i))) ((> i m))
                          (declare (type fixnum i))
                          (setf (aref uq (1- i)) (aref uq i)))
                        (do ((i (1+ fout) (1+ i))) ((> i n))
                          (declare (type fixnum i))
                          (%copy-c-pt-row qw i qw (1- i))))
                      (decf n)
                      (decf m)
                      ;; no more internal knots: finished
                      (when (= n p)
                        (return))
                      ;; update error bounds
                      (let ((k (nl-max (- r p) (1+ p)))
                            (l (nl-min n (- (+ r p) s))))
                        (declare (type fixnum k l))
                        (do ((i k (1+ i))) ((> i l))
                          (declare (type fixnum i))
                          (when (and (/= (aref uq i) (aref uq (1+ i))) (/= (aref br i) +crv-adv--norem+))
                            (setf (aref br i) (crv-remove-knot-max-err curq i (aref sr i))
                                  pending nil)))))))
                 (t
                  ;; knot not removable
                  (setf (aref br r) +crv-adv--norem+))))))))
      ;; compact the output curve: it keeps the knots and control points left,
      ;; as N_CrvSetSizeIndices records them in C
      (when (< n ns)
        (setf (curve-knots curq) (subseq uq 0 (1+ m))
              (curve-cpoints curq) (let ((q (make-cpoint-array (1+ n))))
                                     (dotimes (i (1+ n) q)
                                       (%copy-c-pt-row qw i q i))))
        (crv-compress curq)
        (setf pending nil))
      (when pending
        (error pending))
      curq)))
