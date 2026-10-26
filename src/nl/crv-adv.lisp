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
    (%row-pt dd 1 d1)
    (%row-pt dd 2 d2)
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

;;; ------------------------------------------------------------- L4-g0

(defun %crv-adv-sum-diff-rows (opr pw qw rw n rname)
  "Rows 0..N of cpoint-array RW = PW[i] + QW[i] (OPR +plus+) or PW[i] -
QW[i] (+minus+), N_Sum2CPts / N_Diff2CPts; RW may be PW.  Signals :cal-err
for another OPR (before writing anything)."
  (declare (type fixnum opr n) (type cpoint-array pw qw rw))
  (unless (or (= opr +plus+) (= opr +minus+))
    (nurbs-error :cal-err rname))
  (let ((a (make-cpoint)) (b (make-cpoint)))
    (declare (dynamic-extent a b))
    (dotimes (i (1+ n) rw)
      (%row-c-pt pw i a)
      (%row-c-pt qw i b)
      (if (= opr +plus+)
          (sum-2-c-pts a b a)
          (diff-2-c-pts a b a))
      (%c-pt-row a rw i))))

(defun crv-sum-diff-crv (cur-p cur-q opr)
  "The sum (OPR +plus+) or difference (+minus+) of curves CUR-P and CUR-Q, a
fresh curve.  Rational: N_p/D_p +/- N_q/D_q = (N_p D_q +/- N_q D_p)/(D_p D_q)
by the symbolic product operators, a non-rational input being made rational
in place for the computation and made non-rational again at the end (an
error on the way leaves it rational, as in the C).  Non-rational: copies of
the curves are degree elevated to the higher degree and refined to common
knots, and their control points added or subtracted.  Signals :cal-err for
another OPR, and what the product, elevation and refinement operators
signal.
NLib: N_CrvSumDiffCrv."
  (declare (type curve cur-p cur-q) (type fixnum opr))
  (let ((rat-p (is-crv-rat cur-p))
        (rat-q (is-crv-rat cur-q)))
    (cond
      ((or rat-p rat-q)
       ;; make a non-rational curve rational
       (unless rat-p (crv-non-rat-to-rat cur-p))
       (unless rat-q (crv-non-rat-to-rat cur-q))
       ;; numerators and denominators
       (multiple-value-bind (num-p den-p) (crv-get-num-and-denom cur-p)
         (multiple-value-bind (num-q den-q) (crv-get-num-and-denom cur-q)
           (let* ((den (crv-func-multiply-crv-func den-p den-q))
                  (num (crv-func-multiply-crv den-q num-p))
                  (tmp (crv-func-multiply-crv den-p num-q))
                  (pw (curve-cpoints num))
                  (qw (curve-cpoints tmp))
                  (n (1- (array-dimension qw 0))))
             (declare (type fixnum n))
             ;; add the terms of the numerator, in place
             (%crv-adv-sum-diff-rows opr pw qw pw n 'crv-sum-diff-crv)
             (let ((cur-r (create-crv-from-num-and-denom num den)))
               (unless rat-p (crv-rat-to-non-rat cur-p))
               (unless rat-q (crv-rat-to-non-rat cur-q))
               cur-r)))))
      (t
       (let* ((cur-pw (crv-copy cur-p))
              (cur-qw (crv-copy cur-q))
              (p (curve-degree cur-pw))
              (q (curve-degree cur-qw))
              ;; n and m follow the C: the curve whose data were fetched last
              (n (1- (array-dimension (curve-cpoints cur-qw) 0)))
              (m (1- (length (curve-knots cur-qw)))))
         (declare (type fixnum p q n m))
         ;; degree elevate
         (when (/= p q)
           (if (< p q)
               (progn
                 (setf cur-pw (crv-elevate-degree cur-pw (- q p)))
                 (setf p (curve-degree cur-pw)
                       n (1- (array-dimension (curve-cpoints cur-pw) 0))
                       m (1- (length (curve-knots cur-pw)))))
               (progn
                 (setf cur-qw (crv-elevate-degree cur-qw (- p q)))
                 (setf q (curve-degree cur-qw)
                       n (1- (array-dimension (curve-cpoints cur-qw) 0))
                       m (1- (length (curve-knots cur-qw)))))))
         ;; refine to common knots
         (basis-get-span-count (curve-knots cur-pw) p)
         (basis-get-span-count (curve-knots cur-qw) q)
         (multiple-value-bind (kxp kxq)
             (get-compatible-knot-array-mult (curve-knots cur-pw) (curve-knots cur-qw) p)
           (declare (type real-vector kxp kxq))
           (when (> (length kxp) 0)
             (setf cur-pw (crv-refine cur-pw kxp))
             (setf p (curve-degree cur-pw)
                   n (1- (array-dimension (curve-cpoints cur-pw) 0))
                   m (1- (length (curve-knots cur-pw)))))
           (when (> (length kxq) 0)
             (setf cur-qw (crv-refine cur-qw kxq))
             (setf q (curve-degree cur-qw)
                   n (1- (array-dimension (curve-cpoints cur-qw) 0))
                   m (1- (length (curve-knots cur-qw))))))
         ;; the output curve
         (let ((rw (make-cpoint-array (1+ n)))
               (ur (make-real-vector (1+ m)))
               (up (curve-knots cur-pw)))
           (declare (type cpoint-array rw) (type real-vector ur up))
           (%crv-adv-sum-diff-rows opr (curve-cpoints cur-pw) (curve-cpoints cur-qw) rw n
                                   'crv-sum-diff-crv)
           (dotimes (i (1+ m))
             (setf (aref ur i) (aref up i)))
           (%make-curve p ur rw)))))))

(defun srf-reparam-func (cur-p cfn)
  "Curve CUR-P reparametrized by the B-spline function CFN (f' > 0 on its
domain [c,d], f(c) and f(d) CUR-P's ends), a fresh curve of degree p*q on
CFN's domain: the curve and the function are refined so that their knots
map onto each other, each Bezier segment pair is composed by BEZ-REPARAM,
and the knots between the composed segments are removed as far as the
continuity allows.  Signals :deg-err when p*q exceeds *dmax*, :num-err when
the refined span counts differ, :knt-err for an interior knot of
multiplicity above the degree, and what the knot mapping and refinement
signal.  (The C's in-place case, curP = curQ, computes the same curve.)
NLib: N_SrfReparamFunc."
  (declare (type curve cur-p) (type cfun cfn))
  (let* ((p (curve-degree cur-p))
         (q (cfun-degree cfn))
         (pq (* p q))
         (ptol *ptol*))
    (declare (type fixnum p q pq) (type double-float ptol))
    (when (> pq (the fixnum *dmax*))
      (nurbs-error :deg-err 'srf-reparam-func))
    ;; map the knot vectors onto each other
    (basis-get-span-count (curve-knots cur-p) p)
    (basis-get-span-count (cfun-knots cfn) q)
    (multiple-value-bind (xu xs)
        (map-knots-between-crv-func-and-knot-vector cfn (curve-knots cur-p) p)
      (declare (type real-vector xu xs))
      (let* ((mxu (1- (length xu)))
             (mxs (1- (length xs)))
             ;; refine the curve and the function
             (cur-r (if (>= mxu 0) (crv-refine cur-p xu) cur-p))
             (cfn-r (if (>= mxs 0) (crv-func-refine cfn xs) cfn))
             (pw (curve-cpoints cur-r))
             (up (curve-knots cur-r))
             (mp (1- (length up)))
             (fu (cfun-vals cfn-r))
             (s (cfun-knots cfn-r))
             (ms (1- (length s))))
        (declare (type fixnum mxu mxs mp ms) (type cpoint-array pw) (type real-vector up fu s))
        (setf p (curve-degree cur-r)
              q (cfun-degree cfn-r))
        (let ((i (basis-get-span-count up p))
              (j (basis-get-span-count s q)))
          (declare (type fixnum i j))
          (when (/= i j)
            (nurbs-error :num-err 'srf-reparam-func))
          (let* ((nq (* i pq))
                 (mq (+ nq pq 1))
                 (qw (make-cpoint-array (1+ nq)))
                 (uq (make-real-vector (1+ mq)))
                 ;; local memory
                 (bw (make-cpoint-array (1+ p)))
                 (nw (make-cpoint-array (1+ p)))
                 (rw (make-cpoint-array (1+ pq)))
                 (fb (make-real-vector (1+ q)))
                 (fn (make-real-vector (1+ q)))
                 (alfu (make-real-vector (1+ p)))
                 (omau (make-real-vector (1+ p)))
                 (alfs (make-real-vector (1+ q)))
                 (omas (make-real-vector (1+ q)))
                 ;; initialize
                 (au p) (bu (1+ p))
                 (as q) (bs (1+ q))
                 (ru -1) (rs -1)
                 (ixu 0) (ixs 0)
                 (cind 1)
                 (kind (1+ pq))
                 (rem -1))
            (declare (type fixnum nq mq au bu as bs ru rs ixu ixs cind kind rem)
                     (type cpoint-array qw bw nw rw) (type real-vector uq fb fn alfu omau alfs omas))
            (dotimes (k (1+ pq))
              (setf (aref uq k) (aref s as)))
            (%copy-c-pt-row pw 0 qw 0)
            (dotimes (k (1+ p))
              (%copy-c-pt-row pw k bw k))
            (dotimes (k (1+ q))
              (setf (aref fb k) (aref fu k)))
            ;; loop through the knot vectors: extract the Bezier segments,
            ;; reparametrize each, remove the knot between it and the last
            (do () ((not (and (< bu mp) (< bs ms))))
              (let ((mlu 0) (miu 0) (mls 0) (mis 0) (mi 0) (tt 0) (lbz 0) (rbz 0))
                (declare (type fixnum mlu miu mls mis mi tt lbz rbz))
                ;; multiplicities of the knots
                (let ((i0 bu))
                  (declare (type fixnum i0))
                  (do () ((not (and (< bu mp) (= (aref up bu) (aref up (1+ bu))))))
                    (incf bu))
                  (setf mlu (1+ (- bu i0))
                        miu mlu))
                (let ((j0 bs))
                  (declare (type fixnum j0))
                  (do () ((not (and (< bs ms) (= (aref s bs) (aref s (1+ bs))))))
                    (incf bs))
                  (setf mls (1+ (- bs j0))
                        mis mls))
                ;; multiplicities of internal knots
                (when (and (< bu mp) (< bs ms))
                  (when (or (> mlu p) (> mls q))
                    (nurbs-error :knt-err 'srf-reparam-func)))
                ;; adjust the multiplicities; the output knot's multiplicity
                (when (and (= miu 1) (<= ixu mxu))
                  (when (< (abs (- (aref up bu) (aref xu ixu))) ptol)
                    (setf miu 0)
                    (incf ixu)))
                (when (and (= mis 1) (<= ixs mxs))
                  (when (< (abs (- (aref s bs) (aref xs ixs))) ptol)
                    (setf mis 0)
                    (incf ixs)))
                (setf mi (cond ((= mis 0) (+ (- pq p) miu))
                               ((= miu 0) (+ (- pq q) mis))
                               (t (nl-max (+ (- pq p) miu) (+ (- pq q) mis)))))
                ;; insert knots to get the Bezier segments
                (setf ru (- p mlu)
                      rs (- q mls)
                      tt (- pq mi))
                (setf lbz (if (> rem 0) (truncate (+ rem 2) 2) 1))
                (setf rbz (if (> tt 0) (- pq (truncate (+ tt 1) 2)) pq))
                (when (> ru 0)
                  (let ((num (- (aref up bu) (aref up au))))
                    (declare (type double-float num))
                    (do ((k p (1- k))) ((not (> k mlu)))
                      (declare (type fixnum k))
                      (setf (aref alfu (- k mlu 1)) (/ num (- (aref up (+ au k)) (aref up au))))
                      (setf (aref omau (- k mlu 1)) (- 1.0d0 (aref alfu (- k mlu 1)))))
                    (do ((jj 1 (1+ jj))) ((> jj ru))
                      (declare (type fixnum jj))
                      (let ((save (- ru jj))
                            (s0 (+ mlu jj)))
                        (declare (type fixnum save s0))
                        (do ((k p (1- k))) ((< k s0))
                          (declare (type fixnum k))
                          (%combine-2-c-pt-rows (aref alfu (- k s0)) bw k (aref omau (- k s0)) bw (1- k) bw k))
                        (%copy-c-pt-row bw p nw save)))))
                (when (> rs 0)
                  (let ((num (- (aref s bs) (aref s as))))
                    (declare (type double-float num))
                    (do ((k q (1- k))) ((not (> k mls)))
                      (declare (type fixnum k))
                      (setf (aref alfs (- k mls 1)) (/ num (- (aref s (+ as k)) (aref s as))))
                      (setf (aref omas (- k mls 1)) (- 1.0d0 (aref alfs (- k mls 1)))))
                    (do ((jj 1 (1+ jj))) ((> jj rs))
                      (declare (type fixnum jj))
                      (let ((save (- rs jj))
                            (s0 (+ mls jj)))
                        (declare (type fixnum save s0))
                        (do ((k q (1- k))) ((< k s0))
                          (declare (type fixnum k))
                          (setf (aref fb k) (+ (* (aref alfs (- k s0)) (aref fb k))
                                               (* (aref omas (- k s0)) (aref fb (1- k))))))
                        (setf (aref fn save) (aref fb q))))))
                ;; reparametrize the Bezier segment
                (let* ((f0 (aref fb 0))
                       (fac (/ 1.0d0 (- (aref fb q) (aref fb 0)))))
                  (declare (type double-float f0 fac))
                  (setf (aref fb 0) 0.0d0
                        (aref fb q) 1.0d0)
                  (do ((k 1 (1+ k))) ((>= k q))
                    (declare (type fixnum k))
                    (setf (aref fb k) (* fac (- (aref fb k) f0)))))
                (bez-reparam bw p fb q nil rw)
                ;; remove the knot S[as]
                (when (> rem 1)
                  (let ((frst (- kind 2))
                        (lst kind)
                        (den (- (aref s bs) (aref s as))))
                    (declare (type fixnum frst lst) (type double-float den))
                    (do ((k 1 (1+ k))) ((>= k rem))
                      (declare (type fixnum k))
                      (let* ((ii frst)
                             (jj lst)
                             (l (+ (- jj kind) 1)))
                        (declare (type fixnum ii jj l))
                        (do () ((not (> (- jj ii) k)))
                          (when (< ii cind)
                            (let* ((alf (/ (- (aref s bs) (aref uq ii)) (- (aref s as) (aref uq ii))))
                                   (oma (- 1.0d0 alf)))
                              (declare (type double-float alf oma))
                              (%combine-2-c-pt-rows alf qw ii oma qw (1- ii) qw ii)))
                          (when (>= jj lbz)
                            (let* ((bet (/ (- (aref s bs) (aref uq (- jj k))) den))
                                   (omb (- 1.0d0 bet)))
                              (declare (type double-float bet omb))
                              (%combine-2-c-pt-rows bet rw l omb rw (1+ l) rw l)))
                          (incf ii)
                          (decf jj)
                          (decf l)))
                      (decf frst)
                      (incf lst))))
                ;; load the knot vector and the control points
                (when (/= au p)
                  (dotimes (k (- pq rem))
                    (setf (aref uq kind) (aref s as))
                    (incf kind)))
                (do ((k lbz (1+ k))) ((> k rbz))
                  (declare (type fixnum k))
                  (%copy-c-pt-row rw k qw cind)
                  (incf cind))
                ;; initialize for the next pass
                (cond
                  ((and (< bu mp) (< bs ms))
                   (dotimes (k ru)
                     (%copy-c-pt-row nw k bw k))
                   (do ((k ru (1+ k))) ((> k p))
                     (declare (type fixnum k))
                     (%copy-c-pt-row pw (+ (- bu p) k) bw k))
                   (dotimes (k rs)
                     (setf (aref fb k) (aref fn k)))
                   (do ((k rs (1+ k))) ((> k q))
                     (declare (type fixnum k))
                     (setf (aref fb k) (aref fu (+ (- bs q) k))))
                   (setf au bu)
                   (incf bu)
                   (setf rem (- pq mi))
                   (setf as bs)
                   (incf bs))
                  (t
                   ;; the end knot
                   (dotimes (k (1+ pq))
                     (setf (aref uq kind) (aref s bs))
                     (incf kind))))))
            ;; compact the output curve
            (let ((cur-q (%make-curve pq uq qw)))
              (when (< cind (1+ nq))
                (setf (curve-knots cur-q) (subseq uq 0 kind)
                      (curve-cpoints cur-q) (let ((c (make-cpoint-array cind)))
                                              (dotimes (k cind c)
                                                (%copy-c-pt-row qw k c k))))
                (crv-compress cur-q))
              cur-q)))))))

;;; ------------------------------------------------------------- L5-g1

(defun %crv-adv-become-srf (sur new)
  "Make surface SUR, in place, the surface NEW (sharing NEW's arrays), as the
C's in-place operators (surP = surQ) overwrite their input.  Returns SUR."
  (declare (type surface sur new))
  (setf (surface-degree-u sur) (surface-degree-u new)
        (surface-degree-v sur) (surface-degree-v new)
        (surface-knots-u sur) (surface-knots-u new)
        (surface-knots-v sur) (surface-knots-v new)
        (surface-cpoints sur) (surface-cpoints new))
  sur)

(defun crv-make-compatible-with-srf (cur sur dir ts te)
  "Make curve CUR and surface SUR's DIR direction (+udir+, else the v
direction) compatible, both modified in place: the same rationality (the
non-rational one gets unit weights), both reparametrized onto [TS, TE], the
lower degree elevated to the higher, and the knot vectors merged (knots closer
than 1e-8 moved to their average, in place) and inserted.  Returns CUR.
Signals what the elevation, the merge and the refinement signal (:cal-err
from the surface operators for a DIR that is neither +udir+ nor +vdir+).
NLib: N_CrvMakeCompatibleWithSrf."
  (declare (type curve cur) (type surface sur) (type fixnum dir) (type double-float ts te))
  (let* ((pc (curve-degree cur))
         (uc (curve-knots cur))
         (mc (1- (length uc)))
         (ratc (is-crv-rat cur))
         (rats (is-srf-rat sur)))
    (declare (type fixnum pc mc) (type real-vector uc))
    ;; make the rationality of each compatible
    (when (and ratc (not rats))
      (srf-non-rat-to-rat sur))
    (when (and rats (not ratc))
      (crv-non-rat-to-rat cur))
    ;; the u/v-directional data
    (let* ((ps (if (= dir +udir+) (surface-degree-u sur) (surface-degree-v sur)))
           (us (if (= dir +udir+) (surface-knots-u sur) (surface-knots-v sur)))
           (ms (1- (length us)))
           (pint (make-real-vector 2)))
      (declare (type fixnum ps ms) (type real-vector us pint))
      ;; reparametrize to a common interval
      (create-interval pint ts te)
      (when (or (/= (aref uc 0) ts) (/= (aref uc mc) te))
        (basis-reparam uc pc pint))
      (when (or (/= (aref us 0) ts) (/= (aref us ms) te))
        (basis-reparam us ps pint))
      ;; degree elevate if necessary
      (cond ((< pc ps)
             (%crv-become cur (crv-elevate-degree cur (- ps pc))))
            ((< ps pc)
             (%crv-adv-become-srf sur (srf-elevate-degree sur (- pc ps) dir))))
      ;; merge the knot vectors
      (let* ((knt (vector (curve-knots cur)
                          (if (= dir +udir+) (surface-knots-u sur) (surface-knots-v sur))))
             (knx (get-compatible-knot-array knt 1))
             (kx0 (svref knx 0))
             (kx1 (svref knx 1)))
        (declare (type simple-vector knt knx) (type real-vector kx0 kx1))
        ;; refine the curve and/or the surface
        (when (>= (1- (length kx0)) 0)
          (%crv-become cur (crv-refine cur kx0)))
        (when (>= (1- (length kx1)) 0)
          (%crv-adv-become-srf sur (srf-insert-knots sur kx1 dir)))))
    cur))

(defun crv-get-rat-second-deriv (cur)
  "The second derivative curve of curve CUR, fresh.  A non-rational CUR gives
CRV-NON-RAT-GET-FIRST-DERIV of order 2; a rational one N/D gives, by the
symbolic operators, (D D N'' - 2 D D' N' + N (2 D' D' - D D'')) / D^3, the
numerator elevated by 2 to the denominator's degree.  Signals what the
symbolic operators signal (:inp-err for degree 1, :deg-err when 3p exceeds
*dmax*).
NLib: N_CrvGetRatSecondDeriv."
  (declare (type curve cur))
  (unless (is-crv-rat cur)
    (return-from crv-get-rat-second-deriv (crv-non-rat-get-first-deriv cur 2)))
  ;; numerator and denominator
  (multiple-value-bind (n d) (crv-get-num-and-denom cur)
    (declare (type curve n) (type cfun d))
    ;; the denominator of the output curve
    (let* ((dd (crv-func-multiply-crv-func d d))
           (den (crv-func-multiply-crv-func d dd)))
      (declare (type cfun dd den))
      (let ((knx (basis-increase-knot-mult (cfun-knots den) (cfun-degree den) 2)))
        (declare (type real-vector knx))
        (when (>= (1- (length knx)) 0)
          (setf den (crv-func-refine den knx))))
      ;; the terms of the numerator: D D N''
      (let* ((npp (crv-non-rat-get-first-deriv n 2))
             (ddnpp (crv-func-multiply-crv dd npp))
             ;; 2 D D' N'
             (dp (crv-func-deriv d 1))
             (ddp (crv-func-multiply-crv-func d dp)))
        (declare (type curve npp ddnpp) (type cfun dp ddp))
        (crv-func-multiply-constant 2.0d0 ddp)
        (let* ((np (crv-non-rat-get-first-deriv n 1))
               (ddpnp (crv-func-multiply-crv ddp np)))
          (declare (type curve np ddpnp))
          (let ((knx (basis-increase-knot-mult (curve-knots ddpnp) (curve-degree ddpnp) 1)))
            (declare (type real-vector knx))
            (when (>= (1- (length knx)) 0)
              (setf ddpnp (crv-refine ddpnp knx))))
          ;; N (2 D' D' - D D'')
          (let ((dpdp (crv-func-multiply-crv-func dp dp)))
            (declare (type cfun dpdp))
            (crv-func-multiply-constant 2.0d0 dpdp)
            (let ((knx (basis-increase-knot-mult (cfun-knots dpdp) (cfun-degree dpdp) 1)))
              (declare (type real-vector knx))
              (when (>= (1- (length knx)) 0)
                (setf dpdp (crv-func-refine dpdp knx))))
            (let* ((dpp (crv-func-deriv d 2))
                   (ddpp (crv-func-multiply-crv-func d dpp))
                   (dsum (crv-func-sum-diff-crv-func dpdp ddpp +minus+ +yes+))
                   (ddppn (crv-func-multiply-crv dsum n))
                   ;; add the terms of the numerator
                   (nsum (crv-sum-diff-crv ddnpp ddpnp +minus+))
                   (num (crv-sum-diff-crv nsum ddppn +plus+)))
              (declare (type cfun dpp ddpp dsum) (type curve ddppn nsum num))
              ;; degree elevate the numerator; the output curve
              (create-crv-from-num-and-denom (crv-elevate-degree num 2) den))))))))

(defun crv-rat-get-first-deriv (cur)
  "The first derivative curve of curve CUR, fresh.  A non-rational CUR gives
CRV-NON-RAT-GET-FIRST-DERIV of order 1; a rational one N/D gives, by the
symbolic operators, (N' D - N D') / D^2, the numerator elevated by 1 to the
denominator's degree.  Signals what the symbolic operators signal (:deg-err
when 2p exceeds *dmax*).
NLib: N_CrvRatGetFirstDeriv."
  (declare (type curve cur))
  (unless (is-crv-rat cur)
    (return-from crv-rat-get-first-deriv (crv-non-rat-get-first-deriv cur 1)))
  ;; numerator and denominator
  (multiple-value-bind (n d) (crv-get-num-and-denom cur)
    (declare (type curve n) (type cfun d))
    ;; the denominator of the output curve
    (let ((den (crv-func-multiply-crv-func d d)))
      (declare (type cfun den))
      (let ((knx (basis-increase-knot-mult (cfun-knots den) (cfun-degree den) 1)))
        (declare (type real-vector knx))
        (when (>= (1- (length knx)) 0)
          (setf den (crv-func-refine den knx))))
      ;; the terms of the numerator
      (let* ((np (crv-non-rat-get-first-deriv n 1))
             (dp (crv-func-deriv d 1))
             (npd (crv-func-multiply-crv d np))
             (ndp (crv-func-multiply-crv dp n))
             ;; add the terms of the numerator
             (num (crv-sum-diff-crv npd ndp +minus+)))
        (declare (type curve np npd ndp num) (type cfun dp))
        ;; degree elevate the numerator; the output curve
        (create-crv-from-num-and-denom (crv-elevate-degree num 1) den)))))

;;; ------------------------------------------------------------- L6-g1

(defun crvs-make-compatible-approx (cur kk tol)
  "Make the curves CUR[0..KK] (a simple-vector of curves, each modified in
place) compatible: exactly first (crvs-make-compatible-knot-tol with knot
tolerance *ptol*), then, unless TOL is below *ptol*, by removing from all
curves together every knot whose removal keeps each curve's accumulated
error bound within TOL (scaled by the smallest weight ratio over the curves
when they are rational), the knot of smallest error sum first, the bound
weighted by the maxima of the basis functions involved; knots whose removal
would push a weight of any curve outside [*wmin*, *wmax*] stay.  Returns
CUR.  Signals what crvs-make-compatible-knot-tol signals, and the error of
the last failed basis-find-all-span-maxima when no knot at all was removed
(the C carries on past that failure with fallback maxima, and its error flag
survives only then).
NLib: N_CrvsMakeCompatibleApprox."
  (declare (type (or null simple-vector) cur) (type fixnum kk) (type double-float tol))
  ;; make the curves exactly compatible, then remove unwanted knots
  (crvs-make-compatible-knot-tol cur kk (the double-float *ptol*))
  (when (or (null cur) (< tol (the double-float *ptol*)))
    (return-from crvs-make-compatible-approx cur))
  (let* ((ck (svref cur kk))
         (n (curve-n ck))
         (p (curve-degree ck))
         (m (curve-m ck))
         ;; the C's U: curve CUR[KK]'s knot vector, shifted in place
         (u (curve-knots ck))
         (ns n)
         (rat nil)
         (s 0) (r 0)
         (lam 0d0) (oml 0d0)
         (lto 0d0)
         ;; the C's error flag, left set by a failed basis-find-all-span-maxima
         (pending nil))
    (declare (type curve ck) (type fixnum n p m ns s r) (type real-vector u)
             (type double-float lam oml lto))
    ;; removal tolerance for rational curves
    (when (is-crv-rat (svref cur 0))
      (let ((alf 0d0))
        (declare (type double-float alf))
        (multiple-value-bind (wmin wmax pmin pmax) (crv-get-min-max-weights-and-pts (svref cur 0))
          (declare (type double-float wmin pmax) (ignore wmax pmin))
          (setf alf (/ (* tol wmin) (+ 1d0 pmax))))
        (do ((ll 1 (1+ ll))) ((> ll kk))
          (declare (type fixnum ll))
          (multiple-value-bind (wmin wmax pmin pmax) (crv-get-min-max-weights-and-pts (svref cur ll))
            (declare (type double-float wmin pmax) (ignore wmax pmin))
            (let ((bet (/ (* tol wmin) (+ 1d0 pmax))))
              (declare (type double-float bet))
              (when (< bet alf)
                (setf alf bet)))))
        (setf tol alf
              rat t)))
    (setf lto (* +crv-adv--cto+ (abs (- (aref u m) (aref u 0)))))
    (let ((rw (make-array (1+ kk)))
          (br (make-real-matrix (1+ kk) (1+ m)))
          (b (make-real-vector (1+ kk)))
          (sr (make-index-vector (1+ m)))
          (er (make-real-matrix (1+ kk) (1+ m)))
          (te (make-real-matrix (1+ kk) (1+ m)))
          (minl (make-real-vector (1+ p)))
          (maxl (make-real-vector (1+ p)))
          (minr (make-real-vector (1+ p)))
          (maxr (make-real-vector (1+ p)))
          (mx (make-real-vector (+ p 2)))
          (alfs (make-real-vector (1+ (* 2 p))))
          (omas (make-real-vector (1+ (* 2 p))))
          (bets (make-real-vector (1+ (* 2 p))))
          (ombs (make-real-vector (1+ (* 2 p)))))
      (declare (type simple-vector rw) (type real-matrix br er te)
               (type real-vector b minl maxl minr maxr mx alfs omas bets ombs)
               (type index-vector sr))
      (dotimes (ll (1+ kk))
        (setf (svref rw ll) (make-cpoint-array (1+ (* 2 p)))))
      ;; initialize
      (dotimes (j (1+ m))
        (setf (aref sr j) 0)
        (dotimes (ll (1+ kk))
          (setf (aref br ll j) +bigd+
                (aref er ll j) 0d0)))
      ;; knot removal errors for each distinct knot of each curve
      (setf r (1+ p))
      (do () ((not (<= r n)))
        (let ((i r))
          (declare (type fixnum i))
          (do () ((not (and (<= r n) (= (aref u r) (aref u (1+ r)))))) (incf r))
          (setf (aref sr r) (1+ (- r i)))
          (dotimes (ll (1+ kk))
            (setf (aref br ll r) (crv-remove-knot-max-err (svref cur ll) r (aref sr r))))
          (incf r)))
      (flet ((span-maxima (i mins maxs)
               ;; the C's knt is curve CUR[KK]'s knot vector, whose highest
               ;; index follows the removals; a failure is recorded, not raised
               (handler-case
                   (progn (basis-find-all-span-maxima (subseq u 0 (1+ m)) i p lto mins maxs)
                          (setf pending nil)
                          t)
                 (nurbs-error (c)
                   (setf pending c)
                   nil))))
        ;; try to remove each knot from each curve
        (do () (nil)
          (let ((bsum +bigd+)
                (k 0) (l 0)
                (rmf t))
            (declare (type double-float bsum) (type fixnum k l))
            ;; the knot with the smallest error sum
            (do ((j (1+ p) (1+ j))) ((> j (- m p 1)))
              (declare (type fixnum j))
              (when (and (/= (aref br 0 j) +bigd+) (/= (aref br 0 j) +crv-adv--norem+))
                (let ((alf 0d0))
                  (declare (type double-float alf))
                  (dotimes (ll (1+ kk))
                    (setf alf (+ alf (aref br ll j))))
                  (when (< alf bsum)
                    (setf bsum alf
                          s (aref sr j)
                          r j)
                    (dotimes (ll (1+ kk))
                      (setf (aref b ll) (aref br ll j)))))))
            ;; no more removable knots: finished
            (when (= bsum +bigd+)
              (return))
            ;; error of removal
            (cond
              ((oddp (+ p s))
               ;; maxima of the basis functions over each span
               (setf k (truncate (+ p s 1) 2)
                     l (+ (- r k) p 1))
               (let* ((alf (/ (- (aref u r) (aref u (- r k)))
                              (- (aref u (+ (- r k) p 1)) (aref u (- r k)))))
                      (bet (/ (- (aref u r) (aref u (+ (- r k) 1)))
                              (- (aref u (+ (- r k) p 2)) (aref u (+ (- r k) 1)))))
                      (omb (- 1d0 bet)))
                 (declare (type double-float alf bet omb))
                 (setf lam (/ alf (+ alf bet))
                       oml (- 1d0 lam))
                 (unless (span-maxima (- r k) minl maxl)
                   (dotimes (i (1+ p))
                     (setf (aref minl i) 0d0
                           (aref maxl i) 1d0)))
                 (unless (span-maxima (+ (- r k) 1) minr maxr)
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
               (unless (span-maxima (- r k) minl mx)
                 (dotimes (i (1+ p))
                   (setf (aref mx i) 1d0)))))
            ;; check the error
            (do ((j (- r k) (1+ j))) ((> j l))
              (declare (type fixnum j))
              (when (/= (aref u j) (aref u (1+ j)))
                (dotimes (ll (1+ kk))
                  (setf (aref te ll j) (+ (aref er ll j) (* (aref mx (+ (- j r) k)) (aref b ll))))
                  (when (> (aref te ll j) tol)
                    (setf rmf nil)
                    (return))))
              (unless rmf
                (return)))
            (cond
              (rmf
               ;; error test passed: update the error vector and remove the knot
               (do ((j (- r k) (1+ j))) ((> j l))
                 (declare (type fixnum j))
                 (when (/= (aref u j) (aref u (1+ j)))
                   (dotimes (ll (1+ kk))
                     (setf (aref er ll j) (aref te ll j)))))
               (let* ((fout (truncate (- (* 2 r) s p) 2))
                      (first (- r p))
                      (last (- r s))
                      (off (1- first))
                      (wfl t))
                 (declare (type fixnum fout first last off))
                 ;; save alphas and betas
                 (let ((i first) (j last))
                   (declare (type fixnum i j))
                   (do () ((not (> (- j i) 0)))
                     (setf (aref alfs (- i first)) (/ (- (aref u (+ i p 1)) (aref u i))
                                                      (- (aref u r) (aref u i))))
                     (setf (aref omas (- i first)) (- 1d0 (aref alfs (- i first))))
                     (setf (aref bets (- j first)) (/ (- (aref u (+ j p 1)) (aref u j))
                                                      (- (aref u (+ j p 1)) (aref u r))))
                     (setf (aref ombs (- j first)) (- 1d0 (aref bets (- j first))))
                     (incf i) (decf j)))
                 ;; remove the knot from all curves
                 (dotimes (ll (1+ kk))
                   (let ((pw (curve-cpoints (svref cur ll)))
                         (rl (svref rw ll))
                         (i first) (j last) (ii 1) (jj (- last off)))
                     (declare (type cpoint-array pw rl) (type fixnum i j ii jj))
                     (%copy-c-pt-row pw off rl 0)
                     (%copy-c-pt-row pw (1+ last) rl (- (1+ last) off))
                     ;; new control points for one removal step
                     (fast
                       (do () ((not (> (- j i) 0)))
                         (%combine-2-c-pt-rows (aref alfs (- i first)) pw i
                                               (aref omas (- i first)) rl (1- ii) rl ii)
                         (%combine-2-c-pt-rows (aref bets (- j first)) pw j
                                               (aref ombs (- j first)) rl (1+ jj) rl jj)
                         (incf i) (decf j) (incf ii) (decf jj)))
                     ;; disallowed weights
                     (when rat
                       (let ((i first) (j last) (wmin +bigd+) (wmax +smad+))
                         (declare (type fixnum i j) (type double-float wmin wmax))
                         (do () ((not (> (- j i) 0)))
                           (let ((wi (aref rl (- i off) 3))
                                 (wj (aref rl (- j off) 3)))
                             (declare (type double-float wi wj))
                             (when (< wi wmin) (setf wmin wi))
                             (when (< wj wmin) (setf wmin wj))
                             (when (> wi wmax) (setf wmax wi))
                             (when (> wj wmax) (setf wmax wj))
                             (incf i) (decf j)))
                         (when (or (< wmin (the double-float *wmin*))
                                   (> wmax (the double-float *wmax*)))
                           (setf wfl nil)
                           (return))))
                     (when (oddp (+ p s))
                       (%combine-2-c-pt-rows lam rl (1+ jj) oml rl (1- ii) rl (1+ jj)))))
                 (cond
                   ((not wfl)
                    ;; a weight out of bounds: the knot stays
                    (dotimes (ll (1+ kk))
                      (setf (aref br ll r) +crv-adv--norem+)))
                   (t
                    ;; save the control points
                    (dotimes (ll (1+ kk))
                      (let ((pw (curve-cpoints (svref cur ll)))
                            (rl (svref rw ll))
                            (i first) (j last))
                        (declare (type cpoint-array pw rl) (type fixnum i j))
                        (fast
                          (do () ((not (> (- j i) 0)))
                            (%copy-c-pt-row rl (- i off) pw i)
                            (%copy-c-pt-row rl (- j off) pw j)
                            (incf i) (decf j)))))
                    ;; shift down some parameters
                    (when (= s 1)
                      (dotimes (ll (1+ kk))
                        (setf (aref er ll (1- r)) (nl-max (aref er ll (1- r)) (aref er ll r)))))
                    (when (> s 1)
                      (setf (aref sr (1- r)) (1- (aref sr r))))
                    (do ((j (1+ r) (1+ j))) ((> j m))
                      (declare (type fixnum j))
                      (setf (aref sr (1- j)) (aref sr j))
                      (dotimes (ll (1+ kk))
                        (setf (aref br ll (1- j)) (aref br ll j)
                              (aref er ll (1- j)) (aref er ll j))))
                    ;; shift down knots and control points
                    (dotimes (ll (1+ kk))
                      (let ((pw (curve-cpoints (svref cur ll)))
                            (ul (curve-knots (svref cur ll))))
                        (declare (type cpoint-array pw) (type real-vector ul))
                        (fast
                          (do ((j (1+ r) (1+ j))) ((> j m))
                            (declare (type fixnum j))
                            (setf (aref ul (1- j)) (aref ul j)))
                          (do ((j (1+ fout) (1+ j))) ((> j n))
                            (declare (type fixnum j))
                            (%copy-c-pt-row pw j pw (1- j))))))
                    (decf n)
                    (decf m)
                    ;; no more internal knots: finished
                    (when (= n p)
                      (return))
                    ;; update error bounds
                    (let ((k (nl-max (- r p) (1+ p)))
                          (l (nl-min n (- (+ r p) s))))
                      (declare (type fixnum k l))
                      (do ((j k (1+ j))) ((> j l))
                        (declare (type fixnum j))
                        (when (/= (aref u j) (aref u (1+ j)))
                          (dotimes (ll (1+ kk))
                            (when (/= (aref br ll j) +crv-adv--norem+)
                              (setf (aref br ll j)
                                    (crv-remove-knot-max-err (svref cur ll) j (aref sr j))
                                    pending nil))))))))))
              (t
               ;; knot not removable
               (dotimes (ll (1+ kk))
                 (setf (aref br ll r) +crv-adv--norem+)))))))
      ;; compact the output curves: they keep the knots and control points
      ;; left, as N_CrvSetSizeIndices records them in C
      (when (< n ns)
        (dotimes (ll (1+ kk))
          (let* ((c (svref cur ll))
                 (pw (curve-cpoints c)))
            (declare (type curve c) (type cpoint-array pw))
            (setf (curve-knots c) (subseq (curve-knots c) 0 (1+ m))
                  (curve-cpoints c) (let ((q (make-cpoint-array (1+ n))))
                                      (dotimes (i (1+ n) q)
                                        (%copy-c-pt-row pw i q i))))
            (crv-compress c)))
        (setf pending nil))
      (when pending
        (error pending))
      cur)))

(defun srf-reparam-mult-knots (curp tol)
  "A fresh copy of curve CURP reparametrized roughly by arc length: each knot
of multiplicity at least the degree (the end knots included) is moved to the
arc length (to relative tolerance TOL) from the start to it, the knots
between two such knots are rescaled linearly, and the domain becomes [0,
length].  The control points are copied unchanged.  The C writes into a
caller's curve, or into CURP itself; Batten returns the new curve.  Signals
what crv-arc-length signals.
NLib: N_SrfReparamMultKnots."
  (declare (type curve curp) (type double-float tol))
  (let* ((n (curve-n curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (pw (curve-cpoints curp))
         (uq (make-real-vector (1+ m)))
         (qw (make-cpoint-array (1+ n)))
         ;; the C measures a copy of curP
         (cura (crv-copy curp))
         (len (crv-arc-length cura (aref up 0) (aref up m) tol +relative+))
         (ua (aref up 0))
         (ul 0d0) (ur 0d0)
         (jj p)
         (kk (1+ p))
         (oldmlt -1))
    (declare (type fixnum n p m jj kk oldmlt) (type real-vector up uq)
             (type cpoint-array pw qw) (type double-float len ua ul ur))
    (dotimes (ii (1+ p))
      (setf (aref uq ii) 0d0))
    ;; at each knot of full multiplicity (>= the degree) set the parameter to
    ;; the arc length there, and rescale the knots since the previous one
    (do () ((not (< kk m)))
      (let ((ii kk))
        (declare (type fixnum ii))
        ;; multiplicity of this knot
        (do () ((not (and (< kk m) (= (aref up kk) (aref up (1+ kk))))))
          (incf kk))
        (let ((mlt (1+ (- kk ii))))
          (declare (type fixnum mlt))
          (when (>= mlt p)
            (let* ((ub (aref up kk))
                   (d (crv-arc-length cura ua ub tol +relative+))
                   (fact 0d0))
              (declare (type double-float ub d fact))
              (setf ur (+ ur d))
              (when (= kk m)
                (setf ur len))
              (setf fact (/ (- ur ul) (- ub ua)))
              (do ((i (1+ jj) (1+ i))) ((> i (- kk mlt)))
                (declare (type fixnum i))
                (setf (aref uq i) (+ ul (* fact (- (aref up i) ua)))))
              (when (> oldmlt 0)
                (do ((i (1+ (- jj oldmlt)) (1+ i))) ((> i jj))
                  (declare (type fixnum i))
                  (setf (aref uq i) ul)))
              (setf oldmlt mlt
                    ua ub
                    ul ur
                    jj kk)))))
      (incf kk))
    (dotimes (ii (1+ p))
      (setf (aref uq (+ n ii 1)) len))
    ;; the control points
    (dotimes (ii (1+ n))
      (%copy-c-pt-row pw ii qw ii))
    (%make-curve p uq qw)))

;;; ------------------------------------------------------------- L7-g1

(defconstant +crv-adv--kn-tol+ 1d-6
  "KN_TOL in NL_CrvAdv.cpp: the smallest new knot increment, relative to the domain.")

(defconstant +crv-adv--eps+ 1d-3
  "eps in NL_CrvAdv.cpp (NLIB_UNUSED): N_CrvRemoveAllKnotsArcLen's arc length tolerance.")

(defun %crv-adv-normalize (v)
  "N_VectorNormalize( V, &V, &mag ) with its error flag kept: (values mag
failed)."
  (declare (type point v))
  (handler-case (values (vector-normalize v v) nil)
    (nurbs-error () (values 0d0 t))))

(defun %crv-adv-end-tangent (pw a b c q v)
  "ST_curxcdr's end tangent: Q the Euclidean point of row A of PW, V the unit
vector from row B's point to Q, or from row C's point when that one is zero
length.  Returns the magnitude.  Signals :geo-err when both are zero."
  (declare (type cpoint-array pw) (type fixnum a b c) (type point q v))
  (let ((cp (make-cpoint)) (r (make-point)))
    (declare (dynamic-extent cp r))
    (c-pt-to-pt-euclid (%row-c-pt pw a cp) q)
    (c-pt-to-pt-euclid (%row-c-pt pw b cp) r)
    (vector-dir r q v)
    (multiple-value-bind (mag failed) (%crv-adv-normalize v)
      (declare (type double-float mag))
      (when failed
        (c-pt-to-pt-euclid (%row-c-pt pw c cp) r)
        (vector-dir r q v)
        (multiple-value-setq (mag failed) (%crv-adv-normalize v)))
      (when failed
        (nurbs-error :geo-err 'crv-adv--st-curxcdr))
      mag)))

(defun crv-adv--st-curxcdr (curp dist flg end cont depth)
  "Curve CURP extended at END (+start+ or +end+) by about the arc length DIST
(+absolute+ FLG) or DIST times its arc length (+relative+), as a fresh curve.
CONT +g-1+ adds a straight segment (degree elevated), +g-2+ a mirror image of
the end segment of that length (a longer one in steps, recursively), +cmax+
extrapolates the end polynomial piece (its derivatives at the end), fits its
length by secants, and removes the old end knot.  DEPTH counts the recursion:
past 50 signals :quit.  Signals :inp-err for a curve of zero length (when it is
measured), :geo-err for a zero end tangent, :wei-err for extrapolated weights
outside [*wmin*, *wmax*], and what the arc lengths, knot insertion and removal
signal.  Another CONT leaves the C's output uninitialized; here it is zeros.
NLib: ST_curxcdr."
  (declare (type curve curp) (type double-float dist) (type fixnum flg end cont depth))
  (when (> depth 50)
    (nurbs-error :quit 'crv-adv--st-curxcdr))
  (let* ((n (curve-n curp))
         (pw (curve-cpoints curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (ldst dist)
         (d1 0d0)
         (uu 0d0)
         (span 0) (mult 0) (mq 0)
         (ktol (* +crv-adv--kn-tol+ (- (aref up m) (aref up 0)))))
    (declare (type fixnum n p m span mult mq) (type cpoint-array pw) (type real-vector up)
             (type double-float ldst d1 uu ktol))
    ;; the extension distance
    (when (or (= flg +relative+) (= cont +g-2+))
      (setf d1 (crv-arc-length curp (aref up 0) (aref up m) 0.01d0 +relative+))
      (when (= flg +relative+)
        (setf ldst (* ldst d1)))
      (when (and (<= d1 (the double-float *mtol*)) (<= d1 (the double-float *ptol*)))
        (nurbs-error :inp-err 'crv-adv--st-curxcdr))
      ;; a G2 extension longer than the curve goes in two steps
      (when (and (= cont +g-2+) (> ldst (* 0.9313d0 d1)))
        (let ((dd (* 0.4671d0 d1)))
          (declare (type double-float dd))
          (return-from crv-adv--st-curxcdr
            (crv-adv--st-curxcdr
             (crv-adv--st-curxcdr curp (- ldst dd) +absolute+ end cont (1+ depth))
             dd +absolute+ end cont (1+ depth))))))
    ;; the parameter of the end segment to mirror (G2), the end knot to
    ;; extrapolate (CMAX), and the output's highest knot index
    (cond
      ((= cont +g-2+)
       (let ((du (* (- (aref up m) (aref up 0)) (/ ldst d1)))
             (ul 0d0) (ur 0d0) (dl 0d0) (dr 0d0))
         (declare (type double-float du ul ur dl dr))
         (setf uu (- (aref up 0) 10d0))
         (if (= end +start+)
             (progn
               (setf ur (aref up 0) dr 0d0)
               (do () (nil)
                 (setf ul ur dl dr ur (+ ur du))
                 (when (> ur (aref up m))
                   (setf ur (aref up m)))
                 (setf dr (crv-arc-length curp (aref up 0) ur 0.005d0 +relative+))
                 (when (or (< (abs (- ldst dr)) (* ldst 0.05d0)) (= ur (aref up m)))
                   (setf uu ur)
                   (return))
                 (when (> dr ldst)
                   (return)))
               (when (< uu (aref up 0))
                 (dotimes (ii 4)
                   (setf uu (+ ul (* (- ldst dl) (/ (- ur ul) (- dr dl)))))
                   (let ((dd (crv-arc-length curp (aref up 0) uu 0.005d0 +relative+)))
                     (declare (type double-float dd))
                     (when (< (abs (- ldst dd)) (* ldst 0.05d0))
                       (return))
                     (if (> dd ldst)
                         (setf ur uu dr dd)
                         (setf ul uu dl dd)))))
               (multiple-value-setq (span mult) (basis-find-span-and-mult up p uu +right+))
               (when (and (/= uu (aref up (1+ span)))
                          (< (- uu (aref up span)) ktol)
                          (< (+ (aref up span) ktol) (aref up (1+ span))))
                 (setf uu (+ (aref up span) ktol)))
               (setf mq (+ m span)))
             (progn
               (setf ul (aref up m) dl 0d0)
               (do () (nil)
                 (setf ur ul dr dl ul (- ul du))
                 (when (< ul (aref up 0))
                   (setf ul (aref up 0)))
                 (setf dl (crv-arc-length curp ul (aref up m) 0.005d0 +relative+))
                 (when (or (< (abs (- ldst dl)) (* ldst 0.05d0)) (= ul (aref up 0)))
                   (setf uu ul)
                   (return))
                 (when (> dl ldst)
                   (return)))
               (when (< uu (aref up 0))
                 (dotimes (ii 4)
                   (setf uu (- ur (* (- ldst dr) (/ (- ur ul) (- dl dr)))))
                   (let ((dd (crv-arc-length curp uu (aref up m) 0.005d0 +relative+)))
                     (declare (type double-float dd))
                     (when (< (abs (- ldst dd)) (* ldst 0.05d0))
                       (return))
                     (if (> dd ldst)
                         (setf ul uu dl dd)
                         (setf ur uu dr dd)))))
               (multiple-value-setq (span mult) (basis-find-span-and-mult up p uu +left+))
               (when (and (/= uu (aref up span))
                          (< (- (aref up (1+ span)) uu) ktol)
                          (> (- (aref up (1+ span)) ktol) (aref up span)))
                 (setf uu (- (aref up (1+ span)) ktol)))
               (setf mq (- (* 2 m) span 1))))))
      ((= cont +cmax+)
       (setf mult 1)
       (if (= end +start+)
           (progn
             (setf uu (aref up (1+ p)))
             (do () ((not (< (+ p mult 1) m)))
               (if (= uu (aref up (+ p mult 1)))
                   (incf mult)
                   (return))))
           (progn
             (setf uu (aref up (- m p 1)))
             (do () ((not (> (- m p mult 1) 0)))
               (if (= uu (aref up (- m p mult 1)))
                   (incf mult)
                   (return)))))
       (setf mq (- (+ m (* 2 p)) mult)))
      (t
       (setf mq (+ m p))))
    ;; the new control points and knots
    (let* ((qw (make-cpoint-array (- mq p)))
           (uq (make-real-vector (1+ mq)))
           (curq (%make-curve p uq qw))
           (kk 0) (kq 0) (kp 0)
           (p1 (make-point)) (v (make-point)) (cp (make-cpoint)) (cq (make-cpoint)))
      (declare (type cpoint-array qw) (type real-vector uq) (type curve curq)
               (type fixnum kk kq kp) (type point p1 v) (type cpoint cp cq))
      (flet ((load-rest ()
               ;; the remaining control points and knots
               (dotimes (ii (1+ n))
                 (%copy-c-pt-row pw ii qw kk)
                 (incf kk))
               (do ((ii 1 (1+ ii))) ((> ii m))
                 (declare (type fixnum ii))
                 (setf (aref uq kq) (aref up kp))
                 (incf kq)
                 (incf kp))))
        (cond
          ((= cont +g-1+)
           (let ((dd (float p 1d0)) (mag 0d0) (du 0d0))
             (declare (type double-float dd mag du))
             (if (= end +start+)
                 (progn
                   ;; the new start point on the end tangent
                   (setf mag (%crv-adv-end-tangent pw 0 1 2 p1 v))
                   (vector-blend-pt ldst v p1)
                   (multiple-value-bind (wx wy wz) (pt-to-xyz p1)
                     (declare (type double-float wx wy wz))
                     (let ((w (aref pw 1 3)))
                       (declare (type double-float w))
                       (when (/= w +no-w+)
                         (setf wx (* wx w) wy (* wy w) wz (* wz w)))
                       (when (= (aref pw 0 2) +no-z+)
                         (setf wz +no-z+))
                       (%c-pt-row (c-pt-from-wx-wy-wz wx wy wz w cp) qw 0)))
                   ;; p - 1 new interior control points (degree elevation)
                   (do ((ii 1 (1+ ii))) ((>= ii p))
                     (declare (type fixnum ii))
                     (let ((d (/ (float ii 1d0) dd)))
                       (declare (type double-float d))
                       (%combine-2-c-pt-rows (- 1d0 d) qw 0 d pw 0 qw ii)))
                   ;; the new start knot
                   (setf du (/ (* ldst (- (aref up (1+ p)) (aref up 0))) (* dd mag)))
                   (when (< du ktol)
                     (setf du ktol))
                   (setf uu (- (aref up 0) du))
                   (dotimes (ii (1+ p))
                     (setf (aref uq ii) uu))
                   (setf kk p kq (1+ p) kp 1))
                 (progn
                   ;; the new end point on the end tangent
                   (setf mag (%crv-adv-end-tangent pw n (1- n) (- n 2) p1 v))
                   (vector-blend-pt ldst v p1)
                   (multiple-value-bind (wx wy wz) (pt-to-xyz p1)
                     (declare (type double-float wx wy wz))
                     (let ((w (aref pw (1- n) 3)))
                       (declare (type double-float w))
                       (when (/= w +no-w+)
                         (setf wx (* wx w) wy (* wy w) wz (* wz w)))
                       (when (= (aref pw 0 2) +no-z+)
                         (setf wz +no-z+))
                       (%c-pt-row (c-pt-from-wx-wy-wz wx wy wz w cp) qw (+ n p))))
                   (do ((ii 1 (1+ ii))) ((>= ii p))
                     (declare (type fixnum ii))
                     (let ((d (/ (float ii 1d0) dd)))
                       (declare (type double-float d))
                       (%combine-2-c-pt-rows (- 1d0 d) pw n d qw (+ n p) qw (+ n ii))))
                   ;; the new end knot
                   (setf du (/ (* ldst (- (aref up m) (aref up (- m p 1)))) (* dd mag)))
                   (when (< du ktol)
                     (setf du ktol))
                   (setf uu (+ (aref up m) du))
                   (dotimes (ii (1+ p))
                     (setf (aref uq (+ m ii)) uu))
                   (setf kk 0 kq 0 kp 0))))
           (load-rest))
          ((= cont +g-2+)
           ;; make uu a knot of multiplicity p
           (let* ((curr (if (< mult p)
                            (prog1 (crv-insert-knot curp uu (- p mult))
                              (when (= end +end+)
                                (incf span (- p mult))))
                            curp))
                  (nr (curve-n curr))
                  (rw (curve-cpoints curr))
                  (mr (curve-m curr))
                  (urk (curve-knots curr))
                  (pln (%make-plane (make-point) (make-point))))
             (declare (type curve curr) (type fixnum nr mr) (type cpoint-array rw)
                      (type real-vector urk) (type plane pln))
             (if (= end +start+)
                 (progn
                   ;; reflect through the plane normal to the start tangent
                   (%crv-adv-end-tangent pw 0 1 2 p1 v)
                   (create-plane-pt-normal pln p1 v)
                   (setf kk 0)
                   (do ((ii span (1- ii))) ((<= ii 0))
                     (declare (type fixnum ii))
                     (%c-pt-row (cpt-reflect (%row-c-pt rw ii cp) pln cq) qw kk)
                     (incf kk))
                   (let ((dd (* 2d0 (aref up 0))))
                     (declare (type double-float dd))
                     (dotimes (ii (1+ p))
                       (setf (aref uq ii) (- dd uu)))
                     (setf kq (1+ p))
                     (do ((ii span (1- ii))) ((<= ii p))
                       (declare (type fixnum ii))
                       (setf (aref uq kq) (- dd (aref up ii)))
                       (incf kq)))
                   (setf kp 1))
                 (progn
                   ;; reflect through the plane normal to the end tangent
                   (%crv-adv-end-tangent pw n (1- n) (- n 2) p1 v)
                   (create-plane-pt-normal pln p1 v)
                   (setf kk (1+ n))
                   (do ((ii (1- nr) (1- ii))) ((< ii (- span p)))
                     (declare (type fixnum ii))
                     (%c-pt-row (cpt-reflect (%row-c-pt rw ii cp) pln cq) qw kk)
                     (incf kk))
                   (setf kk 0)
                   (let ((dd (* 2d0 (aref up m))))
                     (declare (type double-float dd))
                     (dotimes (ii (1+ p))
                       (setf (aref uq (+ (- mq p) ii)) (- dd uu)))
                     (setf kq (- mq p 1))
                     (do ((ii (1+ span) (1+ ii))) ((>= ii (- mr p)))
                       (declare (type fixnum ii))
                       (setf (aref uq kq) (- dd (aref urk ii)))
                       (decf kq)))
                   (setf kq 0 kp 0))))
           (load-rest))
          ((= cont +cmax+)
           (let* ((curr (if (< mult p) (crv-insert-knot curp uu (- p mult)) curp))
                  (nr (curve-n curr))
                  (rw (curve-cpoints curr))
                  (mr (curve-m curr))
                  (urk (curve-knots curr))
                  (rat (is-crv-rat curr))
                  (dw (make-cpoint-array (1+ p)))
                  (fp (float p 1d0))
                  (u1 0d0) (u2 0d0) (u3 0d0) (d2 0d0) (dl 0d0) (dr 0d0))
             (declare (type curve curr) (type fixnum nr mr) (type cpoint-array rw dw)
                      (type real-vector urk) (type double-float fp u1 u2 u3 d2 dl dr))
             (flet ((scale-row (alpha arr i)
                      (declare (type double-float alpha) (type cpoint-array arr) (type fixnum i))
                      (%c-pt-row (scale-c-pt alpha (%row-c-pt arr i cp) cp) arr i))
                    (binom (kp i k)
                      (declare (type fixnum kp i k))
                      (* (float kp 1d0) (float (aref +pascal-tri+ i k) 1d0))))
               (if (= end +start+)
                   (progn
                     (dotimes (ii (1+ nr))
                       (%copy-c-pt-row rw ii qw (+ ii p)))
                     (do ((ii 1 (1+ ii))) ((> ii mr))
                       (declare (type fixnum ii))
                       (setf (aref uq (+ ii p)) (aref urk ii)))
                     ;; derivatives from the right
                     (setf dr (/ fp (- (aref urk (1+ p)) (aref urk 0))))
                     (do ((ii 1 (1+ ii))) ((> ii p))
                       (declare (type fixnum ii))
                       (%copy-c-pt-row rw ii dw ii)
                       (let ((kp (if (evenp ii) 1 -1)))
                         (declare (type fixnum kp))
                         (dotimes (k ii)
                           (%vector-blend-c-pt-row (binom kp ii k) rw k dw ii)
                           (setf kp (- kp))))
                       (scale-row dr dw ii)
                       (setf dr (/ (* dr (- fp (float ii 1d0))) (- (aref urk (1+ p)) (aref urk 0)))))
                     ;; the new start segment, until its arc length is right
                     (setf u1 0d0 d1 0d0 u2 (- (aref urk (1+ p)) (aref urk 0)))
                     (dotimes (ii 5)
                       (block iteration
                         (dotimes (jj (1+ p))
                           (setf (aref uq jj) (- (aref urk 0) u2)))
                         (setf dl (/ u2 fp))
                         (do ((jj 1 (1+ jj))) ((> jj p))
                           (declare (type fixnum jj))
                           (let ((kp (if (evenp jj) -1 1)))
                             (declare (type fixnum kp))
                             (%copy-c-pt-row dw jj qw (- p jj))
                             (scale-row (* (float (- kp) 1d0) dl) qw (- p jj))
                             (when (< jj p)
                               (setf dl (/ (* dl u2) (- fp (float jj 1d0)))))
                             (dotimes (k jj)
                               (%vector-blend-c-pt-row (binom kp jj k) qw (- p k) qw (- p jj))
                               (setf kp (- kp)))))
                         ;; no unacceptable weights
                         (when rat
                           (let ((jj 0))
                             (declare (type fixnum jj))
                             (do () ((not (< jj p)))
                               (let ((w (aref qw jj 3)))
                                 (declare (type double-float w))
                                 (when (or (< w (the double-float *wmin*)) (> w (the double-float *wmax*)))
                                   (return)))
                               (incf jj))
                             (when (< jj p)
                               (when (or (= ii 4) (<= u2 u1))
                                 (nurbs-error :wei-err 'crv-adv--st-curxcdr))
                               (setf u2 (* 0.5d0 (+ u1 u2)))
                               (return-from iteration))))
                         ;; the arc length of the new segment
                         (setf d2 (crv-arc-length curq (aref uq 0) (aref uq (1+ p)) 0.005d0 +relative+))
                         (when (< (abs (- ldst d2)) (* ldst 0.03d0))
                           (return))
                         (when (float-op-is-bad (- u2 u1) (- d2 d1) +division+)
                           (return))
                         (setf u3 u2
                               u2 (+ u1 (* (- ldst d1) (/ (- u2 u1) (- d2 d1))))
                               u1 u3
                               d1 d2)))
                     ;; remove the original start knot
                     (setf curq (nth-value 1 (crv-remove-knot curq (aref up 0) p *mtol*))))
                   (progn
                     (dotimes (ii (1+ nr))
                       (%copy-c-pt-row rw ii qw ii))
                     (dotimes (ii mr)
                       (setf (aref uq ii) (aref urk ii)))
                     ;; derivatives from the left
                     (setf dl (/ (float (- p) 1d0) (- (aref urk mr) (aref urk (- mr p 1)))))
                     (do ((ii 1 (1+ ii))) ((> ii p))
                       (declare (type fixnum ii))
                       (%copy-c-pt-row rw (- nr ii) dw ii)
                       (let ((kp (if (evenp ii) 1 -1)))
                         (declare (type fixnum kp))
                         (dotimes (k ii)
                           (%vector-blend-c-pt-row (binom kp ii k) rw (- nr k) dw ii)
                           (setf kp (- kp))))
                       (scale-row dl dw ii)
                       (setf dl (/ (* (- dl) (- fp (float ii 1d0)))
                                   (- (aref urk mr) (aref urk (- mr p 1))))))
                     ;; the new end segment, until its arc length is right
                     (setf u1 0d0 d1 0d0 u2 (- (aref urk mr) (aref urk (- mr p 1))))
                     (dotimes (ii 5)
                       (block iteration
                         (dotimes (jj (1+ p))
                           (setf (aref uq (+ mr jj)) (+ (aref urk mr) u2)))
                         (setf dr (/ u2 fp))
                         (do ((jj 1 (1+ jj))) ((> jj p))
                           (declare (type fixnum jj))
                           (let ((kp (if (evenp jj) -1 1)))
                             (declare (type fixnum kp))
                             (%copy-c-pt-row dw jj qw (+ nr jj))
                             (scale-row dr qw (+ nr jj))
                             (when (< jj p)
                               (setf dr (/ (* dr u2) (- fp (float jj 1d0)))))
                             (dotimes (k jj)
                               (%vector-blend-c-pt-row (binom kp jj k) qw (+ nr k) qw (+ nr jj))
                               (setf kp (- kp)))))
                         ;; no unacceptable weights
                         (when rat
                           (let ((jj 1))
                             (declare (type fixnum jj))
                             (do () ((not (<= jj p)))
                               (let ((w (aref qw (+ nr jj) 3)))
                                 (declare (type double-float w))
                                 (when (or (< w (the double-float *wmin*)) (> w (the double-float *wmax*)))
                                   (return)))
                               (incf jj))
                             (when (<= jj p)
                               (when (or (= ii 4) (<= u2 u1))
                                 (nurbs-error :wei-err 'crv-adv--st-curxcdr))
                               (setf u2 (* 0.5d0 (+ u1 u2)))
                               (return-from iteration))))
                         ;; the arc length of the new segment
                         (setf d2 (crv-arc-length curq (aref uq (- mq p 1)) (aref uq mq) 0.005d0 +relative+))
                         (when (< (abs (- ldst d2)) (* ldst 0.03d0))
                           (return))
                         (when (float-op-is-bad (- u2 u1) (- d2 d1) +division+)
                           (return))
                         (setf u3 u2
                               u2 (+ u1 (* (- ldst d1) (/ (- u2 u1) (- d2 d1))))
                               u1 u3
                               d1 d2)))
                     ;; remove the original end knot
                     (setf curq (nth-value 1 (crv-remove-knot curq (aref up m) p *mtol*)))))
               ;; remove the inserted instances of uu
               (when (< mult p)
                 (setf curq (nth-value 1 (crv-remove-knot curq uu (- p mult) *mtol*)))))))))
      curq)))

(defun crv-extend-by-dist (curp dist flg end cont)
  "Curve CURP extended at END (+start+ or +end+) by about the arc length DIST
(FLG +absolute+) or DIST times CURP's arc length (+relative+), as a fresh
curve, with continuity CONT at the join: +g-1+ (a straight segment), +g-2+ (a
reflected end segment) or +cmax+ (the end piece extrapolated, no knot at the
join).  The C writes into a caller's curve, or CURP itself; Batten returns the
new curve.  Signals what crv-adv--st-curxcdr signals.
NLib: N_CrvExtendByDist."
  (declare (type curve curp) (type double-float dist) (type fixnum flg end cont))
  (crv-adv--st-curxcdr curp dist flg end cont 0))

(defun crv-remove-all-knots-arc-len (curp tol)
  "Curve CURP with every knot removable within TOL removed, after
reparametrizing it roughly by arc length (srf-reparam-mult-knots at relative
tolerance 1e-3), as a fresh curve.  The C writes into a caller's curve, or CURP
itself; Batten returns the new curve.  Signals what the reparametrization and
the removal signal.
NLib: N_CrvRemoveAllKnotsArcLen."
  (declare (type curve curp) (type double-float tol))
  (crv-remove-knots (srf-reparam-mult-knots curp +crv-adv--eps+) tol))
