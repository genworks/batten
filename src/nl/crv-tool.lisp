;;;; SPDX-FileCopyrightText: Copyright (c) 2020-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Ported from NLib's NL_CrvTool.cpp (usd-brep source/SMLib/src).
;;;; Modified: translated to Common Lisp.  Output curves the C sizes or
;;;; allocates on a stack are returned as fresh curves.

(in-package :batten)

;;; ------------------------------------------------- row helpers (not NLib)

;;; N_CopyCPt and N_Combine2CPts on elements of NL_CPOINT arrays, done on rows
;;; of cpoint-arrays with the same arithmetic, every input coordinate read
;;; before any output is written, as C's by-value arguments are.

(declaim (inline %crv-tool-copy-c-row %crv-tool-combine-2-c-rows %crv-tool-row-to-cpoint))

(defun %crv-tool-copy-c-row (a i b j)
  "Row J of cpoint-array B := row I of A, as copy-c-pt."
  (declare (type cpoint-array a b) (type index i j))
  (let ((x (aref a i 0)) (y (aref a i 1)) (z (aref a i 2)) (w (aref a i 3)))
    (setf (aref b j 0) x (aref b j 1) y (aref b j 2) z (aref b j 3) w)))

(defun %crv-tool-combine-2-c-rows (alpha a ia beta b ib c ic)
  "Row IC of C := ALPHA * row IA of A + BETA * row IB of B, as combine-2-c-pts."
  (declare (type double-float alpha beta) (type cpoint-array a b c) (type index ia ib ic))
  (let ((ax (aref a ia 0)) (ay (aref a ia 1)) (az (aref a ia 2)) (aw (aref a ia 3))
        (bx (aref b ib 0)) (by (aref b ib 1)) (bz (aref b ib 2)) (bw (aref b ib 3)))
    (declare (type double-float ax ay az aw bx by bz bw))
    (setf (aref c ic 0) (+ (* alpha ax) (* beta bx))
          (aref c ic 1) (+ (* alpha ay) (* beta by))
          (aref c ic 2) (cond ((and (/= az +no-z+) (/= bz +no-z+)) (+ (* alpha az) (* beta bz)))
                              ((/= az +no-z+) (* alpha az))
                              ((/= bz +no-z+) (* beta bz))
                              (t +no-z+))
          (aref c ic 3) (if (and (/= aw +no-w+) (/= bw +no-w+)) (+ (* alpha aw) (* beta bw)) +no-w+))))

(defun %crv-tool-row-to-cpoint (a i cp)
  "Copy row I of cpoint-array A into cpoint CP and return CP."
  (declare (type cpoint-array a) (type index i) (type cpoint cp))
  (set-cpoint cp (aref a i 0) (aref a i 1) (aref a i 2) (aref a i 3)))

;;; ------------------------------------------------------- knot insertion

(defun crv-insert-knot (cur-p u r)
  "Curve CUR-P with the knot U, strictly inside the domain, inserted R times,
as a fresh curve.  Signals :par-err when U is on or outside an end of the
domain, when R is negative, or when U's multiplicity plus R exceeds the
degree.  NLib: N_CrvInsertKnot (the C's in-place case, curP = curQ, computes
the same curve)."
  (declare (type curve cur-p) (type double-float u) (type fixnum r))
  (let* ((n (curve-n cur-p))
         (pw (curve-cpoints cur-p))
         (p (curve-degree cur-p))
         (m (curve-m cur-p))
         (up (curve-knots cur-p)))
    (declare (type fixnum n p m) (type cpoint-array pw) (type real-vector up))
    ;; check parameter and get knot span
    (knot-vector-is-end-param up u 'crv-insert-knot)
    (multiple-value-bind (spn mlt) (basis-find-span-and-mult up p u +left+)
      (declare (type fixnum spn mlt))
      (when (or (> (+ mlt r) p) (< r 0))
        (nurbs-error :par-err 'crv-insert-knot))
      (let ((k (- spn p))
            (qw (make-cpoint-array (+ n r 1)))
            (uq (make-real-vector (+ m r 1)))
            (rw (make-cpoint-array (1+ p))))
        (declare (type fixnum k) (type cpoint-array qw rw) (type real-vector uq))
        (fast
          ;; auxiliary control points
          (do ((i 0 (1+ i))) ((> i (- p mlt)))
            (declare (type fixnum i))
            (%copy-c-pt-row pw (+ k i) rw i))
          ;; unaltered control points
          (do ((i 0 (1+ i))) ((> i k))
            (declare (type fixnum i))
            (%copy-c-pt-row pw i qw i))
          (do ((i (- spn mlt) (1+ i))) ((> i n))
            (declare (type fixnum i))
            (%copy-c-pt-row pw i qw (+ i r)))
          ;; now insert the knot
          (do ((i 1 (1+ i))) ((> i r))
            (declare (type fixnum i))
            (setf k (+ (- spn p) i))
            (do ((j 0 (1+ j))) ((> j (- p i mlt)))
              (declare (type fixnum j))
              (let* ((alf (/ (- u (aref up (+ k j))) (- (aref up (+ spn j 1)) (aref up (+ k j)))))
                     (oma (- 1d0 alf)))
                (declare (type double-float alf oma))
                (%combine-2-c-pt-rows alf rw (1+ j) oma rw j rw j)))
            (%copy-c-pt-row rw 0 qw k)
            (%copy-c-pt-row rw (- p i mlt) qw (- (+ spn r) i mlt)))
          ;; load the remaining control points
          (do ((i (1+ k) (1+ i))) ((>= i (- spn mlt)))
            (declare (type fixnum i))
            (%copy-c-pt-row rw (- i k) qw i))
          ;; load the knot vector
          (do ((i 0 (1+ i))) ((> i spn))
            (declare (type fixnum i))
            (setf (aref uq i) (aref up i)))
          (do ((i 1 (1+ i))) ((> i r))
            (declare (type fixnum i))
            (setf (aref uq (+ i spn)) u))
          (do ((i (1+ spn) (1+ i))) ((> i m))
            (declare (type fixnum i))
            (setf (aref uq (+ i r)) (aref up i))))
        (%make-curve p uq qw)))))

;;; ------------------------------------------------------- decomposition

(defun %next-bez-cpoints (nqw)
  "NQW, the next Bezier segment's control points, checked: where the C would
write through its still-NULL pointer (an unclamped single-span curve) and
crash, signal instead."
  (if (typep nqw 'cpoint-array)
      nqw
      (nurbs-error :cal-err 'crv-decompose-bez)))

(defun crv-decompose-bez (cur)
  "The Bezier segments of curve CUR, extracted without knot refinement: one
degree-p curve with knots (a ... a b ... b) per non-zero span.  Returns
(values curves k): CURVES a simple-vector of the k+1 fresh segments.
NLib: N_CrvDecomposeBez."
  (declare (type curve cur))
  (let* ((p (curve-degree cur))
         (pw (curve-cpoints cur))
         (up (curve-knots cur))
         (m (curve-m cur))
         (nsp (basis-get-span-count up p))
         (cura (make-array nsp)))
    (declare (type fixnum p m nsp) (type cpoint-array pw) (type real-vector up)
             (type simple-vector cura))
    (dotimes (i nsp)
      (setf (svref cura i) (%make-curve p (make-real-vector (+ p p 2)) (make-cpoint-array (1+ p)))))
    (fast
      (let ((alfs (make-array (1+ p) :element-type 'double-float :initial-element 0d0))
            (omas (make-array (1+ p) :element-type 'double-float :initial-element 0d0))
            (is p)
            (ie (1+ p))
            (iq -1)
            (nqw nil))
        (declare (dynamic-extent alfs omas) (type fixnum is ie iq))
        ;; skip over first knot over-multiplicity (bad data correction)
        (do () ((not (and (< ie m) (= (aref up (1- ie)) (aref up ie)))))
          (incf is)
          (incf ie))
        (let ((qw (curve-cpoints (svref cura 0))))
          (declare (type cpoint-array qw))
          (dotimes (i (1+ p))
            (%copy-c-pt-row pw (- (+ i is) p) qw i)))
        (do () ((not (< ie m)))
          (incf iq)
          (when (>= iq nsp)
            ;; NL_QUIT: the C fails without setting an error code
            (nurbs-error :quit 'crv-decompose-bez))
          (let ((qw (curve-cpoints (svref cura iq)))
                (uq (curve-knots (svref cura iq))))
            (declare (type cpoint-array qw) (type real-vector uq))
            (when (< iq (1- nsp))
              (setf nqw (curve-cpoints (svref cura (1+ iq)))))
            (let ((i ie))
              (declare (type fixnum i))
              (do () ((not (and (< ie m) (= (aref up ie) (aref up (1+ ie))))))
                (incf ie))
              (let* ((mlt (+ (- ie i) 1))
                     (r (- p mlt)))
                (declare (type fixnum mlt r))
                ;; insert the knot on the rhs only
                (when (< mlt p)
                  (let ((num (- (aref up ie) (aref up is))))
                    (declare (type double-float num))
                    (do ((i p (1- i))) ((<= i mlt))
                      (declare (type fixnum i))
                      (setf (aref alfs (- i mlt 1)) (/ num (- (aref up (+ is i)) (aref up is))))
                      (setf (aref omas (- i mlt 1)) (- 1d0 (aref alfs (- i mlt 1))))))
                  (do ((i 1 (1+ i))) ((> i r))
                    (declare (type fixnum i))
                    (let ((s (+ mlt i))
                          (save (- r i)))
                      (declare (type fixnum s save))
                      (do ((j p (1- j))) ((< j s))
                        (declare (type fixnum j))
                        (%combine-2-c-pt-rows (aref alfs (- j s)) qw j
                                              (aref omas (- j s)) qw (1- j)
                                              qw j))
                      (when (< ie m)
                        (%copy-c-pt-row qw p (%next-bez-cpoints nqw) save)))))
                ;; the Bezier segment's knot vector
                (dotimes (i (1+ p))
                  (setf (aref uq i) (aref up is)
                        (aref uq (+ i p 1)) (aref up ie)))
                ;; the next segment's remaining control points
                (when (< ie m)
                  (when (< r 0) (setf r 0))
                  (do ((i r (1+ i))) ((> i p))
                    (declare (type fixnum i))
                    (%copy-c-pt-row pw (+ (- ie p) i) (%next-bez-cpoints nqw) i)))
                (setf is ie
                      ie (1+ ie)))))))
      (values cura (1- nsp)))))

(defun crv-decompose-continuity (curp)
  "Split curve CURP at its interior knots of multiplicity at least its degree,
into pieces each at least C1 inside.  Returns (values curves k): CURVES a
simple-vector of the k+1 fresh pieces (one copy of CURP when it has no such
knot, k = 0).  NLib: N_CrvDecomposeContinuity."
  (declare (type curve curp))
  (let* ((n (curve-n curp))
         (pw (curve-cpoints curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (r (1+ (truncate n p)))
         (ks (make-index-vector (1+ r)))
         (ke (make-index-vector (1+ r)))
         (ns (make-index-vector (1+ r)))
         (ne (make-index-vector (1+ r)))
         (ie (1+ p)))
    (declare (type fixnum n p m r ie) (type cpoint-array pw) (type real-vector up)
             (type index-vector ks ke ns ne))
    (setf r -1
          (aref ns 0) 0
          (aref ks 0) ie)
    (do () ((not (< ie m)))
      (let ((i ie))
        (declare (type fixnum i))
        (do () ((not (and (< ie m) (= (aref up ie) (aref up (1+ ie))))))
          (incf ie))
        (when (>= (+ (- ie i) 1) p)
          (incf r)
          (setf (aref ne r) (1- i)
                (aref ke r) i)
          (when (< ie m)
            (setf (aref ns (1+ r)) (aref ne r)
                  (aref ks (1+ r)) (1+ ie))))
        (incf ie)))
    (let ((cura (make-array (1+ r))))
      (dotimes (i (1+ r))
        (let* ((ni (- (aref ne i) (aref ns i)))
               (mi (+ (- (aref ke i) (aref ks i)) (* 2 p) 1))
               (qw (make-cpoint-array (1+ ni)))
               (uq (make-real-vector (1+ mi)))
               (l 0))
          (declare (type fixnum ni mi l))
          (dotimes (j (1+ ni))
            (%crv-tool-copy-c-row pw (+ (aref ns i) j) qw j))
          (dotimes (j (1+ p))
            (setf (aref uq j) (aref up (1- (aref ks i)))))
          (do ((j (1+ p) (1+ j))) ((> j ni))
            (declare (type fixnum j))
            (setf (aref uq j) (aref up (+ (aref ks i) l)))
            (incf l))
          (do ((j (1+ ni) (1+ j))) ((> j mi))
            (declare (type fixnum j))
            (setf (aref uq j) (aref up (aref ke i))))
          (setf (svref cura i) (%make-curve p uq qw))))
      (values cura r))))

;;; ---------------------------------------------------------- refinement

(defun crv-refine (curp knx)
  "Curve CURP with the knots of KNX (a non-decreasing real-vector strictly
inside CURP's domain) inserted, as a fresh curve.  Signals :inp-err when KNX
is empty and :par-err when its first or last knot is not inside the domain.
NLib: N_CrvRefine."
  (declare (type curve curp) (type real-vector knx))
  (let* ((n (curve-n curp))
         (pw (curve-cpoints curp))
         (p (curve-degree curp))
         (m (curve-m curp))
         (up (curve-knots curp))
         (r (1- (length knx)))
         (x knx))
    (declare (type fixnum n p m r) (type cpoint-array pw) (type real-vector up x))
    (when (< r 0)
      (nurbs-error :inp-err 'crv-refine))
    (knot-vector-is-end-param up (aref x 0) 'crv-refine)
    (knot-vector-is-end-param up (aref x r) 'crv-refine)
    (let* ((qw (make-cpoint-array (+ n r 2)))
           (uq (make-real-vector (+ m r 2)))
           (a (basis-find-span up p (aref x 0) +left+))
           (b (1+ (basis-find-span up p (aref x r) +left+))))
      (declare (type cpoint-array qw) (type real-vector uq) (type fixnum a b))
      (fast
        (do ((j 0 (1+ j))) ((> j a))
          (declare (type fixnum j))
          (setf (aref uq j) (aref up j)))
        (do ((j (+ b p) (1+ j))) ((> j m))
          (declare (type fixnum j))
          (setf (aref uq (+ j r 1)) (aref up j)))
        (do ((j 0 (1+ j))) ((> j (- a p)))
          (declare (type fixnum j))
          (%crv-tool-copy-c-row pw j qw j))
        (do ((j (1- b) (1+ j))) ((> j n))
          (declare (type fixnum j))
          (%crv-tool-copy-c-row pw j qw (+ j r 1)))
        (let ((i (+ b p -1))
              (k (+ b p r)))
          (declare (type fixnum i k))
          (do ((j r (1- j))) ((< j 0))
            (declare (type fixnum j))
            (do () ((not (and (<= (aref x j) (aref up i)) (> i a))))
              (%crv-tool-copy-c-row pw (- i p 1) qw (- k p 1))
              (setf (aref uq k) (aref up i))
              (decf k)
              (decf i))
            (%crv-tool-copy-c-row qw (- k p) qw (- k p 1))
            (do ((l 1 (1+ l))) ((> l p))
              (declare (type fixnum l))
              (let ((tt (+ (- k p) l))
                    (alf (- (aref uq (+ k l)) (aref x j))))
                (declare (type fixnum tt) (type double-float alf))
                (if (<= (abs alf) 0d0)
                    (%crv-tool-copy-c-row qw tt qw (1- tt))
                    (let* ((alf (/ alf (- (aref uq (+ k l)) (aref up (+ (- i p) l)))))
                           (oma (- 1d0 alf)))
                      (declare (type double-float alf oma))
                      (%crv-tool-combine-2-c-rows alf qw (1- tt) oma qw tt qw (1- tt))))))
            (setf (aref uq k) (aref x j))
            (decf k))))
      (%make-curve p uq qw))))

;;; ---------------------------------------------------------- segments

(defun crv-extract-crv-seg (cur-p ul ur)
  "The segment of curve CUR-P over [UL, UR], as a fresh clamped curve.
Signals :inp-err when UR <= UL or the segment has no control points.
NLib: N_CrvExtractCrvSeg (the C's in-place case, curP = curQ, computes the
same curve)."
  (declare (type curve cur-p) (type double-float ul ur))
  (let* ((pw (curve-cpoints cur-p))
         (p (curve-degree cur-p))
         (m (curve-m cur-p))
         (up (curve-knots cur-p)))
    (declare (type fixnum p m) (type cpoint-array pw) (type real-vector up))
    ;; check parameters
    (when (<= ur ul)
      (nurbs-error :inp-err 'crv-extract-crv-seg))
    ;; find knot spans and set new indexes
    (multiple-value-bind (spl mll) (basis-find-span-and-mult up p ul +left+)
      (declare (type fixnum spl mll))
      (multiple-value-bind (spr mlr) (basis-find-span-and-mult up p ur +left+)
        (declare (type fixnum spr mlr))
        (let ((is (- spl p)))
          (declare (type fixnum is))
          (when (= ur (aref up (- m p)))
            (setf spr m
                  mlr (1+ p)))
          (let* ((ie (- spr mlr))
                 (nn (- ie is))
                 (mm (+ (- spr spl mlr) (* 2 p) 1)))
            (declare (type fixnum ie nn mm))
            (when (< nn 0)
              (nurbs-error :inp-err 'crv-extract-crv-seg))
            (let ((qw (make-cpoint-array (1+ nn)))
                  (uq (make-real-vector (1+ mm))))
              (declare (type cpoint-array qw) (type real-vector uq))
              (fast
                ;; initial control points
                (do ((i is (1+ i))) ((> i ie))
                  (declare (type fixnum i))
                  (%copy-c-pt-row pw i qw (- i is)))
                ;; insert the left knot
                (let ((ll (- spl p)))
                  (declare (type fixnum ll))
                  (do ((i 1 (1+ i))) ((> i (- p mll)))
                    (declare (type fixnum i))
                    (do ((j 0 (1+ j))) ((> j (- p i mll)))
                      (declare (type fixnum j))
                      (let* ((left (aref up (+ ll i j)))
                             (alf (/ (- ul left) (- (aref up (+ spl j 1)) left)))
                             (oma (- 1d0 alf)))
                        (declare (type double-float left alf oma))
                        (%combine-2-c-pt-rows alf qw (1+ j) oma qw j qw j)))))
                ;; insert the right knot
                (let ((lr (- spr p))
                      (lk (+ (- nn p) mlr)))
                  (declare (type fixnum lr lk))
                  (do ((i 1 (1+ i))) ((> i (- p mlr)))
                    (declare (type fixnum i))
                    (do ((j (- p i mlr) (1- j))) ((< j 0))
                      (declare (type fixnum j))
                      (let ((k (+ lk i j))
                            (left (aref up (+ lr i j))))
                        (declare (type fixnum k) (type double-float left))
                        (when (< left ul)
                          (setf left ul))
                        (let* ((alf (/ (- ur left) (- (aref up (+ spr j 1)) left)))
                               (oma (- 1d0 alf)))
                          (declare (type double-float alf oma))
                          (%combine-2-c-pt-rows alf qw k oma qw (1- k) qw k))))))
                ;; load the knot vector
                (let ((j -1))
                  (declare (type fixnum j))
                  (dotimes (i (1+ p))
                    (setf (aref uq (incf j)) ul))
                  (do ((i (1+ spl) (1+ i))) ((> i (- spr mlr)))
                    (declare (type fixnum i))
                    (setf (aref uq (incf j)) (aref up i)))
                  (dotimes (i (1+ p))
                    (setf (aref uq (incf j)) ur))))
              (%make-curve p uq qw))))))))

;;; ------------------------------------------------------------ knot removal

(defun crv-remove-knot-max-err (cur r s)
  "The maximum error of removing knot U[R], of multiplicity S, once from curve
CUR.  U[R] must be an interior knot with U[R] /= U[R+1] and 0 < S <= p; this
is not checked.  NLib: N_CrvRemoveKnotMaxErr."
  (declare (type curve cur) (type fixnum r s))
  (let* ((pw (curve-cpoints cur))
         (p (curve-degree cur))
         (u (curve-knots cur))
         (rw (make-array (list (1+ (* 2 p)) 4) :element-type 'double-float :initial-element 0d0))
         (ca (make-cpoint)) (cb (make-cpoint)) (cc (make-cpoint))
         (i (- r p))
         (j (- r s))
         (ii 1)
         (jj (+ (- p s) 1)))
    (declare (type cpoint-array pw rw) (type fixnum p i j ii jj) (type real-vector u)
             (dynamic-extent rw ca cb cc))
    (fast
      (%crv-tool-copy-c-row pw (- r p 1) rw 0)
      (%crv-tool-copy-c-row pw (+ (- r s) 1) rw (+ (- p s) 2))
      (do () ((not (> (- j i) 0)))
        (let* ((alf (/ (- (aref u (+ i p 1)) (aref u i)) (- (aref u r) (aref u i))))
               (oma (- 1d0 alf))
               (bet (/ (- (aref u (+ j p 1)) (aref u j)) (- (aref u (+ j p 1)) (aref u r))))
               (omb (- 1d0 bet)))
          (declare (type double-float alf oma bet omb))
          (%crv-tool-combine-2-c-rows alf pw i oma rw (1- ii) rw ii)
          (%crv-tool-combine-2-c-rows bet pw j omb rw (1+ jj) rw jj)
          (incf i)
          (decf j)
          (incf ii)
          (decf jj))))
    (if (< (- j i) 0)
        (dist-cpt-cpt-homo (%crv-tool-row-to-cpoint rw (1- ii) ca)
                           (%crv-tool-row-to-cpoint rw (1+ jj) cb))
        (let* ((alf (/ (- (aref u r) (aref u i)) (- (aref u (+ i p 1)) (aref u i))))
               (oma (- 1d0 alf)))
          (declare (type double-float alf oma))
          (combine-2-c-pts alf (%crv-tool-row-to-cpoint rw (1+ ii) ca)
                           oma (%crv-tool-row-to-cpoint rw (1- ii) cb) cc)
          (dist-cpt-cpt-homo (%crv-tool-row-to-cpoint pw i ca) cc)))))
