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

;;; ------------------------------------------------------------- L18-g0

(defun %tess-surface ()
  "A surface to tessellate: a random one of low degree (rational or not),
now and then a plane (every control point on one plane), a gently curved
one, or a bilinear parallelogram -- a rectangle half the time, whose corner
edges are perpendicular (ST_TessPolygonizeTrimCrv's special case for degree-1
trimming curves).  Now and then closed in u, in v or both (its last row or
column of control points a copy of its first).  Half the time scaled by
0.01 to 10, which changes the subdivision's size estimates against the
work it does."
  (let ((s (%tess-open-surface)))
    (let* ((net (surface-cpoints s))
           (n (1- (array-dimension net 0))) (m (1- (array-dimension net 1))))
      (when (zerop (rnd 2))
        (let ((f (expt 10d0 (- (rnd 3d0) 2d0))))
          (dotimes (i (1+ n))
            (dotimes (j (1+ m))
              (dotimes (c 3) (setf (aref net i j c) (* f (aref net i j c))))))))
      (when (and (> n 1) (zerop (rnd 6)))
        (dotimes (j (1+ m)) (dotimes (c 4) (setf (aref net n j c) (aref net 0 j c)))))
      (when (and (> m 1) (zerop (rnd 8)))
        (dotimes (i (1+ n)) (dotimes (c 4) (setf (aref net i m c) (aref net i 0 c))))))
    s))

(defun %tess-open-surface ()
  "The surfaces of %tess-surface before any closing."
  (case (rnd 8)
    ((0 1 2) (random-surface :max-degree 3 :max-extra 3))
    ((3 4)
     ;; a plane, or nearly one
     (let* ((s (random-surface :max-degree 3 :max-extra 3))
            (net (surface-cpoints s))
            (a (- (rnd 2d0) 1d0)) (b (- (rnd 2d0) 1d0)) (c (- (rnd 20d0) 10d0))
            (bump (if (zerop (rnd 2)) 0d0 (rnd 5d0))))
       (dotimes (i (array-dimension net 0) s)
         (dotimes (j (array-dimension net 1))
           (let ((w (if (= (aref net i j 3) +no-w+) 1d0 (aref net i j 3))))
             (setf (aref net i j 2)
                   (* w (+ (* a (/ (aref net i j 0) w)) (* b (/ (aref net i j 1) w)) c
                           (- (rnd (* 2 bump)) bump)))))))))
    (5
     ;; a gently curved sheet over a grid
     (let* ((s (random-surface :max-degree 3 :max-extra 3 :rational nil))
            (net (surface-cpoints s))
            (n (array-dimension net 0)) (m (array-dimension net 1))
            (sx (+ 20d0 (rnd 100d0))) (sy (+ 20d0 (rnd 100d0))))
       (dotimes (i n s)
         (dotimes (j m)
           (setf (aref net i j 0) (* sx (/ (float i 1d0) (1- n)))
                 (aref net i j 1) (* sy (/ (float j 1d0) (1- m)))
                 (aref net i j 2) (- (rnd 20d0) 10d0))))))
    (t
     ;; a bilinear parallelogram
     (let* ((net (make-cpoint-net 2 2))
            (o (list (- (rnd 100d0) 50d0) (- (rnd 100d0) 50d0) (- (rnd 100d0) 50d0)))
            (phi (rnd (* 2 pi)))
            (lu (+ 5d0 (rnd 100d0))) (lv (+ 5d0 (rnd 100d0)))
            (skew (if (zerop (rnd 2)) 0d0 (- (rnd 1.2d0) 0.6d0)))
            (du (list (* lu (cos phi)) (* lu (sin phi)) 0d0))
            (dv (list (* lv (- (sin (+ phi skew)))) (* lv (cos (+ phi skew))) (if (zerop (rnd 2)) 0d0 (rnd 10d0)))))
       (dotimes (i 2)
         (dotimes (j 2)
           (dotimes (c 3)
             (setf (aref net i j c) (+ (nth c o) (* i (nth c du)) (* j (nth c dv)))))
           (setf (aref net i j 3) +no-w+)))
       (multiple-value-bind (ulo uhi) (batten-random::random-domain)
         (multiple-value-bind (vlo vhi) (batten-random::random-domain)
           (%make-surface 1 1 (batten-random::random-knots 1 2 :lo ulo :hi uhi)
                          (batten-random::random-knots 1 2 :lo vlo :hi vhi) net)))))))

(defun %tess-curve-loop (cx cy rx ry clockwise)
  "A random closed trimming loop (a simple-vector of curves) around (CX, CY)
within radii RX, RY, clockwise or not."
  (case (rnd 6)
    (0 (ellipse-curves cx cy rx ry :clockwise clockwise))
    (1
     ;; a polygon of lines
     (let* ((k (+ 3 (rnd 5))) (pts (make-array k))
            (dir (if clockwise -1d0 1d0)) (start (rnd (* 2 pi))))
       (dotimes (i k)
         (let ((theta (+ start (* dir 2 pi (/ (+ i (rnd 0.5d0)) k)))))
           (setf (svref pts i) (cons (+ cx (* rx (cos theta))) (+ cy (* ry (sin theta)))))))
       (let ((v (make-array k)))
         (dotimes (i k v)
           (let ((a (svref pts i)) (b (svref pts (mod (1+ i) k))))
             (setf (svref v i) (line-curve (car a) (cdr a) (car b) (cdr b)
                                           :degree (1+ (rnd 2)) :rational (zerop (rnd 3)))))))))
    (t
     (let ((pieces (1+ (rnd 4))) (degree (1+ (rnd 3))))
       (loop-curves cx cy rx ry :clockwise clockwise :pieces pieces :degree degree
                                :points (max (1+ degree) (ceiling 10 pieces))
                                :rational (zerop (rnd 3)) :wobble (rnd 0.1d0)
                                :start (rnd (* 2 pi))
                                :interior (if (zerop (rnd 3)) :random :even))))))

(defun %tess-holes (cx cy rx ry)
  "0, 1 or 2 clockwise inner loops well inside the ellipse at (CX, CY) of
radii RX, RY (the region an outer loop around them surely contains): a
simple-vector of loops."
  (case (rnd 3)
    (0 (vector))
    (1 (vector (%tess-curve-loop (+ cx (* rx (- (rnd 0.4d0) 0.2d0))) (+ cy (* ry (- (rnd 0.4d0) 0.2d0)))
                                 (* rx (+ 0.08d0 (rnd 0.12d0))) (* ry (+ 0.08d0 (rnd 0.12d0))) t)))
    (t (let ((dy (* ry (- (rnd 0.2d0) 0.1d0))))
         (vector (%tess-curve-loop (- cx (* 0.35d0 rx)) (+ cy dy)
                                   (* rx (+ 0.06d0 (rnd 0.09d0))) (* ry (+ 0.06d0 (rnd 0.12d0))) t)
                 (%tess-curve-loop (+ cx (* 0.35d0 rx)) (- cy dy)
                                   (* rx (+ 0.06d0 (rnd 0.09d0))) (* ry (+ 0.06d0 (rnd 0.12d0))) t))))))

(defun %tess-box-loop (ul ur vb vt &key (cuts 0))
  "The parameter domain's boundary as an outer loop: its sides as straight
curves of degree 1..3, rational or not, counterclockwise, each now and then
split in two or three, starting at a random one of the pieces; with CUTS
corners (0..4) cut off by a curve (a line or a quadratic arc, bulging in or
out) from a point on the side before the corner to one on the side after it."
  (let* ((corners (vector (cons ul vb) (cons ur vb) (cons ur vt) (cons ul vt)))
         (cut (let ((v (make-array 4 :initial-element nil)) (k 0))
                (do () ((>= k (min cuts 4)) v)
                  (let ((c (rnd 4))) (unless (svref v c) (setf (svref v c) t) (incf k))))))
         (pieces nil))                  ; (:line a b) or (:cut a b corner), in order
    (flet ((along (a b s) (cons (+ (car a) (* s (- (car b) (car a)))) (+ (cdr a) (* s (- (cdr b) (cdr a)))))))
      ;; side k runs from corner k to corner k+1; a cut corner's ends are on
      ;; its sides at fraction 0.1..0.3 of their length from it
      (let ((ends (make-array 4)) (starts (make-array 4)))
        (dotimes (k 4)
          (let ((a (svref corners k)) (b (svref corners (mod (1+ k) 4))))
            (setf (svref starts k) (if (svref cut k) (along a b (+ 0.1d0 (rnd 0.2d0))) a)
                  (svref ends k) (if (svref cut (mod (1+ k) 4)) (along b a (+ 0.1d0 (rnd 0.2d0))) b))))
        (dotimes (k 4)
          (let ((a (svref starts k)) (b (svref ends k)))
            (case (if (zerop (rnd 4)) (+ 2 (rnd 2)) 1)
              (1 (push (list :line a b) pieces))
              (2 (let ((m (along a b (+ 0.2d0 (rnd 0.6d0)))))
                   (push (list :line a m) pieces) (push (list :line m b) pieces)))
              (t (let* ((s (+ 0.2d0 (rnd 0.3d0)))
                        (m1 (along a b s)) (m2 (along a b (+ s 0.1d0 (rnd (- 0.8d0 s))))))
                   (push (list :line a m1) pieces) (push (list :line m1 m2) pieces)
                   (push (list :line m2 b) pieces))))
            (let ((k1 (mod (1+ k) 4)))
              (when (svref cut k1)
                (push (list :cut b (svref starts k1) (svref corners k1)) pieces))))))
      (let* ((pieces (coerce (nreverse pieces) 'simple-vector))
             (np (length pieces)) (r (rnd np))
             (curves (make-array np)))
        (dotimes (i np curves)
          (destructuring-bind (kind a b &optional corner) (svref pieces (mod (+ i r) np))
            (setf (svref curves i)
                  (if (or (eq kind :line) (zerop (rnd 3)))
                      (line-curve (car a) (cdr a) (car b) (cdr b)
                                  :degree (if (zerop (rnd 2)) 1 (1+ (rnd 3)))
                                  :rational (zerop (rnd 4))
                                  :interior (if (zerop (rnd 4)) :random :even))
                      (let* ((mid (along a b 0.5d0))
                             (ctl (along mid corner (- (rnd 1.3d0) 0.5d0)))
                             (cp (make-cpoint-array 3)))
                        (setf (aref cp 0 0) (car a) (aref cp 0 1) (cdr a)
                              (aref cp 1 0) (car ctl) (aref cp 1 1) (cdr ctl)
                              (aref cp 2 0) (car b) (aref cp 2 1) (cdr b))
                        (dotimes (j 3) (setf (aref cp j 2) +no-z+ (aref cp j 3) +no-w+))
                        (%make-curve 2 (make-array 6 :element-type 'double-float
                                                     :initial-contents '(0d0 0d0 0d0 1d0 1d0 1d0))
                                     cp))))))))))

(defun %tess-polygon-loop (pts)
  "The closed polygon through the points PTS (conses (x . y), in order) as a
loop of straight curves, of degree 1 mostly."
  (let* ((k (length pts)) (v (make-array k)))
    (dotimes (i k v)
      (let ((a (nth i pts)) (b (nth (mod (1+ i) k) pts)))
        (setf (svref v i) (line-curve (car a) (cdr a) (car b) (cdr b)
                                      :degree (if (zerop (rnd 4)) 2 1)))))))

(defun %tess-grid-holes (ul ur vb vt fx0 fx1 fy0 fy1)
  "One or two clockwise polygonal holes in the part [FX0, FX1] x [FY0, FY1]
(fractions) of the domain [UL, UR] x [VB, VT], their vertices on a dyadic
grid of the domain, so on or near the subdivision's lines: rectangles and
diamonds, whose vertices and edges touch the subdivision rectangles' sides."
  (let ((g (nth (rnd 3) '(4 8 16))))
    (flet ((at (fx fy) (cons (+ ul (* (- ur ul) fx)) (+ vb (* (- vt vb) fy))))
           (snap (f) (/ (float (round (* f g)) 1d0) g)))
      (let ((holes nil) (count (1+ (rnd 2))))
        (dotimes (k count (coerce holes 'simple-vector))
          (let* ((x0 (+ fx0 (* (- fx1 fx0) (/ (float k 1d0) count))))
                 (x1 (+ fx0 (* (- fx1 fx0) (/ (float (1+ k) 1d0) count))))
                 ;; cell of the hole, inside [x0, x1] x [fy0, fy1] with a margin
                 (a (snap (+ x0 (* 0.25d0 (- x1 x0))))) (b (snap (- x1 (* 0.25d0 (- x1 x0)))))
                 (c (snap (+ fy0 (* 0.25d0 (- fy1 fy0))))) (d (snap (- fy1 (* 0.25d0 (- fy1 fy0))))))
            (when (and (< a b) (< c d) (> a x0) (< b x1) (> c fy0) (< d fy1))
              (let ((mx (* 0.5d0 (+ a b))) (my (* 0.5d0 (+ c d))))
                (push (%tess-polygon-loop
                       (case (rnd 4)
                         (0 (list (at a c) (at a d) (at b d) (at b c)))          ; rectangle, clockwise
                         (1 (list (at mx c) (at a my) (at mx d) (at b my)))      ; diamond, clockwise
                         (t
                          ;; a rectangle off the grid with a notch whose apex
                          ;; is on it: the apex turns on a rectangle's side
                          ;; (MX, MY on the grid; the apex's other coordinate off it)
                          (let* ((a (+ a (* (- mx a) (+ 0.05d0 (rnd 0.3d0)))))
                                 (b (- b (* (- b mx) (+ 0.05d0 (rnd 0.3d0)))))
                                 (c (+ c (* (- my c) (+ 0.05d0 (rnd 0.3d0)))))
                                 (d (- d (* (- d my) (+ 0.05d0 (rnd 0.3d0)))))
                                 (rx (+ a (* (- b a) (+ 0.2d0 (rnd 0.6d0)))))
                                 (ry (+ c (* (- d c) (+ 0.2d0 (rnd 0.6d0))))))
                            (case (rnd 4)
                              (0 (list (at a c) (at a d) (at b d) (at b c) (at rx my)))
                              (1 (list (at a c) (at a d) (at rx my) (at b d) (at b c)))
                              (2 (list (at a c) (at mx ry) (at a d) (at b d) (at b c)))
                              (t (list (at a c) (at a d) (at b d) (at mx ry) (at b c))))))))
                      holes)))))))))

(defun %tess-loops (ul ur vb vt)
  "Trimming loops for a surface on [UL, UR] x [VB, VT], inside N_TessTrimmedSrf's
domain (disjoint, simple, outer loops counterclockwise, inner clockwise):
(values cuo cui m ho hi hs)."
  (let* ((w (- ur ul)) (h (- vt vb))
         (cx (+ ul (* 0.5d0 w))) (cy (+ vb (* 0.5d0 h)))
         (outer nil) (inner nil))
    (case (rnd 11)
      ((0 1 2)
       ;; untrimmed: the domain's boundary, or with corners cut off
       (push (%tess-box-loop ul ur vb vt :cuts (if (zerop (rnd 2)) 0 (1+ (rnd 4)))) outer)
       (push (%tess-holes cx cy (* 0.4d0 w) (* 0.4d0 h)) inner))
      ((8 9)
       ;; the domain's boundary with holes on the subdivision's grid
       (push (%tess-box-loop ul ur vb vt) outer)
       (push (%tess-grid-holes ul ur vb vt 0d0 1d0 0d0 1d0) inner))
      (10
       ;; an outer polygon on the grid, with holes on it
       (let ((f (nth (rnd 3) '(0.125d0 0.25d0 0.0625d0))))
         (flet ((at (fx fy) (cons (+ ul (* w fx)) (+ vb (* h fy)))))
           (push (%tess-polygon-loop
                  (if (zerop (rnd 2))
                      (list (at f f) (at (- 1 f) f) (at (- 1 f) (- 1 f)) (at f (- 1 f)))
                      (list (at 0.5d0 f) (at (- 1 f) 0.5d0) (at 0.5d0 (- 1 f)) (at f 0.5d0))))
                 outer))
         (push (%tess-grid-holes ul ur vb vt 0.375d0 0.625d0 0.375d0 0.625d0) inner)))
      ((3 4)
       ;; a curved outer loop inside the domain
       (let ((rx (* w (+ 0.3d0 (rnd 0.18d0)))) (ry (* h (+ 0.3d0 (rnd 0.18d0)))))
         (push (%tess-curve-loop cx cy rx ry nil) outer)
         (push (%tess-holes cx cy (* 0.6d0 rx) (* 0.6d0 ry)) inner)))
      (5
       ;; the domain with its upper right corner cut off by a curve
       (let* ((um (+ ul (* w (+ 0.3d0 (rnd 0.5d0))))) (vm (+ vb (* h (+ 0.3d0 (rnd 0.5d0)))))
              (mid (cons (+ um (* (- ur um) (+ 0.3d0 (rnd 0.6d0))))
                         (+ vm (* (- vt vm) (+ 0.3d0 (rnd 0.6d0))))))
              (cut (if (zerop (rnd 3))
                       (line-curve ur vm um vt :degree (1+ (rnd 2)))
                       (%make-curve 2 (make-array 6 :element-type 'double-float
                                                    :initial-contents '(0d0 0d0 0d0 1d0 1d0 1d0))
                                    (let ((cp (make-cpoint-array 3)))
                                      (setf (aref cp 0 0) ur (aref cp 0 1) vm
                                            (aref cp 1 0) (car mid) (aref cp 1 1) (cdr mid)
                                            (aref cp 2 0) um (aref cp 2 1) vt)
                                      (dotimes (i 3 cp)
                                        (setf (aref cp i 2) +no-z+ (aref cp i 3) +no-w+)))))))
         (push (vector (line-curve ul vb ur vb) (line-curve ur vb ur vm) cut
                       (line-curve um vt ul vt) (line-curve ul vt ul vb))
               outer)
         (push (if (zerop (rnd 2))
                   (vector)
                   (vector (%tess-curve-loop (+ ul (* 0.25d0 w)) (+ vb (* 0.25d0 h))
                                             (* 0.1d0 w) (* 0.1d0 h) t)))
               inner)))
      (t
       ;; two outer loops side by side
       (dolist (fx '(0.25d0 0.75d0))
         (let ((ox (+ ul (* fx w))) (rx (* w (+ 0.15d0 (rnd 0.08d0)))) (ry (* h (+ 0.3d0 (rnd 0.15d0)))))
           (push (%tess-curve-loop ox cy rx ry nil) outer)
           (push (%tess-holes ox cy (* 0.6d0 rx) (* 0.6d0 ry)) inner)))))
    (let* ((cuo (coerce (nreverse outer) 'simple-vector))
           (cui (coerce (nreverse inner) 'simple-vector))
           (m (1- (length cuo)))
           (ho (make-array (1+ m) :element-type 'fixnum))
           (hi (make-array (1+ m) :element-type 'fixnum))
           (hs (make-array (1+ m))))
      (dotimes (i (1+ m))
        (setf (aref ho i) (1- (length (svref cuo i)))
              (aref hi i) (1- (length (svref cui i)))
              (svref hs i) (map '(simple-array fixnum (*)) (lambda (l) (1- (length l))) (svref cui i))))
      (values cuo cui m ho hi hs))))

(defun %tess-trim-case ()
  "A random input to N_TessTrimmedSrf: (list sur cuo cui m ho hi hs epc eps tol),
EPS and EPC between about a thousandth and a tenth of the surface's size."
  (let* ((sur (%tess-surface))
         (size (sqrt (max 1d-6 (values (batten::tess-srf-area sur)))))
         (eps (* size (expt 10d0 (- (+ 1d0 (rnd 2d0))))))
         (epc (* eps (+ 0.3d0 (rnd 2d0))))
         (big (* eps (+ 0.001d0 (rnd 0.01d0))))
         (tol (expt 10d0 (- (+ 6d0 (rnd 3d0))))))
    (multiple-value-bind (ul ur vb vt) (batten::srf-get-parameter-bounds sur)
      ;; TOL about 1e-7, small against eps; now and then larger than 0.001
      ;; eps, which the C lowers to that -- when that is still small: the C
      ;; scales it by max(lu, lv) against a domain rescaled to [0, lu] x [0,
      ;; lv], and a tolerance near the domain's size collapses the trimming
      ;; polygons (outside the C's domain)
      (multiple-value-bind (a lu lv) (batten::tess-srf-area sur)
        (declare (ignore a))
        (when (and (zerop (rnd 6))
                   (<= (* 0.001d0 eps (max lu lv)) (* 1d-2 (min lu lv))))
          (setf tol big)))
      (multiple-value-bind (cuo cui m ho hi hs) (%tess-loops ul ur vb vt)
        (list sur cuo cui m ho hi hs epc eps tol)))))

(defun %tess-rows (dt hd n)
  "The rows 0..N of a triangulation DT, each cut to its HD[i]+1 entries."
  (let ((out (make-array (1+ n))))
    (dotimes (i (1+ n) out)
      (setf (svref out i) (subseq (svref dt i) 0 (1+ (aref hd i)))))))

(defmacro %with-c-tess-trimmed-srf ((err sp u v n dt hd sg) case &body body)
  "Run N_TessTrimmedSrf on CASE (a %tess-trim-case list); bind ERR (the C
failed), SP (the foreign surface, rescaled by it), U, V, N, DT, HD (the
output cells) and SG (their stack) for BODY."
  (let ((c (gensym)))
    `(let ((,c ,case))
       (destructuring-bind (sur cuo cui m ho hi hs epc eps tol) ,c
         (with-foreign ((,sp sur :surface) (co cuo :curve-tree) (ci cui :curve-tree)
                        (hop ho :index-vector) (hip hi :index-vector) (hsp hs :index-lists))
           (with-stacks (,sg)
             (with-scalars ((,u :pointer) (,v :pointer) (,n :long) (,dt :pointer) (,hd :pointer))
               (let ((,err (c-flag-error-p
                            (c-tess-trimmed-srf ,sp co ci m hop hip hsp epc eps tol ,u ,v ,n ,dt ,hd ,sg))))
                 ,@body))))))))

(defun %lisp-tess-trimmed-srf (case)
  "Run tess-trimmed-srf on a copy of CASE's surface: (list surface u v n dt hd),
U, V cut to N+1 and DT's rows to their HD[i]+1 entries."
  (destructuring-bind (sur cuo cui m ho hi hs epc eps tol) case
    (let ((s (batten::srf-copy sur)))
      (multiple-value-bind (u v n dt hd) (batten::tess-trimmed-srf s cuo cui m ho hi hs epc eps tol)
        (list s (subseq u 0 (1+ n)) (subseq v 0 (1+ n)) n (%tess-rows dt hd n) hd)))))

(define-oracle-test tess-trimmed-srf (:file "tessellate" :trials 15000)
  (let ((case (%tess-trim-case)))
    (agree case
           (%lisp-tess-trimmed-srf case)
           (%with-c-tess-trimmed-srf (err sp u v n dt hd sg) case
             (values err
                     (unless err
                       (let* ((n (cffi:mem-ref n :long))
                              (hdv (from-foreign (cffi:mem-ref hd :pointer) :index-vector (1+ n))))
                         (list (from-foreign sp :surface)
                               (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ n))
                               (from-foreign (cffi:mem-ref v :pointer) :real-vector (1+ n))
                               n
                               (from-foreign (cffi:mem-ref dt :pointer) :index-lists
                                             (map 'list #'1+ hdv))
                               hdv))))))))

(define-oracle-test tess-trimmed-srf-triangles (:file "tessellate" :trials 8000)
  ;; tess-trimmed-srf's triangulation handed to tess-get-triangles, on each side
  (let ((case (%tess-trim-case)))
    (agree case
           (destructuring-bind (s u v n dt hd) (%lisp-tess-trimmed-srf case)
             (declare (ignore s u v))
             (multiple-value-list (batten::tess-get-triangles dt hd n)))
           (%with-c-tess-trimmed-srf (err sp u v n dt hd sg) case
             (if err
                 (values t nil)
                 (with-scalars ((alf :pointer) (bet :pointer) (gam :pointer) (k :long))
                   (let ((err (c-flag-error-p
                               (c-tess-get-triangles (cffi:mem-ref dt :pointer) (cffi:mem-ref hd :pointer)
                                                     (cffi:mem-ref n :long) alf bet gam k sg))))
                     (values err
                             (unless err
                               (let ((k (cffi:mem-ref k :long)))
                                 (list (from-foreign (cffi:mem-ref alf :pointer) :index-vector (1+ k))
                                       (from-foreign (cffi:mem-ref bet :pointer) :index-vector (1+ k))
                                       (from-foreign (cffi:mem-ref gam :pointer) :index-vector (1+ k))
                                       k)))))))))))
