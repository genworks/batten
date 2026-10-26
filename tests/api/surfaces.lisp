;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Tests of the API in src/api/surfaces.lisp: each function against the C routines
;;;; it wraps (or the port routines, where the C has no single counterpart).

(in-package :batten-tests)

;;; ---------------------------------------------------------------- helpers

(defun %api-srf-snapshot (x)
  "A copy of X (a surface, or a point or other array) sharing nothing with it."
  (if (surface-p x)
      (%copy-surface-deep x)
      (let ((c (make-array (array-dimensions x) :element-type (array-element-type x))))
        (dotimes (i (array-total-size x) c)
          (setf (row-major-aref c i) (row-major-aref x i))))))

(defun %api-srf-walk (x fn)
  "Call FN on every surface, curve and array inside X (lists, arrays of
objects, surfaces and curves with their arrays)."
  (cond ((consp x) (%api-srf-walk (car x) fn) (%api-srf-walk (cdr x) fn))
        ((surface-p x)
         (funcall fn x)
         (mapc fn (list (surface-knots-u x) (surface-knots-v x) (surface-cpoints x))))
        ((curve-p x)
         (funcall fn x)
         (mapc fn (list (curve-knots x) (curve-cpoints x))))
        ((and (arrayp x) (not (stringp x)) (eq (array-element-type x) t))
         (dotimes (i (array-total-size x)) (%api-srf-walk (row-major-aref x i) fn)))
        ((arrayp x) (funcall fn x))))

(defun %api-srf-shares-p (results arg)
  "True when RESULTS holds ARG or one of its arrays."
  (let ((mine nil))
    (%api-srf-walk arg (lambda (o) (push o mine)))
    (when (vectorp arg) (push arg mine))
    (%api-srf-walk results (lambda (o) (when (member o mine :test #'eq) (return-from %api-srf-shares-p t))))
    nil))

(defmacro %api-srf-pure ((&rest args) form)
  "FORM (the API's results, a list) with :unchanged appended when every one of
ARGS is as it was before FORM and shares nothing with the results, else
:modified."
  (let ((snaps (gensym "SNAPS")) (results (gensym "RESULTS")) (vals (gensym "ARGS")))
    `(let* ((,vals (list ,@args))
            (,snaps (mapcar #'%api-srf-snapshot ,vals))
            (,results ,form))
       (append ,results
               (list (if (and (every (lambda (a s) (values (same a s))) ,vals ,snaps)
                              (notany (lambda (a) (%api-srf-shares-p ,results a)) ,vals))
                         :unchanged
                         :modified))))))

(defmacro %api-srf-c (form)
  "FORM's (values failed results) with :unchanged appended to the results,
matching %API-SRF-PURE."
  (let ((f (gensym)) (r (gensym)))
    `(multiple-value-bind (,f ,r) ,form
       (values ,f (unless ,f (append ,r (list :unchanged)))))))

(defun %api-srf-dir ()
  "(values :u-or-:v NLib-flag)."
  (if (zerop (rnd 2)) (values :u +udir+) (values :v +vdir+)))

(defun %api-srf-sequence (v)
  "Real-vector V as itself, a list or a general vector, at random."
  (case (rnd 3) (0 v) (1 (coerce v 'list)) (t (coerce v 'simple-vector))))

(defun %api-srf-side ()
  "(values :left-or-:right NLib-flag)."
  (if (zerop (rnd 2)) (values :left +left+) (values :right +right+)))

;;; --------------------------------------------------------- generic functions

(define-oracle-test api-domain-surface (:file "api-surfaces" :trials 500)
  (let ((s (%some-srf)))
    (agree s
           (%api-srf-pure (s) (multiple-value-list (domain s)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-scalars ((ul :double) (ur :double) (vb :double) (vt :double))
                (c-srf-get-parameter-bounds sp ul ur vb vt)
                (values nil (list (cffi:mem-ref ul :double) (cffi:mem-ref ur :double)
                                  (cffi:mem-ref vb :double) (cffi:mem-ref vt :double)))))))))

(define-oracle-test api-bounding-box-surface (:file "api-surfaces" :trials 500)
  (let ((s (%some-srf)))
    (agree s
           (%api-srf-pure (s) (multiple-value-list (bounding-box s)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-foreign-out ((box :box))
                (if (c-flag-error-p (c-srf-get-b-box sp box))
                    (values t nil)
                    (let ((b (from-foreign box :box)))
                      (values nil (list (make-point (aref b 0) (aref b 2) (aref b 4))
                                        (make-point (aref b 1) (aref b 3) (aref b 5))))))))))))

(define-oracle-test api-copy-surface (:file "api-surfaces" :trials 500)
  (let ((s (%some-srf)))
    (agree s
           (%api-srf-pure (s) (list (copy s)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-stacks (st)
                (with-foreign-out ((sq :surface))
                  (if (c-flag-error-p (c-srf-copy sp sq st))
                      (values t nil)
                      (values nil (list (from-foreign sq :surface)))))))))))

(define-oracle-test api-transform-surface (:file "api-surfaces" :trials 500)
  (let* ((s (%some-srf))
         (rm (%random-transform))
         ;; now and then a general array of reals rather than a real-matrix
         (matrix (if (zerop (rnd 3))
                     (let ((a (make-array '(4 4))))
                       (dotimes (i 4 a) (dotimes (j 4) (setf (aref a i j) (aref rm i j)))))
                     rm)))
    (agree (list :surface s :matrix rm)
           (%api-srf-pure (s matrix) (list (transform s matrix)))
           (%api-srf-c (%api-srf-c-transform s rm)))))

(defun %api-srf-c-transform (s rm)
  "Surface S with every control point transformed by N_TransformCPt with the
4x4 matrix RM: (values nil (list surface))."
  (let* ((r (%copy-surface-deep s)) (pw (surface-cpoints r)))
    (with-foreign ((rmp rm :real-matrix))
      (let ((frm (%foreign-rmatrix 3 3 rmp 0 4)))
        (dotimes (i (array-dimension pw 0))
          (dotimes (j (array-dimension pw 1))
            (let ((cp (make-cpoint (aref pw i j 0) (aref pw i j 1) (aref pw i j 2) (aref pw i j 3))))
              (with-foreign ((a cp :cpoint))
                (with-foreign-out ((q :cpoint))
                  (c-transform-c-pt a frm q)
                  (let ((out (from-foreign q :cpoint)))
                    (dotimes (k 4) (setf (aref pw i j k) (aref out k)))))))))))
    (values nil (list r))))

(define-oracle-test api-translate-surface (:file "api-surfaces" :trials 500)
  (let* ((s (%some-srf)) (tv (random-point))
         (vector (if (zerop (rnd 2)) tv (coerce tv 'list))))
    (agree (list :surface s :vector tv)
           (%api-srf-pure (s tv) (list (translate s vector)))
           (%api-srf-c
            (with-foreign ((sp s :surface) (tp tv :point))
              (c-srf-translate sp tp)
              (values nil (list (from-foreign sp :surface))))))))

(define-oracle-test api-insert-knots-surface (:file "api-surfaces" :trials 1000)
  (multiple-value-bind (direction dir) (%api-srf-dir)
    (let* ((s (%some-srf))
           (knx (if (= dir +vdir+)
                    (%random-insert-knots (surface-knots-v s) (surface-degree-v s))
                    (%random-insert-knots (surface-knots-u s) (surface-degree-u s))))
           ;; the API sorts: hand it the knots reversed now and then
           (knots (%api-srf-sequence (if (zerop (rnd 2)) (reverse knx) knx))))
      (agree (list :surface s :knots knx :direction direction)
             (%api-srf-pure (s) (list (insert-knots s knots :direction direction)))
             (%api-srf-c
              (with-foreign ((sp s :surface))
                (with-foreign-out ((sq :surface))
                  (with-stacks (stp)
                    (with-stacks (stq)
                      (if (zerop (length knx))
                          (if (c-flag-error-p (c-srf-copy sp sq stp))
                              (values t nil)
                              (values nil (list (from-foreign sq :surface))))
                          (with-foreign ((kx knx :knots))
                            (if (c-flag-error-p (c-srf-insert-knots sp kx dir sq stp stq))
                                (values t nil)
                                (values nil (list (from-foreign sq :surface)))))))))))))))

(define-oracle-test api-elevate-degree-surface (:file "api-surfaces" :trials 1000)
  (multiple-value-bind (direction dir) (%api-srf-dir)
    (let ((s (%some-srf)) (by (case (rnd 30) (0 -1) (t (rnd 4)))))
      (agree (list :surface s :by by :direction direction)
             (%api-srf-pure (s) (list (elevate-degree s by :direction direction)))
             (%api-srf-c
              (with-foreign ((sp s :surface))
                (with-foreign-out ((sq :surface))
                  (with-stacks (stp)
                    (with-stacks (stq)
                      (if (c-flag-error-p (c-srf-elevate-degree sp by dir sq stp stq))
                          (values t nil)
                          (values nil (list (from-foreign sq :surface)))))))))))))

(define-oracle-test api-reparameterize-surface (:file "api-surfaces" :trials 500)
  (let* ((s (%some-srf))
         (ku (surface-knots-u s)) (kv (surface-knots-v s)))
    (flet ((ends (a b)
             (case (rnd 4)
               (0 (values a b))
               (1 (values a (+ a 0.5d0 (rnd 10d0))))
               (t (let ((lo (- (rnd 20d0) 10d0))) (values lo (+ lo 0.5d0 (rnd 10d0))))))))
      (multiple-value-bind (us ue) (ends (aref ku 0) (aref ku (1- (length ku))))
        (multiple-value-bind (vs ve) (ends (aref kv 0) (aref kv (1- (length kv))))
          (agree (list :surface s :us us :ue ue :vs vs :ve ve)
                 (%api-srf-pure (s) (list (reparameterize s us ue vs ve)))
                 (%api-srf-c
                  (with-foreign ((sp s :surface))
                    (c-srf-reparam sp us ue vs ve)
                    (values nil (list (from-foreign sp :surface)))))))))))

(define-oracle-test api-split-surface (:file "api-surfaces" :trials 1000)
  (multiple-value-bind (direction dir) (%api-srf-dir)
    (let* ((s (%some-srf))
           (at (%srf-tool-param (if (= dir +vdir+) (surface-knots-v s) (surface-knots-u s)))))
      (agree (list :surface s :at at :direction direction)
             (%api-srf-pure (s) (multiple-value-list (split s at :direction direction)))
             (%api-srf-c
              (with-foreign ((sp s :surface))
                (with-foreign-out ((sl :surface) (sr :surface))
                  (with-stacks (sg)
                    (if (c-flag-error-p (c-srf-split sp at dir sl sr sg))
                        (values t nil)
                        (values nil (list (from-foreign sl :surface) (from-foreign sr :surface))))))))))))

;;; ------------------------------------------------------------------ queries

(define-oracle-test api-surface-normal (:file "api-surfaces" :trials 1000)
  (multiple-value-bind (u-side ufl) (%api-srf-side)
    (multiple-value-bind (v-side vfl) (%api-srf-side)
      (let* ((s (random-surface))
             (u (random-parameter (surface-knots-u s))) (v (random-parameter (surface-knots-v s))))
        ;; now and then a pole: the first row of control points collapsed
        (when (zerop (rnd 6))
          (let ((pw (surface-cpoints s)))
            (dotimes (j (array-dimension pw 1))
              (dotimes (k 4) (setf (aref pw 0 j k) (aref pw 0 0 k))))
            (setf u (aref (surface-knots-u s) 0))))
        (agree (list :surface s :u u :v v :u-side u-side :v-side v-side)
               (%api-srf-pure (s) (multiple-value-list (batten::surface-normal s u v :u-side u-side :v-side v-side)))
               (%api-srf-c
                (with-foreign ((sp s :surface))
                  (with-foreign-out ((p :point) (su :point) (sv :point) (n :point))
                    (values (c-flag-error-p (c-srf-eval-pt-pt-deriv-normal sp u v ufl vfl p su sv n))
                            (list (from-foreign p :point) (from-foreign n :point)
                                  (from-foreign su :point) (from-foreign sv :point)))))))))))

(define-oracle-test api-surface-grid (:file "api-surfaces" :trials 1000)
  (multiple-value-bind (u-side ufl) (%api-srf-side)
    (multiple-value-bind (v-side vfl) (%api-srf-side)
      (let* ((s (%some-srf))
             (nn (rnd 9)) (mm (rnd 9))
             (u (%srf-adv-grid (surface-knots-u s) nn))
             (v (%srf-adv-grid (surface-knots-v s) mm)))
        (agree (list :surface s :us u :vs v :u-side u-side :v-side v-side)
               (%api-srf-pure (s) (list (batten::surface-grid s (%api-srf-sequence u) (%api-srf-sequence v)
                                                              :u-side u-side :v-side v-side)))
               (%api-srf-c
                (with-foreign ((sp s :surface) (cu u :real-vector) (cv v :real-vector))
                  (with-foreign-out ((pts :point-net (1+ nn) (1+ mm)))
                    (if (c-flag-error-p (c-srf-eval-pt-grid sp cu cv nn mm ufl vfl pts))
                        (values t nil)
                        (values nil (list (from-foreign pts :point-net (1+ nn) (1+ mm)))))))))))))

(defun %api-srf-c-closest-start (s p)
  "The start surface-closest-point takes for point P on surface S, from the
C's grid of the same samples: (values u v)."
  (let* ((ku (surface-knots-u s)) (kv (surface-knots-v s))
         (k (* 2 (1+ (surface-n s)))) (l (* 2 (1+ (surface-m s))))
         (us (batten::%surface-sample-parameters (aref ku 0) (aref ku (1- (length ku))) k))
         (vs (batten::%surface-sample-parameters (aref kv 0) (aref kv (1- (length kv))) l))
         (grid (with-foreign ((sp s :surface) (cu us :real-vector) (cv vs :real-vector))
                 (with-foreign-out ((pts :point-net (1+ k) (1+ l)))
                   (when (c-flag-error-p (c-srf-eval-pt-grid sp cu cv k l +left+ +left+ pts))
                     (error "N_SrfEvalPtGrid failed on the closest-point samples"))
                   (from-foreign pts :point-net (1+ k) (1+ l)))))
         (best most-positive-double-float) (bu 0d0) (bv 0d0))
    (dotimes (i (1+ k))
      (dotimes (j (1+ l))
        (let ((d (batten::dist-pt-pt (make-point (aref grid i j 0) (aref grid i j 1) (aref grid i j 2)) p)))
          (when (< d best) (setf best d bu (aref us i) bv (aref vs j))))))
    (values bu bv)))

(define-oracle-test api-surface-closest-point (:file "api-surfaces" :trials 500)
  (let* ((s (random-surface :max-degree 4 :max-extra 5 :two-d (zerop (rnd 6))))
         (s (progn (when (zerop (rnd 4)) (%close-surface s +udir+))
                   (when (zerop (rnd 4)) (%close-surface s +vdir+))
                   s))
         (ku (surface-knots-u s)) (kv (surface-knots-v s))
         (on (batten::srf-eval-pt s (random-parameter ku) (random-parameter kv) +left+ +left+ (make-point)))
         (p (case (rnd 4)
              (0 on)
              (1 (%pt+ on 1d0 (random-point 5d0)))
              (2 (%pt+ on 1d0 (random-point 50d0)))
              (t (random-point))))
         (guess (zerop (rnd 3)))
         (u0 (when guess (random-parameter ku)))
         (v0 (when guess (random-parameter kv)))
         (top (case (rnd 3) (0 batten::*mtol*) (1 1d-4) (t 1d-2)))
         (toc (case (rnd 3) (0 1d-5) (1 1d-3) (t 1d-9)))
         (error-p (zerop (rnd 2))))
    (agree (list :surface s :p p :u0 u0 :v0 v0 :top top :toc toc :error-p error-p)
           (%api-srf-pure (s p)
             (multiple-value-list
              (if guess
                  (batten::surface-closest-point s (coerce p 'list) :u u0 :v v0 :tolerance top
                                                 :cosine-tolerance toc :error-p error-p)
                  (batten::surface-closest-point s p :tolerance top :cosine-tolerance toc
                                                 :error-p error-p))))
           (%api-srf-c
            (multiple-value-bind (su sv) (if guess (values u0 v0) (%api-srf-c-closest-start s p))
              (with-foreign ((sp s :surface) (pp p :point))
                (with-foreign-out ((q :point) (d :point-net 3 3))
                  (with-scalars ((u :double) (v :double) (dist :double))
                    ;; the C reports the best point found even when it fails
                    (let ((failed (c-flag-error-p (c-get-closest-pt-on-srf sp pp su sv top toc u v q d))))
                      (if (and failed error-p)
                          (values t nil)
                          (progn
                            (c-dist-pt-pt q pp dist)
                            (values nil (list (cffi:mem-ref u :double) (cffi:mem-ref v :double)
                                              (from-foreign q :point) (cffi:mem-ref dist :double)
                                              (not failed))))))))))))))

(define-oracle-test api-isocurve (:file "api-surfaces" :trials 1000)
  (let* ((s (%some-srf))
         (at-u (zerop (rnd 2)))
         (dir (if at-u +vdir+ +udir+))
         (tt (%srf-tool-param (if at-u (surface-knots-u s) (surface-knots-v s)))))
    (agree (list :surface s (if at-u :u :v) tt)
           (%api-srf-pure (s) (list (if at-u (batten::isocurve s :u tt) (batten::isocurve s :v tt))))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-foreign-out ((cur :curve))
                (with-stacks (sc)
                  (if (c-flag-error-p (c-srf-extract-iso-crv sp tt dir cur sc))
                      (values t nil)
                      (values nil (list (from-foreign cur :curve)))))))))))

(define-oracle-test api-subsurface (:file "api-surfaces" :trials 1000)
  (let ((s (%some-srf)))
    (multiple-value-bind (ul ur) (%patch-interval (surface-knots-u s))
      (multiple-value-bind (vl vr) (%patch-interval (surface-knots-v s))
        (agree (list :surface s :ul ul :ur ur :vl vl :vr vr)
               (%api-srf-pure (s) (list (batten::subsurface s ul ur vl vr)))
               (%api-srf-c
                (with-foreign ((sp s :surface))
                  (with-foreign-out ((sq :surface))
                    (with-stacks (stp)
                      (with-stacks (stq)
                        (if (c-flag-error-p (c-srf-extract-patch sp ul ur vl vr sq stp stq))
                            (values t nil)
                            (values nil (list (from-foreign sq :surface))))))))))))))

(define-oracle-test api-bezier-patches (:file "api-surfaces" :trials 500)
  (let ((s (%some-srf)))
    (agree s
           (%api-srf-pure (s) (list (batten::bezier-patches s)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-stacks (sq)
                (cffi:with-foreign-object (cell :pointer)
                  (setf (cffi:mem-ref cell :pointer) (cffi:null-pointer))
                  (with-scalars ((kk :long) (ll :long))
                    (if (c-flag-error-p (c-srf-decompose-to-bez sp cell kk ll sq))
                        (values t nil)
                        (let* ((k (cffi:mem-ref kk :long)) (l (cffi:mem-ref ll :long))
                               (a (cffi:mem-ref cell :pointer))
                               (out (make-array (list (1+ k) (1+ l)))))
                          (dotimes (i (1+ k))
                            (let ((row (cffi:mem-aref a :pointer i)))
                              (dotimes (j (1+ l))
                                (setf (aref out i j) (from-foreign (cffi:mem-aref row :pointer j) :surface)))))
                          (values nil (list out))))))))))))

(define-oracle-test api-offset-surface (:file "api-surfaces" :trials 200)
  (let* ((s (%srf-offset-input))
         (d (* (if (zerop (rnd 2)) 1d0 -1d0) (rnd-range 0.1d0 4d0)))
         (pp (1+ (rnd 3))) (qq (1+ (rnd 3)))
         (parameterization (case (rnd 6) ((0 1 2) :inherited) (3 :uniform) (4 :chord-length) (t :centripetal)))
         (par (ecase parameterization
                (:inherited +inherited+) (:uniform +uniform+)
                (:chord-length +chordlength+) (:centripetal +centripetal+)))
         (tol (case (rnd 3) (0 0.05d0) (1 0.5d0) (t 2d0)))
         (eps (case (rnd 4) (0 1d-6) (1 1d-4) (2 1d-3) (t (* 0.1d0 tol)))))
    (agree (list :surface s :d d :pp pp :qq qq :par par :eps eps :tol tol)
           (%api-srf-pure (s)
             (list (batten::offset-surface s d :tolerance tol :type-tolerance eps
                                               :degree-u pp :degree-v qq
                                               :parameterization parameterization)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-stacks (sg)
                (with-foreign-out ((sq :surface))
                  (let ((failed (c-flag-error-p (c-srf-offset sp d pp qq par eps tol sq sg))))
                    (values failed (unless failed (list (from-foreign sq :surface))))))))))))

(define-oracle-test api-surface-swap-uv (:file "api-surfaces" :trials 500)
  (let ((s (%some-srf)))
    (agree s
           (%api-srf-pure (s) (list (batten::surface-swap-uv s)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-stacks (st)
                (values (c-flag-error-p (c-swap-uv sp st))
                        (list (from-foreign sp :surface)))))))))

(define-oracle-test api-surface-closed-p (:file "api-surfaces" :trials 1000)
  (let* ((s (%some-srf))
         (direction (case (rnd 3) (0 nil) (1 :u) (t :v))))
    (when (zerop (rnd 3)) (%close-surface s +udir+))
    (when (zerop (rnd 3)) (%close-surface s +vdir+))
    (agree (list :surface s :direction direction)
           (%api-srf-pure (s) (list (if (batten::surface-closed-p s :direction direction) 1 0)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (values nil
                      (list (ecase direction
                              (:u (c-srf-is-closed sp +udir+))
                              (:v (c-srf-is-closed sp +vdir+))
                              ((nil) (if (or (= 1 (c-srf-is-closed sp +udir+))
                                             (= 1 (c-srf-is-closed sp +vdir+)))
                                         1 0))))))))))

(define-oracle-test api-surface-flat-p (:file "api-surfaces" :trials 1000)
  (let* ((s (%sb-flattish-surface))
         (eps (case (rnd 4) (0 0d0) (1 1d-5) (2 5d-2) (t (rnd 50d0)))))
    (agree (list :surface s :tolerance eps)
           (%api-srf-pure (s) (list (if (batten::surface-flat-p s :tolerance eps) 1 0)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-scalars ((flt :short) (dir :short))
                (let ((rc (c-srf-is-flat-cheap sp eps flt dir)))
                  (values (c-flag-error-p rc)
                          (list (if (= (cffi:mem-ref flt :short) +yes+) 1 0))))))))))

(define-oracle-test api-surface-type (:file "api-surfaces" :trials 500)
  (let ((s (%st-surface))
        (tol (case (rnd 4) (0 1d-3) (1 1d-4) (t 1d-6))))
    (agree (list :surface s :tolerance tol)
           (%api-srf-pure (s) (multiple-value-list (batten::surface-type s :tolerance tol)))
           (%api-srf-c
            (with-foreign ((sp s :surface))
              (with-stacks (sg)
                (with-foreign-out ((faa :point) (fbb :point) (fcc :point)
                                   (fdd :point) (fxx :point) (fyy :point))
                  (with-scalars ((ca :pointer) (cb :pointer) (r1 :double) (r2 :double)
                                 (ac :double) (al :double) (h :double) (stp :short) (dir :short))
                    (let ((failed (c-flag-error-p
                                   (c-srf-type sp tol ca cb faa fbb fcc fdd fxx fyy
                                               r1 r2 ac al h stp dir sg))))
                      (flet ((crv (cell)
                               (let ((pc (cffi:mem-ref cell :pointer)))
                                 (if (cffi:null-pointer-p pc) nil (from-foreign pc :curve)))))
                        (values failed
                                (unless failed
                                  (multiple-value-list
                                   (batten::%surface-type-properties
                                    (cffi:mem-ref stp :short) (cffi:mem-ref dir :short)
                                    (from-foreign faa :point) (from-foreign fbb :point)
                                    (from-foreign fcc :point) (from-foreign fdd :point)
                                    (from-foreign fxx :point) (from-foreign fyy :point)
                                    (cffi:mem-ref r1 :double) (cffi:mem-ref r2 :double)
                                    (cffi:mem-ref ac :double) (cffi:mem-ref al :double)
                                    (cffi:mem-ref h :double) (crv ca) (crv cb)))))))))))))))
