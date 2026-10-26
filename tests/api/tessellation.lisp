;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Tests of the API in src/api/tessellation.lisp: tessellate-surface against
;;;; the C routines it wraps -- N_TessTrimmedSrf, N_TessGetTriangles, and
;;;; N_SrfEvalPt on the surface the tessellator rescaled.

(in-package :batten-tests)

(defun %api-tess-case ()
  "A tessellator case (as %tess-trim-case) and the same loops in
tessellate-surface's form: (values case outer-loops holes), OUTER-LOOPS nil
(and the case's loops the domain's boundary) a third of the time."
  (let ((case (%tess-trim-case)))
    (if (zerop (rnd 3))
        (destructuring-bind (sur cuo cui m ho hi hs epc eps tol) case
          (declare (ignore cuo cui m ho hi hs))
          (values (list sur (vector (batten::%domain-loop sur)) (vector (vector)) 0
                        (make-array 1 :element-type 'fixnum :initial-element 3)
                        (make-array 1 :element-type 'fixnum :initial-element -1)
                        (vector (make-array 0 :element-type 'fixnum))
                        epc eps tol)
                  nil nil))
        (destructuring-bind (sur cuo cui m &rest rest) case
          (declare (ignore sur rest))
          (let ((outer (map 'list #'(lambda (loop) (coerce loop 'list)) cuo))
                (inner (map 'list #'(lambda (loops) (map 'list #'(lambda (loop) (coerce loop 'list)) loops)) cui)))
            (values case outer (if (zerop m) (first inner) inner)))))))

(define-oracle-test api-tessellate-surface (:file "api-tessellation" :trials 1500)
  (multiple-value-bind (case outer holes) (%api-tess-case)
    (destructuring-bind (sur cuo cui m ho hi hs epc eps tol) case
      (declare (ignore cuo cui m ho hi hs))
      (agree (list :case case :outer outer :holes holes)
             (let ((copy (batten::srf-copy sur)))
               (prog1 (multiple-value-list
                       (tessellate-surface sur :tolerance eps :trim-tolerance epc :parameter-tolerance tol
                                               :outer-loops outer :holes holes))
                 (unless (equalp copy sur) (error "tessellate-surface changed its surface"))))
             (%with-c-tess-trimmed-srf (err sp u v n dt hd sg) case
               (if err
                   (values t nil)
                   (with-scalars ((alf :pointer) (bet :pointer) (gam :pointer) (k :long))
                     (let ((terr (c-flag-error-p
                                  (c-tess-get-triangles (cffi:mem-ref dt :pointer) (cffi:mem-ref hd :pointer)
                                                        (cffi:mem-ref n :long) alf bet gam k sg))))
                       (if terr
                           (values t nil)
                           (let* ((n (cffi:mem-ref n :long))
                                  (k (cffi:mem-ref k :long))
                                  (us (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ n)))
                                  (vs (from-foreign (cffi:mem-ref v :pointer) :real-vector (1+ n)))
                                  (a (from-foreign (cffi:mem-ref alf :pointer) :index-vector (1+ k)))
                                  (b (from-foreign (cffi:mem-ref bet :pointer) :index-vector (1+ k)))
                                  (g (from-foreign (cffi:mem-ref gam :pointer) :index-vector (1+ k)))
                                  (points (make-point-array (1+ n)))
                                  (parameters (make-array (list (1+ n) 2) :element-type 'double-float))
                                  (triangles (make-array (list (1+ k) 3) :element-type 'fixnum)))
                             (dotimes (i (1+ n))
                               (with-foreign-out ((pt :point))
                                 (c-srf-eval-pt sp (aref us i) (aref vs i) +left+ +left+ pt)
                                 (let ((p (from-foreign pt :point)))
                                   (setf (aref points i 0) (aref p 0) (aref points i 1) (aref p 1)
                                         (aref points i 2) (aref p 2))))
                               (setf (aref parameters i 0) (aref us i) (aref parameters i 1) (aref vs i)))
                             (dotimes (j (1+ k))
                               (setf (aref triangles j 0) (aref a j) (aref triangles j 1) (aref b j)
                                     (aref triangles j 2) (aref g j)))
                             (values nil (list points triangles parameters))))))))))))
