;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/funcs-basic.lisp against NLib.
;;;;
;;;; Most of these routines make or read pointer assignments in NLib's
;;;; NL_CFUN / NL_SFUN structs; the tests build the C structs, call the C, and
;;;; read back what its outputs point at.

(in-package :batten-tests)

(defmacro %fb-slot (ptr struct slot)
  "A slot of one of NLib's structs as batten-oracle defines them."
  `(cffi:foreign-slot-value ,ptr '(:struct ,(intern (string struct) :batten-oracle))
                            ',(intern (string slot) :batten-oracle)))

(defun %fb-cvalue (cvl)
  "The control values an NL_CVALUE pointer holds, as a real-vector."
  (from-foreign (%fb-slot cvl nl-cvalue fu) :real-vector (1+ (%fb-slot cvl nl-cvalue n))))

(defun %fb-svalue (svl)
  "The control values an NL_SVALUE pointer holds, as a real-matrix."
  (from-foreign (%fb-slot svl nl-svalue fuv) :real-matrix
                (1+ (%fb-slot svl nl-svalue n)) (1+ (%fb-slot svl nl-svalue m))))

(defun %fb-near (k) "K, or one either side of it." (+ k (- (rnd 3) 1)))

(defun %fb-empty-cfun () (batten::%make-cfun 0 (make-real-vector 0) (make-real-vector 0)))

(defun %fb-empty-sfun ()
  (batten::%make-sfun 0 0 (make-real-vector 0) (make-real-vector 0) (make-real-matrix 0 0)))

(define-oracle-test crv-is-func-sized (:file "funcs-basic" :trials 500)
  (let* ((f (random-cfun))
         (nf (%fb-near (1- (length (cfun-vals f)))))
         (mk (%fb-near (1- (length (cfun-knots f))))))
    (agree (list :cfun f :nf nf :mk mk)
           (list (batten::crv-is-func-sized f nf mk 'crv-is-func-sized))
           (with-foreign ((fp f :cfun))
             (values (c-flag-error-p (c-crv-is-func-sized fp nf mk (cffi:null-pointer)))
                     (list nil))))))

(define-oracle-test c-value-from-array (:file "funcs-basic" :trials 500)
  (let* ((len (1+ (rnd 15)))
         (fu (random-real-vector len))
         (n (if (zerop (rnd 2)) (1- len) (rnd len))))
    (agree (list :fu fu :n n)
           (let ((cvl (batten::c-value-from-array fu n)))
             (list (1- (length cvl)) cvl))
           (with-foreign ((f fu :real-vector))
             (cffi:with-foreign-object (cvl '(:struct batten-oracle::nl-cvalue))
               (c-c-value-from-array cvl f n)
               (values nil (list (%fb-slot cvl nl-cvalue n) (%fb-cvalue cvl))))))))

(define-oracle-test c-func-from-knot-vector (:file "funcs-basic" :trials 500)
  (let ((f (random-cfun)) (p (1+ (rnd 7))))
    (agree (list :cfun f :p p)
           (list (batten::c-func-from-knot-vector (%fb-empty-cfun) (cfun-vals f) p (cfun-knots f)))
           (with-foreign ((fp f :cfun))
             (with-foreign-out ((g :cfun))
               (c-c-func-from-knot-vector g (%fb-slot fp nl-cfun cvl) p (%fb-slot fp nl-cfun knt))
               (values nil (list (from-foreign g :cfun))))))))

(define-oracle-test crv-func-cntrl-val (:file "funcs-basic" :trials 500)
  (let ((f (random-cfun)))
    (agree (list :cfun f)
           (multiple-value-list (batten::crv-func-cntrl-val f))
           (with-foreign ((fp f :cfun))
             (cffi:with-foreign-objects ((n :long) (fu :pointer))
               (c-crv-func-cntrl-val fp n fu)
               (let ((nn (cffi:mem-ref n :long)))
                 (values nil (list nn (from-foreign (cffi:mem-ref fu :pointer) :real-vector (1+ nn))))))))))

(define-oracle-test crv-func-cntrl-val-knot-vector (:file "funcs-basic" :trials 500)
  (let ((f (random-cfun)))
    (agree (list :cfun f)
           (multiple-value-list (batten::crv-func-cntrl-val-knot-vector f))
           (with-foreign ((fp f :cfun))
             (cffi:with-foreign-objects ((fu :pointer) (knt :pointer) (u :pointer))
               (c-crv-func-cntrl-val-knot-vector fp fu knt u)
               (let* ((n (%fb-slot (%fb-slot fp nl-cfun cvl) nl-cvalue n))
                      (k (cffi:mem-ref knt :pointer))
                      (m (%fb-slot k nl-knotvector m)))
                 (values nil (list (from-foreign (cffi:mem-ref fu :pointer) :real-vector (1+ n))
                                   (from-foreign k :knots)
                                   (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ m))))))))))

(define-oracle-test crv-func-cntrl-val-knots (:file "funcs-basic" :trials 500)
  (let ((f (random-cfun)))
    (agree (list :cfun f)
           (multiple-value-list (batten::crv-func-cntrl-val-knots f))
           (with-foreign ((fp f :cfun))
             (cffi:with-foreign-objects ((fu :pointer) (u :pointer))
               (c-crv-func-cntrl-val-knots fp fu u)
               (let ((n (%fb-slot (%fb-slot fp nl-cfun cvl) nl-cvalue n))
                     (m (%fb-slot (%fb-slot fp nl-cfun knt) nl-knotvector m)))
                 (values nil (list (from-foreign (cffi:mem-ref fu :pointer) :real-vector (1+ n))
                                   (from-foreign (cffi:mem-ref u :pointer) :real-vector (1+ m))))))))))

(define-oracle-test crv-func-set-ptrs (:file "funcs-basic" :trials 500)
  ;; the C sets only the pointers, keeping the old highest indexes, so the
  ;; new arrays are as long as the old ones
  (let* ((f1 (random-cfun))
         (p (cfun-degree f1))
         (n+1 (length (cfun-vals f1)))
         (f2 (batten::%make-cfun p (random-knots p n+1 :lo (- (rnd 4d0) 2d0) :hi (+ 3d0 (rnd 4d0)))
                                 (random-real-vector n+1))))
    (agree (list :f1 f1 :f2 f2)
           (list (batten::crv-func-set-ptrs (batten::copy-cfun f1) (cfun-vals f2) (cfun-knots f2)))
           (with-foreign ((a f1 :cfun) (b f2 :cfun))
             (c-crv-func-set-ptrs a (%fb-slot (%fb-slot b nl-cfun cvl) nl-cvalue fu)
                                  (%fb-slot (%fb-slot b nl-cfun knt) nl-knotvector u))
             (values nil (list (from-foreign a :cfun)))))))

(define-oracle-test srf-is-func-sized (:file "funcs-basic" :trials 500)
  (let* ((f (random-sfun))
         (vals (sfun-vals f))
         (nf (%fb-near (1- (array-dimension vals 0))))
         (mf (%fb-near (1- (array-dimension vals 1))))
         (rk (%fb-near (1- (length (sfun-knots-u f)))))
         (sk (%fb-near (1- (length (sfun-knots-v f))))))
    (agree (list :sfun f :nf nf :mf mf :rk rk :sk sk)
           (list (batten::srf-is-func-sized f nf mf rk sk 'srf-is-func-sized))
           (with-foreign ((fp f :sfun))
             (values (c-flag-error-p (c-srf-is-func-sized fp nf mf rk sk (cffi:null-pointer)))
                     (list nil))))))

(define-oracle-test s-value-from-array (:file "funcs-basic" :trials 500)
  (let* ((rows (1+ (rnd 8))) (cols (1+ (rnd 8)))
         (fuv (random-real-matrix rows cols))
         (whole (zerop (rnd 2)))
         (n (if whole (1- rows) (rnd rows)))
         (m (if whole (1- cols) (rnd cols))))
    (agree (list :fuv fuv :n n :m m)
           (let ((svl (batten::s-value-from-array fuv n m)))
             (list (1- (array-dimension svl 0)) (1- (array-dimension svl 1)) svl))
           (with-foreign ((f fuv :real-matrix))
             (cffi:with-foreign-object (svl '(:struct batten-oracle::nl-svalue))
               (c-s-value-from-array svl f n m)
               (values nil (list (%fb-slot svl nl-svalue n) (%fb-slot svl nl-svalue m)
                                 (%fb-svalue svl))))))))

(define-oracle-test s-func-from-knot-vectors (:file "funcs-basic" :trials 500)
  (let ((f (random-sfun)) (p (1+ (rnd 5))) (q (1+ (rnd 5))))
    (agree (list :sfun f :p p :q q)
           (list (batten::s-func-from-knot-vectors (%fb-empty-sfun) (sfun-vals f) p q
                                                   (sfun-knots-u f) (sfun-knots-v f)))
           (with-foreign ((fp f :sfun))
             (with-foreign-out ((g :sfun))
               (c-s-func-from-knot-vectors g (%fb-slot fp nl-sfun svl) p q
                                           (%fb-slot fp nl-sfun knu) (%fb-slot fp nl-sfun knv))
               (values nil (list (from-foreign g :sfun))))))))

(define-oracle-test srf-func-cntrl-val (:file "funcs-basic" :trials 500)
  (let ((f (random-sfun)))
    (agree (list :sfun f)
           (multiple-value-list (batten::srf-func-cntrl-val f))
           (with-foreign ((fp f :sfun))
             (cffi:with-foreign-objects ((n :long) (m :long) (fuv :pointer))
               (c-srf-func-cntrl-val fp n m fuv)
               (let ((nn (cffi:mem-ref n :long)) (mm (cffi:mem-ref m :long)))
                 (values nil (list nn mm (from-foreign (cffi:mem-ref fuv :pointer) :real-matrix
                                                       (1+ nn) (1+ mm))))))))))

(define-oracle-test s-func-get-knot-vectors (:file "funcs-basic" :trials 500)
  (let ((f (random-sfun)))
    (agree (list :sfun f)
           (multiple-value-list (batten::s-func-get-knot-vectors f))
           (with-foreign ((fp f :sfun))
             (cffi:with-foreign-objects ((knu :pointer) (knv :pointer))
               (c-s-func-get-knot-vectors fp knu knv)
               (values nil (list (from-foreign (cffi:mem-ref knu :pointer) :knots)
                                 (from-foreign (cffi:mem-ref knv :pointer) :knots))))))))

(define-oracle-test s-func-get-array-and-knot-vectors (:file "funcs-basic" :trials 500)
  (let ((f (random-sfun)))
    (agree (list :sfun f)
           (multiple-value-list (batten::s-func-get-array-and-knot-vectors f))
           (with-foreign ((fp f :sfun))
             (cffi:with-foreign-objects ((svl :pointer) (p :short) (q :short) (knu :pointer) (knv :pointer))
               (c-s-func-get-array-and-knot-vectors fp svl p q knu knv)
               (values nil (list (%fb-svalue (cffi:mem-ref svl :pointer))
                                 (cffi:mem-ref p :short) (cffi:mem-ref q :short)
                                 (from-foreign (cffi:mem-ref knu :pointer) :knots)
                                 (from-foreign (cffi:mem-ref knv :pointer) :knots))))))))
