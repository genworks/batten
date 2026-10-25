;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/geom.lisp against NLib.

(in-package :batten-tests)

(defun some-vector ()
  "A random vector, now and then zero or tiny."
  (case (rnd 12)
    (0 (make-point 0d0 0d0 0d0))
    (1 (make-point 1d-310 0d0 0d0))
    (t (random-point))))

(defmacro vector-to-vector-test (name lisp-fn c-fn)
  "Test a routine (A B C) writing C from vectors A and B."
  `(define-oracle-test ,name (:file "geom" :trials 3000)
     (let ((a (some-vector)) (b (some-vector)))
       (agree (list a b)
              (list (,lisp-fn a b (make-point)))
              (with-foreign ((ap a :point) (bp b :point))
                (with-foreign-out ((c :point))
                  (,c-fn ap bp c)
                  (values nil (list (from-foreign c :point)))))))))

(vector-to-vector-test vector-cross batten::vector-cross c-vector-cross)
(vector-to-vector-test vector-diff batten::vector-diff c-vector-diff)
(vector-to-vector-test vector-sum batten::vector-sum c-vector-sum)
(vector-to-vector-test vector-dir batten::vector-dir c-vector-dir)

(define-oracle-test vector-dot (:file "geom" :trials 3000)
  (let ((a (some-vector)) (b (some-vector)))
    (agree (list a b)
           (list (batten::vector-dot a b))
           (with-foreign ((ap a :point) (bp b :point))
             (with-scalars ((d :double))
               (c-vector-dot ap bp d)
               (values nil (list (cffi:mem-ref d :double))))))))

(define-oracle-test vector-magnitude (:file "geom" :trials 3000)
  (let ((a (some-vector)))
    (agree a
           (list (batten::vector-magnitude a))
           (with-foreign ((ap a :point))
             (with-scalars ((d :double))
               (c-vector-magnitude ap d)
               (values nil (list (cffi:mem-ref d :double))))))))

(define-oracle-test vector-scale (:file "geom" :trials 3000)
  (let ((a (some-vector)) (alpha (rnd-range -3d0 3d0)))
    (agree (list a alpha)
           (list (batten::vector-scale a alpha (make-point)))
           (with-foreign ((ap a :point))
             (with-foreign-out ((b :point))
               (c-vector-scale ap alpha b)
               (values nil (list (from-foreign b :point))))))))

(define-oracle-test vector-combine (:file "geom" :trials 3000)
  (let ((alf (rnd-range -3d0 3d0)) (a (some-vector)) (bet (rnd-range -3d0 3d0)) (b (some-vector)))
    (agree (list alf a bet b)
           (list (batten::vector-combine alf a bet b (make-point)))
           (with-foreign ((ap a :point) (bp b :point))
             (with-foreign-out ((c :point))
               (c-vector-combine alf ap bet bp c)
               (values nil (list (from-foreign c :point))))))))

(define-oracle-test vector-pt-along-vector (:file "geom" :trials 3000)
  (let ((p (random-point)) (alpha (rnd-range -3d0 3d0)) (v (some-vector)))
    (agree (list p alpha v)
           (list (batten::vector-pt-along-vector p alpha v (make-point)))
           (with-foreign ((pp p :point) (vp v :point))
             (with-foreign-out ((q :point))
               (c-vector-pt-along-vector pp alpha vp q)
               (values nil (list (from-foreign q :point))))))))

(define-oracle-test dist-pt-pt (:file "geom" :trials 3000)
  (let ((p (random-point)) (q (random-point)))
    (agree (list p q)
           (list (batten::dist-pt-pt p q))
           (with-foreign ((pp p :point) (qp q :point))
             (with-scalars ((d :double))
               (c-dist-pt-pt pp qp d)
               (values nil (list (cffi:mem-ref d :double))))))))

(define-oracle-test dist-cpt-cpt-homo (:file "geom" :trials 3000)
  (let ((p (random-cpoint :two-d (zerop (rnd 4)))) (q (random-cpoint :two-d (zerop (rnd 4)))))
    (agree (list p q)
           (list (batten::dist-cpt-cpt-homo p q))
           (with-foreign ((pp p :cpoint) (qp q :cpoint))
             (with-scalars ((d :double))
               (c-dist-cpt-cpt-homo pp qp d)
               (values nil (list (cffi:mem-ref d :double))))))))

(define-oracle-test vector-normalize (:file "geom" :trials 3000)
  (let ((a (some-vector)))
    (agree a
           (let ((b (make-point))) (list b (batten::vector-normalize a b)))
           (with-foreign ((ap a :point))
             (with-foreign-out ((b :point))
               (with-scalars ((mag :double))
                 (values (c-flag-error-p (c-vector-normalize ap b mag))
                         (list (from-foreign b :point) (cffi:mem-ref mag :double)))))))))

(define-oracle-test vector-normalize-ref (:file "geom" :trials 3000)
  (let ((a (some-vector)))
    (agree a
           (list (batten::vector-normalize-ref (copy-seq a)))
           (with-foreign ((ap a :point))
             (values (c-flag-error-p (c-vector-normalize-ref ap))
                     (list (from-foreign ap :point)))))))
