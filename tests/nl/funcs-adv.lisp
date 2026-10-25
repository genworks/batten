;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/funcs-adv.lisp against NLib.

(in-package :batten-tests)

(defun %fa-copy-cfun (f)
  "A copy of curve function F with its own arrays."
  (%make-cfun (cfun-degree f) (copy-seq (cfun-knots f)) (copy-seq (cfun-vals f))))

(defun %fa-insert-knots (knots p)
  "Knots to insert into KNOTS (degree P): sorted, strictly inside the domain,
some existing interior knots and some repeated, every resulting multiplicity
at most P.  Now and then an end knot or one outside (:par-err), or none at
all (:inp-err)."
  (let* ((lo (aref knots 0)) (hi (aref knots (1- (length knots))))
         (interior (remove-if (lambda (x) (or (= x lo) (= x hi))) (coerce knots 'list)))
         (cands nil))
    (dotimes (i (1+ (rnd 7)))
      (push (case (rnd 6)
              (0 (if interior (nth (rnd (length interior)) interior) (rnd-range lo hi)))
              (1 (if cands (first cands) (rnd-range lo hi)))
              (t (rnd-range lo hi)))
            cands))
    (setf cands (sort (remove-if (lambda (x) (or (<= x lo) (>= x hi))) cands) #'<))
    (let ((accepted nil))
      (dolist (x cands)
        (when (<= (+ (count x knots :test #'=) (count x accepted :test #'=) 1) p)
          (push x accepted)))
      (setf accepted (nreverse accepted))
      (case (rnd 25)
        (0 (setf accepted nil))
        (1 (push lo accepted))
        (2 (setf accepted (append accepted (list hi))))
        (3 (setf accepted (append accepted (list (+ hi 1d0))))))
      (coerce accepted 'real-vector))))

(define-oracle-test crv-func-refine (:file "funcs-adv" :trials 4000)
  (let* ((f (random-cfun))
         (knx (%fa-insert-knots (cfun-knots f) (cfun-degree f)))
         (in-place (zerop (rnd 2))))
    (agree (list :cfun f :knx knx :in-place in-place)
           (list (batten::crv-func-refine f knx))
           (with-foreign ((fp f :cfun) (kx knx :knots))
             (with-stacks (sp)
               (if in-place
                   (if (c-flag-error-p (c-crv-func-refine fp kx fp sp sp))
                       (values t nil)
                       (values nil (list (from-foreign fp :cfun))))
                   (with-stacks (sq)
                     (with-foreign-out ((fq :cfun))
                       (if (c-flag-error-p (c-crv-func-refine fp kx fq sp sq))
                           (values t nil)
                           (values nil (list (from-foreign fq :cfun))))))))))))

(define-oracle-test crv-func-deriv (:file "funcs-adv" :trials 4000)
  (let* ((f (random-cfun))
         (p (cfun-degree f))
         (d (case (rnd 20) (0 -1) (1 (1+ p)) (t (rnd (1+ p))))))
    (agree (list :cfun f :d d)
           (list (batten::crv-func-deriv f d))
           (with-foreign ((fp f :cfun))
             (with-stacks (sg)
               (with-foreign-out ((der :cfun))
                 (if (c-flag-error-p (c-crv-func-deriv fp d der sg))
                     (values t nil)
                     (values nil (list (from-foreign der :cfun))))))))))

(define-oracle-test crv-func-multiply-constant (:file "funcs-adv" :trials 2000)
  (let ((f (random-cfun))
        (alpha (case (rnd 8) (0 0d0) (1 -1d0) (t (rnd-range -10d0 10d0)))))
    (agree (list :cfun f :alpha alpha)
           (list (batten::crv-func-multiply-constant alpha (%fa-copy-cfun f)))
           (with-foreign ((fp f :cfun))
             (c-crv-func-multiply-constant alpha fp)
             (values nil (list (from-foreign fp :cfun)))))))
