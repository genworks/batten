;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; The differential test harness.  A test is one random trial, run TRIALS
;;;; times; inside it, AGREE runs the Lisp port and the C original on the same
;;;; input and records whether they agree: bit-identical results, or both
;;;; failing (the Lisp signals NURBS-ERROR exactly when the C returns NL_YES).
;;;;
;;;;   (define-oracle-test crv-eval (:file "crv-nurb" :trials 2000)
;;;;     (let* ((c (random-curve)) (u (random-parameter (curve-knots c))) (flg (random-side)))
;;;;       (agree (list :curve c :u u :flg flg)
;;;;              (list (batten::crv-eval c u flg (make-point)))
;;;;              (with-foreign ((cur c :curve))
;;;;                (with-foreign-out ((pt :point))
;;;;                  (values (c-flag-error-p (c-crv-eval cur u flg pt))
;;;;                          (list (from-foreign pt :point))))))))
;;;;
;;;; The Lisp form returns a list of results; the C form returns
;;;; (values error-p list-of-results).  Results compare element by element with
;;;; SAME: numbers by EQL (bits, not tolerance), arrays element-wise, Batten
;;;; structures slot by slot, lists recursively.

(defpackage :batten-tests
  (:use :common-lisp :batten :batten-oracle :batten-random)
  (:export #:run-all #:define-oracle-test #:agree #:same #:coverage-report
           #:random-side #:*verbose*))

(in-package :batten-tests)

;;; The generated C-<lisp name> bindings live in BATTEN-ORACLE unexported (there
;;; are hundreds); every test file sees them through this import.
(eval-when (:compile-toplevel :load-toplevel :execute)
  (do-symbols (s :batten-oracle)
    (when (and (eq (symbol-package s) (find-package :batten-oracle))
               (> (length (symbol-name s)) 2)
               (string= "C-" (symbol-name s) :end2 2)
               (fboundp s))
      (import s :batten-tests))))

(defvar *tests* (make-hash-table :test 'equal) "name -> (file trials function)")
(defvar *tally* nil "The tally of the test now running.")
(defvar *verbose* nil)

(defstruct tally (name "") (file "") (trials 0) (identical 0) (both-failed 0)
  (mismatch 0) (lisp-only-failed 0) (c-only-failed 0) (max-rel 0d0) (first-mismatch nil)
  (crashed nil))

(defmacro define-oracle-test (name (&key file (trials 2000)) &body body)
  "Register a differential test NAME for NLib source FILE (a stem such as
\"crv-nurb\"); BODY is one random trial."
  `(setf (gethash ,(string-downcase (string name)) *tests*)
         (list ,file ,trials (lambda () ,@body))))

(defun rel-diff (a b)
  (if (and (realp a) (realp b))
      (let ((d (abs (- a b))) (s (max (abs a) (abs b))))
        (if (zerop s) 0d0 (float (/ d s) 1d0)))
      1d0))

(defun same (a b)
  "(values identical-p max-relative-difference) for two results."
  (cond
    ((and (floatp a) (floatp b))
     (if (or (eql a b) (and (/= a a) (/= b b))) (values t 0d0) (values nil (rel-diff a b))))
    ((and (numberp a) (numberp b)) (values (= a b) (if (= a b) 0d0 (rel-diff a b))))
    ((and (arrayp a) (arrayp b) (not (stringp a)))
     (if (not (equal (array-dimensions a) (array-dimensions b)))
         (values nil 1d0)
         (let ((ok t) (worst 0d0))
           (dotimes (i (array-total-size a) (values ok worst))
             (multiple-value-bind (s r) (same (row-major-aref a i) (row-major-aref b i))
               (unless s (setf ok nil worst (max worst r))))))))
    ((and (consp a) (consp b))
     (multiple-value-bind (s1 r1) (same (car a) (car b))
       (multiple-value-bind (s2 r2) (same (cdr a) (cdr b))
         (values (and s1 s2) (max r1 r2)))))
    ((and (typep a 'structure-object) (typep b 'structure-object) (eq (type-of a) (type-of b)))
     (same (structure-slots a) (structure-slots b)))
    (t (values (equal a b) (if (equal a b) 0d0 1d0)))))

(defun structure-slots (s)
  "The slot values of a Batten structure, in order."
  (etypecase s
    (curve (list (curve-degree s) (curve-knots s) (curve-cpoints s)))
    (surface (list (surface-degree-u s) (surface-degree-v s) (surface-knots-u s)
                   (surface-knots-v s) (surface-cpoints s)))
    (cfun (list (cfun-degree s) (cfun-knots s) (cfun-vals s)))
    (sfun (list (sfun-degree-u s) (sfun-degree-v s) (sfun-knots-u s) (sfun-knots-v s) (sfun-vals s)))
    (batten::lineseg (list (batten::lineseg-p s) (batten::lineseg-v s) (batten::lineseg-bounded s)))
    (batten::plane (list (batten::plane-p s) (batten::plane-n s)))))

(defmacro agree (case lisp-form c-form)
  "Run LISP-FORM (a list of results; may signal NURBS-ERROR) and C-FORM
((values error-p results)) and record their agreement in the running tally.
CASE describes the input, kept for the first mismatch."
  `(let ((lisp-results nil) (lisp-failed nil))
     (handler-case (setf lisp-results ,lisp-form)
       (nurbs-error () (setf lisp-failed t)))
     (multiple-value-bind (c-failed c-results) ,c-form
       (record-agreement lisp-failed lisp-results c-failed c-results (lambda () ,case)))))

(defun record-agreement (lisp-failed lisp-results c-failed c-results case-thunk)
  (let ((tl *tally*))
    (incf (tally-trials tl))
    (flet ((note-mismatch (rel)
             (incf (tally-mismatch tl))
             (setf (tally-max-rel tl) (max (tally-max-rel tl) rel))
             (unless (tally-first-mismatch tl)
               (setf (tally-first-mismatch tl)
                     (list :case (funcall case-thunk) :lisp-failed lisp-failed :c-failed c-failed
                           :lisp lisp-results :c c-results)))))
      (cond ((and lisp-failed c-failed) (incf (tally-both-failed tl)))
            (lisp-failed (incf (tally-lisp-only-failed tl)) (note-mismatch 1d0))
            (c-failed (incf (tally-c-only-failed tl)) (note-mismatch 1d0))
            (t (multiple-value-bind (ok rel) (same lisp-results c-results)
                 (if ok (incf (tally-identical tl)) (note-mismatch rel))))))))

(defun random-side () (if (zerop (rnd 2)) +left+ +right+))

(defun run-test (name &key trials)
  (destructuring-bind (file default-trials fn) (gethash name *tests*)
    (let ((*tally* (make-tally :name name :file file)))
      (handler-case
          (dotimes (i (or trials default-trials))
            (funcall fn))
        (error (e) (setf (tally-crashed *tally*) (princ-to-string e))))
      *tally*)))

(defun report (tl)
  (format t "~&~12a ~36a ~6d trials ~6d same ~4d both-fail ~4d MISMATCH~@[ (max rel ~,2e)~]~@[  CRASHED: ~a~]~%"
          (tally-file tl) (tally-name tl) (tally-trials tl) (tally-identical tl)
          (tally-both-failed tl) (tally-mismatch tl)
          (when (plusp (tally-mismatch tl)) (tally-max-rel tl))
          (tally-crashed tl))
  (when (and *verbose* (tally-first-mismatch tl))
    (let ((*print-length* 12) (*print-level* 4))
      (format t "   first mismatch: ~s~%" (tally-first-mismatch tl))))
  (and (zerop (tally-mismatch tl)) (not (tally-crashed tl)) (plusp (tally-trials tl))))

(defun run-all (&key file name trials (seed 20261025))
  "Run the registered tests -- all, those of NLib source FILE, or the one NAME
-- and print a table.  True when every test agreed with NLib on every trial."
  (let ((names (sort (let (ns)
                       (maphash (lambda (k v) (when (and (or (null file) (equal (first v) file))
                                                         (or (null name) (equal k (string-downcase (string name)))))
                                                (push k ns)))
                                *tests*)
                       ns)
                     #'string< :key (lambda (n) (concatenate 'string (first (gethash n *tests*)) "/" n)))))
    (let ((results (mapcar (lambda (n) (reseed (+ seed (sxhash n))) (run-test n :trials trials)) names)))
      (let ((ok (every #'identity (mapcar #'report results))))
        (format t "~&~d tests, ~:[SOME DISAGREE~;all agree with NLib~]~%" (length results) ok)
        (values ok results)))))

;;; ------------------------------------------------------------------ coverage

(defun read-manifest (&optional (path (asdf:system-relative-pathname "batten" "tools/manifest.sexp")))
  (with-open-file (s path) (let ((*read-eval* nil)) (read s))))

(defun coverage-report (&key file (stream t))
  "For each manifest routine (optionally of one NLib FILE stem), whether its Lisp
name is defined and whether a test of that name is registered.  Returns
(values defined-count tested-count to-port-count missing-list)."
  (let ((missing nil) (defined 0) (tested 0) (to-port 0))
    (dolist (e (read-manifest))
      (let* ((status (getf e :status))
             (lisp (getf e :lisp))
             (sym (find-symbol (string-upcase lisp) :batten)))
        (unless (eq status :vanish)
          (when (or (null file) (equal file (manifest-file-stem e)))
            (incf to-port)
            (let ((is-defined (and sym (fboundp sym)))
                  (is-tested (gethash lisp *tests*)))
              (when is-defined (incf defined))
              (when is-tested (incf tested))
              (unless is-defined (push (list (getf e :c) lisp (getf e :wave)) missing)))))))
    (format stream "~&~d of ~d routines defined, ~d with a test of the same name~%" defined to-port tested)
    (values defined tested to-port (nreverse missing))))

(defun manifest-file-stem (e)
  "\"NL_CrvNurb.cpp\" -> \"crv-nurb\"."
  (let* ((f (getf e :file))
         (base (subseq f 3 (position #\. f)))
         (out (make-string-output-stream)))
    (dotimes (i (length base))
      (let ((c (char base i)))
        (when (and (> i 0) (upper-case-p c)
                   (or (lower-case-p (char base (1- i))) (digit-char-p (char base (1- i)))))
          (write-char #\- out))
        (write-char (char-downcase c) out)))
    (get-output-stream-string out)))
