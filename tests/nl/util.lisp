;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Differential tests for src/nl/util.lisp against NLib.

(in-package :batten-tests)

(define-oracle-test xyz-to-2d-c-pt-array (:file "util" :trials 1000)
  (let* ((n (rnd 8)) (m (rnd 8))
         (wx (random-real-matrix (1+ n) (1+ m)))
         (wy (random-real-matrix (1+ n) (1+ m)))
         (wz (if (zerop (rnd 5))
                 (make-array (list (1+ n) (1+ m)) :element-type 'double-float :initial-element +no-z+)
                 (random-real-matrix (1+ n) (1+ m))))
         (w (if (zerop (rnd 2))
                (make-array (list (1+ n) (1+ m)) :element-type 'double-float :initial-element +no-w+)
                (random-real-matrix (1+ n) (1+ m) 5d0))))
    (agree (list :wx wx :wy wy :wz wz :w w :n n :m m)
           (list (batten::xyz-to-2d-c-pt-array wx wy wz w n m))
           (with-stacks (s)
             (with-foreign ((fx wx :real-matrix) (fy wy :real-matrix) (fz wz :real-matrix)
                            (fw w :real-matrix))
               (let ((pw (c-xyz-to-2d-c-pt-array fx fy fz fw n m s)))
                 (values (cffi:null-pointer-p pw)
                         (unless (cffi:null-pointer-p pw)
                           (list (from-foreign pw :cpoint-net (1+ n) (1+ m)))))))))))

(defun %file-length-or-nil (path)
  (with-open-file (in path :direction :input :if-does-not-exist nil
                           :element-type '(unsigned-byte 8))
    (and in (file-length in))))

;;; N_FileClose closes a stream N_FileOpen opened: both sides open a fresh
;;; file for writing, close it, and report that the file exists, empty.
(define-oracle-test file-close (:file "util" :trials 50)
  (let ((lisp-path (format nil "/tmp/batten-file-close-lisp-~d" (rnd 1000000000)))
        (c-path (format nil "/tmp/batten-file-close-c-~d" (rnd 1000000000))))
    (unwind-protect
         (agree (list :lisp-path lisp-path :c-path c-path)
                (let ((s (open lisp-path :direction :output :if-exists :supersede)))
                  (batten::file-close s)
                  (list (not (open-stream-p s)) (%file-length-or-nil lisp-path)))
                (cffi:with-foreign-strings ((fname c-path) (mode "w"))
                  (let ((fp (c-file-open fname mode)))
                    (if (cffi:null-pointer-p fp)
                        (values t nil)
                        (progn (c-file-close fp)
                               (values nil (list t (%file-length-or-nil c-path))))))))
      (when (probe-file lisp-path) (delete-file lisp-path))
      (when (probe-file c-path) (delete-file c-path)))))

;;; N_ErrClear resets NLib's error record; Batten keeps none (errors are
;;; conditions), so the Lisp side has nothing pending after ERR-CLEAR, and
;;; the C side must report no error after N_ErrSet then N_ErrClear.
(define-oracle-test err-clear (:file "util" :trials 50)
  (let ((code (1+ (rnd 20))))
    (agree (list :code code)
           (list (batten::err-clear) 0)
           (cffi:with-foreign-string (rname "err-clear-test")
             (c-err-set code rname)
             (c-err-clear)
             (values nil (list nil (c-err-get-type)))))))
