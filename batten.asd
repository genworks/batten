;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0

(defsystem "batten"
  :description "NURBS curves and surfaces in portable Common Lisp: a port of NLib from NVIDIA's open-source SMLib."
  :version "0.0.1"
  :author "Genworks International"
  :license "Apache-2.0"
  :serial t
  :components ((:module "src"
                :serial t
                :components ((:file "package")
                             (:file "types")
                             (:file "math")
                             (:file "basis")
                             (:file "curve"))))
  :in-order-to ((test-op (test-op "batten/tests"))))

;;; Random curves from a seeded generator of its own, so every implementation
;;; draws the same ones.  Portable; no dependencies.
(defsystem "batten/random"
  :depends-on ("batten")
  :components ((:module "tests" :serial t :components ((:file "random") (:file "digest")))))

;;; The oracle: NLib itself, the C functions each port came from, loaded from
;;; an SMLib 10 build through CFFI.  Development and testing only.
(defsystem "batten/oracle"
  :depends-on ("batten" "cffi")
  :components ((:module "tests" :components ((:file "oracle")))))

(defsystem "batten/tests"
  :depends-on ("batten/random" "batten/oracle")
  :components ((:module "tests" :components ((:file "differential"))))
  :perform (test-op (o c) (uiop:symbol-call :batten-tests :run-all)))

;;; Timings of the Lisp routines alone: runs on any implementation, and the
;;; numbers compare across implementations because the curves are the same.
(defsystem "batten/bench"
  :depends-on ("batten/random")
  :components ((:module "bench" :components ((:file "bench")))))

;;; The same timings beside NLib's C, through the oracle.
(defsystem "batten/bench-c"
  :depends-on ("batten/bench" "batten/oracle")
  :components ((:module "bench" :components ((:file "bench-c")))))
