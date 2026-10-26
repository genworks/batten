;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0

(defun batten-around-compile (thunk)
  "Compile Batten's own files without the per-function records CCL and
Allegro keep by default -- source notes and local-symbol tables on CCL,
source-file and cross-reference information on Allegro -- which grew a saved
image by a third again (CCL +7.7 to +5.2 MB, Allegro +10.5 to +6.6 MB).  Bound
here, so a user's own compilations keep their settings.  Batten's functions
are then not found by M-. on Allegro; the source is plain to read."
  (let (#+ccl (ccl:*save-source-locations* nil)
        #+ccl (ccl:*save-local-symbols* nil)
        #+ccl (ccl:*fasl-save-local-symbols* nil)
        #+allegro (excl:*record-source-file-info* nil)
        #+allegro (excl:*record-xref-info* nil))
    (funcall thunk)))

(defsystem "batten"
  :description "NURBS curves and surfaces in portable Common Lisp: a port of NLib from NVIDIA's open-source SMLib."
  :version "0.0.1"
  :author "Genworks International"
  :license "Apache-2.0"
  :serial t
  :around-compile batten-around-compile
  :components ((:module "src"
                :serial t
                :components ((:file "package")
                             (:file "constants")
                             (:file "types")
                             ;; one file per NLib source file, infrastructure first,
                             ;; then the dependency waves of tools/manifest.sexp
                             (:module "nl"
                              :serial t
                              :components ((:file "math")
                                           (:file "geom")
                                           (:file "util")
                                           (:file "basis-basic")
                                           (:file "basis-adv")
                                           (:file "crv-basic")
                                           (:file "srf-basic")
                                           (:file "funcs-basic")
                                           (:file "frame-basic")
                                           (:file "crv-nurb")
                                           (:file "srf-nurb")
                                           (:file "bezier")
                                           ;; wave 1
                                           (:file "crv-tool")
                                           (:file "srf-tool")
                                           ;; wave 2
                                           (:file "funcs-adv")
                                           (:file "crv-common")
                                           (:file "srf-geom")
                                           (:file "crv-geom")
                                           ;; wave 3
                                           (:file "crv-fit")
                                           (:file "crv-adv")
                                           (:file "srf-common")
                                           ;; wave 4
                                           (:file "srf-fit")
                                           (:file "crv-approx")
                                           (:file "crv-offset")
                                           (:file "spiral")
                                           ;; wave 5
                                           (:file "srf-adv")
                                           (:file "srf-approx")
                                           ;; wave 6
                                           (:file "srf-offset")
                                           ;; level 14: the tessellator
                                           (:file "tessellate")))
                             ;; the readable API over the port (API.md)
                             (:module "api"
                              :serial t
                              :components ((:file "generics")
                                           (:file "evaluation")
                                           (:file "construction")
                                           (:file "curve-queries")
                                           (:file "curve-operations")
                                           (:file "surfaces")
                                           (:file "geometry")
                                           (:file "tessellation"))))))
  :in-order-to ((test-op (test-op "batten/tests"))))

;;; Random curves and surfaces from a seeded generator of its own, so every
;;; implementation draws the same ones.  Portable; no dependencies.
(defsystem "batten/random"
  :depends-on ("batten")
  :components ((:module "tests" :serial t :components ((:file "random") (:file "digest")))))

;;; The oracle: NLib itself (tools/build-oracle.sh), through CFFI, with
;;; bindings generated for every routine in tools/manifest.sexp.
(defsystem "batten/oracle"
  :depends-on ("batten" "cffi")
  :serial t
  :components ((:module "tests"
                :serial t
                :components ((:file "oracle")
                             (:file "oracle-bindings")
                             (:file "marshal")))))

;;; Differential tests: one file per NLib source file, registered with
;;; DEFINE-ORACLE-TEST and run by RUN-ALL.
(defsystem "batten/tests"
  :depends-on ("batten/random" "batten/oracle")
  :serial t
  :components ((:module "tests"
                :serial t
                :components ((:file "harness")
                             (:module "nl"
                              :serial t
                              :components ((:file "math") (:file "geom") (:file "util")
                                           (:file "basis-basic") (:file "basis-adv")
                                           (:file "crv-basic") (:file "srf-basic")
                                           (:file "funcs-basic") (:file "frame-basic")
                                           (:file "crv-nurb") (:file "srf-nurb") (:file "bezier")
                                           (:file "crv-tool") (:file "srf-tool")
                                           (:file "funcs-adv") (:file "crv-common")
                                           (:file "srf-geom") (:file "crv-geom")
                                           (:file "crv-fit") (:file "crv-adv") (:file "srf-common")
                                           (:file "srf-fit") (:file "crv-approx")
                                           (:file "crv-offset") (:file "spiral")
                                           (:file "srf-adv") (:file "srf-approx")
                                           (:file "srf-offset") (:file "tessellate")))
                             (:module "api"
                              :serial t
                              :components ((:file "construction") (:file "curve-queries")
                                           (:file "curve-operations") (:file "surfaces")
                                           (:file "geometry") (:file "tessellation"))))))
  :perform (test-op (o c) (uiop:symbol-call :batten-tests :run-all)))

;;; Timings of the Lisp routines alone: any implementation.
(defsystem "batten/bench"
  :depends-on ("batten/random")
  :components ((:module "bench" :components ((:file "bench")))))

;;; The same timings beside NLib's C, through the oracle.
(defsystem "batten/bench-c"
  :depends-on ("batten/bench" "batten/oracle")
  :components ((:module "bench" :components ((:file "bench-c")))))
