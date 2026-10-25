#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
# SPDX-License-Identifier: Apache-2.0
#
# run-tests.sh [FILE-STEM|all] [TRIALS] [TEST-NAME]
#
# Run Batten's differential tests in a fresh, throwaway SBCL against the NLib
# oracle, so that any number of these can run side by side without sharing a
# Lisp image.  FILE-STEM picks one NLib file's tests ("crv-nurb"); TRIALS
# overrides each test's trial count; TEST-NAME runs one test.  Exit status 0
# when every test agrees with NLib.  Also prints the manifest coverage for
# FILE-STEM.
#
# BATTEN_ORACLE_DIR  where tools/build-oracle.sh put the oracle libraries
# BATTEN_TEST_IMAGE  an SBCL image with CFFI (default: the Gendl SBCL image)
set -euo pipefail
FILE=${1:-all}
TRIALS=${2:-nil}
NAME=${3:-}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
ORACLE=${BATTEN_ORACLE_DIR:-/home/dcooper8/projects/tmp-shots/batten-scratch/oracle}
IMAGE=${BATTEN_TEST_IMAGE:-gornskew/gendl:devo-sbcl}
if [ "$FILE" = all ]; then FILEARG=nil; else FILEARG="\"$FILE\""; fi
if [ -n "$NAME" ]; then NAMEARG="\"$NAME\""; else NAMEARG=nil; fi

docker run --rm -v "$ROOT":/batten:ro -v "$ORACLE":/oracle:ro \
  -e BATTEN_ORACLE_DIR=/oracle/ \
  --entrypoint /opt/gendl/gdl/program/gdl-sbcl "$IMAGE" --non-interactive --eval "
(progn
  (asdf:register-preloaded-system \"cffi\")
  (push #p\"/batten/\" asdf:*central-registry*)
  (let ((warnings nil))
    (handler-bind ((warning (lambda (w)
                              (unless (typep w 'style-warning) (push (princ-to-string w) warnings))
                              (muffle-warning w))))
      (asdf:load-system \"batten/tests\"))
    (when warnings
      (format t \"~&COMPILER WARNINGS (~d):~%~{  ~a~%~}\" (length warnings) warnings)))
  (uiop:symbol-call :batten-oracle :load-oracle)
  (format t \"~&=== batten tests ===~%\")
  (let ((*print-length* 8) (*print-level* 4))
    (setf (symbol-value (find-symbol \"*VERBOSE*\" :batten-tests)) t)
    (let ((ok (uiop:symbol-call :batten-tests :run-all :file $FILEARG :trials $TRIALS :name $NAMEARG)))
      (uiop:symbol-call :batten-tests :coverage-report :file $FILEARG)
      (sb-ext:exit :code (if ok 0 1)))))" 2>&1 |
  # drop the image's banner and the compiler's chatter up to the marker, but
  # keep the lines around any error before it (a failed load)
  awk '/^=== batten tests ===/ {on=1; next}
       on {print; next}
       /COMPILER WARNINGS|caught (ERROR|WARNING)|^Unhandled|READ error/ {show=14}
       show > 0 {print; show--}'
exit "${PIPESTATUS[0]}"
