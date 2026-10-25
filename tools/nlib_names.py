# SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
# SPDX-License-Identifier: Apache-2.0
"""The mechanical NLib -> Lisp naming rule, shared by the oracle generator.
The same rule is BATTEN::NLIB-LISP-NAME in src/naming.lisp; the two must agree.

  N_CrvEval              -> crv-eval
  N_Combine2CPts         -> combine-2-c-pts
  N_CPtToPtEuclid        -> c-pt-to-pt-euclid
  N_RealMatrixLuDecompose-> real-matrix-lu-decompose
  N_XYZTo2dCPtArray      -> xyz-to-2d-c-pt-array
  ST_findStep in NL_CrvNurb.cpp -> crv-nurb--st-find-step   (file-private)
"""
import re


def hyphenate(s):
    out = []
    for i, c in enumerate(s):
        prev = s[i - 1] if i else ""
        nxt = s[i + 1] if i + 1 < len(s) else ""
        if i and c == "_":
            out.append("-")
            continue
        if i and c.isupper() and (prev.islower() or prev.isdigit()):
            out.append("-")
        elif i and c.isupper() and prev.isupper() and nxt.islower():
            out.append("-")
        elif i and c.isdigit() and prev.isalpha():
            out.append("-")
        out.append(c.lower())
    return re.sub(r"-+", "-", "".join(out)).strip("-")


def file_stem(nlib_file):
    """NL_CrvNurb.cpp -> crv-nurb"""
    base = nlib_file.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    return hyphenate(re.sub(r"^NL_", "", base))


def lisp_name(c_name, nlib_file=None, static=False):
    if c_name.startswith("ST_") or static:
        rest = c_name[3:] if c_name.startswith("ST_") else re.sub(r"^N_", "", c_name)
        return f"{file_stem(nlib_file)}--st-{hyphenate(rest)}"
    for prefix in ("NI_", "N_"):
        if c_name.startswith(prefix):
            return hyphenate(c_name[len(prefix):])
    return hyphenate(c_name)


if __name__ == "__main__":
    for n in ["N_CrvEval", "N_Combine2CPts", "N_CPtToPtEuclid", "N_RealMatrixLuDecompose",
              "N_XYZTo2dCPtArray", "N_SrfGetBBox", "N_Crv2dTo3d", "N_IsCrvRat", "N_BasisFindSpan"]:
        print(f"{n:28} {lisp_name(n)}")
    print("ST_findStep NL_CrvNurb.cpp ", lisp_name("ST_findStep", "NL_CrvNurb.cpp"))
