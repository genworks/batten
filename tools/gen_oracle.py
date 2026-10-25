# SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
# SPDX-License-Identifier: Apache-2.0
"""gen_oracle.py SMLIB_SRC MANIFEST ORACLE_SO CPP_SYMBOLS OUT_SHIMS_CPP OUT_BINDINGS_LISP

For every non-static routine in the manifest that the oracle library exports,
write a shim with C linkage that takes pointers where the routine takes
structs by value (CFFI cannot pass structs by value without libffi), and a
CFFI binding for the shim: BATTEN-ORACLE::C-<lisp-name>.  Each binding
carries the routine's C signature as a comment.
"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from nlib_names import hyphenate  # noqa: E402

src, manifest, oracle_so, cpp_symbols, out_cpp, out_lisp = sys.argv[1:7]

entries = []
for m in re.finditer(r'\(:c "(\w+)" :lisp "([^"]+)" :file "([^"]+)" :line (\d+)(.*?)\n', open(manifest).read()):
    entries.append((m.group(1), m.group(2), m.group(3), int(m.group(4)), ":static t" in m.group(5)))

exported = set(l.split()[-1] for l in subprocess.run(
    ["nm", "-D", "--defined-only", oracle_so], capture_output=True, text=True).stdout.splitlines()
    if " T " in l)
cpp_linkage = dict(re.findall(r'\("(\w+)" \. "(\w+)"\)', open(cpp_symbols).read()))

BYVAL = {"NL_POINT", "NL_VECTOR", "NL_CPOINT", "NL_CVECTOR", "NL_INTERVAL", "NL_RECTANGLE",
         "NL_MINMAXBOX", "NL_LINESEG", "NL_PLANE"}
SCALAR = {"NL_REAL": ":double", "NL_PARAMETER": ":double", "double": ":double",
          "NL_INDEX": ":long", "NL_INTEGER": ":long", "NL_BOOLEAN": ":long", "long": ":long",
          "NL_FLAG": ":short", "NL_DEGREE": ":short", "short": ":short",
          "int": ":int", "ULONG": ":unsigned-long", "unsigned long": ":unsigned-long",
          "NL_MATRIXTYPE": ":int", "float": ":float"}
RETURN = dict(SCALAR, NL_VOID=":void", void=":void")


def strip_comments(s):
    s = re.sub(r'/\*.*?\*/', ' ', s, flags=re.S)
    return re.sub(r'//[^\n]*', ' ', s)


def signature(path, line):
    L = open(path, encoding="latin-1").read().split("\n")
    text = strip_comments("\n".join(L[line - 1: line + 60]))
    start = text.index("(")
    depth, i = 0, start
    while True:
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    head, params = text[:start], text[start + 1:i]
    head = re.sub(r'\b(static|NL_INLINE|NL_NOINLINE|GW_EXPORT|inline)\b', ' ', head)
    head = " ".join(head.split())
    ret_type, name = head.rsplit(" ", 1) if " " in head else ("int", head)
    if name.startswith("*"):
        ret_type, name = ret_type + " *", name[1:]
    ps = []
    depth, cur = 0, ""
    for ch in params:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            ps.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        ps.append(cur)
    out = []
    for p in ps:
        p = " ".join(p.split())
        if p in ("void", "NL_VOID", ""):
            continue
        m = re.match(r'^(.*?)([A-Za-z_]\w*)\s*((?:\[[^\]]*\])*)$', p)
        ptype, pname, arr = m.group(1).strip(), m.group(2), m.group(3)
        if not ptype:            # unnamed parameter
            ptype, pname = pname, f"arg{len(out)}"
        out.append((ptype + (" *" if arr else ""), pname, p))
    return ret_type.strip(), name, out


def base(t):
    return " ".join(w for w in t.replace("*", " ").replace("&", " ").split() if w != "const")


def is_ptr(t):
    return "*" in t or "&" in t


shims, bindings, skipped = [], [], []
for c_name, lisp, fname, line, static in entries:
    if static or (c_name not in exported and c_name not in cpp_linkage):
        skipped.append(c_name)
        continue
    try:
        ret, name, params = signature(os.path.join(src, fname), line)
    except Exception as e:  # noqa: BLE001
        skipped.append(f"{c_name}(parse: {e})")
        continue
    if name != c_name:
        skipped.append(f"{c_name}(found {name})")
        continue
    rbase = base(ret)
    ret_byval = (not is_ptr(ret)) and rbase in BYVAL
    if not is_ptr(ret) and not ret_byval and rbase not in RETURN:
        skipped.append(f"{c_name}(return {ret})")
        continue
    sparams, call_args, lparams, bad = [], [], [], None
    for i, (ptype, pname, orig) in enumerate(params):
        b = base(ptype)
        if is_ptr(ptype):
            sparams.append(f"void * a{i}")
            # references and pointers alike arrive as void *; the adapters let
            # the compiler deduce the exact parameter type (arrays included)
            call_args.append(f"oracle_ref{{a{i}}}" if "&" in ptype else f"oracle_ptr{{a{i}}}")
            lparams.append((pname, ":pointer"))
        elif b in ("NL_STRING", "TCHAR *"):
            # NLib's routine-name argument: a string pointer, passed through
            sparams.append(f"NL_STRING a{i}")
            call_args.append(f"a{i}")
            lparams.append((pname, ":pointer"))
        elif b in BYVAL:
            sparams.append(f"{b} * a{i}")
            call_args.append(f"*a{i}")
            lparams.append((pname, ":pointer"))
        elif b in SCALAR:
            sparams.append(f"{b} a{i}")
            call_args.append(f"a{i}")
            lparams.append((pname, SCALAR[b]))
        else:
            bad = f"param {orig}"
            break
    if bad:
        skipped.append(f"{c_name}({bad})")
        continue
    proto = ""
    if c_name in cpp_linkage:
        proto = f"{ret} {c_name}({', '.join(p[2] for p in params)});\n"
    if ret_byval:
        sparams.append(f"{rbase} * result")
        body = f"*result = {c_name}({', '.join(call_args)});"
        sret, lret = "void", ":void"
        lparams.append(("result", ":pointer"))
    else:
        sret = "void *" if is_ptr(ret) else ret
        body = (f"return (void *) {c_name}({', '.join(call_args)});" if is_ptr(ret)
                else (f"{c_name}({', '.join(call_args)});" if RETURN.get(rbase) == ":void"
                      else f"return {c_name}({', '.join(call_args)});"))
        lret = ":pointer" if is_ptr(ret) else RETURN[rbase]
    shims.append(f'{proto}extern "C" {sret} shim_{c_name}({", ".join(sparams) or "void"}) {{ {body} }}\n')
    seen, lp = {}, []
    for pname, ltype in lparams:
        s = hyphenate(pname)
        if s in ("t", "nil", "pi"):
            s += "*"
        while s in seen:
            s += "*"
        seen[s] = 1
        lp.append(f"({s} {ltype})")
    csig = f"{ret} {c_name}({', '.join(p[2] for p in params)})"
    bindings.append(f';; {csig}\n(cffi:defcfun ("shim_{c_name}" c-{lisp}) {lret}' +
                    (f'\n  {" ".join(lp)})\n' if lp else ")\n"))

with open(out_cpp, "w") as o:
    o.write("// Generated by tools/gen_oracle.py: C-linkage, pointer-only shims over NLib.\n")
    o.write("#include <nurbs.h>\n#include <NL_Globals.h>\n\n")
    o.write("struct oracle_ptr { void *p; template <class T> operator T *() const { return (T *) p; } };\n")
    o.write("struct oracle_ref { void *p; template <class T> operator T &() const { return *(T *) p; } };\n\n")
    o.writelines(shims)
with open(out_lisp, "w") as o:
    o.write(";;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International\n")
    o.write(";;;; SPDX-License-Identifier: Apache-2.0\n;;;;\n")
    o.write(";;;; CFFI bindings to NLib through the oracle's shims, one per routine in\n")
    o.write(";;;; tools/manifest.sexp: C-<lisp name>.  Structs passed by value in C are\n")
    o.write(";;;; passed by pointer here.  Generated by tools/gen_oracle.py.\n\n")
    o.write("(in-package :batten-oracle)\n\n")
    o.write("\n".join(bindings))
print(f"{len(shims)} shims/bindings; skipped {len(skipped)}")
for s in skipped:
    if "(" in s:
        print("  skipped", s)
