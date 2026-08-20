#!/usr/bin/env python3
"""check_numbers.py - is every number printed in the SU(7) paper greppable in an archived run?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

The paper's Data availability section CLAIMS every displayed number regenerates from the ancillary
scripts and is literally present in the archived stdout beside them.  This is what makes that claim
true rather than decorative.  Same design as Part V's verifier, and the same two lessons burnt into
it: keep the decimal point and match on token boundaries (a gate that matches substrings passed
`1.816` on the strength of `816`), and treat an ESCAPED percent as content, because the naive
r"%.*" silently swallows every number to the right of a printed percentage.

Run:  python check_numbers.py     (from part_vi/paper/)
"""
import os
import re

TEXS = ["su7_proton_row.tex", "su7_proton_row_es.tex"]
# The paper quotes Part V for the instrument's external validation, so the archive
# it is checked against has to include Part V's runs too -- a gate that cannot see
# the run behind a number reports that number as unbacked when it is not.
OUT = ["../outputs", "../../part_v/outputs"]

ALLOW = {
    # --- drawing parameters of fig:twoloop, not measurements: TikZ opacities
    "0.45",
    # --- an equation number of a DIFFERENT paper: GuoDu's eq. (5.13)
    "5.13",
    # --- their own equation numbers, quoted throughout
    "11", "12", "13", "31", "37", "40", "41", "42", "43", "44", "46", "47",
    "53", "54", "55", "56", "57", "58", "62", "63", "64", "65", "66", "67",
    "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79",
    "80", "81", "82", "19",
    # --- identifiers: arXiv, journals, volumes, pages, years
    "2503.04090", "2312.08608", "1001.0768", "0304220", "1903.08359",
    "1804.06012", "1105.0541", "1704.04840", "0108049", "0204037",
    "0012092", "240",
    "1979", "1983", "1989", "1992", "2008", "2018", "1566", "1571", "1221",
    "1566", "62", "368", "43", "39", "693", "065", "015022", "158", "126",
    "190", "233", "309", "141", "1959", "07", "98", "3",
    # --- quoted from THEIR paper or from the PDG, not measured by us
    "125", "127", "80.4", "0.63", "2.0", "7.5", "3.8", "6.1",
    "126.8", "125.5", "125.1", "126.4", "126.2",
    "0.043", "0.081", "0.021", "0.026",
    # --- structural counts stated in words (dimensions, index labels, ranks)
    "7", "28", "48", "84", "21", "35", "56", "44", "24", "10", "15",
    "16", "20", "26", "2", "4", "5", "6", "8", "9", "14", "118", "300",
    "1", "0",
}

# Exact rationals are PRINTED AS RATIONALS by the runs (su7_vacuum.py and
# su7_repair_space.sage keep them in QQ), so they are checked as the literal
# string "a/b" rather than through a decimal expansion that no run ever emits.
FRACTIONS = ["27/8", "5/4", "9/8", "15/8", "29/8", "19/8", "11/8", "5/8",
             "1/8", "27/46", "27/26", "15/4", "3/4", "1/4", "5/2", "1/2",
             "3/2", "1/6", "1/3", "1/9", "19/46", "1/26", "1/27", "1/18"]


def tokens(s):
    return set(re.findall(r"(?<![\w.])\d+(?:\.\d+)?(?![\w.])", s))


def audit(tex, have, HAVE_RAW):
    s = open(tex, encoding="utf-8").read()
    # (?<!\\): an ESCAPED percent is content, not a comment
    s = re.sub(r"(?<!\\)%.*", "", s)
    s = re.sub(r"\\begin\{thebibliography\}.*", "", s, flags=re.S)
    # layout is not data
    s = re.sub(r"\\definecolor\{[^}]*\}\{RGB\}\{[^}]*\}", "", s)
    s = re.sub(r"[\d.]+\\textwidth", "", s)
    s = re.sub(r"\\(geometry|includegraphics|usepackage|documentclass)\[[^\]]*\]", "", s)
    s = re.sub(r"p\{[\d.]+\\textwidth\}|\[[\d.]+pt\]|\{[\d.]+cm\}", "", s)
    s = re.sub(r"\\setlength\{\\[a-zA-Z]+\}\{[^}]*\}", "", s)
    s = re.sub(r"\\rowcolor\{[^}]*\}|\\cellcolor\{[^}]*\}", "", s)
    s = re.sub(r"[a-zA-Z]+!\d+", "", s)                    # colour mixes: hdrblue!12
    s = re.sub(r"\\\\\[\d+pt\]", "", s)                    # \\[7pt]
    s = re.sub(r"\\fboxrule|\\fboxsep", "", s)
    # a thin space inside a numeral is typesetting: 18\,648 IS the number 18648
    s = re.sub(r"(?<=\d)\\,(?=\d)", "", s)

    found, missing, allowed, frac = [], [], [], []
    for t in sorted(tokens(s), key=lambda x: (-len(x), x)):
        if len(t) < 2 and t not in ALLOW:
            continue
        (allowed if t in ALLOW else found if t in have else missing).append(t)
    # exact rationals, checked as the literal string the runs print
    for f in FRACTIONS:
        if f in s:
            (frac if f in HAVE_RAW else missing).append("the rational %s" % f)

    print("%s: every printed number against %s" % (tex, ", ".join(d + "/*" for d in OUT)))
    print("  greppable in an archived run : %d" % len(found))
    print("  exact rationals, as printed  : %d" % len(frac))
    print("  declared non-measurements    : %d" % len(allowed))
    print("  NOT FOUND                    : %d" % len(missing))
    for t in missing:
        print("     %-14s <-- archive its run, or remove it from the paper" % t)
    return missing


def main():
    corpus = ""
    for d in OUT:
        for fn in sorted(os.listdir(d)):
            p = os.path.join(d, fn)
            if os.path.isfile(p):
                corpus += open(p, encoding="utf-8", errors="ignore").read() + "\n"
    have = tokens(corpus)
    # a run that printed 2.73e+00 DID archive 2.73
    have |= set(re.findall(r"(?<![\w.])(\d+(?:\.\d+)?)e[+-]?\d+", corpus))
    bad = []
    for tex in TEXS:
        bad += audit(tex, have, corpus)
        print()
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
