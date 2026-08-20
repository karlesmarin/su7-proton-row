#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

THE SIXTH ROW: PRE-REGISTERING THE TEST THAT WOULD SEPARATE THE TWO READINGS.

ANCHOR_SECTION_31.md section 5 states the confound and stops there: the two rows
whose alpha_min we fail to reproduce, (3) and (4), are exactly the two containing
a 48 -- AND exactly the two with the smallest alpha_min.  With five rows the two
readings cannot be separated, and no sixth published row exists.

That is a reason to write the prediction down NOW, not a reason to drop it.  A
content with a 48 and a LARGE alpha_min, or one with no 48 and a SMALL one,
separates them, and the two hypotheses predict OPPOSITE things for it:

    "the 48 is the locus"     -> the discrepancy tracks the 48 content
    "small alpha is the locus" -> the discrepancy tracks alpha_min

So this searches their own allowed content -- up to six fermions from 7, 28, 48
and 84, their own sentence before Table 1 -- for the two discriminating kinds,
and prints what OUR instrument predicts for each.  Whoever computes such a row
next, they or anyone, can then read the answer straight off.

[[the-conclusion-line-is-a-prediction]]
"""
import itertools
import json
import math
import os
from fractions import Fraction as F

import numpy as np

P = lambda *a: print(*a, flush=True)

src = open("su7_vacuum.py", encoding="utf-8").read().split("\n")
# cut at a NAMED line, never at a number: su7_vacuum.py grew and 102 fell inside
# minimise(), which then compiled without its return and gave None to everyone.
cut = next(i for i, l in enumerate(src) if l.startswith("REC = {"))
g = {"__name__": "probe"}
exec(compile("\n".join(src[:cut]), "su7_vacuum_head", "exec"), g)
V, minimise = g["V"], g["minimise"]
FINE = g["GRID"]
COARSE = np.linspace(0.0, 1.0, 1501)

TYPES = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
         ("48", 1, 1), ("84", 1, 1)]
DVAL = {("7", 1, 1): F(-3, 4), ("7", 1, -1): F(1), ("28", 1, 1): F(2),
        ("28", 1, -1): F(1, 4), ("48", 1, 1): F(0), ("84", 1, 1): F(5, 4)}
GAUGE = F(-27, 8)
MW, G4 = 80.4, 0.63
KC = 2.245624                       # K = m_h a / sqrt(F'') = KC * g4

# their five rows, as the reference scale
THEIRS = {"(1)": 0.043, "(2)": 0.081, "(3)": 0.021, "(4)": 0.026, "(5)": 0.043}
OURS = {"(1)": 0.0553, "(2)": 0.0831, "(3)": 0.0436, "(4)": 0.0469, "(5)": 0.0553}
N48 = {"(1)": 0, "(2)": 0, "(3)": 3, "(4)": 2, "(5)": 0}

P("=" * 78)
P("A -- THE CONFOUND, AS THEIR OWN FIVE ROWS PRESENT IT")
P("=" * 78)
P("  %-6s %-8s %-9s %-8s %s" % ("case", "n(48)", "a theirs", "a ours", "ratio"))
for cs in THEIRS:
    P("  %-6s %-8d %-9.4f %-8.4f %.2f"
      % (cs, N48[cs], THEIRS[cs], OURS[cs], OURS[cs] / THEIRS[cs]))
P("")
P("  n(48) and a_theirs are perfectly rank-correlated across the five, so the")
P("  ratio column cannot tell them apart.  Every row with a 48 is a row with a")
P("  small alpha.  THAT is the confound, and it needs a sixth row to break.")

# --------------------------------------------------------------- the search
P("")
P("=" * 78)
P("B -- SEARCHING THEIR OWN ALLOWED CONTENT FOR THE TWO DISCRIMINATING KINDS")
P("=" * 78)
P("  Up to six fermions from 7, 28, 48, 84 -- their own sentence before Table 1.")
P("  Keep only contents that break electroweak symmetry, D > 0, which is exact.")
P("")
cands = []
for n in range(1, 7):
    for combo in itertools.combinations_with_replacement(TYPES, n):
        D = GAUGE + sum(DVAL[t] for t in combo)
        if D <= 0:
            continue
        cands.append((combo, D))
P("  contents enumerated : %d" % sum(1 for n in range(1, 7)
                                     for _ in itertools.combinations_with_replacement(TYPES, n)))
P("  of those, D > 0     : %d" % len(cands))

g["GRID"] = COARSE
rows = []
for combo, D in cands:
    content = []
    for t in set(combo):
        content.append((t[0], t[1], t[2], combo.count(t)))
    a = minimise(content)
    if not (0.005 < a < 0.30):
        continue
    n48 = sum(1 for t in combo if t[0] == "48")
    rows.append({"combo": combo, "D": D, "a": a, "n48": n48})
P("  with a minimum strictly inside (0.005, 0.30) : %d" % len(rows))

# A sixth row is only useful if it is a row somebody would PUBLISH.  Their own
# selection was m_h in 125-127 GeV, so a candidate outside that band by a factor
# is idle: nobody will ever compute it and the prediction would never be read.
for r in rows:
    content = [(t[0], t[1], t[2], r["combo"].count(t)) for t in set(r["combo"])]
    fpp = float(V(content, r["a"], d=2)[0])
    r["fpp"] = fpp
    r["mh"] = KC * G4 * math.sqrt(max(fpp, 0.0)) / r["a"]
band = [r for r in rows if 120 <= r["mh"] <= 132]
P("")
P("  and of those, with m_h in 120-132 GeV -- a band around their own")
P("  selection of 125-127, which is what makes a row publishable : %d" % len(band))
P("")
P("  THE QUESTION THIS ACTUALLY ANSWERS: inside that band, is n(48) still")
P("  locked to a small alpha?  If it is, no sixth row can EVER separate the two")
P("  readings within their content, and the confound is structural.")
P("")
P("     %-8s %-10s %-10s %s" % ("n(48)", "contents", "alpha range", "m_h range"))
for k in sorted({r["n48"] for r in band}):
    sub = [r for r in band if r["n48"] == k]
    P("     %-8d %-10d %-10s %s"
      % (k, len(sub), "%.4f-%.4f" % (min(r["a"] for r in sub), max(r["a"] for r in sub)),
         "%.0f-%.0f" % (min(r["mh"] for r in sub), max(r["mh"] for r in sub))))
with48 = [r["a"] for r in band if r["n48"] >= 1]
no48 = [r["a"] for r in band if r["n48"] == 0]
P("")
if with48 and no48:
    overlap = min(max(with48), max(no48)) - max(min(with48), min(no48))
    P("     alpha with a 48 : %.4f - %.4f      alpha with none : %.4f - %.4f"
      % (min(with48), max(with48), min(no48), max(no48)))
    P("     the two ranges OVERLAP by %.4f : %s" % (max(overlap, 0.0),
                                                    "yes" if overlap > 0 else "NO"))
    P("     >> so a publishable row DOES exist that breaks the lock." if overlap > 0
      else "     >> so inside their content the confound is STRUCTURAL.")
else:
    P("     one of the two kinds is EMPTY in the publishable band: %d with a 48, "
      "%d without" % (len(with48), len(no48)))
    P("     >> inside their content the confound cannot be broken at all.")

g["GRID"] = FINE
lo_a = max(OURS[c] for c in ("(3)", "(4)"))        # 0.0469: "small" ends here
hi_a = min(OURS[c] for c in ("(1)", "(5)"))        # 0.0553: "large" starts here
rows = band or rows


def refine(r):
    content = [(t[0], t[1], t[2], r["combo"].count(t)) for t in set(r["combo"])]
    a = minimise(content)
    # V carries its own derivative order, so F'' is analytic here, not a difference
    fpp = float(V(content, a, d=2)[0])
    return a, fpp, KC * G4 * math.sqrt(max(fpp, 0.0)) / a


def show(title, sel, want):
    P("")
    P("  %s" % title)
    if not sel:
        P("     none found in their allowed content.")
        return None
    sel = sorted(sel, key=lambda r: want(r))[:3]
    P("     %-42s %-6s %-9s %-9s %s"
      % ("content", "n(48)", "a ours", "F''", "m_h at g4=0.63"))
    best = None
    for r in sel:
        a, fpp, mh = refine(r)
        lab = " + ".join("%d x %s(%+d,%+d)" % (r["combo"].count(t), t[0], t[1], t[2])
                         for t in sorted(set(r["combo"])))
        P("     %-42s %-6d %-9.4f %-9.2f %.1f" % (lab, r["n48"], a, fpp, mh))
        if best is None:
            best = dict(r, a=a, fpp=fpp, mh=mh, label=lab)
    return best


A = show("(i) WITH a 48 and a LARGE alpha_min -- the kind their table has not got:",
         [r for r in rows if r["n48"] >= 1 and r["a"] >= hi_a],
         lambda r: -r["a"])
B = show("(ii) with NO 48 and a SMALL alpha_min -- the other discriminator:",
         [r for r in rows if r["n48"] == 0 and r["a"] <= lo_a],
         lambda r: r["a"])

# ------------------------------------------------------------ the prediction
P("")
P("=" * 78)
P("C -- THE PREDICTION, WRITTEN BEFORE ANY SUCH ROW EXISTS")
P("=" * 78)
rat_small = sum(OURS[c] / THEIRS[c] for c in ("(3)", "(4)")) / 2
rat_large = sum(OURS[c] / THEIRS[c] for c in ("(1)", "(2)", "(5)")) / 3
P("  On their five rows the ratio a_ours/a_theirs is %.2f on the two rows with a"
  % rat_small)
P("  48 and %.2f on the three without.  For a sixth row the two readings split:" % rat_large)
P("")
for tag, r, s48, salpha in (("(i)  48-rich, large alpha", A, rat_small, rat_large),
                            ("(ii) no 48, small alpha", B, rat_large, rat_small)):
    if r is None:
        continue
    P("  %s" % tag)
    P("     content            %s" % r["label"])
    P("     ours               a_min = %.4f,  m_h = %.1f GeV at g4 = 0.63"
      % (r["a"], r["mh"]))
    P("     if the 48 is the locus       -> published a_min near %.4f  (ratio %.2f)"
      % (r["a"] / s48, s48))
    P("     if small alpha is the locus  -> published a_min near %.4f  (ratio %.2f)"
      % (r["a"] / salpha, salpha))
    P("     the two differ by a factor %.2f, which two significant figures resolve."
      % (max(s48, salpha) / min(s48, salpha)))
    P("")
P("  Neither reading is adopted here.  What is fixed is the number each of them")
P("  commits to, on a content that does not yet exist in print, so that whoever")
P("  computes one can read the answer off instead of re-opening the argument.")

REC = {"source": "su7_sixth_row.py", "ratio_with48": rat_small,
       "ratio_no48": rat_large,
       "predictions": [{"kind": k, "content": r["label"], "a_ours": r["a"],
                        "mh_ours": r["mh"], "n48": r["n48"]}
                       for k, r in (("48-rich, large alpha", A),
                                    ("no 48, small alpha", B)) if r]}
_d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper_data")
with open(os.path.join(_d, "su7_sixth_row.json"), "w") as fh:
    json.dump(REC, fh, indent=1)
P("  [paper_data/su7_sixth_row.json written]")
