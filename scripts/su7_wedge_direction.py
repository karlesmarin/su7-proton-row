#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

WHICH SIDE OF THE WEDGE DOES THE ANCHOR DATA POINT TO?

su7_loop_order_margin.py measures the wedge as a tolerance on ONE number, the
fermion-to-gauge weight ratio: the verdict holds iff it lies in (27/46, 27/26),
i.e. -41.304 % / +3.846 % about 1.  The upper margin is thin, and the two-loop
pure-gauge factor is larger than it.  Stated alone, that reads as an alarm.

It is not, and the paper's own anchor data is what defuses it.  Two facts nobody
had put together:

  1. alpha_min is MONOTONE INCREASING in the fermion weight.
  2. Our alpha_min is TOO LARGE on every row.

So the reweighting the anchor column asks for is a SUPPRESSION of the fermion
sector relative to the gauge sector -- ratio < 1 -- which is the direction AWAY
from the thin ceiling.  A correction that moves the ratio up towards 27/26 does
not merely risk the verdict: it makes the discrepancy it was invoked to explain
WORSE.

This measures both, per row, on their own five contents.
"""
from fractions import Fraction as F
import numpy as np

P = lambda *a: print(*a, flush=True)

src = open("su7_vacuum.py", encoding="utf-8").read().split("\n")
# cut at a NAMED line, never at a number: su7_vacuum.py grew and 102 fell inside
# minimise(), which then compiled without its return and gave None to everyone.
cut = next(i for i, l in enumerate(src) if l.startswith("REC = {"))
g = {"__name__": "probe"}
exec(compile("\n".join(src[:cut]), "su7_vacuum_head", "exec"), g)
# a coarser entry grid: the bracket is refined 25 times afterwards either way,
# and lam is only wanted to about a per cent.  su7_vacuum.py's own 40001-point
# grid is kept there; this is a speed choice local to this script, and the
# control below is that rows (1),(2) must reproduce ANCHOR_SECTION_31.md's
# published lam to three decimals.
g["GRID"] = np.linspace(0.0, 1.0, 2001)
minimise, V = g["minimise"], g["V"]
PUBLISHED = {"(1)": 0.848, "(2)": 0.962, "(3)": 0.810, "(4)": 0.800, "(5)": 0.848}

ROWS = {
    "(1)": ([("28", 1, -1, 1), ("84", 1, 1, 4)], 0.043),
    "(2)": ([("28", 1, 1, 1), ("84", 1, 1, 4)], 0.081),
    "(3)": ([("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], 0.021),
    "(4)": ([("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], 0.026),
    "(5)": ([("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], 0.043),
}
LO, HI = F(27, 46), F(27, 26)

P("=" * 78)
P("A -- alpha_min IS MONOTONE INCREASING IN THE FERMION WEIGHT")
P("=" * 78)
P("  %-6s %s" % ("case", "  ".join("lam=%.2f" % l for l in (0.6, 0.8, 1.0, 1.2, 1.4))))
mono = True
for cs, (c, _) in ROWS.items():
    xs = [minimise(c, lam=l) for l in (0.6, 0.8, 1.0, 1.2, 1.4)]
    mono &= all(b >= a for a, b in zip(xs, xs[1:]))
    P("  %-6s %s" % (cs, "  ".join("%7.4f" % x for x in xs)))
P("")
P("  monotone increasing on every row : %s" % mono)
assert mono
P("  >> so raising the fermion weight RAISES alpha_min, and lowering it lowers it.")

P("")
P("=" * 78)
P("B -- THE RATIO EACH ROW ASKS FOR, AND WHERE IT SITS IN THE WEDGE")
P("=" * 78)
P("  Solve  argmin F(alpha; lam) = alpha_theirs  for lam, per row.  lam is the")
P("  fermion weight with the gauge part held at 1, i.e. exactly the ratio the")
P("  wedge constrains.")
P("")
P("  %-6s %-10s %-10s %-9s %-10s %s"
  % ("case", "theirs", "ours(1)", "lam", "in wedge?", "side"))
lams = []
for cs, (c, at) in ROWS.items():
    lo, hi = 0.3, 2.0
    for _ in range(28):                       # bisection on a monotone function
        mid = 0.5 * (lo + hi)
        if minimise(c, lam=mid) < at:
            lo = mid
        else:
            hi = mid
    lam = 0.5 * (lo + hi)
    lams.append(lam)
    inw = float(LO) < lam < float(HI)
    P("  %-6s %-10.4f %-10.4f %-9.4f %-10s %-9s published %.3f %s"
      % (cs, at, minimise(c), lam, "yes" if inw else "NO",
         "below 1" if lam < 1 else "ABOVE 1", PUBLISHED[cs],
         "OK" if abs(lam - PUBLISHED[cs]) < 2e-3 else "*** MISMATCH ***"))
    assert abs(lam - PUBLISHED[cs]) < 2e-3, (cs, lam)
P("")
P("  wedge = (%s, %s) = (%.4f, %.4f)" % (LO, HI, float(LO), float(HI)))
P("  every row inside the wedge          : %s" % all(float(LO) < l < float(HI) for l in lams))
P("  every row BELOW 1                   : %s" % all(l < 1 for l in lams))
P("  the largest of them                 : %.4f, against the ceiling %.4f"
  % (max(lams), float(HI)))
assert all(float(LO) < l < float(HI) for l in lams)
assert all(l < 1 for l in lams)

P("")
P("=" * 78)
P("C -- SO THE TWO EXPOSURES ARE MUTUALLY EXCLUSIVE")
P("=" * 78)
P("  The verdict breaks only ABOVE %.4f.  The anchor column is reproduced only"
  % float(HI))
P("  BELOW 1, at %.4f to %.4f.  The two live on opposite sides of 1."
  % (min(lams), max(lams)))
P("")
P("  distance from 1 to the ceiling        : %.4f" % (float(HI) - 1))
P("  distance from 1 down to the fits      : %.4f to %.4f"
  % (1 - max(lams), 1 - min(lams)))
P("  distance from the fits to the ceiling : %.4f to %.4f"
  % (float(HI) - max(lams), float(HI) - min(lams)))
P("")
P("  >> A correction that pushes the ratio UP towards 27/26 does not merely")
P("     endanger the verdict.  It makes alpha_min LARGER, and ours is already")
P("     too large on all five rows -- so it makes the discrepancy it would be")
P("     invoked to explain WORSE.  The two-loop pure-gauge factor is exactly")
P("     such a push: suppressing the gauge sector by delta is a ratio 1/(1-delta)")
P("     above 1.")
P("")
P("  >> Therefore, if loop order IS the anchor residual, its net effect on this")
P("     ratio has to be a SUPPRESSION of the fermion sector -- the fermionic")
P("     two-loop term outrunning the pure-gauge one and reversing its sign --")
P("     and that direction moves the verdict AWAY from breaking, not towards it.")
P("     The alarming direction and the explanatory direction are opposite.")
P("")
P("  This does not make the verdict two-loop safe: nobody has computed the term")
P("  and the ratio could still land above 27/26 while the anchor stays unexplained.")
P("  What it removes is the reading that the leading candidate for the residual")
P("  is ALSO the thing most likely to break the headline.  It is the opposite.")

# ---------------------------------------------------------------- for the figure
import json
import math
import os

REC = {
    "source": "su7_wedge_direction.py",
    "wedge": {"lo": float(LO), "hi": float(HI), "lo_exact": "27/46",
              "hi_exact": "27/26"},
    "rows": [{"case": cs, "a_theirs": ROWS[cs][1], "lam": lam}
             for cs, lam in zip(ROWS, lams)],
    "twoloop": [{"g4": g,
                 "delta": 5 * g * g * 7 / (16 * math.pi ** 2),
                 "ratio": 1 / (1 - 5 * g * g * 7 / (16 * math.pi ** 2))}
                for g in (0.55, 0.63, 0.70, 0.80)],
}
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper_data")
with open(os.path.join(d, "su7_wedge_direction.json"), "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("  [paper_data/su7_wedge_direction.json written for the figure]")
