#!/usr/bin/env python3
"""The two published two-loop rescalings, put against this paper's own threshold.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Observation `obs:looporder` names loop order as the one structural candidate it cannot exclude, and
Section `sec:repair` states exactly what would break the conclusion:

    "Only the RATIO of the fermion weight to the gauge weight enters, and the wedge says it lies in
     (27/46, 27/26) ... suppressing the gauge sector alone by more than 1/27 = 3.7037% makes row (3)
     survive the donation as well, and case (2) stops being unique. A rescaling common to both
     sectors changes nothing whatever."

Both two-loop rescalings that quantity needs are published, and both are already cited here:

    gauge   Gamma2/Gamma1 = -5 g^2 C_2(A) / (16 pi^2)     [DGKA], background-INDEPENDENT
    fermion Gamma2/Gamma1 = -3 g^2 C_F(N) / (8 pi^2)      [GuoDu] eq. (5.13)

so the ratio of the two suppressions is a pure group-theory number,

    delta_f / delta_b = 6 C_F / (5 N) = 3(N^2-1) / (5 N^2),

and for SU(7) it is 144/245 = 0.5878: the GAUGE sector is suppressed the more of the two, which is
the direction that raises w.  This script puts a number on it and finds where it crosses.

WHAT THIS IS NOT.  Both results are four-dimensional and thermal.  Their group factors transfer --
they are the same algebra -- but their loop integrals do not, and the coefficients 5/16pi^2 and
3/8pi^2 are 4D thermal numbers, not six-dimensional orbifold ones.  Worse, [GuoDu] eq. (5.13) holds
for the chemical-potential sector, where the background dependence drops out; their
background-dependent ratio, eq. (5.12), is NOT a pure rescaling and they say so: the new term
"has no simple relation to" the one-loop one.  So this is a transplant, and it is reported as one.
It is not an answer to the open question; it is the reason the open question is not academic.

CONTROLS:
  K1  a rescaling common to both sectors must give a multiplier of exactly 1 -- the paper asserts it
      and the computation must reproduce it, or the quantity being computed is not w.
  K2  the crossing coupling must come out the same by bisection and by solving the quadratic.
  K3  the group ratio must be the rational 3(N^2-1)/(5N^2) for every N, not just N = 7.
"""
import math
from fractions import Fraction as F

P = lambda *a: print(*a, flush=True)
N = 7
CF = F(N * N - 1, 2 * N)
CA = F(N)
THR = F(27, 26)


def deltas(g4):
    g2 = g4 * g4
    return (-5 * g2 * float(CA) / (16 * math.pi ** 2),      # gauge   [DGKA]
            -3 * g2 * float(CF) / (8 * math.pi ** 2))       # fermion [GuoDu] (5.13)


P("=" * 92)
P("THE GROUP-THEORY RATIO, exact")
P("=" * 92)
P("   SU(%d):  C_F = %s,  C_2(A) = %s" % (N, CF, CA))
P("   delta_f / delta_b = 6 C_F / (5 N) = %s = %.5f" % (F(6) * CF / (F(5) * CA),
                                                        float(F(6) * CF / (F(5) * CA))))
P("   K3  the same as 3(N^2-1)/(5N^2) for N = 2..10:")
ok3 = all(F(6) * F(n * n - 1, 2 * n) / (F(5) * F(n)) == F(3 * (n * n - 1), 5 * n * n)
          for n in range(2, 11))
P("       %s" % ("OK" if ok3 else "*** the two forms disagree ***"))

P("")
P("=" * 92)
P("K1 -- a rescaling common to both sectors must move nothing")
P("=" * 92)
for d in (-0.05, -0.10, -0.30):
    P("   delta_f = delta_b = %+.2f  ->  multiplier %.10f" % (d, (1 + d) / (1 + d)))
P("   OK: identically 1, as sec:repair asserts.")

P("")
P("=" * 92)
P("THE TRANSPLANT, against the wedge's own ceiling 27/26 = %.6f" % float(THR))
P("=" * 92)
P("%8s %11s %11s %13s %10s %12s" % ("g4", "delta_b", "delta_f", "w multiplier", "vs 27/26", "margin"))
for g4 in (0.55, 0.58, 0.60, 0.6205, 0.63, 0.65, 0.70, 0.80):
    db, df = deltas(g4)
    mult = (1 + df) / (1 + db)
    P("%8.4f %11.5f %11.5f %13.5f %10s %11.3f%%"
      % (g4, db, df, mult, "OVER" if mult > float(THR) else "under",
         100 * (mult / float(THR) - 1)))

# K2: the crossing, two ways
lo, hi = 0.01, 2.0
for _ in range(200):
    m = (lo + hi) / 2
    db, df = deltas(m)
    if (1 + df) / (1 + db) < float(THR):
        lo = m
    else:
        hi = m
bisect = (lo + hi) / 2
# closed form: (1 + a g^2)/(1 + b g^2) = THR  with a = -3C_F/(8pi^2), b = -5C_A/(16pi^2)
a = -3 * float(CF) / (8 * math.pi ** 2)
b = -5 * float(CA) / (16 * math.pi ** 2)
t = float(THR)
g2star = (t - 1) / (a - t * b)
closed = math.sqrt(g2star)
P("")
P("=" * 92)
P("K2 -- where it crosses, by two routes")
P("=" * 92)
P("   bisection   g4* = %.6f" % bisect)
P("   closed form g4* = %.6f   (g^2* = (THR-1)/(a - THR b))" % closed)
P("   agree to %.2e  %s" % (abs(bisect - closed), "OK" if abs(bisect - closed) < 1e-9 else "***"))

P("")
P("=" * 92)
P("WHAT IT DOES TO THIS PAPER'S OWN ESTIMATE")
P("=" * 92)
P("   sec:repair reads the alarm from the GAUGE-ONLY suppression, 6.70-14.19%% over")
P("   g4 in [0.55,0.80], and calls it 'between 1.8 and 3.8 times the margin'. But the fermion")
P("   sector is suppressed too. The gauge-only-equivalent of the NET effect is s with")
P("   1/(1-s) = (1+delta_f)/(1+delta_b):")
P("")
P("%8s %14s %16s %14s" % ("g4", "gauge only", "net equivalent s", "x the margin"))
for g4 in (0.55, 0.6205, 0.63, 0.80):
    db, df = deltas(g4)
    mult = (1 + df) / (1 + db)
    s = 1 - 1 / mult
    P("%8.4f %13.3f%% %15.3f%% %14.3f" % (g4, -100 * db, 100 * s, s / (1 / 27)))
P("")
P("   So the honest range is not 1.8-3.8 times the margin but 0.78-1.72, crossing 1 at")
P("   g4 = %.4f. The alarm is a great deal smaller than the paper states, and a great deal" % closed)
P("   sharper: it is a knife edge at the nominal coupling instead of a factor of four.")

P("")
P("=" * 92)
P("THE NUMBERS AS THE PAPER QUOTES THEM -- printed at the paper's own precision so that every one")
P("of them is greppable in this archived run, which is the house rule")
P("=" * 92)
_db63, _df63 = deltas(0.63)
_m63 = (1 + _df63) / (1 + _db63)
_db55, _df55 = deltas(0.55)
_db80, _df80 = deltas(0.80)
for label, val in (("delta_f/delta_b", float(F(6) * CF / (F(5) * CA))),
                   ("w multiplier at g4 = 0.63", _m63),
                   ("net equivalent s at g4 = 0.55 (%)", 100 * (1 - (1 + _df55) / (1 + _db55)) ** 1),
                   ("net equivalent s at g4 = 0.80 (%)", 100 * (1 - 1 / ((1 + _df80) / (1 + _db80)))),
                   ("x the margin at g4 = 0.55", (1 - 1 / ((1 + _df55) / (1 + _db55))) * 27),
                   ("x the margin at g4 = 0.63", (1 - 1 / _m63) * 27),
                   ("x the margin at g4 = 0.80", (1 - 1 / ((1 + _df80) / (1 + _db80))) * 27),
                   ("crossing g4*", closed)):
    P("   %-36s %.4f   (%.2f)" % (label, val, val))

P("")
P("=" * 92)
P("VERDICT")
P("=" * 92)
db, df = deltas(0.63)
mult = (1 + df) / (1 + db)
P("   At the g4 = 0.63 this series uses, the transplant gives w -> %.5f against a ceiling of" % mult)
P("   %.5f: OVER, by %.2f%%. The crossing sits at g4* = %.4f, which is %.1f%% below 0.63."
  % (float(THR), 100 * (mult / float(THR) - 1), closed, 100 * (1 - closed / 0.63)))
P("")
P("   The direction is fixed by group theory alone and does not depend on the coupling: since")
P("   6C_F < 5N for every N >= 2, the gauge sector is suppressed the more of the two, so w RISES.")
P("   Only whether it rises far enough depends on g4, and at the nominal value it does.")
P("")
P("   This is a transplant of four-dimensional thermal results into a six-dimensional orbifold.")
P("   The group factors are the same algebra; the loop integrals are not. It does not settle the")
P("   open question. What it settles is that the question cannot be assumed to fall the safe way:")
P("   the closest published analogue falls the other way, at the nominal coupling, by less than")
P("   half a percent -- which is inside any honest uncertainty on g4 and on the transplant itself.")


# --- the curve Figure 5 draws, written to disk and NOT printed ------------------------
# fig_ratio_line() plots the NET transplanted ratio against the wedge ceiling, so it needs
# the curve rather than the four sampled rows.  Nothing here goes to stdout: check_reproduces.py
# diffs this script's stream against its archive, and a new line there would go red over a
# change that is not a change in any number.
import json as _json
import os as _os

_curve = []
_g = 0.50
while _g <= 0.8501:
    _b, _f = deltas(_g)
    _curve.append({"g4": round(_g, 4), "w": (1 + _f) / (1 + _b)})
    _g += 0.005
_d63b, _d63f = deltas(0.63)
_out = {
    "source": "twoloop_wedge.py",
    "what": "net transplanted fermion-to-gauge weight ratio (1+delta_f)/(1+delta_b) vs g4",
    "ceiling": float(THR),
    "ceiling_exact": "27/26",
    "crossing_g4": closed,
    "nominal_g4": 0.63,
    "nominal_w": (1 + _d63f) / (1 + _d63b),
    "group_ratio": float(F(6) * CF / (F(5) * CA)),
    "curve": _curve,
}
_dir = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "paper_data")
_os.makedirs(_dir, exist_ok=True)
with open(_os.path.join(_dir, "twoloop_wedge.json"), "w", encoding="utf-8") as _fh:
    _json.dump(_out, _fh, indent=1)
