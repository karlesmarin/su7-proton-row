#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

CAN A BRANE FIELD FINISH THE HALF-GENERATION THE 48 CARRIES?

This is question 5 of the paper, aimed at the paper's own section 4, and it is
asked before publishing rather than after. If it opens, section 4's
classification is incomplete and the headline goes with it.

FIRST, A CORRECTION TO OUR OWN ARITHMETIC. su7_can_the_adjoint_host.py checked
the partner charge with `want = extra(L) - 1`. That is wrong. The ladder has
L = (5,6,7^k) and e_R = (6,7^{k+1}), so with +1/2 per index 6 and -1/2 per index 7

    extra(L) = (1-k)/2 ,   extra(e_R) = -k/2 ,   so   extra(e_R) = extra(L) - 1/2

and the hypercharges differ by the Higgs's own, Y(e_R) = Y(L) - 1/2. Those are
exactly what <A_5> does to a component (index 5 -> index 7 changes Y by -1/2 and
extra by -1/2). So the two "independent routes" that script claimed are ONE
route, and one of them was run with the wrong target. Redone here from scratch.

The conclusion of that script survives -- there is no partner in the 48 -- but it
survives on one argument, not two, and the reason is worth having right, because
the whole question below is about what happens when the partner is NOT in the 48.

THE STAKES, which are the reason this runs before publication and not after:
donating a 48 costs D = 0 exactly. If the 48 can carry a lepton generation with a
brane partner, no row pays anything for the escape, every row keeps its published
D, and "case (2) is the unique survivor" is false.
"""
import re
from fractions import Fraction as F

import numpy as np

P = lambda *a: print(*a, flush=True)
H = F(1, 2)
N = 7

YIND = [F(0), F(0), F(0), H, H, F(-1), F(0)]        # their eq. (41), per index
EIND = [F(0), F(0), F(0), F(0), F(0), H, -H]        # the U(1)' charge, per index
COL, WEAK = {0, 1, 2}, {3, 4}

_src = open("su7_vacuum.py", encoding="utf-8").read()
_fn = re.search(r"\ndef terms\(.*?\n(?=\n\ndef |\n\nGRID|\Z)", _src, re.S)
_ns = {}
exec(_fn.group(0), _ns)
terms = _ns["terms"]


def qn(i, j):
    return dict(su3=(i in COL) - (j in COL), su2=(i in WEAK) - (j in WEAK),
                Y=YIND[i] - YIND[j], extra=EIND[i] - EIND[j])


P("=" * 78)
P("A -- THE LADDER, AND WHAT A PARTNER MUST CARRY")
P("=" * 78)
P("  %-6s %-14s %-14s %-14s %s" % ("rung k", "extra(L)", "extra(e_R)", "boxes L", "boxes e_R"))
for k in range(0, 4):
    P("  %-6d %-14s %-14s %-14d %d" % (k, F(1 - k, 2), F(-k, 2), 2 + k, 2 + k))
P("")
P("  <A_5> is index 5 -> index 7, so it moves a component by")
P("     dY = Y_7 - Y_5 = %s     d(extra) = e_7 - e_5 = %s"
  % (YIND[6] - YIND[4], EIND[6] - EIND[4]))
P("  which is exactly (Y, extra)(e_R) - (Y, extra)(L). The Yukawa condition and")
P("  the charge condition are THE SAME CONDITION -- not two.")
assert YIND[6] - YIND[4] == -H and EIND[6] - EIND[4] == -H

P("")
P("=" * 78)
P("B -- WHAT THE 48 ACTUALLY CARRIES, COUNTED HONESTLY")
P("=" * 78)
P("  A lepton doublet, all-left-handed, is an SU(3) singlet SU(2) doublet of")
P("  Y = -1/2. Y = +1/2 is its CONJUGATE, which is the same physical state in a")
P("  real representation and not a second generation.")
P("")
Ls, conj = [], []
for i in range(N):
    for j in range(N):
        if i == j:
            continue
        q = qn(i, j)
        if q["su3"] == 0 and abs(q["su2"]) == 1 and abs(q["Y"]) == H:
            (Ls if q["Y"] == -H else conj).append((i, j, q))
for i, j, q in Ls + conj:
    P("  (%d,%dbar)  Y = %-5s extra = %-5s rung = %-4s %s"
      % (i + 1, j + 1, q["Y"], q["extra"], 1 - 2 * q["extra"],
         "lepton doublet" if q["Y"] == -H else "its conjugate"))
P("")
K48 = 1 - 2 * Ls[0][2]["extra"]
P("  >> %d states of Y = -1/2 and %d conjugates. So the 48 carries ONE lepton"
  % (len(Ls), len(conj)))
P("     doublet -- an SU(2) pair -- at rung %s. NOT four doublets at two rungs,"
  % K48)
P("     which is how su7_can_the_adjoint_host.py reported it and how the paper")
P("     now repeats it: the rung-0 pair IS the rung-2 pair conjugated, and a")
P("     conjugate is not a second generation.")
assert len(Ls) == 2 and K48 == 2

P("")
P("=" * 78)
P("C -- ITS PARTNER: NOWHERE IN THEIR BULK CONTENT")
P("=" * 78)
wantY, wantE = Ls[0][2]["Y"] - H, Ls[0][2]["extra"] - H
P("  the partner must be an SU(3) x SU(2) singlet with Y = %s, extra = %s"
  % (wantY, wantE))
P("")
hits = [(i, j) for i in range(N) for j in range(N) if i != j
        and qn(i, j)["su3"] == 0 and qn(i, j)["su2"] == 0
        and qn(i, j)["Y"] == wantY and qn(i, j)["extra"] == wantE]
P("  in the 48                      : %s" % (hits or "ABSENT"))
P("  as a bulk tensor it is (6,7^%d), i.e. %d boxes -- and their content stops"
  % (K48 + 1, K48 + 2))
P("  at three. That IS the four-box exclusion, and it applies to the PARTNER.")
assert not hits

P("")
P("=" * 78)
P("C2 -- THE LAST THING THAT COULD KILL THIS: DOES THE DOUBLET EVEN SURVIVE?")
P("=" * 78)
P("  A component of a bulk multiplet gives a 4D zero mode only if its orbifold")
P("  parity matches. su7_anomaly_channels.py's rule, from their eqs. (37)-(40):")
P("  a weight survives in a multiplet of parities (eta, eta') iff")
P("  P5 P5' of that component equals eta*eta'. If the 48's doublet fails this in")
P("  a 48(+,+), there is no generation to host and the whole question dies here.")
P("")
PI = [1, 1, 1, -1, -1, 1, -1]      # P5P5' per index: diag(1,1,1,1,1,-1,-1) x diag(1,1,1,-1,-1,-1,1)
for i, j, q in Ls:
    p = PI[i] * PI[j]
    P("  (%d,%dbar)  P5P5' = %+d   survives in a 48(+,+) (eta*eta' = +1): %s"
      % (i + 1, j + 1, p, "YES" if p == 1 else "no"))
sur = [1 for i, j, _ in Ls if PI[i] * PI[j] == 1]
P("")
P("  >> %d of %d survive. The doublet is NOT projected out." % (len(sur), len(Ls)))
P("     (Control that had to be able to fire: of the 48's ten states at")
P("     |q| = 1/2, only two carry P5P5' = +1 -- su7_gauge_from_group.py's")
P("     table -- so this test rejects 8 of 10 and is not vacuous.)")
assert len(sur) == len(Ls) == 2

P("")
P("=" * 78)
P("D -- SO PUT IT ON A BRANE. WHAT DOES THAT COST THE POTENTIAL?")
P("=" * 78)
W = F(3, 4)


def D_of(rep, eta, etap):
    A = B = F(0)
    for m, s, c in terms(rep, eta, etap):
        (A, B) = (A + m * c * c, B) if s == 1 else (A, B + m * c * c)
    return A - W * B


P("  %-10s %s" % ("multiplet", "D lost when it is donated to host a generation"))
for rep in ("7", "28", "48", "84"):
    P("  %-10s %s" % (rep + "(+,+)", D_of(rep, 1, 1)))
P("")
P("  >> the 48 is the ONLY one that is free. A generation hosted in a 48 with a")
P("     brane partner takes NOTHING out of the Higgs potential.")
assert D_of("48", 1, 1) == 0

GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]
T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)]),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)]),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)]),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)]),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)])]


def D_row(content, donate84=0):
    tot = F(0)
    for m, s, c in GAUGE:
        tot += F(m).limit_denominator() * c * c * (1 if s == 1 else -W)
    for rep, e, ep, mult in content:
        tot += mult * D_of(rep, e, ep)
    return tot - donate84 * D_of("84", 1, 1)


P("")
P("  %-6s %-16s %-16s %-16s %s"
  % ("case", "D published", "after an 84", "after a 48", "has a 48?"))
for tag, c in T1:
    n48 = sum(m for r, _, _, m in c if r == "48")
    P("  %-6s %-16s %-16s %-16s %s"
      % (tag, D_row(c), D_row(c, 1), D_row(c), "yes, %d" % n48 if n48 else "no"))

P("")
P("=" * 78)
P("E -- THE CONSEQUENCE, STATED BEFORE IT IS ARGUED AWAY")
P("=" * 78)
alive84 = [t for t, c in T1 if D_row(c, 1) > 0]
alive48 = [t for t, c in T1 if D_row(c) > 0 and any(r == "48" for r, _, _, _ in c)]
P("  rows surviving the vacuum condition after donating an 84 : %s" % " ".join(alive84))
P("  rows that COULD instead donate a 48, and survive at no cost: %s" % " ".join(alive48))
P("")
P("  The nu_R condition admits only rows (2) and (3) (sect. 3). Intersecting:")
P("    with the 84 route  : %s   <- the paper's headline"
  % " ".join(t for t in alive84 if t in ("(2)", "(3)")))
P("    with the 48 route  : %s"
  % " ".join(t for t in alive48 if t in ("(2)", "(3)")))
P("")
P("  Case (3) is the row with three 48s. If the 48 route is open it pays")
P("  nothing, keeps 9/8, and survives -- and the headline is false.")
P("")
P("  What is NOT an obstruction, checked rather than assumed:")
P("   - anomalies. The assignment that needs rung 2 is l = (1/2, 1/2, -1/2) at")
P("     X_Q = -1/18, and it is one of the FOURTEEN that cancel all six channels")
P("     (su7_family_u1.py). It was rejected by sect. 4 on rungs alone.")
P("   - the charge lattice. extra = %s is on the bulk lattice (1/2)Z anyway, so"
  % wantE)
P("     not even the permissive reading is needed.")
P("   - the doublet. It is there, in the 48, with the right Y and the right")
P("     extra, as section B shows.")
P("")
P("  So the exclusion of rung 2 does NOT rest on anomalies, charges or content.")
P("  It rests on the partner having to be reachable by <A_5>, i.e. on the Yukawa")
P("  being GAUGE-GENERATED. A brane partner is not reachable by <A_5>: its")
P("  Yukawa would be a brane-localised coupling with a free coefficient. That is")
P("  a real price -- it is the calculability that gauge-Higgs unification exists")
P("  to buy, and Maru's own review calls this obstruction generic -- but it is a")
P("  price, not a prohibition, and this paper cannot claim otherwise.")
