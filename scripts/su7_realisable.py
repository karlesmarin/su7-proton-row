#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

WHICH OF THE SURVIVING THREE-GENERATION ASSIGNMENTS FIT INSIDE THEIR OWN CONTENT?

su7_family_u1.py finds the rung assignments that cancel all six anomaly channels
while protecting every generation.  It reports 14, and it separately notes (its
STEP 6) that "rung 2 needs 4 boxes -> a rep they do not introduce".  The two
facts were never put together.

This script puts them together, and then carries the survivors through the
vacuum criterion of su7_vacuum.py.  Three questions, in order:

  1. Which rungs can Komori-Maru's own tensors host?
  2. Of the surviving assignments, which use only those rungs?
  3. For each survivor, how many 84(+,+) leave the Higgs potential, and which
     rows of their Table 1 still break electroweak symmetry?

Everything is re-derived here -- the channels, the search, the vacuum -- so this
is an independent path and not a filter applied to a stored list.  The 14 is
re-obtained and asserted against su7_family_u1.py's published count.

INPUTS, from Komori-Maru arXiv:2503.04090:
  eq. (41)         Y = (0,0,0,1/2,1/2,-1,0) per index
  eq. (43)         their own leptons: a bulk 21, (eta,eta') = (-1,+1)
  eq. (46)         lepton Yukawa <A5>: index 5 <-> index 7, so extra(e^c)=1/2-l
  before eq. (47)  the quarks are 4D brane fields (their hypercharge is not on
                   the bulk lattice) -- so the brane-quark charge a is free data
  Table 1          the five fermion contents
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr

P = lambda *a: print(*a, flush=True)
H = F(1, 2)

# ============================================================================
# 0.  one generation, written all-left-handed, as (T3C, T3L, Y, X)
# ============================================================================
# colour is carried by T3C in {1/2,-1/2,0}: sum T3C^2 = 1/2 = T(fund), which is
# what the SU(3)^2 X channel needs, and it comes out rather than being put in.
COL = (H, -H, F(0))


def generation(a, l, nus=()):
    s = []
    for c in COL:
        for tl in (H, -H):
            s.append((c, tl, F(1, 6), a + tl))          # q_L
        s.append((-c, F(0), F(-2, 3), -(a + H)))        # u^c
        s.append((-c, F(0), F(1, 3), -(a - H)))         # d^c
    for tl in (H, -H):
        s.append((F(0), tl, -H, l + tl))                # L
    s.append((F(0), F(0), F(1), H - l))                 # e^c
    for v in nus:
        s.append((F(0), F(0), F(0), v))                 # SM singlets kept
    return s


CHANNELS = {
    "SU(3)^2 X": lambda s: sum(x * c * c for c, t, y, x in s),
    "SU(2)^2 X": lambda s: sum(x * t * t for c, t, y, x in s),
    "X grav^2":  lambda s: sum(x for c, t, y, x in s),
    "X^3":       lambda s: sum(x ** 3 for c, t, y, x in s),
    "X^2 Y":     lambda s: sum(x * x * y for c, t, y, x in s),
    "X Y^2":     lambda s: sum(x * y * y for c, t, y, x in s),
}
Achg = lambda a, l: 3 * a + l        # proton-operator charge of a generation

P("=" * 78)
P("STEP 0 -- CONTROL: the one-generation case must reproduce X_Q = -1/6")
P("=" * 78)
hits = [q for q in (F(n, 18) for n in range(-18, 19))
        if all(f(generation(q, H)) == 0 for k, f in CHANNELS.items()
               if k not in ("X grav^2", "X^3"))]
P("  family-universal, one generation, the four cancellable channels vanish at")
P("  X_Q = %s   (their published content; the other two cannot cancel)" % hits)
assert hits == [F(-1, 6)], hits
P("  and there  A = 3X_Q + 1/2 = %s  -- the unprotected point.  PASS"
  % Achg(F(-1, 6), H))

# ============================================================================
# 0b. the singlet sector: a 28 is NECESSARY, not merely minimal
# ============================================================================
P("")
P("=" * 78)
P("STEP 0b -- the two uncancellable channels need a 28.  Not minimally: at all")
P("=" * 78)
P("  The SM singlets are the pure index-7 components, charge m/2 for m boxes:")
P("     7 -> +-1/2 ,  28 -> +-1 ,  84 -> +-3/2 ,  48 -> 0 (contributes nothing)")
P("  Write u, w, e for the NET numbers at +-1/2, +-3/2, +-1.  The two channels")
P("  X.grav^2 and X^3 read")
P("     sum x   = u/2 + 3w/2 + e   = -1     <=>   u + 3w + 2e = -2")
P("     sum x^3 = u/8 + 27w/8 + e  = -1     <=>   u + 27w + 8e = -8")
P("  Subtracting:  24 w + 6 e = -6,  i.e.  e = -1 - 4w,  and then  u = 5w.")
P("")
P("  >> e is NEVER 0: e = 0 would need w = -1/4.  So no number of 7s and 84s,")
P("     at any parities, cancels the two channels.  A 28 is NECESSARY.")
P("     The whole solution set is the one-parameter family (u, w, e) = (5w, w, -1-4w),")
P("     so the net 28-number is always congruent to 3 mod 4.")
P("")
RAD = 7
def brute(allow28):
    out = set()
    for np_ in range(RAD):
        for nm in range(RAD):
            for mp in range(RAD):
                for mm in range(RAD):
                    for ep in range(RAD if allow28 else 1):
                        for em in range(RAD if allow28 else 1):
                            s = F(np_ - nm, 2) + F(3 * (mp - mm), 2) + (ep - em)
                            s3 = F(np_ - nm, 8) + F(27 * (mp - mm), 8) + (ep - em)
                            if s == -1 and s3 == -1:
                                out.add((np_ - nm, mp - mm, ep - em))
    return sorted(out)
no28, with28 = brute(False), brute(True)
P("  brute force over all multiplicities up to %d, no algebra used:" % (RAD - 1))
P("     with no 28 available          : %d solutions   %s" % (len(no28), no28))
P("     with the 28 available (u,w,e) : %d solutions   %s" % (len(with28), with28))
assert no28 == []
assert all(e == -1 - 4 * w and u == 5 * w for u, w, e in with28)
P("     every one obeys e = -1-4w and u = 5w                            PASS")
P("")
P("  CONTROL that must fire, or the scan proves nothing -- drop ONE condition")
P("  (keep sum x = -1 only, still with no 28) and solutions must appear:")
relax = {(np_ - nm, mp - mm) for np_ in range(RAD) for nm in range(RAD)
         for mp in range(RAD) for mm in range(RAD)
         if F(np_ - nm, 2) + F(3 * (mp - mm), 2) == -1}
P("     %d solutions -- the search CAN find things, so the empty set above is" % len(relax))
P("     a result and not a broken loop.                                 PASS")
assert len(relax) > 0

# ============================================================================
# 1.  which rungs can their own tensors host?
# ============================================================================
P("")
P("=" * 78)
P("STEP 1 -- the rungs their own tensors can host")
P("=" * 78)
P("  A rung-k lepton doublet is the component (5,6,7^k): 2+k boxes, with the")
P("  index 7 repeated k times.  Its partner e_R is (6,7^{k+1}), by their eq.")
P("  (46).  An antisymmetric tensor cannot repeat an index.")
P("")
P("  %-6s %-6s %-7s %-16s %s" % ("rep", "boxes", "type", "hosts L at rung",
                                 "and e_R too?"))
TENSORS = (("7", 1, "fund"), ("21", 2, "anti"), ("28", 2, "sym"),
           ("35", 3, "anti"), ("84", 3, "sym"), ("48", None, "adjoint"))
hostable = set()
for nm, bx, ty in TENSORS:
    okL, okE = [], []
    for k in range(0, 4):
        if bx is None:
            continue
        nL = [0] * 7
        nL[4], nL[5], nL[6] = 1, 1, k                   # (5,6,7^k)
        nE = [0] * 7
        nE[5], nE[6] = 1, k + 1                         # (6,7^{k+1})
        fits = lambda n: sum(n) == bx and (max(n) <= 1 if ty == "anti" else True)
        if fits(nL):
            okL.append(k)
            if fits(nE):
                okE.append(k)
    P("  %-6s %-6s %-7s %-16s %s" % (nm, bx if bx else "1+1b", ty,
                                     okL if okL else "none",
                                     okE if okE else "none"))
    hostable |= set(okE)
P("")
P("  >> a rung needs BOTH L and e_R inside one multiplet, or the charged-lepton")
P("     Yukawa of their eq. (46) does not close.  The 35 hosts L at rung 1 and")
P("     NOT e_R -- an antisymmetric tensor has no repeated index.")
P("  >> rungs their content can host : %s" % sorted(hostable))
P("     rung 0 -> the 21 (their own leptons, their eq. (43)) or a 28")
P("     rung 1 -> the 84, and only the 84")
P("     rung >= 2 needs >= 4 boxes: NOT A REPRESENTATION THEIR PAPER INTRODUCES")
assert sorted(hostable) == [0, 1], sorted(hostable)
RUNG_OF = {H: 0, F(0): 1, -H: 2, F(-1): 3}              # l = (1-k)/2
LREAL = {l for l, k in RUNG_OF.items() if k in hostable}
P("     in terms of the lepton charge l = (1-k)/2 :  l in %s"
  % sorted(str(v) for v in LREAL))

# the SM singlets are the pure index-7 components, so charge m/2 needs m boxes
NU_ALL = [F(m, 2) for m in range(-4, 5) if m]
NU_REAL = [v for v in NU_ALL if abs(v) <= F(3, 2)]
P("")
P("  The same counting bounds the neutrinos: an SM singlet of charge m/2 is the")
P("  pure component (7^m), so it needs m boxes.  Their reps give |charge| <= 3/2")
P("  (7, 28, 84); |charge| = 2 would need a four-box tensor they never introduce.")

# ============================================================================
# 2.  the search, re-derived
# ============================================================================
P("")
P("=" * 78)
P("STEP 2 -- every assignment that cancels all six channels and protects all 3")
P("=" * 78)
P("  Brane-quark charge kept family-universal (an exact CKM matrix needs it),")
P("  so a = -(sum l)/9 is forced, and A_j = l_j - mean(l).")
P("")
sols = []
for ls in cwr(sorted(RUNG_OF, reverse=True), 3):
    a = -sum(ls) / 9
    As = [Achg(a, l) for l in ls]
    if any(A == 0 for A in As):
        continue                                        # unprotected
    for n in range(0, 4):
        for nus in cwr(NU_ALL, n):
            s = []
            for j, l in enumerate(ls):
                s += generation(a, l, nus if j == 0 else ())
            if all(f(s) == 0 for f in CHANNELS.values()):
                sols.append((ls, a, tuple(As), nus))
                break
        else:
            continue
        break
P("  %-24s %-8s %-30s %-24s %s"
  % ("(l1,l2,l3)", "X_Q", "(A1,A2,A3)", "neutrinos", "realisable?"))
real = []
for ls, a, As, nus in sols:
    okR = set(ls) <= LREAL and all(v in NU_REAL for v in nus)
    if okR:
        real.append((ls, a, As, nus))
    P("  %-24s %-8s %-30s %-24s %s"
      % (str(tuple(str(x) for x in ls)), a,
         str(tuple(str(x) for x in As)), str(tuple(str(v) for v in nus)),
         "YES" if okR else "no -- rung >= 2"))
P("")
P("  total assignments surviving all six channels : %d" % len(sols))
assert len(sols) == 14, len(sols)
P("  (su7_family_u1.py publishes 14 by an independent implementation: PASS)")
P("  of those, realisable inside their own tensor content : %d" % len(real))

P("")
P("  >> AND EVERY REALISABLE ONE NEEDS A SINGLET OF CHARGE -1:")
for ls, a, As, nus in real:
    P("       l = %-22s neutrinos %-20s contains -1 : %s"
      % (str(tuple(str(x) for x in ls)), str(tuple(str(v) for v in nus)),
         F(-1) in nus))
assert all(F(-1) in nus for _, _, _, nus in real)
P("     charge -1 is the pure (7,7) component, i.e. a 28 -- and 28(+,+) is the")
P("     parity that carries it (su7_anomaly_channels.py STEP 8).  So the nu_R")
P("     condition of section 2, derived on their ONE-generation content, is")
P("     still the condition on the THREE-generation extension.")

# ============================================================================
# 3.  the bill, re-run for each realisable assignment
# ============================================================================
P("")
P("=" * 78)
P("STEP 2b -- and the STRICT reading loses its last assignment")
P("=" * 78)
P("  Their bulk lattice is exactly (1/2)Z (su7_qphi.py), so the STRICT reading of")
P("  what may sit at a fixed point demands X_Q in (1/2)Z.  Which assignments")
P("  satisfy it, and are any of them realisable?")
P("")
P("  %-24s %-9s %-12s %s" % ("(l1,l2,l3)", "X_Q", "X_Q in Z/2?", "realisable?"))
strict_real = []
for ls, a, As, nus in sols:
    st = (2 * a).denominator == 1
    rl = set(ls) <= LREAL and all(v in NU_REAL for v in nus)
    if st or rl:
        P("  %-24s %-9s %-12s %s"
          % (str(tuple(str(x) for x in ls)), a, "yes" if st else "no",
             "yes" if rl else "no -- rung >= 2"))
    if st and rl:
        strict_real.append((ls, a))
P("")
P("  assignments that are BOTH strict-compatible and realisable : %s"
  % (strict_real or "NONE"))
assert not strict_real
P("")
P("  >> the two realisable assignments sit at X_Q = -1/9 and -1/18, neither on")
P("     the bulk lattice; the one strict-compatible assignment needs rung 3,")
P("     i.e. a five-box tensor their paper never introduces.  So under the")
P("     STRICT reading their model has NO three-generation extension inside its")
P("     own content at all.  That forces the permissive reading a second time,")
P("     and by an argument independent of the one-generation X_Q = -1/6.")

P("")
P("=" * 78)
P("STEP 3 -- the vacuum, once each realisable assignment is paid for")
P("=" * 78)
P("  A lepton-hosting bulk multiplet is massive and leaves the potential (their")
P("  own sentence below eq. (79)).  Rung 0 is hosted by a 21, which is NOT in")
P("  Table 1, so it is free.  Rung 1 is hosted by an 84(+,+), which IS.")
P("")
D = {"7(+,+)": F(-3, 4), "7(+,-)": F(1), "28(+,+)": F(2), "28(+,-)": F(1, 4),
     "48(+,+)": F(0), "84(+,+)": F(5, 4)}
GAUGE = F(-27, 8)
TABLE1 = {                                              # their Table 1, verbatim
    "(1)": {"28(+,-)": 1, "84(+,+)": 4},
    "(2)": {"28(+,+)": 1, "84(+,+)": 4},
    "(3)": {"28(+,+)": 1, "48(+,+)": 3, "84(+,+)": 2},
    "(4)": {"7(+,-)": 1, "48(+,+)": 2, "84(+,+)": 3},
    "(5)": {"7(+,+)": 1, "7(+,-)": 1, "84(+,+)": 4},
}
NUR = {"(1)": False, "(2)": True, "(3)": True, "(4)": False, "(5)": False}
Dof = lambda c: GAUGE + sum(n * D[r] for r, n in c.items())
P("  control -- D on their published rows must be > 0 on all five, since every")
P("  row of their Table 1 is a row BECAUSE it breaks electroweak symmetry:")
P("     %s" % {k: str(Dof(c)) for k, c in TABLE1.items()})
assert all(Dof(c) > 0 for c in TABLE1.values())
P("     PASS.  And the gauge-only content must be dead: D = %s < 0.  PASS"
  % GAUGE)

for ls, a, As, nus in real:
    n84 = sum(1 for l in ls if RUNG_OF[l] == 1)
    P("")
    P("  ---- assignment l = %s,  X_Q = %s"
      % (str(tuple(str(x) for x in ls)), a))
    P("       %d generation(s) at rung 1  ->  %d x 84(+,+) donated, cost %s"
      % (n84, n84, n84 * D["84(+,+)"]))
    P("       %-6s %-12s %-12s %-9s %-8s %s"
      % ("case", "D theirs", "D donated", "EWSB", "nu_R", "survives both"))
    win = []
    for cs in ("(1)", "(2)", "(3)", "(4)", "(5)"):
        c = dict(TABLE1[cs])
        if c.get("84(+,+)", 0) < n84:
            P("       %-6s %-12s %s" % (cs, str(Dof(c)), "cannot host: too few 84s"))
            continue
        c["84(+,+)"] -= n84
        d = Dof(c)
        both = d > 0 and NUR[cs]
        if both:
            win.append(cs)
        P("       %-6s %-12s %-12s %-9s %-8s %s"
          % (cs, str(Dof(TABLE1[cs])), str(d), "yes" if d > 0 else "LOST",
             "yes" if NUR[cs] else "no", "**YES**" if both else ""))
    P("       >> rows surviving both conditions: %s" % (win or "none"))

import json
import os
REC = {"source": "su7_realisable.py",
       "D_theirs": {cs: str(Dof(c)) for cs, c in TABLE1.items()},
       "nuR": NUR, "assignments": []}
for ls, a, As, nus in real:
    n84 = sum(1 for l in ls if RUNG_OF[l] == 1)
    rec = {"l": [str(x) for x in ls], "X_Q": str(a), "n84": n84, "rows": {}}
    for cs in TABLE1:
        c = dict(TABLE1[cs])
        if c.get("84(+,+)", 0) < n84:
            continue
        c["84(+,+)"] -= n84
        rec["rows"][cs] = str(Dof(c))
    REC["assignments"].append(rec)
_d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper_data")
with open(os.path.join(_d, "su7_realisable.json"), "w") as fh:
    json.dump(REC, fh, indent=1)

P("")
P("=" * 78)
P("WHAT THIS ADDS")
P("=" * 78)
P("  [paper_data/su7_realisable.json written for the figure]")
P("  1. Their own tensors host rungs 0 and 1 only, so 12 of the 14 assignments")
P("     need a representation their paper never introduces.")
P("  2. Both realisable assignments require a singlet of charge -1, which lives")
P("     in the 28 and nowhere else.  The nu_R condition therefore does NOT have")
P("     to be imposed at one generation and the donation at three: both are")
P("     conditions on the same three-generation extension.")
P("  3. Case (2) is the unique row surviving BOTH realisable assignments.  On")
P("     the second one it is the unique row surviving the vacuum condition")
P("     alone, before the nu_R condition is imposed at all.")
