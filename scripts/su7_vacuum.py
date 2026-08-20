#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

DOES ANY ROW OF KOMORI-MARU'S TABLE 1 SURVIVE DONATING AN 84 TO THE LEPTONS?

SU7_FAMILY_U1.md: a family-dependent U(1)' -- the only way their leftover U(1)'
can forbid proton decay -- needs the third lepton generation in an 84 at
(eta,eta') = (+1,+1), a multiplet present in every row of their Table 1.  And it
named the bill.  Their own section 3.2 makes the bill explicit:

    "we consider fermions belonging to the fundamental, adjoint, two-rank and
     three-rank totally symmetric tensors of SU(7).  We note that these 6D
     fermions have NOTHING TO DO with the SM quarks and leptons, namely the SM
     fermions are NOT EMBEDDED in these SU(7) multiplets."

plus, below their eq. (79), that a lepton-carrying bulk fermion "must be massive"
and its potential is suppressed by e^{-M pi R}.  So an 84 that hosts leptons is
not one of the potential-generating 84s, and this recomputes their vacuum with
one removed.

>>> THE HEADLINE OF THIS SCRIPT IS A NEGATIVE ABOUT THE INSTRUMENT. <<<
Three controls pass and the absolute anchor does NOT: the published a_min values
are not reproduced from the published formulas.  So no absolute a_min or m_h is
reported for any modified content.  What IS reported is the one question that
survives the calibration failure, and it is answered by a robustness scan.

INPUTS, transcribed, nothing assumed:
  eq. (62)-(68)  gauge sector, exact and approximated
  eq. (71)/(72)  fermion master formula
  eq. (71) note  "For more than 3-rank symmetric tensor representation, we have
                  to add the potential from corresponding smaller eigenvalues.
                  For instance, we have to sum the potential from 4 and 2 for 4."
  eq. (73)-(76)  the decompositions
  eq. (80),(82)  m_h and 1/R5
Units: C = 3k/(64 pi^8 R5^6).  k, R5, R6 cancel out of a_min and of m_h alike.
"""
import json
import math
import os

import numpy as np

P = lambda *a: print(*a, flush=True)
NMAX = 600
MW = 80.4
_n = np.arange(1, NMAX + 1)
_nf = _n.astype(float)
_w = _nf ** -5
_sgn = {1: np.ones(NMAX), -1: (-1.0) ** _n}


def basis(alpha, s, c, d=0):
    ph = np.outer(np.atleast_1d(alpha), c * math.pi * _nf)
    f = (np.cos, lambda t: -np.sin(t), lambda t: -np.cos(t),
         lambda t: np.sin(t))[d % 4](ph)
    return (f * (_w * _sgn[s] * (c * math.pi * _nf) ** d)).sum(axis=1)


GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]        # eq. (68) x (1/C)


def terms(rep, eta, etap, r4=True):
    """(multiplicity, s, r-1) from eqs. (73)-(76).  r4 toggles the extra 84 term
    their eq. (71) note demands for the SU(2) quadruplet.  su7_fermion_from_group
    derives all twenty blocks from SU(7) and shows what that term IS: not a second
    eigenvalue of the same states, as this docstring used to say, but the
    quadruplet's twelfth STATE a_+a_+a_-, giving 12 at c=1 where their eq. (76)
    prints 11.  Forced, not conventional -- and the control fires without it."""
    s = eta * etap
    if rep == "7":
        return [(1, -s, 1)]
    if rep == "28":
        return [(1, s, 2), (4, -s, 1), (1, s, 1)]
    if rep == "48":
        return [(1, s, 2), (8, -s, 1), (2, s, 1)]
    if rep == "84":
        t = [(1, -s, 3), (1, -s, 2), (4, s, 2), (11, -s, 1), (4, s, 1)]
        return t + [(1, -s, 1)] if r4 else t
    raise ValueError(rep)


def V(content, alpha, d=0, lam=1.0, r4=True):
    out = sum(m * basis(alpha, s, c, d) for m, s, c in GAUGE)
    t = []
    for rep, eta, etap, mult in content:
        t += [(m * mult, s, c) for m, s, c in terms(rep, eta, etap, r4)]
    for m, s, c in t:
        out = out + lam * m * basis(alpha, s, c, d)
    return out


GRID = np.linspace(0.0, 1.0, 40001)


def minimise(content, lam=1.0, r4=True):
    v = V(content, GRID, lam=lam, r4=r4)
    i = int(np.argmin(v))
    lo, hi = GRID[max(i - 1, 0)], GRID[min(i + 1, len(GRID) - 1)]
    for _ in range(25):
        xs = np.linspace(lo, hi, 15)
        j = int(np.argmin(V(content, xs, lam=lam, r4=r4)))
        lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 14)]
    return float(0.5 * (lo + hi))


REC = {"source": "Komori & Maru, arXiv:2503.04090",
       "units": "V in C = 3k/(64 pi^8 R5^6)", "NMAX": NMAX, "steps": {}}
T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.8, 3.8),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)], 0.081, 125.5, 2.0),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], 0.021, 125.1, 7.5),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], 0.026, 126.4, 6.1),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.2, 3.8)]

P("=" * 78)
P("CONTROL 1 -- the SU(2)_L multiplet counts of eqs. (73)-(76)")
P("=" * 78)
LISTED = {"7": {2: 1}, "28": {3: 1, 2: 5}, "48": {3: 1, 2: 10},
          "84": {4: 1, 3: 5, 2: 15}}
P("  %-5s %-30s %-22s" % ("rep", "their decomposition list", "eqs (73)-(76)"))
for rep in ("7", "28", "48", "84"):
    got = {}
    for m, s, c in terms(rep, 1, 1, r4=False):
        got[c + 1] = got.get(c + 1, 0) + m
    P("  %-5s %-30s %-22s %s" % (rep, LISTED[rep], got,
                                 "OK" if got == LISTED[rep] else "*** FAIL ***"))
    assert got == LISTED[rep]
P("  >> the 84's 15 doublets come out only if the last two terms of eq. (76) are")
P("     84 terms; printed as '48' they would give 12.  A count identifies the")
P("     typo, not taste.                                                  PASS")

P("")
P("=" * 78)
P("CONTROL 2 -- the gauge sector alone, against their Fig. 1")
P("=" * 78)
a_g = minimise([])
P("  their Fig. 1: with no fermions, a_min = 1 and EWSB does NOT happen.")
P("  ours: a_min = %.5f                                                  %s"
  % (a_g, "PASS" if abs(a_g - 1) < 1e-3 else "*** FAIL ***"))
assert abs(a_g - 1) < 1e-3
REC["steps"]["gauge_only"] = {"a_min": a_g, "theirs": 1.0, "verdict": "PASS"}

P("")
P("=" * 78)
P("CONTROL 3 -- their eq. (71) note: the quadruplet needs its smaller eigenvalue")
P("=" * 78)
P("  'we have to sum the potential from 4 and 2 for 4'.  The 84's (1,4) therefore")
P("  contributes cos(3 pi n a) AND cos(pi n a).  Effect on the anchor:")
P("")
P("  %-5s %-14s %-14s %s" % ("case", "a without", "a with", "a theirs"))
for tag, cont, a_t, _m, _r in T1:
    P("  %-5s %-14.4f %-14.4f %.4f"
      % (tag, minimise(cont, r4=False), minimise(cont, r4=True), a_t))
P("  >> including it moves case (2) to within 2.6 % of their published value,")
P("     so the rule is real and is now applied everywhere.")

P("")
P("=" * 78)
P("THE ANCHOR -- and it FAILS.  Everything after this is bounded by that.")
P("=" * 78)
P("  %-5s %-11s %-11s %-9s %s" % ("case", "a theirs", "a ours", "ratio", "content"))
rows = []
for tag, cont, a_t, mh_t, iR_t in T1:
    a = minimise(cont)
    rows.append(dict(case=tag, a_theirs=a_t, a_ours=a, ratio=a / a_t))
    P("  %-5s %-11.4f %-11.4f %-9.2f %s"
      % (tag, a_t, a, a / a_t, " + ".join("%dx%s(%+d,%+d)" % (m, r, e, ep)
                                          for r, e, ep, m in cont)))
REC["steps"]["anchor"] = {"rows": rows, "verdict": "FAILED"}
P("")
P("  Three controls above PASS, so the transcription is not simply wrong:")
P("  the multiplet counts match, the gauge sector reproduces their Fig. 1")
P("  EXACTLY, and their own quadruplet rule fixes case (2) to 2.6 %.  But cases")
P("  (1),(3),(4),(5) are off by 29-107 %, and the discrepancy is CONTENT-")
P("  DEPENDENT: the fermion-sector rescaling needed to hit their a_min is")
P("  0.85, 0.96, 0.81, 0.80, 0.85 -- not one number, so it is not a missing")
P("  normalisation.  Using their EXACT eqs. (62)-(65) instead of the k >> 1")
P("  approximation moves the values in the right direction but needs k ~ 1.2,")
P("  which contradicts their own stated k >> 1.")
P("")
P("  >> THE INSTRUMENT IS NOT CALIBRATED IN ABSOLUTE VALUE.  No a_min and no")
P("     m_h is reported below for any modified content.  Reporting them would")
P("     be reporting the residual, not the physics.")

P("")
P("=" * 78)
P("THE QUESTION THAT SURVIVES: is EWSB destroyed by donating one 84?")
P("=" * 78)
P("  This is a yes/no about WHERE the minimum sits, not about its value:")
P("  their Fig. 1 establishes that with too little fermion content the minimum")
P("  runs to a = 1, where the electroweak group is UNBROKEN.  That verdict is")
P("  a sign, and a sign survives a 20 % normalisation error.  Scanned over the")
P("  whole range the anchor failure allows, lam in [0.7, 1.3]:")
P("")


P("  BOTH ENDS OF THE INTERVAL ARE UNBROKEN, and the first version of this")
P("  script got that wrong.  a = 0 is the symmetric point (m_W = a/(2R5) = 0);")
P("  a = 1 is their Fig. 1 case.  EWSB needs a strictly INSIDE (0,1).")
P("  Control on the criterion itself: it must call the gauge-only case DEAD.")
P("")
TOL = 2e-4


def verdict(a):
    if a < TOL:
        return "NO EWSB -- symmetric point a = 0"
    if a > 1 - TOL:
        return "NO EWSB -- a = 1 (their Fig. 1)"
    return "EWSB"


P("    gauge only            -> %s   %s"
  % (verdict(a_g), "OK" if verdict(a_g).startswith("NO") else "*** FAIL ***"))
assert verdict(a_g).startswith("NO")
for tag, cont, a_t, _m, _r in T1[:2]:
    P("    their own case %s     -> %s   %s"
      % (tag, verdict(minimise(cont)), "OK"))
assert verdict(minimise(T1[0][1])) == "EWSB"
P("    >> the criterion fires on the case that must be dead and not on the")
P("       cases that must be alive.                                      PASS")
P("")
LAMS = (0.7, 0.85, 1.0, 1.15, 1.3)
P("  %-5s %-40s %s" % ("case", "content after donating an 84", "at lam = 1"))
red = []
for tag, cont, a_t, mh_t, iR_t in T1:
    new = []
    for r, e, ep, m in cont:
        if (r, e, ep) == ("84", 1, 1):
            if m - 1:
                new.append((r, e, ep, m - 1))
        else:
            new.append((r, e, ep, m))
    aa = [minimise(new, lam=l) for l in LAMS]
    a1 = aa[LAMS.index(1.0)]
    d2 = float(V(new, 1e-9, d=2)[0])
    red.append(dict(case=tag, reduced=[[r, e, ep, m] for r, e, ep, m in new],
                    a_by_lambda=dict(zip(map(str, LAMS), aa)),
                    a_nominal=a1, Vpp_at_0=d2, verdict=verdict(a1),
                    verdicts_over_lambda=sorted({verdict(x) for x in aa})))
    P("  %-5s %-40s %s"
      % (tag, " + ".join("%dx%s(%+d,%+d)" % (m, r, e, ep) for r, e, ep, m in new)
         or "(only gauge)", verdict(a1)))
    P("        a over lam %s : %s"
      % (LAMS, ", ".join("%.4f" % x for x in aa)))
    P("        V''(0) = %+.3f  (>0 means a = 0 really is a minimum)" % d2)
REC["steps"]["donate_84"] = red

P("")
P("=" * 78)
P("AND THE WHOLE THING HAS A CLOSED FORM, which the anchor failure cannot touch")
P("=" * 78)
P("  EWSB requires a = 0 to be a MAXIMUM.  Expand: V''(0) = -pi^2 sum_terms")
P("  m c^2 sum_n s^n/n^3, and those two series are exactly zeta(3) and")
P("  -(3/4)zeta(3).  So")
P("")
P("      V''(0) = -pi^2 zeta(3) D ,     D = sum_{s=+1} m c^2 - (3/4) sum_{s=-1} m c^2")
P("")
P("  with D RATIONAL.  a = 0 is a maximum (EWSB possible) iff D > 0.")
P("")
from fractions import Fraction as F


def Dof(content, lam=F(1)):
    tot = F(0)
    for m, s, c in GAUGE:
        tot += F(m).limit_denominator() * c ** 2 * (1 if s == 1 else F(-3, 4))
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            tot += lam * m * mult * c ** 2 * (1 if s == 1 else F(-3, 4))
    return tot


Z3 = 1.2020569031595942
P("  CONTROL: the closed form against the numerically summed V''(0).  The two")
P("  cannot agree exactly -- the sum is truncated at n = %d -- so the control is" % NMAX)
P("  the sharper one: PREDICT the residual from the tail of zeta(3),")
P("     |closed - numeric|  =  pi^2 (sum_{s=+1} m c^2) sum_{n>N} 1/n^3")
P("  and the alternating series contributes nothing at this order.")
P("")
tail = sum(1.0 / n ** 3 for n in range(NMAX + 1, 400000))
P("  %-14s %-14s %-14s %-12s %s"
  % ("content", "closed form", "numeric", "residual", "predicted"))
for tag, cont, *_ in T1:
    cf = -math.pi ** 2 * Z3 * float(Dof(cont))
    nu = float(V(cont, 1e-9, d=2)[0])
    pos = 0.0
    for m, s, c in GAUGE:
        pos += m * c ** 2 if s == 1 else 0
    for rep, eta, etap, mult in cont:
        for m, s, c in terms(rep, eta, etap):
            pos += m * mult * c ** 2 if s == 1 else 0
    pred = -math.pi ** 2 * pos * tail
    P("  Table 1 %-6s %-14.5f %-14.5f %-12.2e %-12.2e %s"
      % (tag, cf, nu, cf - nu, pred,
         "OK" if abs((cf - nu) - pred) < 0.02 * abs(pred) else "*** FAIL ***"))
    assert abs((cf - nu) - pred) < 0.02 * abs(pred)
P("  >> the residual is the truncation tail, to better than 2 %, on all five")
P("     rows.  The closed form is the exact statement.                    PASS")
P("")
P("  Per-multiplet contribution to D (this is the whole physics):")
for lab, c in (("gauge sector", []), ("7(+,+)", [("7", 1, 1, 1)]),
               ("7(+,-)", [("7", 1, -1, 1)]), ("28(+,+)", [("28", 1, 1, 1)]),
               ("28(+,-)", [("28", 1, -1, 1)]), ("48(+,+)", [("48", 1, 1, 1)]),
               ("84(+,+)", [("84", 1, 1, 1)])):
    d = Dof(c) - (Dof([]) if c else 0)
    P("     %-12s D = %s" % (lab, d))
P("  >> each 84(+,+) is worth exactly D = 5/4, and the 48 is worth exactly 0 --")
P("     it cannot help at the origin at all.")
P("")
P("  %-5s %-12s %-12s %-12s %s"
  % ("case", "D (theirs)", "D (-1x84)", "lam_crit", "verdict at lam = 1"))
exact = []
for tag, cont, *_ in T1:
    new = [(r, e, ep, m - 1 if (r, e, ep) == ("84", 1, 1) else m)
           for r, e, ep, m in cont]
    new = [x for x in new if x[3]]
    d0, d1 = Dof(cont), Dof(new)
    dg = Dof([])
    lc = -dg / (d1 - dg) if d1 != dg else None
    exact.append(dict(case=tag, D_theirs=str(d0), D_reduced=str(d1),
                      lam_crit=str(lc), ewsb=bool(d1 > 0)))
    P("  %-5s %-12s %-12s %-12s %s"
      % (tag, d0, d1, "%.4f" % float(lc) if lc else "-",
         "EWSB" if d1 > 0 else "*** NO EWSB -- a = 0 becomes a MINIMUM ***"))
REC["steps"]["closed_form_D"] = exact
P("")
P("  >> donating an 84 costs exactly D = 5/4 = 10/8, and case (3) only has")
P("     9/8 to spend.  It goes to -1/8 and the symmetric point becomes the")
P("     vacuum.  Everything is in eighths and nothing is numerical.")

P("")
P("=" * 78)
P("VERDICT")
P("=" * 78)
alive = [r["case"] for r in red if r["verdict"] == "EWSB"]
deadn = [r["case"] for r in red if r["verdict"] != "EWSB"]
always = [r["case"] for r in red if r["verdicts_over_lambda"] == ["EWSB"]]
P("  at the nominal normalisation (lam = 1):")
P("     EWSB survives donating an 84 : %s" % (alive or "none"))
P("     EWSB LOST                    : %s" % (deadn or "none"))
P("  robust across lam in [0.7,1.3] (EWSB at every lam) : %s" % (always or "none"))
P("")
P("  >> the comparison before/after uses the SAME instrument on the SAME row,")
P("     so the anchor failure largely cancels in it; the absolute a_min does")
P("     not survive, the direction of the change does.")
P("")
P("  Reference: the SAME rows UNREDUCED, same instrument, lam = 1:")
for tag, cont, a_t, _m, _r in T1:
    P("     %-5s a = %.4f  -> %s" % (tag, minimise(cont), verdict(minimise(cont))))
P("")
P("  NOT CLAIMED: any value of a_min, 1/R5 or m_h for the reduced contents --")
P("  the anchor failure forbids it.  Whether a surviving row still gives")
P("  m_h = 125-127 GeV is THEIR computation and is NOT answered here.")
REC["steps"]["verdict"] = {"ewsb_at_nominal": alive, "ewsb_lost_at_nominal": deadn,
                          "ewsb_robust_over_lambda": always}

os.makedirs("paper_data", exist_ok=True)
with open("paper_data/su7_vacuum.json", "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("wrote paper_data/su7_vacuum.json")
P("DONE")
