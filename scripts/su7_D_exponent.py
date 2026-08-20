#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

THE ONE NUMBER THE WHOLE SELECTION HANGS ON, AND WHETHER IT IS SAFE.

Ask it Socratically, to the root.

  Q. What selects case (2) as the unique survivor?
  A. Two conditions. The nu_R condition is anomalies -- rational, no potential.
     The vacuum condition is D > 0, with the escape costing D(84) = 5/4 and case
     (3) having only 9/8 to spend.  "Everything is in eighths."

  Q. Then the selection never touches alpha_min, the column we fail on?
  A. It does not.  That is why the paper's conclusion survives the anchor
     failure.

  Q. But D is not pure group theory either. Where does its 3/4 come from?
  A. V''(0) = -pi^2 sum m c^2 sum_n s^n/n^3.  The two series are zeta(3) and
     -eta(3), and eta(3)/zeta(3) = 1 - 2^{1-3} = 3/4.

  Q. And where does the 3 come from?
  A. From the potential going as 1/n^5: two alpha-derivatives bring down n^2.

  Q. So the 3/4 is a consequence of the EXPONENT of the potential -- which is a
     property of the same instrument whose normalisation we cannot match. If the
     sum went as 1/n^6, the coefficient would be 1 - 2^{1-4} = 7/8, and the
     margin the selection rests on is one eighth. Does case (2) survive that?

That last question is the root, and it is the one this script answers. It is not
rhetorical: Part V's own AHMN gate had a live |k|^-6 reading alongside the |k|^-5
one, and both reproduced that paper's vacuum to three digits. So the exponent is
exactly the kind of thing a vacuum position does not pin down.

The test sweeps the coefficient over every exponent p that could plausibly sit in
such a potential, recomputes the five rows' D and the escape cost, and asks
whether "case (2) is the unique survivor" is a statement about SU(7) or a
statement about our choice of 5.
"""
import math
import re
from fractions import Fraction as F

P = lambda *a: print(*a, flush=True)

_src = open("su7_vacuum.py", encoding="utf-8").read()
_fn = re.search(r"\ndef terms\(.*?\n(?=\n\ndef |\n\nGRID|\Z)", _src, re.S)
_ns = {}
exec(_fn.group(0), _ns)
terms = _ns["terms"]

GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]        # their eq. (68) x (1/C)
T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)]),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)]),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)]),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)]),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)])]


def eta_over_zeta(p):
    """the alternating series relative to the plain one at exponent p:
    eta(p)/zeta(p) = 1 - 2^{1-p}.  Rational for every integer p."""
    return 1 - F(1, 2 ** (p - 1))


def Dof(content, w):
    """D at alternating-weight w.  w = 3/4 is the potential going as 1/n^5."""
    tot = F(0)
    for m, s, c in GAUGE:
        tot += F(m).limit_denominator() * c ** 2 * (1 if s == 1 else -w)
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            tot += m * mult * c ** 2 * (1 if s == 1 else -w)
    return tot


COST = lambda w: Dof([("84", 1, 1, 1)], w) - Dof([], w)     # donating one 84(+,+)

P("=" * 78)
P("THE EXPONENT, AND WHAT IT DOES TO THE SELECTION")
P("=" * 78)
P("  potential ~ 1/n^P  =>  V''(0) ~ 1/n^(P-2)  =>  w = eta(P-2)/zeta(P-2)")
P("")
P("  %-4s %-8s %-9s %-46s %s" % ("P", "w", "cost", "D of the five rows after donating one 84", "survivors"))
rows = {}
for Pw in range(4, 10):
    p = Pw - 2
    if p < 2:
        continue
    w = eta_over_zeta(p)
    cost = COST(w)
    after = [Dof(c, w) - cost for _, c in T1]
    surv = [t for (t, _), d in zip(T1, after) if d > 0]
    rows[Pw] = (w, cost, after, surv)
    P("  %-4d %-8s %-9s %-46s %s"
      % (Pw, w, cost, "  ".join("%s" % d for d in after), " ".join(surv) or "NONE"))

P("")
P("  the published reading is P = 5 (w = 3/4), the first row above with p = 3.")
P("")
P("=" * 78)
P("BOTH CONDITIONS, NOT ONE -- AND THE WINDOW IN w")
P("=" * 78)
P("  The paper's claim is not 'case (2) is the only row with D > 0'. It is that")
P("  case (2) is the only row satisfying BOTH necessary conditions. The nu_R")
P("  condition of the anomaly section admits only rows (2) and (3), and carries")
P("  no w at all: it is anomalies and tensor content. So the question is")
P("  narrower than the sweep above --")
P("")
P("      for which w does the vacuum condition keep (2) and kill (3)?")
P("")
P("  D_i(w) - cost(w) is LINEAR in w, so this is two inequalities and the answer")
P("  is an interval with exact rational endpoints.")
P("")


def AB(content):
    """(A, B) with D = A - w B."""
    A = B = F(0)
    for m, s, c in GAUGE:
        (A if s == 1 else B).__class__          # keep it explicit below
        if s == 1:
            A += F(m).limit_denominator() * c ** 2
        else:
            B += F(m).limit_denominator() * c ** 2
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            if s == 1:
                A += m * mult * c ** 2
            else:
                B += m * mult * c ** 2
    return A, B


A84, B84 = AB([("84", 1, 1, 1)])
A0, B0 = AB([])
dA, dB = A84 - A0, B84 - B0                     # the escape cost, as A - w B


def net(content):
    """D_i - cost = (A_i - dA) - w (B_i - dB), returned as the pair."""
    A, B = AB(content)
    return A - dA, B - dB


n2, n3 = net(T1[1][1]), net(T1[2][1])
P("  case (2) after donating:  D = %s - %s w" % (n2[0], n2[1]))
P("  case (3) after donating:  D = %s - %s w" % (n3[0], n3[1]))
P("")
# (2) lives: A2 - w B2 > 0.  (3) dies: A3 - w B3 <= 0.
lo2 = F(n2[0], n2[1]) if n2[1] > 0 else None    # (2) alive for w < lo2
lo3 = F(n3[0], n3[1]) if n3[1] > 0 else None    # (3) dead  for w >= lo3
P("  (2) survives while w %s %s" % ("<" if n2[1] > 0 else ">", lo2))
P("  (3) is killed once w %s %s" % (">=" if n3[1] > 0 else "<=", lo3))
P("")
lo, hi = lo3, lo2
P("  >> the headline holds exactly on   w in [%s, %s)   =   [%.4f, %.4f)"
  % (lo, hi, float(lo), float(hi)))
W5 = eta_over_zeta(3)
inside = lo <= W5 < hi
P("     the published reading sits at w = %s = %.4f, which is inside: %s"
  % (W5, float(W5), inside))
assert inside
P("     margin to the lower edge %.4f, to the upper edge %.4f"
  % (float(W5 - lo), float(hi - W5)))
P("")
P("  And now the exponents, read off that interval rather than re-swept:")
for Pw in range(4, 9):
    w = eta_over_zeta(Pw - 2)
    P("     P = %d   w = %-6s  %s"
      % (Pw, w, "HEADLINE HOLDS" if lo <= w < hi else
         ("case (3) also survives" if w < lo else "no row survives at all")))

P("")
P("=" * 78)
P("VERDICT")
P("=" * 78)
P("  1. The selection is NOT exponent-free. It holds on w in [%s, %s), and of"
  % (lo, hi))
P("     the admissible exponents only P = 5 lands inside: at P = 4 case (3)")
P("     survives too and the row is not unique; at P >= 6 NO row survives at all")
P("     and the minimal extension dies outright. So their eq. (67) -- which we")
P("     prove rather than assume, the brace closing to 3k/4(pi n)^5 with Sage")
P("     simplifying the difference to exactly 0 -- is load-bearing for the")
P("     SELECTION and not only for the instrument.")
P("     Section 7 files that proof under 'what the anchor failure cannot be'.")
P("     It belongs in section 5's dependency list too, and it is not there.")
P("")
P("  2. And the tight-looking lower margin of %.4f is NOT a fragility, which is"
  % float(W5 - lo))
P("     worth stating because it would be easy to sell it as one. w is not a")
P("     free parameter that could drift into the gap: w = 1 - 2^(1-p) takes only")
P("     the values 1/2, 3/4, 7/8, ..., and the nearest admissible neighbours of")
P("     3/4 are 1/2 below and 7/8 above -- both far outside the window.")
P("     The window is narrow in a continuum that the")
P("     physics does not live in. What CAN move the conclusion is a change to")
P("     (A, B) -- new terms, two loops -- and that is the repair wedge of")
P("     section 8, which is measured there and is a different quantity.")
P("")
P("  So: one genuine gap in the ledger, and one alarming-looking number that is")
P("  an artifact of sweeping a quantised constant continuously. Both had to be")
P("  computed to tell them apart.")
P("")
P("  What this does NOT rescue: alpha_min itself. The exponent is exactly the")
P("  kind of structural difference that could produce a row-dependent ratio,")
P("  and this script does not test that -- it tests only the SELECTION. The")
P("  anchor question stays open and stays where section 7 leaves it.")
