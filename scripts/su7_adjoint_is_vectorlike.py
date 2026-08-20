#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

THE BRANE ROUTE IS CLOSED, AND NOT BY THE YUKAWA.

su7_brane_partner.py left the paper's headline resting on a price rather than a
prohibition: the 48 carries a lepton doublet at rung 2, donating a 48 costs D = 0
exactly, the assignment that needs rung 2 is anomaly-free, and their own quarks
are already brane fields -- so a brane-localised charged-lepton partner would let
case (3) survive at no cost and the headline would be false. The only thing
standing in the way was that such a Yukawa would not be gauge-generated.

That was the wrong place to look. The obstruction is not the partner. It is the
DOUBLET.

  The 48 is the adjoint: 7 x 7bar - 1, a REAL representation. Its component
  (i, jbar) has conjugate (j, ibar), which is in the SAME multiplet. Under the
  orbifold, a component survives iff P5P5' = eta*eta' and carries 4D chirality
  -eta*P5. For a component and its conjugate,

      P5P5'(i,jbar) = p_i p_j = P5P5'(j,ibar)      identical
      P5  (i,jbar) = P5_i P5_j = P5  (j,ibar)      identical

  so they survive together, at the SAME chirality. Written all-left-handed, a
  state and its conjugate at the same chirality is precisely a vector-like Dirac
  pair, and every gauge quantum number of the pair sums to zero -- so a Dirac
  mass is allowed by everything, and their own prescription (a conjugate brane
  fermion against every unwanted zero mode) supplies one.

A vector-like pair is not a generation. No brane field repairs that, because the
brane field would be the partner and the defect is upstream of the partner.

This is checked here as a statement about the whole multiplet, not about one
component, and against the complex representations as the control that must
distinguish them -- if the 28 and the 84 came out vector-like too, the test would
be measuring the orbifold and not the reality of the representation, and their
model could not have chiral matter at all.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement as multiset

P = lambda *a: print(*a, flush=True)
H = F(1, 2)
N = 7

P5 = [1, 1, 1, 1, 1, -1, -1]                    # their eq. (37)-(40), per index
P5p = [1, 1, 1, -1, -1, -1, 1]
YIND = [F(0), F(0), F(0), H, H, F(-1), F(0)]    # their eq. (41)
EIND = [F(0), F(0), F(0), F(0), F(0), H, -H]


def adjoint():
    """components (i,jbar) of 7 x 7bar, with the Cartan left aside (q = Y = 0)."""
    return [(i, j) for i in range(N) for j in range(N) if i != j]


def sym(d):
    """components of Sym^d(7) as multisets of indices."""
    return list(multiset(range(N), d))


def qn_adj(c):
    i, j = c
    return dict(p5=P5[i] * P5[j], p5p=P5p[i] * P5p[j],
                Y=YIND[i] - YIND[j], extra=EIND[i] - EIND[j])


def qn_sym(c):
    p5 = p5p = 1
    for i in c:
        p5 *= P5[i]
        p5p *= P5p[i]
    return dict(p5=p5, p5p=p5p, Y=sum(YIND[i] for i in c),
                extra=sum(EIND[i] for i in c))


def survivors(comps, qn, eta, etap):
    """(component, 4D chirality) for what the orbifold keeps."""
    out = []
    for c in comps:
        q = qn(c)
        if q["p5"] * q["p5p"] == eta * etap:
            out.append((c, -eta * q["p5"], q))
    return out


P("=" * 78)
P("A -- IS THE SURVIVING 4D SPECTRUM OF A 48 CHIRAL?")
P("=" * 78)
P("  For every surviving component, look for a survivor with EXACTLY opposite")
P("  Y and extra at the SAME 4D chirality. That pair can be given a gauge-")
P("  invariant Dirac mass, so neither member is chiral matter.")
P("")
for eta, etap in ((1, 1), (1, -1)):
    S = survivors(adjoint(), qn_adj, eta, etap)
    paired = 0
    for c, ch, q in S:
        if any(ch2 == ch and q2["Y"] == -q["Y"] and q2["extra"] == -q["extra"]
               for c2, ch2, q2 in S if c2 != c):
            paired += 1
    P("  48(%+d,%+d): %3d survivors, %3d of them have a same-chirality conjugate"
      % (eta, etap, len(S), paired))
    P("               UNPAIRED (i.e. chiral): %d" % (len(S) - paired))

P("")
P("  And the reason it can never be otherwise, which is one line:")
P("  the conjugate of (i,jbar) is (j,ibar), and")
P("      P5P5'(j,ibar) = p_j p_i = P5P5'(i,jbar)   -> survives with it")
P("      P5  (j,ibar) = P5_j P5_i = P5  (i,jbar)   -> at the same chirality")
P("  It is not a property of this Wilson line or of these parities. It is that")
P("  the adjoint is its own conjugate.")
ok = True
for i, j in adjoint():
    ok &= (qn_adj((i, j))["p5"] == qn_adj((j, i))["p5"]
           and qn_adj((i, j))["p5p"] == qn_adj((j, i))["p5p"])
P("  verified on all %d components: %s" % (len(adjoint()), ok))
assert ok

P("")
P("=" * 78)
P("B -- THE CONTROL: THE COMPLEX REPRESENTATIONS MUST COME OUT CHIRAL")
P("=" * 78)
P("  If Sym^2 and Sym^3 were vector-like too, this test would be measuring the")
P("  orbifold rather than the reality of the representation -- and their model")
P("  could not have chiral matter at all. It must distinguish them.")
P("")
P("  %-10s %-8s %-12s %-14s %s" % ("rep", "real?", "survivors", "same-chir conj", "chiral?"))
for nm, comps, qn, real in (("48", adjoint(), qn_adj, True),
                            ("28=Sym2", sym(2), qn_sym, False),
                            ("84=Sym3", sym(3), qn_sym, False)):
    S = survivors(comps, qn, 1, 1)
    paired = sum(1 for c, ch, q in S
                 if any(ch2 == ch and q2["Y"] == -q["Y"] and q2["extra"] == -q["extra"]
                        for c2, ch2, q2 in S if c2 != c))
    P("  %-10s %-8s %-12d %-14d %s"
      % (nm, "yes" if real else "no", len(S), paired,
         "NO -- vector-like" if paired == len(S) else "yes, %d chiral states" % (len(S) - paired)))

P("")
P("=" * 78)
P("C -- WHAT THIS DOES TO THE BRANE QUESTION")
P("=" * 78)
P("  su7_brane_partner.py asked whether a brane-localised charged-lepton partner")
P("  could complete the 48's rung-2 doublet. The answer is that there is nothing")
P("  to complete. The doublet comes with its own conjugate, at the same 4D")
P("  chirality, and the pair is massive. A brane field would be the partner; the")
P("  defect is one step earlier, in the doublet itself.")
P("")
P("  So the exclusion of the 48 in sect. 4 does NOT rest on the Yukawa being")
P("  gauge-generated, and the headline does not need that precondition. It rests")
P("  on the 48 being a real representation, which no model-building choice can")
P("  change and no localised field can repair.")
P("")
P("  It also explains their own sentence, which we had taken as a stipulation:")
P("  the 7, 28, 48 and 84 are introduced with the Higgs potential as their only")
P("  role. For the 48 that is not a choice -- it is the only thing a real")
P("  representation CAN do here.")
P("")
P("  And it generalises past this paper: on this orbifold no real representation")
P("  can furnish a chiral generation, so the rung ladder of sect. 4 never has to")
P("  consider one. The box count and the reality condition between them cover")
P("  every multiplet their model contains.")
