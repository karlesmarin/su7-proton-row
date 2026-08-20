#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

ADVERSARIAL AUDIT of su7_qphi.py, run BEFORE the note exists.  Five attacks, each
aimed at an assumption that carries weight rather than at the arithmetic.

A1  is "A/q_phi in Z" the right survival condition at all?
A2  the strict reading was killed by ANOMALY CANCELLATION.  Green-Schwarz is
    available (su7_channels_socratic.py A4).  Turn it on and the strict reading
    comes back -- does the verdict then invert?
A3  the whole result is stated in charges read off ONE normalisation of U(1)'
    (their eq. 79).  Does it survive rescaling the generator?
A4  ONE breaking scalar was assumed.  What does a second one do?
A5  is the residual group really exact, or only a leading-order selection rule?
"""
from fractions import Fraction as F
from math import gcd as igcd

P = lambda *a: print(*a, flush=True)
H = F(1, 2)

fails_at = lambda As, q: [j for j, A in enumerate(As) if (A / q).denominator == 1]


def fgcd(a, b):
    a, b = abs(F(a)), abs(F(b))
    if a == 0 or b == 0:
        return a or b
    return F(igcd(a.numerator * b.denominator, b.numerator * a.denominator),
             a.denominator * b.denominator)


P("=" * 78)
P("A1 -- is the survival condition right?")
P("=" * 78)
P("  An operator of charge A is made invariant by dressing with n insertions of")
P("  <phi> (charge +q_phi) or <phi*> (charge -q_phi):  A + n q_phi = 0 for some")
P("  integer n of either sign  <=>  A/q_phi in Z.  Two ways it could be wrong:")
P("")
P("  (i)  phi REAL.  A real scalar cannot carry a U(1) charge at all -- a real")
P("       representation of U(1) is the trivial one -- so a U(1)'-breaking")
P("       scalar is necessarily complex and both signs are available.")
P("  (ii) the residual group.  U(1) acts on charge Q by exp(i theta Q).  <phi>")
P("       survives iff theta q_phi in 2 pi Z, so the unbroken set is")
P("       theta in (2 pi / q_phi) Z, which acts on charge A by exp(2 pi i A/q_phi).")
P("       Invariant for every element iff A/q_phi in Z.  Same condition.")
P("")
P("  >> A1: the condition is the definition, not a model.  SURVIVES.")

P("")
P("=" * 78)
P("A2 -- turn Green-Schwarz ON and give the strict reading back its life")
P("=" * 78)
P("  su7_qphi.py STEP 3 kills the strict reading because anomaly cancellation")
P("  forces X_Q = -1/6, which is off the bulk lattice.  But su7_channels_")
P("  socratic.py A4 states that mixed U(1) anomalies are the EVADABLE class.")
P("  If a Green-Schwarz two-form absorbs them, X_Q is free again -- and under")
P("  the strict reading it is then a HALF-INTEGER.  What happens to the proton?")
P("")
P("  A = 3 X_Q + 1/2 (su7_XQ.py).  For X_Q = m/2 that is A = (3m+1)/2.")
P("")
P("  %-10s %-10s %-14s %s" % ("X_Q", "A", "A = 0?", "q_phi that FAIL, of {1/2, 1, 3/2}"))
never_zero = True
always84 = True
for m in range(-4, 5):
    xq = F(m, 2)
    A = 3 * xq + H
    never_zero &= (A != 0)
    bad = [q for q in (H, F(1), F(3, 2)) if fails_at([A], q)]
    always84 &= (F(3, 2) not in bad)
    P("  %-10s %-10s %-14s %s"
      % (xq, A, "yes -- UNPROTECTED" if A == 0 else "no",
         ", ".join(str(q) for q in bad) if bad else "none"))
P("")
assert never_zero
P("  >> A = (3m+1)/2 is NEVER zero, because 3m+1 is never zero.  So under the")
P("     strict reading the protection cannot fail the way it fails in their")
P("     model as written: the dangerous point is simply not on the lattice.")
P("")
assert always84
P("  >> AND THE ATTACK PRODUCED A THEOREM.  A/q_phi with q_phi = 3/2 is")
P("     (3m+1)/3, and 3m+1 is never divisible by 3.  So:")
P("")
P("        q_phi = 3/2 forbids all four dimension-6 |dB|=1 operators to all")
P("        orders, for EVERY half-integer brane-quark charge, with no family")
P("        dependence and no assumption about the anomalies.")
P("")
P("     3/2 is the charge of the (777) component of the 84 -- a multiplet every")
P("     row of their Table 1 already contains.  Checked to |m| <= 4 above and")
P("     proved for all m by the mod-3 argument.")
P("")
P("  CONTROL that had to fire: the same test at q_phi = 1/2 and 1 must NOT be")
P("  universal, or the statement would be empty.")
b12 = [m for m in range(-4, 5) if fails_at([3 * F(m, 2) + H], H)]
b1 = [m for m in range(-4, 5) if fails_at([3 * F(m, 2) + H], F(1))]
P("     q_phi = 1/2 fails at m = %s" % b12)
P("     q_phi = 1   fails at m = %s" % b1)
assert b12 and b1
P("     both non-empty: the 3/2 result is a real discrimination.  SURVIVES,")
P("     and it strengthens rather than weakens the note.")

P("")
P("=" * 78)
P("A3 -- the verdict must not depend on how U(1)' is normalised")
P("=" * 78)
P("  Their eq. (79) fixes a normalisation.  Any other differs by an overall")
P("  factor c, which multiplies A_j and q_phi alike.  A/q_phi is then invariant.")
P("")
CASES = {"minimal l=(1/2,1/2,0)": [F(1, 6), F(1, 6), F(-1, 3)],
         "strict  l=(1/2,1/2,-1)": [H, H, F(-1)]}
P("  %-24s %-10s %-10s %s" % ("case", "c", "q_phi", "generations dressable"))
inv = True
for name, As in CASES.items():
    base = None
    for c in (F(1), F(2), F(1, 3), F(7, 5), F(-1)):
        row = [len(fails_at([c * A for A in As], c * q)) for q in (H, F(1), F(3, 2))]
        if base is None:
            base = row
        inv &= (row == base)
        P("  %-24s %-10s %-10s %s" % (name, c, "c*{1/2,1,3/2}", row))
assert inv
P("")
P("  >> every verdict is unchanged under every rescaling, including c < 0.")
P("     A3: SURVIVES.  The result is about ratios, which is why it can be")
P("     stated at all without knowing their coupling normalisation.")

P("")
P("=" * 78)
P("A4 -- a SECOND breaking scalar, which the paper also does not exclude")
P("=" * 78)
P("  Two VEVs of charges q1, q2 leave the subgroup that fixes BOTH, i.e. the")
P("  one generated by gcd(q1, q2).  So protection needs A_j / gcd not in Z --")
P("  strictly harder.  The worst case for the strict-compatible solution:")
P("")
As = CASES["strict  l=(1/2,1/2,-1)"]
P("  %-22s %-10s %s" % ("charges present", "gcd", "generations dressable"))
for pair in [(F(3, 2),), (F(3, 2), F(1)), (F(3, 2), H), (F(3, 2), F(1), H)]:
    g = F(0)
    for q in pair:
        g = fgcd(g, q)
    P("  %-22s %-10s %s"
      % (", ".join(str(q) for q in pair), g, fails_at(As, g) or "PROTECTED"))
assert fails_at(As, fgcd(F(3, 2), H))
P("")
P("  >> adding the 7's scalar to the 84's destroys the protection the 84 gave.")
P("     So the claim is not 'the 84 protects' but 'the 84 protects IF it is the")
P("     only U(1)'-breaking VEV'.  A4 does not kill the result; it names the")
P("     hypothesis the result needs, and su7_qphi.py already states it under")
P("     NOT CLAIMED.  SURVIVES, with the scope written in.")

P("")
P("=" * 78)
P("A5 -- exact, or only a leading-order selection rule?")
P("=" * 78)
P("  The residual subgroup of a spontaneously broken GAUGE U(1) is a discrete")
P("  GAUGE symmetry (Krauss-Wilczek 1989), not an accidental global one.  It is")
P("  therefore not violated by gravity or by instantons, and it holds at every")
P("  order in <phi>/M.  This is the whole reason the mechanism is used in the")
P("  literature (JHEP 07 (2008) 065; hep-ph/0012092) and the reason the answer")
P("  to their question can be 'forbidden' rather than 'suppressed'.")
P("")
P("  What it does NOT cover, and the note must say so: operators of dimension")
P("  higher than 6, whose charges are different and were not enumerated here;")
P("  and any charge assignment for fields beyond the ones their eqs. (43)-(47)")
P("  and (76) put in the theory.")
P("")
P("  >> A5: SURVIVES with scope.")

P("")
P("=" * 78)
P("VERDICT")
P("=" * 78)
P("""
  Five attacks, none fatal, and A2 produced the strongest single statement of
  the line:

     q_phi = 3/2 -- the (777) of an 84, which every row of their Table 1 already
     contains -- forbids all four dimension-6 |dB|=1 operators to all orders for
     every half-integer brane-quark charge, provided it is the only U(1)'-
     breaking VEV.

  The dependencies, all named: one breaking scalar (A4), the survival condition
  as defined (A1), dimension 6 (A5).  The verdict is normalisation-independent
  (A3).
""")
P("DONE")
