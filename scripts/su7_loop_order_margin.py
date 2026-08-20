#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

HOW MUCH LOOP ORDER THE REPAIR WEDGE CAN ABSORB.

REPAIR_SPACE.md defends the headline over per-representation reweightings of the
FERMION sector, holding the gauge part at its exact one-loop value -27/8.  That
is right at one loop: the gauge part is not fitted, it is rebuilt from their
eq. (62) with their eq. (67) proved.

But CROSS_PAPER.md section 3 names loop order as the leading remaining candidate
for the anchor residual, and a two-loop term does NOT respect that split: in pure
gauge it multiplies the gauge sector by 1 - 5 g^2 C2(A)/(16 pi^2) (Dumitru-Guo-
Korthals Altes, eq. 66), while with fermions the proportionality is known to fail
(Guo-Du).  So the two sectors are rescaled by different factors, and only their
RATIO matters to the verdict.

This computes what the wedge tolerates in that direction, exactly, and puts it
next to the only size estimate anybody has.  It changes no claim of the paper; it
measures one.
"""
from fractions import Fraction as F
import math

P = lambda *a: print(*a, flush=True)

# ---------------------------------------------------------------- the wedge
# after donating one 84(+,+), with w_g on the gauge part and w_f on the fermions:
#    row (3) must die  :  -(27/8) w_g + (2 + 5/4)  w_f < 0
#    row (2) must live :  -(27/8) w_g + (2 + 15/4) w_f > 0
GAUGE = F(27, 8)
ROW3 = F(2) + F(5, 4)                 # one 28(+,+) + one 84(+,+) left
ROW2 = F(2) + F(15, 4)                # one 28(+,+) + three 84(+,+) left

P("=" * 78)
P("A -- THE WEDGE IS A STATEMENT ABOUT ONE RATIO")
P("=" * 78)
P("  row (3) dies  <=>  w_f/w_g <  %s / %s = %s" % (GAUGE, ROW3, GAUGE / ROW3))
P("  row (2) lives <=>  w_f/w_g >  %s / %s = %s" % (GAUGE, ROW2, GAUGE / ROW2))
lo, hi = GAUGE / ROW2, GAUGE / ROW3
P("")
P("  >> the verdict holds iff   w_f/w_g  in  (%s, %s) = (%.5f, %.5f)"
  % (lo, hi, lo, hi))
assert (lo, hi) == (F(27, 46), F(27, 26)), (lo, hi)
P("     which is REPAIR_SPACE.md's diagonal interval, as it must be: on the")
P("     diagonal w_g = 1 and w_f = w.  PASS")
P("")
P("  the two margins from the unrescaled point w_f/w_g = 1 are NOT symmetric:")
P("     upward   %s - 1 = %s = %.3f %%" % (hi, hi - 1, float(hi - 1) * 100))
P("     downward 1 - %s = %s = %.3f %%" % (lo, 1 - lo, float(1 - lo) * 100))
assert hi - 1 == F(1, 26) and 1 - lo == F(19, 46)

P("")
P("  CONTROL that had to fire: a COMMON rescaling of both sectors must leave")
P("  the verdict untouched, since it is the sign of a homogeneous expression.")
for c in (F(1, 10), F(1, 2), F(1), F(3), F(100)):
    d3 = -GAUGE * c + ROW3 * c
    d2 = -GAUGE * c + ROW2 * c
    assert d3 < 0 < d2
P("     rescaled by 1/10, 1/2, 1, 3, 100: row (3) dead and row (2) alive in all")
P("     five.  PASS -- so nothing below is about an overall normalisation.")

# ------------------------------------------------- what a gauge-only shift costs
P("")
P("=" * 78)
P("B -- A GAUGE-ONLY TWO-LOOP FACTOR, AND WHERE IT BREAKS THE VERDICT")
P("=" * 78)
P("  Suppress the gauge sector alone by delta, i.e. w_g = 1 - delta, w_f = 1.")
P("  The ratio becomes 1/(1-delta), and it leaves the wedge when")
P("     1/(1-delta) = %s   <=>   delta = 1 - %s = %s" % (hi, 1 / hi, 1 - 1 / hi))
dbreak = 1 - 1 / hi
assert dbreak == F(1, 27), dbreak
P("")
P("  >> THE VERDICT TOLERATES A GAUGE-ONLY SUPPRESSION OF AT MOST 1/27 = %.4f %%"
  % (float(dbreak) * 100))
P("     Beyond it row (3) survives the donation too, and case (2) stops being")
P("     unique.  (It is the SAME 1/26 margin as above, written the other way.)")

P("")
P("  The only estimate of that factor anybody has is the pure-gauge two-loop")
P("  ratio of Dumitru-Guo-Korthals Altes, 5 g^2 C2(A)/(16 pi^2), with C2(A) = 7")
P("  for SU(7).  Against the threshold:")
P("")
P("  %-8s %-14s %-14s %s" % ("g_4", "delta", "ratio 1/(1-d)", "verdict"))
for g in (0.55, 0.63, 0.70, 0.80):
    d = 5 * g * g * 7 / (16 * math.pi ** 2)
    r = 1 / (1 - d)
    P("  %-8.2f %-14.4f %-14.4f %s"
      % (g, d, r, "INSIDE the wedge" if r < float(hi) else "OUTSIDE -- verdict flips"))
P("")
P("  >> the estimate is 6.70 % to 14.19 % over the plausible coupling range,")
P("     against a threshold of 3.70 %.  Between 1.8 and 3.8 times the margin,")
P("     at every value of g_4.")

P("")
P("=" * 78)
P("WHAT THIS DOES AND DOES NOT SAY")
P("=" * 78)
P("  It does NOT say the verdict is wrong.  The fermion sector receives its own")
P("  two-loop correction, Guo-Du show it is not the same multiple, and nobody")
P("  has computed either one for this model.  If the two track each other to")
P("  within 4 % the ratio never moves and the wedge never notices.")
P("")
P("  What it says is that the wedge of REPAIR_SPACE.md is a ONE-LOOP robustness")
P("  statement.  It absorbs any reweighting of the fermion representations, and")
P("  it absorbs any common rescaling; it does not absorb a relative shift")
P("  between the fermion and gauge sectors larger than 1/27.  Loop order is the")
P("  one candidate CROSS_PAPER.md section 3 leaves open, and it is exactly the")
P("  direction the wedge is thinnest in.")
P("")
P("  The bounded question that would close it is already stated in the paper:")
P("  whether the two-loop fermionic potential depends on the representation only")
P("  through the multiset of its charges.  TWOLOOP_WEIGHT.md answers that in the")
P("  negative, which is why the ratio cannot be assumed to be 1.")
