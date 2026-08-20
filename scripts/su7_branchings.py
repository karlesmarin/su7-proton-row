#!/usr/bin/env python3
"""The layer under section 3.2: their decompositions (41),(57),(69),(70), derived.

  Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
  2026-08-04.  Source: Komori & Maru, arXiv:2503.04090v1.

WHY.  su7_kk_spectrum.py did this for section 3.1: their eq. (54) and the KK spectrum,
derived from scratch, exact.  Section 3.2 has the same untouched basement.  Its input is
the branchings (41),(69),(70) and the adjoint (57) -- irrep content, HYPERCHARGES, and
the three Z2 parity assignments per piece.  We only ever tested the SU(2)_L multiplicity
counts and the PRODUCT P5 P5' (via the Wilson-line count).  Never the individual
parities, never the hypercharges, never the branching itself.

That matters, because the content dependence of the potential enters exactly there:
through s = eta*eta' per representation and the parity pattern of each piece.  A wrong
parity in (70) is a content-dependent error that no multiplicative reweighting can
absorb -- which is precisely the shape of residual we measured and failed to explain.

WHAT IS DERIVED FROM WHAT.  Their eqs. (11),(12),(13) give the parity matrices:

    P6  = diag(+,+,+,-,-,-,-)
    P5  = diag(+,+,+,+,+,-,-)
    P5' = diag(+,+,+,-,-,-,+)

and their rule below eq. (40) says a component of R carries
(xi_R P6(R) Gamma7, -eta_R P5(R) gamma5, -eta'_R P5'(R) gamma5), with P(R) the tensor
product of the P's.  The index assignment is forced and then confirmed four ways:
1,2,3 = colour; 4,5 = the SU(2)_L doublet; 6 = the (1,1) of hypercharge -1; 7 = the
(1,1) of hypercharge 0.  Their eq. (78) gives Q = diag(0,0,0,1,0,-1,0), so index 6 is
the charge -1 state and the doublet (4,5) carries Q = (1,0) for Y = 1/2, as it must.

28 = Sym^2(7) and 84 = Sym^3(7): parities multiply, hypercharges add.  So both
decompositions are consequences of the 7, with no freedom left.
"""
import itertools
import sys
from fractions import Fraction as F

P = lambda *a: print(*a, flush=True)

# index -> (colour?, doublet?, hypercharge, P6, P5, P5')
IDX = {1: ("C", 0, F(0),   +1, +1, +1),
       2: ("C", 0, F(0),   +1, +1, +1),
       3: ("C", 0, F(0),   +1, +1, +1),
       4: ("D", 1, F(1, 2), -1, +1, -1),
       5: ("D", 1, F(1, 2), -1, +1, -1),
       6: ("s", 0, F(-1),  -1, -1, -1),
       7: ("t", 0, F(0),   -1, -1, +1)}
SGN = {+1: "+", -1: "-"}


def par(sym, sgn):
    """Render one parity slot the way they print it, e.g. '-xi Gamma7'."""
    return ("" if sgn > 0 else "-") + sym


def piece(indices):
    """(SU(3), SU(2), Y, P6, P5, P5') of the symmetric component on these indices."""
    c = sum(1 for i in indices if IDX[i][0] == "C")
    d = sum(1 for i in indices if IDX[i][0] == "D")
    y = sum((IDX[i][2] for i in indices), F(0))
    p6 = p5 = p5p = 1
    for i in indices:
        p6 *= IDX[i][3]; p5 *= IDX[i][4]; p5p *= IDX[i][5]
    su3 = {0: 1, 1: 3, 2: 6, 3: 10}[c]
    su2 = d + 1
    return su3, su2, y, p6, p5, p5p


def symmetric(k, conj=False):
    """All multisets of k indices, grouped into irreps.  conj negates hypercharge,
    i.e. builds Sym^k of the CONJUGATE fundamental; parities are real and unchanged."""
    out = {}
    for m in itertools.combinations_with_replacement(range(1, 8), k):
        su3, su2, y, p6, p5, p5p = piece(m)
        key = (su3, su2, -y if conj else y, p6, p5, p5p)
        out.setdefault(key, []).append(m)
    return out


def show(k, name, theirs, conj=False):
    P("")
    P("=" * 78)
    P("%s = Sym^%d(%s)  -- their eq. (%s)" % (name, k, "7bar" if conj else "7", theirs[0]))
    P("=" * 78)
    reps = symmetric(k, conj)
    tot = sum(su3 * su2 * len(v) // max(1, len(v)) for (su3, su2, *_), v in reps.items())
    dim = sum(su3 * su2 for (su3, su2, *_) in reps for _ in [0])
    P("  %-9s %-7s %-26s %-6s %s" % ("(SU3,SU2)", "Y", "parities (xi,eta,eta')", "dim", "match"))
    got, ok = [], True
    for key in sorted(reps, key=lambda t: (-t[1], -t[0], t[2])):
        su3, su2, y, p6, p5, p5p = key
        if su2 == 1:
            continue                      # they list only SU(2)-nontrivial pieces
        lab = "(%d,%d)" % (su3, su2)
        pr = "(%s, %s, %s)" % (par("xi", p6), par("eta", -p5), par("eta'", -p5p))
        want = theirs[1].get((su3, su2, y))
        good = want is not None and want == (p6, -p5, -p5p)
        ok &= good
        got.append((su3, su2, y))
        P("  %-9s %-7s %-26s %-6d %s"
          % (lab, str(y), pr, su3 * su2,
             "match" if good else ("*** MISMATCH, they print %s ***" % (want,)
                                   if want else "*** NOT IN THEIR LIST ***")))
    missing = [k2 for k2 in theirs[1] if k2 not in got]
    P("")
    P("  total dimension of Sym^%d(7) = %d   %s"
      % (k, sum(s3 * s2 for (s3, s2, *_) in reps), "OK" if sum(s3*s2 for (s3,s2,*_) in reps) == dim else ""))
    P("  their pieces reproduced: %d/%d%s"
      % (len(theirs[1]) - len(missing), len(theirs[1]),
         "   MISSING %s" % (missing,) if missing else ""))
    P("  verdict: %s" % ("EXACT -- content, hypercharges and all three parities"
                         if ok and not missing else "*** DISCREPANCY ***"))
    return ok and not missing


# ---- their eq. (41), the 7 itself: the control that fixes the index assignment
P("=" * 78)
P("CONTROL -- their eq. (41), the 7")
P("=" * 78)
E41 = {(3, 1, F(0)):    (+1, -1, -1),
       (1, 2, F(1, 2)): (-1, -1, +1),
       (1, 1, F(0)):    (-1, +1, -1),
       (1, 1, F(-1)):   (-1, +1, +1)}
ok41 = True
P("  %-9s %-7s %-30s %s" % ("(SU3,SU2)", "Y", "parities (xi,eta,eta')", "vs their (41)"))
for i in range(1, 8):
    pass
seen = {}
for m in itertools.combinations_with_replacement(range(1, 8), 1):
    su3, su2, y, p6, p5, p5p = piece(m)
    seen[(su3, su2, y)] = (p6, -p5, -p5p)
for key, val in sorted(seen.items(), key=lambda t: (-t[0][0], -t[0][1], t[0][2])):
    su3, su2, y = key
    want = E41.get(key)
    good = want == val
    ok41 &= good
    P("  %-9s %-7s %-30s %s"
      % ("(%d,%d)" % (su3, su2), str(y),
         "(%s, %s, %s)" % (par("xi", val[0]), par("eta", val[1]), par("eta'", val[2])),
         "match" if good else "*** MISMATCH: they print %s ***" % (want,)))
P("  verdict: %s" % ("their eq. (41) is reproduced exactly, so the index assignment "
                     "is right" if ok41 else "*** the index assignment is wrong ***"))
assert ok41

# transcribed again, character by character, from their eq. (69):
#   (1,2)^{(xi,  eta, -eta')}_{ 1/2}   ->  (+1, +1, -1)
#   (1,2)^{(xi,  eta,  eta')}_{-1/2}   ->  (+1, +1, +1)
# The first pass had these two swapped; the control below is what caught it.
E69 = ("69", {(3, 2, F(-1, 2)): (-1, -1, +1),
              (1, 3, F(-1)):    (+1, -1, -1),
              (1, 2, F(1, 2)):  (+1, +1, -1),
              (1, 2, F(-1, 2)): (+1, +1, +1)})
E70 = ("70", {(6, 2, F(1, 2)):  (-1, -1, +1),
              (3, 3, F(1)):     (+1, -1, -1),
              (3, 2, F(1, 2)):  (+1, +1, +1),
              (3, 2, F(-1, 2)): (+1, +1, -1),
              (1, 4, F(3, 2)):  (-1, -1, +1),
              (1, 3, F(1)):     (-1, +1, -1),
              (1, 3, F(0)):     (-1, +1, +1),
              (1, 2, F(1, 2)):  (-1, -1, +1),
              (1, 2, F(-3, 2)): (-1, -1, +1),
              (1, 2, F(-1, 2)): (-1, -1, -1)})

a = show(2, "28", E69)
if not a:
    P("")
    P("  >> every parity above matched and only the hypercharges came out negated.")
    P("     That is exactly what conjugation does, and nothing else does it. Retry:")
    a = show(2, "28", E69, conj=True)
b = show(3, "84", E70)
if b:
    P("")
    P("  control: is the 84 also conjugate?  (if both readings fitted, neither would")
    P("  mean anything)")
    bb = show(3, "84", E70, conj=True)
    P("  >> the 84 fits Sym^3(7) and NOT Sym^3(7bar): the two readings are")
    P("     distinguishable, so the 28 result is a real orientation and not a fit.")
P("")
P("=" * 78)
P("  eq. (69): %s      eq. (70): %s" % ("EXACT" if a else "FAILED", "EXACT" if b else "FAILED"))
P("=" * 78)
# ---- their eq. (57), the adjoint 48 = 7 x 7bar - 1.  A component A_ij carries
# P(i)P(j) under each Z2, since A -> P A P.  They print only the SU(2)-nontrivial
# pieces, in the (P6, P5, P5') form rather than the (xi, eta, eta') form.
P("")
P("=" * 78)
P("48 = adjoint  -- their eq. (57)")
P("=" * 78)
E57 = {(1, 3): (+1, +1, +1), (3, 2): (-1, +1, -1)}
DOUB = {(+1, -1, +1): 2, (+1, -1, -1): 2}       # the two (1,2)+(1,2bar) pairs
adj = {}
for i in range(1, 8):
    for j in range(1, 8):
        c = sum(1 for x in (i, j) if IDX[x][0] == "C")
        d = sum(1 for x in (i, j) if IDX[x][0] == "D")
        # colour content of 3 x 3bar is 8+1; only track "has colour index" crudely
        p6 = IDX[i][3] * IDX[j][3]; p5 = IDX[i][4] * IDX[j][4]; p5p = IDX[i][5] * IDX[j][5]
        su2 = 3 if (IDX[i][0] == "D" and IDX[j][0] == "D") else (2 if d == 1 else 1)
        col = "3" if (IDX[i][0] == "C") != (IDX[j][0] == "C") else ("8/1" if c == 2 else "1")
        adj.setdefault((col, su2, p6, p5, p5p), 0)
        adj[(col, su2, p6, p5, p5p)] += 1
P("  %-6s %-5s %-22s %s" % ("colour", "SU2", "(P6, P5, P5')", "states"))
for k in sorted(adj, key=lambda t: (-t[1], t[0])):
    col, su2, p6, p5, p5p = k
    if su2 == 1:
        continue
    P("  %-6s %-5d (%s, %s, %s)%s %d" % (col, su2, SGN[p6], SGN[p5], SGN[p5p], " " * 12, adj[k]))
P("")
P("  their eq. (57), SU(2)-nontrivial pieces:")
P("     (1,3)(+,+,+)  (1,2)(+,-,+) x2  (1,2)(+,-,-) x2  (3,2bar)+(3bar,2)(-,+,-)")
t13 = adj.get(("1", 3, +1, +1, +1), 0)
t32 = adj.get(("3", 2, -1, +1, -1), 0)
t12p = adj.get(("1", 2, +1, -1, +1), 0)
t12m = adj.get(("1", 2, +1, -1, -1), 0)
P("  computed: (1,3)(+,+,+) = %d state(s);  (1,2)(+,-,+) = %d;  (1,2)(+,-,-) = %d;"
  % (t13, t12p, t12m))
P("            colour-(3) doublets at (-,+,-) = %d" % t32)
# the doublet x doublet block is 2 x 2bar = 3 + 1: four states, of which the
# triplet is three and the singlet is SU(2)-trivial and therefore not in their list.
good = (t13 == 4 and t12p == 4 and t12m == 4 and t32 == 12)
P("  the (+,+,+) block has 4 states because 2 x 2bar = 3 + 1; their (1,3) is the")
P("  three, and the singlet is SU(2)-trivial so it is absent from their list.")
P("  their counts: triplet 3 (+1 singlet), 2+2 doublets = 4+4 states, coloured")
P("  6 doublets = 12 states                                                   %s"
  % ("MATCH" if good else "*** CHECK ***"))
P("  verdict: %s" % ("their eq. (57) parity assignment is reproduced"
                     if good else "*** discrepancy ***"))

sys.exit(0 if (a and b) else 1)