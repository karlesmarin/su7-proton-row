#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

ONE AFFIRMATION, NOT TWO: the fate of the proton in Komori-Maru's SU(7) GGHU is a
FUNCTION of a single datum their paper never states -- what U(1)' charges may live
at a fixed point -- and this computes the function completely.

Two places in this programme found the same hole from opposite ends:

  su7_socratic.py (Q2)   whether X_Q = -1/6 (the one point where the U(1)' stops
                         forbidding the proton) is REACHABLE depends on whether
                         brane fields must carry charges of the bulk lattice.
  su7_family_u1.py (9)   the selection rule that forbids the proton to all orders
                         needs A_j / q_phi not in Z, and their paper does not
                         state q_phi, the charge of the U(1)'-breaking scalar.

Both are the SAME unstated global choice.  Written as one statement it becomes a
theorem with a closed form rather than a caveat repeated twice.

  STRICT      fixed-point fields must be liftable to SU(7) reps
              => every charge in the BULK lattice (1/2)Z, and q_phi is one of the
                 finitely many values their own representations supply.
  PERMISSIVE  arbitrary reps of the group unbroken at the fixed point
              => charges free.  Standard orbifold-GUT practice, and the reading
                 their own brane quarks already use.

INPUTS, all from the paper:
  eq. (27)/(31)  Y = (0,0,0,1/2,1/2,-1,0) per index
  eq. (79)       U(1)' = (1/2) diag(0,0,0,1,-1,1,-1);  extra := U(1)' - T3L
  eq. (46)       the lepton Yukawa is <A5> connecting index 5 <-> index 7
  eq. (47)       brane quark Yukawas => extra(u^c) = -(a+1/2), extra(d^c) = -(a-1/2)
  Table 1        the reps in play: 7, 21, 28, 35, 48, 84
  above eq. (80) U(1)' "has to be spontaneously broken by introducing a U(1)'
                 charged 4D scalar field on the fixed point"
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr, combinations
from math import gcd as igcd

P = lambda *a: print(*a, flush=True)
H = F(1, 2)

YIND = (F(0), F(0), F(0), H, H, F(-1), F(0))
COLOUR, WEAK, EXTRA = (0, 1, 2), (3, 4), (5, 6)
t3l = lambda n: F(n[3] - n[4], 2)
ex = lambda n: F(n[5] - n[6], 2)
col = lambda n: n[0] + n[1] + n[2]
Ydot = lambda n: sum(c * k for c, k in zip(YIND, n))


def fgcd(a, b):
    """gcd of two non-negative rationals: the generator of <a, b> as a subgroup of Q."""
    a, b = abs(F(a)), abs(F(b))
    if a == 0:
        return b
    if b == 0:
        return a
    num = igcd(a.numerator * b.denominator, b.numerator * a.denominator)
    return F(num, a.denominator * b.denominator)


def lattice(charges):
    g = F(0)
    for c in charges:
        g = fgcd(g, c)
    return g


def counts(idx, r=7):
    n = [0] * r
    for i in idx:
        n[i] += 1
    return tuple(n)


def components(rep):
    """Index content of every component of the reps Komori-Maru use."""
    if rep == "7":
        return [counts((i,)) for i in range(7)]
    if rep == "21":                                   # Lambda^2
        return [counts(c) for c in combinations(range(7), 2)]
    if rep == "28":                                   # Sym^2
        return [counts(c) for c in cwr(range(7), 2)]
    if rep == "35":                                   # Lambda^3
        return [counts(c) for c in combinations(range(7), 3)]
    if rep == "84":                                   # Sym^3
        return [counts(c) for c in cwr(range(7), 3)]
    if rep == "48":                                   # adjoint: e_ij (i!=j) + Cartan
        out = []
        for i in range(7):
            for j in range(7):
                if i != j:
                    n = [0] * 7
                    n[i], n[j] = 1, -1
                    out.append(tuple(n))
        out += [tuple([0] * 7)] * 6
        return out
    raise ValueError(rep)


REPS = ["7", "21", "28", "35", "48", "84"]
DIM = {"7": 7, "21": 21, "28": 28, "35": 35, "48": 48, "84": 84}

P("=" * 78)
P("STEP 0 -- CONTROL: the reps are built right before anything is concluded")
P("=" * 78)
P("  %-6s %-8s %s" % ("rep", "built", "matches its name?"))
for r in REPS:
    m = len(components(r))
    P("  %-6s %-8d %s" % (r, m, "OK" if m == DIM[r] else "*** MISMATCH ***"))
    assert m == DIM[r]

P("")
P("=" * 78)
P("STEP 1 -- the BULK lattice, derived and not assumed")
P("=" * 78)
P("  extra(n) = (n_6 - n_7)/2 for every component of every bulk multiplet, so")
P("  every bulk charge is a half-integer.  Check it, and check that 1/2 is")
P("  actually attained (a lattice claim needs both).")
P("")
allbulk = []
for r in REPS:
    ch = [ex(n) for n in components(r)]
    allbulk += ch
    P("  %-6s charges in (1/2)Z: %-6s   min %-6s  max %-6s  lattice %s"
      % (r, "yes" if all(2 * c == int(2 * c) for c in ch) else "NO",
         min(ch), max(ch), lattice(ch) or "0 (all neutral)"))
assert all((2 * c).denominator == 1 for c in allbulk)
BULK = lattice(allbulk)
P("")
P("  >> the bulk lattice is exactly (%s)Z.  This is what su7_socratic.py's Q2" % BULK)
P("     called the STRICT reading, and it is the only thing the bulk shows.")
assert BULK == H

P("")
P("=" * 78)
P("STEP 2 -- what q_phi their OWN representations can supply")
P("=" * 78)
P("  Their sentence above eq. (80) needs a 4D scalar that is U(1)'-charged and")
P("  breaks nothing else: a total SM singlet (no colour index, no weak index)")
P("  with Y = 0.  Enumerate those over the reps of their Table 1.")
P("")
P("  %-6s %-26s %-8s %-8s %s" % ("rep", "component", "Y", "q_phi", "usable?"))
supply = {}
for r in REPS:
    for n in components(r):
        if col(n) != 0 or n[3] != 0 or n[4] != 0:
            continue
        if any(k < 0 for k in n):                     # adjoint lowering pieces
            continue
        y, q = Ydot(n), ex(n)
        if y != 0:
            continue
        lab = "(" + ",".join(str(i + 1) * k for i, k in enumerate(n) if k) + ")"
        lab = lab if lab != "()" else "(Cartan)"
        P("  %-6s %-26s %-8s %-8s %s"
          % (r, lab, y, q, "YES" if q != 0 else "no -- uncharged"))
        if q != 0:
            supply.setdefault(abs(q), []).append(r)
P("")
STRICT_QPHI = sorted(supply)
P("  >> under the STRICT reading  |q_phi| in %s" % [str(v) for v in STRICT_QPHI])
P("     supplied by             %s"
  % {str(k): sorted(set(v)) for k, v in sorted(supply.items())})
P("     and the 35 supplies NONE: a pure index-7 component needs a repeated")
P("     index, which an antisymmetric tensor has not.  Same one-line reason")
P("     that kept rung 1 out of the 35 in SU7_FAMILY_U1.md.")
P("")
P("  >> under the PERMISSIVE reading q_phi is a free rational.")

# ------------------------------------------------------------------ the content
def gen(a, l, nus=()):
    """One generation, all-left-handed.  (T3C, T3L, Y, X).  su7_family_u1.py."""
    out = []
    for tc in (H, -H, F(0)):
        for tl in (H, -H):
            out.append((tc, tl, F(1, 6), a + tl))
        out.append((-tc, F(0), F(-2, 3), -(a + H)))
        out.append((-tc, F(0), F(1, 3), -(a - H)))
    for tl in (H, -H):
        out.append((F(0), tl, -H, l + tl))
    out.append((F(0), F(0), F(1), H - l))
    for v in nus:
        out.append((F(0), F(0), F(0), v))
    return out


CH = {
    "SU(3)^2 X": lambda s: sum(x * tc ** 2 for tc, tl, y, x in s),
    "SU(2)^2 X": lambda s: sum(x * tl ** 2 for tc, tl, y, x in s),
    "X grav^2":  lambda s: sum(x for tc, tl, y, x in s),
    "X^3":       lambda s: sum(x ** 3 for tc, tl, y, x in s),
    "X^2 Y":     lambda s: sum(x ** 2 * y for tc, tl, y, x in s),
    "X Y^2":     lambda s: sum(x * y ** 2 for tc, tl, y, x in s),
}
Achg = lambda a, l: 3 * a + l

P("")
P("=" * 78)
P("STEP 3 -- N = 1, THEIR MODEL AS WRITTEN: the strict reading is not available")
P("=" * 78)
P("  SU7_ANOMALY_CHANNELS.md: three channels each force X_Q = -1/6 alone.")
XQ1 = F(-1, 6)
s1 = gen(XQ1, H, (F(-1),))
P("  content = brane quarks at X_Q + the 21's leptons + one nu_R of charge -1")
for k, f in CH.items():
    P("    %-12s = %s" % (k, f(s1)))
assert all(f(s1) == 0 for f in CH.values())
P("  CONTROL: every channel cancels at X_Q = -1/6                        : OK")
P("")
P("  Is that value on the bulk lattice (%s)Z?  %s / %s = %s"
  % (BULK, XQ1, BULK, XQ1 / BULK))
on = (XQ1 / BULK).denominator == 1
P("    -> %s" % ("YES" if on else "NO -- it is not a multiple of 1/2"))
assert not on
P("")
P("  >> THE STRICT READING AND ANOMALY CANCELLATION ARE INCOMPATIBLE for their")
P("     model as written.  Absent Green-Schwarz, their own consistency FORCES")
P("     the permissive reading -- and the permissive reading is exactly what")
P("     puts X_Q = -1/6 within reach.  su7_socratic.py's Q2 is settled here by")
P("     THEIR arithmetic, not by our convention, and on the permissive side,")
P("     which is su7_channels_socratic.py's A5 reached independently.")
P("")
P("  And at A = 3 X_Q + 1/2 = %s no q_phi can help: 0/q is an integer for every"
  % Achg(XQ1, H))
P("  q, so the operators are dressable at every breaking scale.")
assert Achg(XQ1, H) == 0

P("")
P("=" * 78)
P("STEP 4 -- N = 3: recompute the surviving assignments (independent of STEP 8")
P("          of su7_family_u1.py -- the count is the control)")
P("=" * 78)
RUNGS = [H, F(0), -H, F(-1)]
NU = [F(m, 2) for m in range(-4, 5) if m]
sols = []
seen = set()
for ls in cwr(RUNGS, 3):
    a = -sum(ls) / 9
    As = [Achg(a, l) for l in ls]
    if not all(v != 0 for v in As):
        continue
    for nn in range(0, 4):
        for nus in cwr(NU, nn):
            s = []
            for j, l in enumerate(ls):
                s += gen(a, l, nus if j == 0 else ())
            if all(f(s) == 0 for f in CH.values()):
                if (ls, a) not in seen:
                    seen.add((ls, a))
                    sols.append((ls, a, As, nus))
P("  distinct (rung set, X_Q) solutions : %d   (SU7_FAMILY_U1.md says 14)" % len(sols))
assert len(sols) == 14
P("  CONTROL PASSES")

P("")
P("=" * 78)
P("STEP 5 -- which solutions are STRICT-compatible at all?")
P("=" * 78)
P("  A solution lives on the bulk lattice iff every charge it introduces does.")
P("  The rungs and the nu_R charges are half-integers by construction, so the")
P("  question is entirely X_Q: the quark charges are X_Q, X_Q +- 1/2.")
P("")
P("  %-24s %-8s %-28s %-10s %s"
  % ("(l1,l2,l3)", "X_Q", "(A1,A2,A3)", "lattice", "strict?"))
strict_ok = []
for ls, a, As, nus in sols:
    lat = lattice([BULK] + [a])
    ok = (a / BULK).denominator == 1
    if ok:
        strict_ok.append((ls, a, As, nus))
    P("  %-24s %-8s %-28s %-10s %s"
      % (str(tuple(str(x) for x in ls)), a, str(tuple(str(v) for v in As)),
         "(%s)Z" % lat, "YES" if ok else "no"))
P("")
P("  >> %d of %d survive the strict reading, and it is %s"
  % (len(strict_ok), len(sols),
     ", ".join(str(tuple(str(x) for x in ls)) for ls, a, As, nus in strict_ok)))

# ------------------------------------------------------------------ selection rule
def failing_qphi(As):
    """The COMPLETE set of positive q_phi at which the selection rule fails.

    An operator of charge A is dressable by <phi>^n or <phi*>^n iff A/q_phi is an
    integer.  So protection fails at q_phi iff q_phi divides some A_j, i.e. iff
    q_phi = |A_j|/n for a positive integer n.  That set is countable and closed
    form -- no scan is involved.
    """
    out = {}
    for j, A in enumerate(As):
        if A == 0:
            return None                  # every q_phi fails
        out.setdefault(abs(A), []).append(j)
    return out


def fails_at(As, q):
    return [j for j, A in enumerate(As) if (A / q).denominator == 1]


P("")
P("=" * 78)
P("STEP 6 -- THE SELECTION RULE, in closed form instead of a scan")
P("=" * 78)
P("  <phi> of charge q_phi breaks U(1)' to a discrete subgroup; an operator of")
P("  charge A survives dressing iff A/q_phi is an INTEGER.  So")
P("")
P("     protection fails at q_phi  <=>  q_phi = |A_j| / n  for some j, n >= 1.")
P("")
P("  That is the whole answer, and it is finite to state.  For each solution:")
P("")
P("  %-24s %-26s %s" % ("(l1,l2,l3)", "(A1,A2,A3)", "q_phi at which it FAILS"))
closed = {}
for ls, a, As, nus in sols:
    d = sorted({abs(v) for v in As})
    # the union of {A/n} over j, written as the coarsest generators
    desc = " U ".join("%s/n" % v for v in d)
    closed[(ls, a)] = d
    P("  %-24s %-26s %s"
      % (str(tuple(str(x) for x in ls)), str(tuple(str(v) for v in As)), desc))

P("")
P("  CONTROL -- the closed form against a brute scan.  Every rational p/q with")
P("  q <= 36, p/q in (0, 2], for every solution: does the scan ever disagree?")
bad = 0
tested = 0
for ls, a, As, nus in sols:
    d = {abs(v) for v in As}
    for den in range(1, 37):
        for num in range(1, 2 * den + 1):
            q = F(num, den)
            scan = len(fails_at(As, q)) > 0
            closedf = any((A / q).denominator == 1 for A in d)
            tested += 1
            if scan != closedf:
                bad += 1
P("  %d values tested, %d disagreements" % (tested, bad))
assert bad == 0
P("  CONTROL PASSES")

P("")
P("=" * 78)
P("STEP 7 -- and now the two readings give two different, decidable answers")
P("=" * 78)
MIN = [(ls, a, As, nus) for ls, a, As, nus in sols if ls == (H, H, F(0))][0]
P("  (a) THE MINIMAL ASSIGNMENT, l = (1/2, 1/2, 0), X_Q = %s -- permissive only." % MIN[1])
P("      A_j = %s" % str(tuple(str(v) for v in MIN[2])))
P("      fails exactly at q_phi in {1/6/n} U {1/3/n} = {1/(3n) : n >= 1}:")
u = sorted({F(1, 3) / n for n in range(1, 7)} | {F(1, 6) / n for n in range(1, 7)})
P("        %s ..." % ", ".join(str(v) for v in u[-8:][::-1]))
chk = all(any((A / q).denominator == 1 for A in MIN[2]) for q in u)
P("      CONTROL: every member of that set really does fail : %s" % ("OK" if chk else "NO"))
assert chk
P("      largest failing q_phi = %s" % max(u))
P("")
P("      >> THE BOUND, and it is exact: the failing set has a MAXIMUM, so")
P("         every q_phi > 1/3 forbids all four operators for all three")
P("         generations to all orders.  Protection is not a fine-tuning; it is")
P("         a half-line.")
assert max(abs(v) for v in MIN[2]) == F(1, 3)
P("")
P("      Against what their reps supply, %s:" % [str(v) for v in STRICT_QPHI])
for q in STRICT_QPHI:
    f = fails_at(MIN[2], q)
    P("        q_phi = %-5s  -> %d of 3 generations dressable  %s"
      % (q, len(f), "PROTECTED" if not f else "proton allowed via %s" % f))
assert all(not fails_at(MIN[2], q) for q in STRICT_QPHI)
P("      >> EVERY q_phi their own representations can supply PROTECTS.  Failure")
P("         needs q_phi <= 1/3, finer than anything in their spectrum.")
P("")
S = strict_ok[0]
P("  (b) THE ONE STRICT-COMPATIBLE SOLUTION, l = %s, X_Q = %s."
  % (str(tuple(str(x) for x in S[0])), S[1]))
P("      A_j = %s -- all half-integers, so the failing set is coarse:"
  % str(tuple(str(v) for v in S[2])))
for q in STRICT_QPHI:
    f = fails_at(S[2], q)
    P("        q_phi = %-5s  (%s)  -> %s"
      % (q, "+".join(sorted(set(supply[q]))),
         "PROTECTED" if not f else "proton allowed via generation(s) %s" % f))
P("      the bound here is q_phi > 1, since max |A_j| = %s"
  % max(abs(v) for v in S[2]))
assert max(abs(v) for v in S[2]) == F(1)
P("      >> here the two smallest breaking scalars FAIL, and protection needs")
P("         q_phi = 3/2 -- available only in the 84, the same multiplet that")
P("         SU7_FAMILY_U1.md needed to host rung 1 and that every row of their")
P("         Table 1 already contains.")
assert fails_at(S[2], F(1, 2)) and fails_at(S[2], F(1)) and not fails_at(S[2], F(3, 2))

P("")
P("=" * 78)
P("STEP 8 -- CONTROL that must fire: the family-universal case is unprotectable")
P("=" * 78)
uni = [Achg(F(-1, 6), H)] * 3
P("  l = (1/2,1/2,1/2) -> X_Q = -1/6, A_j = %s" % str(tuple(str(v) for v in uni)))
P("  failing_qphi -> %s" % ("EVERY q_phi" if failing_qphi(uni) is None else "?"))
assert failing_qphi(uni) is None
P("  >> A_j = 0 and 0/q is an integer for every q: no choice of breaking scalar")
P("     protects anything.  The control fires, as it must.")

P("")
P("=" * 78)
P("THE AFFIRMATION -- one statement, two readings, both decided")
P("=" * 78)
P("""
  Komori-Maru's paper does not state what U(1)' charges may live at a fixed
  point.  That single omission decides the proton, and it decides it twice in
  opposite directions:

    * Their model AS WRITTEN (one generation) is anomaly-consistent ONLY under
      the permissive reading -- X_Q = -1/6 is forced by three channels and is
      not on the bulk lattice (1/2)Z.  And under that same permissive reading
      the protection is gone: A = 0, dressable at every q_phi.

    * In the three-generation extension the choice becomes a live variable.
      13 of the 14 anomaly-free assignments exist only permissively; exactly one,
      l = (1/2,1/2,-1) with X_Q = 0, is strict-compatible.  For it the two
      smallest breaking scalars their reps supply (q_phi = 1/2 from the 7, 1 from
      the 28) leave the proton dressable, and only the 84's q_phi = 3/2 protects.
      For the minimal permissive assignment l = (1/2,1/2,0) the failing set is
      exactly {1/(3n)}, and no representation of theirs reaches it.

  So this is not two caveats.  It is one function -- (reading, q_phi) -> verdict --
  computed here in closed form, whose argument their paper leaves unset.  It is a
  well-posed question FOR THE AUTHORS, and it is the same question in both places
  su7_socratic.py and su7_family_u1.py met it.

  AND THE FUNCTION IS KINDER THAN A SCAN SUGGESTS.  The failing set has a
  maximum in every case, so the answer is a half-line and not a fine-tuning:

      minimal permissive assignment   protected for every  q_phi > 1/3
      the strict-compatible one       protected for every  q_phi > 1

  What their paper leaves open is therefore not "is the proton stable?" but
  "is the U(1)'-breaking scalar heavier-charged than 1/3?", which is a single
  number they could state in one line.
""")
P("NOT CLAIMED:")
P("  - that only one scalar breaks U(1)'.  With several, the residual group is")
P("    fixed by gcd(q_1, q_2, ...), which can only make protection HARDER; the")
P("    single-scalar case computed here is therefore the optimistic bound.")
P("  - any VEV, scale, coefficient or lifetime.  None is computed.")
P("  - novelty of the mechanism.  A residual discrete gauge symmetry from a")
P("    broken extra U(1) forbidding the proton operators is Krauss-Wilczek /")
P("    JHEP 07 (2008) 065 / hep-ph/0012092.  Ours is the specialisation.")
P("  - that the anomalies must cancel on localised matter: Green-Schwarz remains")
P("    available, and it is the escape su7_channels_socratic.py's A4 bounds.")
P("")
P("DONE")
