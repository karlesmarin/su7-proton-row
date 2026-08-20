#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

CAN KOMORI-MARU'S MODEL CARRY A FAMILY-DEPENDENT U(1)?

SU7_ANOMALY_CHANNELS.md closed the family-universal case: anomaly cancellation
forces U(1)' = T3L + Y - (B-L), and every dimension-6 |dB|=1 operator is neutral
under it.  The literature's own cure for exactly that disease (arXiv:1001.0768)
is to make the U(1) FAMILY-DEPENDENT -- B - x_i L with sum x_i = 3 and x_i != 1.

In their model the leptons come from a BULK 21, and U(1)' is a fixed generator of
SU(7) (their eq. 79), so a lepton's charge is NOT free data: it is read off the
SU(7) weight.  Family dependence therefore exists only if the three generations
can sit at INEQUIVALENT weights.  That is what this computes.

INPUTS, from the paper:
  eq. (27)/(31)  T3L, Qem  =>  Y = (0,0,0,1/2,1/2,-1,0) per index
  eq. (79)       U(1)' = (1/2) diag(0,0,0,1,-1,1,-1);  extra := U(1)' - T3L
  eq. (46)       the lepton Yukawa is <A5> connecting index 5 <-> index 7
  eq. (47)       the brane quark Yukawas  =>  extra(u^c) = -(a+1/2),
                 extra(d^c) = -(a-1/2),  a := extra(Q)
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement as cwr, product

P = lambda *a: print(*a, flush=True)
H = F(1, 2)

YIND = (F(0), F(0), F(0), H, H, F(-1), F(0))
UP1 = (F(0), F(0), F(0), H, -H, H, -H)
dot = lambda v, n: sum(c * k for c, k in zip(v, n))
t3l = lambda n: F(n[3] - n[4], 2)
ex = lambda n: F(n[5] - n[6], 2)
col = lambda n: n[0] + n[1] + n[2]

P("=" * 78)
P("STEP 0 -- the ladder: what U(1)' charge can an SM lepton doublet have?")
P("=" * 78)
P("  A component of any SU(7) tensor is fixed by its index content n.  Demand")
P("  it be a colour singlet, an SU(2)_L doublet member, and Y = -1/2:")
P("     colour indices 1-3 : none;   exactly one index in {4,5};   Y = -1/2.")
P("  Y per index is (0,0,0,1/2,1/2,-1,0), so the non-weak part must carry")
P("  Y = -1, i.e. exactly one index 6 and any number k of index 7.")
P("")
P("  %-6s %-24s %-10s %-10s %s" % ("boxes", "content of L", "T3L", "Y", "extra(L)"))
ladder = {}
for k in range(0, 5):
    n = [0] * 7
    n[4] = 1        # the lower component, T3L = -1/2
    n[5] = 1
    n[6] = k
    n = tuple(n)
    assert col(n) == 0 and dot(YIND, n) == -H and t3l(n) == -H
    ladder[k] = ex(n)
    P("  %-6d %-24s %-10s %-10s %s"
      % (2 + k, "(5,6" + ",7" * k + ")", t3l(n), dot(YIND, n), ex(n)))
P("")
P("  >> extra(L) = (1 - k)/2 :  a LADDER, spaced by 1/2, indexed by how many")
P("     index-7 boxes the multiplet carries.  It is NOT a free parameter and it")
P("     is NOT unique -- which is exactly what family dependence needs.")
P("")
P("  CONTROL  k = 0 must reproduce their own assignment, L = 21^{56}")
assert ladder[0] == H
P("           extra(L) = 1/2                                          : OK")

P("")
P("=" * 78)
P("STEP 1 -- the Yukawa fixes e^c, and it closes for EVERY rung")
P("=" * 78)
P("  Their eq. (46): <A5> connects index 5 <-> index 7.  So e_R is L with the")
P("  index 5 replaced by a 7, inside the SAME multiplet.")
P("")
P("  %-6s %-18s %-12s %-18s %-12s %s"
  % ("boxes", "L", "extra(L)", "e_R", "extra(e_R)", "extra(L)+extra(e^c)"))
for k in range(0, 4):
    nL = [0] * 7
    nL[4], nL[5], nL[6] = 1, 1, k
    nE = [0] * 7
    nE[5], nE[6] = 1, k + 1            # 5 -> 7
    nL, nE = tuple(nL), tuple(nE)
    assert col(nE) == 0 and t3l(nE) == 0 and dot(YIND, nE) == F(-1)
    lc, ec = ex(nL), -ex(nE)           # e^c is the conjugate of e_R
    P("  %-6d %-18s %-12s %-18s %-12s %s"
      % (2 + k, "(5,6" + ",7" * k + ")", ex(nL),
         "(6" + ",7" * (k + 1) + ")", ex(nE), lc + ec))
    assert lc + ec == H
P("")
P("  >> extra(L) + extra(e^c) = 1/2 on EVERY rung -- the Higgs's own extra")
P("     charge.  The charged-lepton Yukawa is automatic at every rung, because")
P("     it comes from the 6D gauge coupling inside one multiplet.")
P("     e_R is a singlet of the right Y at every rung: verified above.")

# ------------------------------------------------------------------ the content
def gen(a, l, nus=()):
    """One generation, all-left-handed.  (T3C, T3L, Y, X)."""
    out = []
    for tc in (H, -H, F(0)):
        for tl in (H, -H):
            out.append((tc, tl, F(1, 6), a + tl))
        out.append((-tc, F(0), F(-2, 3), -(a + H)))
        out.append((-tc, F(0), F(1, 3), -(a - H)))
    for tl in (H, -H):
        out.append((F(0), tl, -H, l + tl))
    out.append((F(0), F(0), F(1), H - l))       # e^c, from STEP 1
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
A = lambda a, l: 3 * a + l          # the proton-operator charge, su7_XQ.py

P("")
P("=" * 78)
P("STEP 2 -- THE IDENTITY that decides the whole question")
P("=" * 78)
P("  The proton-operator charge of the generation with lepton (a_j, l_j) is")
P("  A_j = 3 a_j + l_j  (su7_XQ.py; and Theorem 2 of halves_theorem.py makes the")
P("  e^c half exactly -A_j).  Compute the SU(2)_L anomaly of N generations:")
P("")
P("  %-42s %s" % ("content", "SU(2)^2 X  vs  (1/2) sum A_j"))
for trial in ([(F(0), H)], [(F(0), H), (F(1, 3), F(0))],
              [(F(-1, 9), H), (F(-1, 9), H), (F(-1, 9), F(0))]):
    s = [st for a, l in trial for st in gen(a, l)]
    lhs = CH["SU(2)^2 X"](s)
    rhs = sum(A(a, l) for a, l in trial) / 2
    P("  %-42s %-12s %-12s %s"
      % (str([(str(a), str(l)) for a, l in trial]), lhs, rhs,
         "OK" if lhs == rhs else "*** MISMATCH ***"))
    assert lhs == rhs
P("")
P("  >>  SU(2)^2 X  =  (1/2) sum_j A_j   exactly.")
P("")
P("  ANOMALY CANCELLATION  demands   sum_j A_j = 0.")
P("  PROTON PROTECTION     demands   A_j != 0  for every j.")
P("")
P("  For N = 1 those are the SAME equation, and that is the whole content of")
P("  SU7_ANOMALY_CHANNELS.md: A = 0 is forced, the proton is unprotected.")
P("  For N >= 2 they are NOT: a sum can vanish with no term vanishing.")
P("  >> THIS IS arXiv:1001.0768's 'sum x_i = 3 with x_i != 1', in their")
P("     variables.  The escape is real, and it is only open for N > 1.")

P("")
P("=" * 78)
P("STEP 3 -- but the quark charges are not free either: CKM")
P("=" * 78)
P("  A family-dependent extra charge on the brane quarks forbids the OFF-")
P("  DIAGONAL quark Yukawas (their eq. 47 is diagonal in flavour), so the CKM")
P("  matrix would have to come from <phi> insertions.  Two regimes:")
P("")
P("   (i)  a_j UNIVERSAL  -- exact CKM at renormalisable level.  Then")
P("        sum_j A_j = 3 N a + sum_j l_j = 0  =>  a = -(sum l)/(3N), and")
P("        A_j = l_j - lbar   with  lbar the MEAN of the l_j.")
P("        Protection  <=>  no l_j equals the mean.")
P("   (ii) a_j FAMILY-DEPENDENT -- more freedom, CKM suppressed by <phi>/M.")
P("")
P("  Regime (i), N = 3, over the rungs l in {1/2, 0, -1/2, -1}:")
P("")
P("  %-22s %-10s %-30s %s" % ("(l1,l2,l3)", "a", "(A1,A2,A3)", "all != 0?"))
RUNGS = [H, F(0), -H, F(-1)]
good_i = []
for ls in cwr(RUNGS, 3):
    a = -sum(ls) / 9
    As = [A(a, l) for l in ls]
    ok = all(v != 0 for v in As)
    if ok:
        good_i.append((ls, a, As))
    P("  %-22s %-10s %-30s %s"
      % (str(tuple(str(x) for x in ls)), a,
         str(tuple(str(v) for v in As)), "YES" if ok else "no"))
P("")
P("  >> %d of %d rung assignments protect every generation."
  % (len(good_i), len(list(cwr(RUNGS, 3)))))
P("     The ones that fail are exactly those where some l_j is the mean --")
P("     in particular ALL-EQUAL rungs, which is the family-universal case.")

P("")
P("=" * 78)
P("STEP 4 -- now impose EVERY anomaly channel, not just the SU(2) one")
P("=" * 78)
P("  The other channels are not linear in A_j, so surviving STEP 3 is necessary")
P("  and not sufficient.  Add right-handed neutrinos: the SM singlets with Y = 0")
P("  are the pure index-7 components, extra = -m/2, so their available charges")
P("  are (+-1/2, +-1, +-3/2, ...) -- the same ladder, signed by chirality.")
P("")
NU = [F(m, 2) for m in range(-4, 5) if m]
P("  searching: rungs (l1,l2,l3) x universal a x up to 3 neutrinos from %s"
  % [str(v) for v in NU])
P("")
sols = []
for ls, a, As in good_i:
    for nn in range(0, 4):
        for nus in cwr(NU, nn):
            s = []
            for j, l in enumerate(ls):
                s += gen(a, l, nus if j == 0 else ())
            if all(f(s) == 0 for f in CH.values()):
                sols.append((ls, a, As, nus))
P("  %-22s %-9s %-28s %-22s" % ("(l1,l2,l3)", "a", "(A1,A2,A3)", "neutrinos"))
seen = set()
for ls, a, As, nus in sols:
    key = (ls, a)
    if key in seen:
        continue
    seen.add(key)
    P("  %-22s %-9s %-28s %-22s"
      % (str(tuple(str(x) for x in ls)), a, str(tuple(str(v) for v in As)),
         str(tuple(str(v) for v in nus))))
if not sols:
    P("  NONE -- no assignment survives every channel.")
P("")
P("  total solutions found (rung set, a) : %d" % len(seen))

P("")
P("=" * 78)
P("STEP 5 -- CONTROL: the family-universal case must come back UNPROTECTED")
P("=" * 78)
P("  Set all three rungs to their own value l = 1/2 (three copies of the 21):")
ls = (H, H, H)
a = -sum(ls) / 9
P("     a = %s,  A_j = %s" % (a, [str(A(a, l)) for l in ls]))
assert all(A(a, l) == 0 for l in ls)
P("     >> every A_j = 0.  The proton is unprotected, which is exactly")
P("        SU7_ANOMALY_CHANNELS.md's result, recovered here as a special case.")
P("     CONTROL PASSES -- the new machinery reproduces the old verdict.")
P("  And a = -1/6, the value that note derived:  %s" % (a == F(-1, 6)))
assert a == F(-1, 6)

P("")
P("=" * 78)
P("STEP 6 -- can their construction actually put a generation on rung k >= 1?")
P("=" * 78)
P("  Rung k needs L at index content (4or5, 6, 7^k), i.e. a multiplet with")
P("  2 + k boxes.  Their reps: 7 (1 box), 21 and 28 (2), 84 (3), 48 (adjoint).")
P("  Which of them CONTAIN such a component at all?")
P("")


def contains(rep_boxes, sym, k):
    """Does a k-rung lepton doublet live in the (anti)symmetric rep?"""
    n = [0] * 7
    n[4], n[5], n[6] = 1, 1, k
    if sum(n) != rep_boxes:
        return False
    return max(n) <= 1 if sym == "anti" else True


P("  %-6s %-8s %-8s %s" % ("rep", "boxes", "type", "rungs it can host"))
for nm, bx, ty in (("7", 1, "fund"), ("21", 2, "anti"), ("28", 2, "sym"),
                   ("84", 3, "sym"), ("35", 3, "anti"), ("112", 4, "mixed")):
    rr = [k for k in range(0, 4) if contains(bx, ty, k)]
    P("  %-6s %-8d %-8s %s" % (nm, bx, ty, rr if rr else "none"))
P("")
P("  >> rung 0 needs 2 boxes  -> the 21 (antisymmetric: indices 5,6 distinct)")
P("     rung 1 needs 3 boxes  -> the 84 (symmetric, indices 5,6,7 distinct, so")
P("                             the ANTIsymmetric 35 hosts it too)")
P("     rung 2 needs 4 boxes  -> a rep they do not introduce")
P("")
P("  >> so the assignment closest to their own model is  l = (1/2, 1/2, 0):")
P("     TWO generations in 21s and ONE in an 84 or a 35.  The 84 is already in")
P("     every row of their Table 1, introduced for the potential.")
ls = (H, H, F(0))
a = -sum(ls) / 9
P("     a = %s, A_j = %s  -- all non-zero"
  % (a, [str(A(a, l)) for l in ls]))
assert all(A(a, l) != 0 for l in ls)
P("     proton suppression per generation, |A_j| in units of <phi>/M:")
for j, l in enumerate(ls, 1):
    P("        generation %d : |A| = %s" % (j, abs(A(a, l))))

P("")
P("=" * 78)
P("STEP 7 -- does the orbifold KEEP a rung-1 lepton generation?  (the real gap)")
P("=" * 78)
P("  Their eqs. (37)-(40), transcribed in su7_signed_census.py, give:")
P("     zero mode  <=>  p5(n) p5'(n) = eta eta' ,  4D chirality  g = -eta p5(n)")
P("  Demand of a candidate rep: L = (5,6,7^k) is a zero mode and LEFT-handed,")
P("  and e_R = (6,7^{k+1}) is a zero mode and RIGHT-handed.  Nothing else fixed.")
P("")
P5V = (1, 1, 1, 1, 1, -1, -1)
P5P = (1, 1, 1, -1, -1, -1, 1)


def par(M, n):
    s = 1
    for i, k in enumerate(n):
        if k % 2:
            s *= M[i]
    return s


def mode(n, eta, etap):
    if par(P5V, n) * par(P5P, n) != eta * etap:
        return None
    return -eta * par(P5V, n)            # -1 = left, +1 = right


P("  %-6s %-6s %-12s %-22s %-10s %s"
  % ("rep", "rung", "(eta,eta')", "L = (5,6,7^k)", "e_R", "verdict"))
hosts = {}
for nm, bx, ty in (("21", 2, "anti"), ("28", 2, "sym"),
                   ("35", 3, "anti"), ("84", 3, "sym")):
    k = bx - 2
    nL = [0] * 7
    nL[4], nL[5], nL[6] = 1, 1, k
    nE = [0] * 7
    nE[5], nE[6] = 1, k + 1
    if ty == "anti" and (max(nL) > 1 or max(nE) > 1):
        P("  %-6s %-6d %-12s %-22s %-10s %s"
          % (nm, k, "-", "-", "-", "e_R needs a repeated index: NOT in an "
             "antisymmetric rep"))
        continue
    nL, nE = tuple(nL), tuple(nE)
    for eta in (1, -1):
        for etap in (1, -1):
            gL, gE = mode(nL, eta, etap), mode(nE, eta, etap)
            ok = (gL == -1 and gE == +1)
            if ok:
                hosts.setdefault(k, []).append((nm, eta, etap))
            P("  %-6s %-6d %-12s %-22s %-10s %s"
              % (nm, k, "(%+d,%+d)" % (eta, etap),
                 {-1: "L", 1: "R", None: "projected out"}[gL],
                 {-1: "L", 1: "R", None: "out"}[gE],
                 "**HOSTS IT**" if ok else "no"))
P("")
P("  >> rung 0 hosted by : %s" % hosts.get(0))
P("  >> rung 1 hosted by : %s" % hosts.get(1))
P("")
P("  CONTROL  rung 0 must come out at the 21 with eta = -1, eta' = +1 --")
P("           Komori-Maru's own eq. (43), which nothing here put in")
assert ("21", -1, 1) in hosts.get(0, [])
P("           it does                                                : OK")
P("")
P("  >> THE 35 IS RULED OUT for rung 1 and the reason is one line: e_R needs")
P("     the index content (6,7,7), a repeated index, which an ANTIsymmetric")
P("     tensor does not have.  The charged lepton, not the doublet, decides it.")
P("  >> the 84 hosts rung 1, and at (eta,eta') = (+1,+1).")
P("")
T1 = {"(1)": [("28", 1, -1, 1), ("84", 1, 1, 4)],
      "(2)": [("28", 1, 1, 1), ("84", 1, 1, 4)],
      "(3)": [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)],
      "(4)": [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)],
      "(5)": [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)]}
P("  AND IT IS ALREADY THERE.  Their Table 1, rows containing an 84(+,+):")
for k in sorted(T1):
    m = [(r, a, b, mm) for r, a, b, mm in T1[k] if (r, a, b) == ("84", 1, 1)]
    P("     %-5s %s" % (k, "%d x 84(+,+)" % m[0][3] if m else "none"))
P("")
P("  >> every row of their Table 1 already contains 84(+,+) multiplets,")
P("     introduced for the Higgs potential and for nothing else.  The multiplet")
P("     that would host a rung-1 lepton generation is not a new ingredient.")

P("")
P("  CONTROL  an SU(2)_L doublet is only a doublet if BOTH components survive")
P("           with the SAME chirality.  Check the upper component (4,6,7^k)")
P("           against the lower one used above, on every hosting rep:")
for k, hs in sorted(hosts.items()):
    for nm, eta, etap in hs:
        up = [0] * 7
        up[3], up[5], up[6] = 1, 1, k
        lo = [0] * 7
        lo[4], lo[5], lo[6] = 1, 1, k
        gu, gl = mode(tuple(up), eta, etap), mode(tuple(lo), eta, etap)
        P("           %-4s rung %d (%+d,%+d):  upper %s / lower %s   %s"
          % (nm, k, eta, etap, {-1: "L", 1: "R", None: "out"}[gu],
             {-1: "L", 1: "R", None: "out"}[gl],
             "OK" if gu == gl == -1 else "*** BROKEN DOUBLET ***"))
        assert gu == gl == -1
P("           the doublet is never split -- indices 4 and 5 carry identical")
P("           parities in both P5 and P5', so this could not have failed here;")
P("           it is recorded because a rep where it DID fail would be excluded.")

P("")
P("=" * 78)
P("STEP 8 -- and what the neutrinos then cost, rep by rep")
P("=" * 78)
P("  From su7_anomaly_channels.py the signed SM-singlet charges available are")
P("  fixed by (rep, parity).  Map each STEP 4 solution onto them:")
SING = {F(1, 2): "7(+,-)", F(-1, 2): "7(-,+)", F(-1): "28(+,+)",
        F(1): "28(-,-)", F(3, 2): "84(+,-)", F(-3, 2): "84(-,+)"}
P("")
P("  %-22s %-9s %-26s %s" % ("(l1,l2,l3)", "a", "neutrinos", "as (rep,parity)"))
best = None
for ls, a, As, nus in sols:
    if not all(v in SING for v in nus):
        continue
    tag = " + ".join(SING[v] for v in nus)
    if ls == (H, H, F(0)) and best is None:
        best = (ls, a, As, nus, tag)
    P("  %-22s %-9s %-26s %s"
      % (str(tuple(str(x) for x in ls)), a,
         str(tuple(str(v) for v in nus)), tag))
P("")
if best:
    ls, a, As, nus, tag = best
    P("  >> the minimal-change assignment, l = (1/2, 1/2, 0):")
    P("        two lepton generations in 21s, the third in an 84(+,+)")
    P("        brane-quark charge  X_Q = %s  (universal, so CKM is unspoiled)" % a)
    P("        neutrinos           %s" % tag)
    P("        proton operators    forbidden for ALL THREE generations,")
    P("                            |A_j| = %s"
      % ", ".join(str(abs(v)) for v in As))
    P("     and 28(+,+) is present in their cases (2) and (3), while 7(-,+) is")
    P("     present in none of the five -- so this costs ONE new multiplet.")

P("")
P("=" * 78)
P("STEP 9 -- and then it is not a suppression, it is a SELECTION RULE")
P("=" * 78)
P("  Their model requires U(1)' to be broken by a 4D scalar phi on the fixed")
P("  point (above their eq. 80).  They do NOT state its charge; su7_u1prime.py")
P("  identified a candidate, the (1,1)_0 of the 7, with q_phi = 1/2.")
P("")
P("  A VEV of charge q_phi does not destroy the symmetry, it breaks")
P("  U(1)' -> a discrete subgroup: charges are conserved MODULO q_phi.  An")
P("  operator of charge A can be dressed with n insertions of <phi> or <phi>*")
P("  only if  A / q_phi  is an INTEGER.  Otherwise it is forbidden to ALL")
P("  ORDERS -- which is what their question asks for ('by some symmetries').")
P("")
QPHI = [F(1, 2), F(1), F(3, 2), F(2), F(1, 6), F(1, 3)]
P("  %-22s %-26s %s" % ("(l1,l2,l3)", "(A1,A2,A3)", "q_phi -> generations still"))
P("  %-22s %-26s %s" % ("", "", "forbidden to ALL orders"))
for ls, a, As, nus in sols:
    if not all(v in SING for v in nus):
        continue
    cells = []
    for q in QPHI:
        n = sum(1 for v in As if (v / q).denominator != 1)
        cells.append("%s:%d/3" % (q, n))
    P("  %-22s %-26s %s"
      % (str(tuple(str(x) for x in ls)), str(tuple(str(v) for v in As)),
         "  ".join(cells)))
P("")
P("  CONTROL that MUST fire: the family-universal case has A_j = 0, and 0/q is")
P("  an integer for every q, so every generation is ALLOWED -- no selection")
P("  rule exists there, at any q_phi.")
for q in QPHI:
    assert (F(0) / q).denominator == 1
P("           verified at every q_phi above                          : OK")
P("")
lsb, ab, Asb, nub, tagb = best
P("  >> AT THE MINIMAL ASSIGNMENT l = (1/2, 1/2, 0), X_Q = %s:" % ab)
for q in (F(1, 2), F(1)):
    r = [str(v / q) for v in Asb]
    n = sum(1 for v in Asb if (v / q).denominator != 1)
    P("        q_phi = %-5s  A_j/q_phi = %-26s  forbidden: %d of 3"
      % (q, ", ".join(r), n))
P("")
P("     with their own natural q_phi = 1/2, A_j/q_phi = (1/3, 1/3, -2/3):")
P("     NONE is an integer, so ALL FOUR dimension-6 operators are forbidden")
P("     for ALL THREE generations, EXACTLY, to all orders in <phi>/M --")
P("     by a residual DISCRETE GAUGE symmetry, which is a Krauss-Wilczek")
P("     object and needs no supersymmetry.")
assert all((v / F(1, 2)).denominator != 1 for v in Asb)
P("     verified                                                     : OK")
P("")
P("  >> and that is a direct answer to the sentence they wrote:")
P("     'How the dangerous proton decay processes can be forbidden, for")
P("      instance, by some symmetries ... are left for our future work.'")

P("")
P("=" * 78)
P("STEP 10 -- THE BILL, and it is charged against their own Table 1")
P("=" * 78)
P("  Their own words, just below their eq. (79):")
P("     'For the leptons to be embedded into the bulk fermion, the bulk fermion")
P("      must be massive to generate lepton masses.  The potential from such a")
P("      massive bulk fermion has an extra factor e^{-M pi R} ~ e^{-O(10)} and")
P("      is very suppressed ... we do not consider the potential from the SM")
P("      fermion contributions in this paper.'")
P("")
P("  So a multiplet that HOSTS leptons is massive, and a massive multiplet does")
P("  NOT contribute to the Higgs potential.  Putting the third generation in an")
P("  84(+,+) therefore REMOVES one massless 84 from the potential of that row:")
P("")
P("  %-6s %-26s %-24s %s"
  % ("case", "their content", "after donating one 84", "still their row?"))
for k in sorted(T1):
    txt = " + ".join("%d x %s(%+d,%+d)" % (m, r, a, b) for r, a, b, m in T1[k])
    new = []
    ok = False
    for r, a, b, m in T1[k]:
        if (r, a, b) == ("84", 1, 1):
            ok = True
            if m - 1:
                new.append("%d x 84(+,+)" % (m - 1))
        else:
            new.append("%d x %s(%+d,%+d)" % (m, r, a, b))
    P("  %-6s %-26s %-24s %s"
      % (k, txt, " + ".join(new) if new else "(nothing left)",
         "NO -- not a row of their Table 1" if ok else "n/a"))
P("")
P("  >> the donated 84 is exactly the multiplet their potential analysis")
P("     counted as massless.  None of the five reduced contents is a row of")
P("     their Table 1, so m_h = 125-127 GeV and 1/R5 would have to be re-")
P("     derived.  That is THEIR computation (their eqs. 68-76), not one that")
P("     can be done from the quantum numbers here.")
P("  >> the honest shape: the family-dependent escape is CONSISTENT at the")
P("     level of anomalies, charges and zero modes, and it is UNTESTED at the")
P("     level of the vacuum.  It is a well-posed question for the authors, not")
P("     a defect in their paper.")

P("")
P("=" * 78)
P("WHAT IS AND IS NOT ESTABLISHED")
P("=" * 78)
P("  ESTABLISHED:")
P("    - extra(L) is a LADDER (1-k)/2 set by the number of index-7 boxes, so")
P("      family dependence is available in this model without new gauge")
P("      structure -- only by putting generations in different-size multiplets.")
P("    - the Yukawa closes on every rung: extra(L) + extra(e^c) = 1/2 always.")
P("    - SU(2)^2 X = (1/2) sum_j A_j EXACTLY, so anomaly cancellation demands")
P("      sum A_j = 0 while protection demands every A_j != 0: incompatible at")
P("      N = 1, compatible at N = 3.  That IS arXiv:1001.0768's structure.")
P("    - with universal brane-quark charge (exact CKM), A_j = l_j - lbar.")
P("    - rung 1 IS kept by the orbifold, in the 84 at (eta,eta') = (+1,+1),")
P("      with L left-handed and e_R right-handed and the doublet unsplit -- and")
P("      the control that rung 0 comes out at the 21 with their own eq. (43)")
P("      parities passes, with nothing put in by hand.")
P("    - 14 (rung set, X_Q) assignments cancel EVERY anomaly channel with all")
P("      A_j != 0, and each names the (rep, parity) of the neutrinos it needs.")
P("    - at q_phi = 1/2 the residual DISCRETE symmetry forbids all four")
P("      operators for all three generations to ALL ORDERS -- a selection rule,")
P("      not a suppression.")
P("  NOT ESTABLISHED, and each of these is a real gap:")
P("    - q_phi.  Their paper does not state the charge of the U(1)'-breaking")
P("      scalar.  1/2 is su7_u1prime.py's candidate; at q_phi = 1/6 or finer")
P("      the selection rule fails (printed in STEP 9).  This is the SAME")
P("      unstated-global-choice dependence su7_socratic.py already found.")
P("    - the potential.  STEP 10: donating an 84 to the leptons removes it from")
P("      the potential by their own argument, and no reduced content is a row")
P("      of their Table 1.  The vacuum must be re-derived.")
P("    - the lepton mass hierarchy, which a rung difference will affect and")
P("      which is not computed.")
P("    - that the mixed U(1) anomalies must cancel at all (Green-Schwarz).")
P("    - the exotic burden: the 84(+,+) has 44 zero modes, of which 3 would be")
P("      kept; the rest need conjugate brane partners, as their own")
P("      prescription already requires for every other row.")
P("")
P("DONE")
