#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

ALL the anomaly channels of Komori-Maru's U(1)', not just the SU(2)_L one --
and the paper's own pairing prescription, which su7_signed_census.py did not use.

su7_signed_census.py closed a chain:  Table 1 row -> S(added) -> X_Q -> the
exponent of the proton-decay suppression.  It rested on TWO things it flagged
itself:

  (i)  only the U(1)' x [SU(2)_L]^2 channel was computed;
  (ii) every bulk zero mode was counted.

(ii) is contradicted by the paper.  Immediately after their eq. (76), verbatim:

    "Introducing the fermions like above, many unwanted massless fermions are
     remained as zero modes.  In order to make these massless fermions heavy,
     we have to introduce the 4D fermion CONJUGATE to the representations and
     charges of the massless fermions and their Dirac mass terms on the fixed
     points."

A conjugate 4D fermion is a localised fermion and contributes to the SAME
localised anomaly with the OPPOSITE sign.  Every unwanted zero mode is therefore
paired away in every channel, and the anomaly of the 4D theory is the anomaly of
the SURVIVING chiral content: the brane quarks and the leptons.  Nothing else.

So this script does (i) and (ii) together.

INPUTS, all transcribed from the paper, none assumed:
  eq. (27)  T15 = (1/2) diag(0,0,0,1,-1,0,0)              = T3L
  eq. (31)  Qem = T15 + sqrt(3) T24 = diag(0,0,0,1,0,-1,0)
            =>  Y = Qem - T3L = (0,0,0,1/2,1/2,-1,0)   [derived, then controlled
                against their eq. (41) hypercharges 0, 1/2, 0, -1 on the 7]
  eq. (79)  U(1)' = (1/2) diag(0,0,0,1,-1,1,-1)
  eq. (43)  the 21's parities: eta = -1, eta' = +1
  eq. (46)  lepton Yukawa  =>  L = 21^{56}, e_R = 21^{67}  (su7_centre.py)
  eq. (47)  brane quark Yukawas  =>  extra(u_R) = a+1/2, extra(d_R) = a-1/2,
            a := extra(q_L) = X_Q free                     (su7_XQ.py)

Anomaly coefficients are computed state by state with Cartan generators, so that
the nonabelian index needs no representation theory:
    sum over the states of a fundamental of  T3^2  =  1/2  =  T(fund),
which is checked as a control before anything else is summed.
"""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement as cwr

P = lambda *a: print(*a, flush=True)
H = F(1, 2)

# ----------------------------------------------------------------- their input
YIND = (F(0), F(0), F(0), H, H, F(-1), F(0))          # from eq. (31)
UP1 = (F(0), F(0), F(0), H, -H, H, -H)                # eq. (79)
QEM = (F(0), F(0), F(0), F(1), F(0), F(-1), F(0))     # eq. (31)
P5 = (1, 1, 1, 1, 1, -1, -1)
P5p = (1, 1, 1, -1, -1, -1, 1)

t3l = lambda n: F(n[3] - n[4], 2)
t3c = lambda n: F(n[0] - n[1], 2)
dot = lambda v, n: sum(c * k for c, k in zip(v, n))
e = lambda i: tuple(1 if k == i else 0 for k in range(7))
add = lambda *v: tuple(sum(c) for c in zip(*v))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))

REPS = {
    "7":  [e(i) for i in range(7)],
    "21": [add(e(i), e(j)) for i, j in combinations(range(7), 2)],
    "28": [add(e(i), e(j)) for i, j in cwr(range(7), 2)],
    "84": [add(e(i), e(j), e(k)) for i, j, k in cwr(range(7), 3)],
    # adjoint = 42 roots + the Cartan, which is 6 states of weight 0.  Those
    # carry no charge at all, so they change no anomaly sum -- but leaving them
    # out makes the zero-mode COUNT wrong by 6, so they go in.
    "48": ([sub(e(i), e(j)) for i in range(7) for j in range(7) if i != j]
           + [(0,) * 7] * 6),
}
assert [len(REPS[k]) for k in ("7", "21", "28", "48", "84")] == [7, 21, 28, 48, 84]


def par(M, n):
    s = 1
    for i, k in enumerate(n):
        if k % 2:
            s *= M[i]
    return s


def zero_modes(wts, eta, etap):
    """su7_signed_census.py's rule, from their eqs. (37)-(40)."""
    return [(n, -eta * par(P5, n)) for n in wts
            if par(P5, n) * par(P5p, n) == eta * etap]


P("=" * 78)
P("STEP 0 -- the transcribed generators, and the controls on them")
P("=" * 78)
P("  their eq. (31):  Qem = diag%s" % (tuple(map(str, QEM)),))
P("  T3L  = (1/2)(n4 - n5)                             their eq. (27)")
P("  Y    = Qem - T3L  =  diag%s        derived" % (tuple(map(str, YIND)),))
P("  U(1)'= diag%s      their eq. (79)" % (tuple(map(str, UP1)),))
P("")
tr = lambda v: sum(v)
P("  CONTROL 0a  Y traceless (a Cartan generator of SU(7))   : %s   (sum = %s)"
  % ("OK" if tr(YIND) == 0 else "FAILED", tr(YIND)))
assert tr(YIND) == 0 and tr(UP1) == 0 and tr(QEM) == 0
P("  CONTROL 0b  their eq. (41): the 7's hypercharges under SU(3)xSU(2)xU(1)_Y")
sev = sorted({dot(YIND, n) for n in REPS["7"]})
P("              computed from Qem - T3L : %s" % sev)
P("              their eq. (41) states   : 0, 1/2, 0, -1")
assert sev == [F(-1), F(0), H]
P("              %s" % "OK -- the derived Y reproduces their stated hypercharges")
P("")
P("  CONTROL 0c  Cartan normalisation, so no index has to be assumed:")
for nm, gen, wts in (("SU(2)_L doublet", t3l, [(0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0)]),
                     ("SU(3)_C triplet", t3c, [e(0), e(1), e(2)])):
    s = sum(gen(n) ** 2 for n in wts)
    P("              sum of T3^2 over a %s = %s  (must be 1/2)" % (nm, s))
    assert s == H
P("              OK -- T(fund) = 1/2 comes out, it is not put in")

# ------------------------------------------------------- the surviving content
P("")
P("=" * 78)
P("STEP 1 -- the SURVIVING 4D chiral content, per the paper's own prescription")
P("=" * 78)
P("  After their eq. (76): every unwanted bulk zero mode is given a Dirac mass")
P("  with a 4D fermion CONJUGATE to it, at a fixed point.  A conjugate localised")
P("  fermion cancels that zero mode in EVERY anomaly channel.  What is left is")
P("  the SM chiral content: brane quarks (eq. 47) + leptons from the 21 (eq. 44).")
P("")


def content(a):
    """One generation, all left-handed Weyl.  (label, T3C, T3L, Y, X)."""
    out = []
    for i, tc in enumerate((H, -H, F(0))):                 # colour 3
        for tl in (H, -H):
            out.append(("Q", tc, tl, F(1, 6), a + tl))     # X = extra + T3L
        out.append(("uc", -tc, F(0), F(-2, 3), -(a + H)))  # 3bar
        out.append(("dc", -tc, F(0), F(1, 3), -(a - H)))
    for tl in (H, -H):
        out.append(("L", F(0), tl, -H, H + tl))
    out.append(("ec", F(0), F(0), F(1), F(0)))
    return out


a0 = F(0)
P("  %-5s %-8s %-8s %-8s %s" % ("field", "T3C", "T3L", "Y", "U(1)'"))
seen = set()
for lab, tc, tl, y, x in content(F(0)):
    if lab in seen:
        continue
    seen.add(lab)
    xs = {"Q": "a + T3L", "uc": "-(a + 1/2)", "dc": "-(a - 1/2)",
          "L": "1/2 + T3L", "ec": "0"}[lab]
    P("  %-5s %-8s %-8s %-8s %s" % (lab, "+-1/2,0" if tc else "0",
                                    "+-1/2" if tl else "0", y, xs))

# ----------------------------------------------------------------- the channels
CH = [
    ("SU(3)^2 x Y",   lambda s: sum(y * tc ** 2 for _, tc, tl, y, x in s)),
    ("SU(2)^2 x Y",   lambda s: sum(y * tl ** 2 for _, tc, tl, y, x in s)),
    ("Y^3",           lambda s: sum(y ** 3 for _, tc, tl, y, x in s)),
    ("Y x grav^2",    lambda s: sum(y for _, tc, tl, y, x in s)),
    ("SU(3)^2 x X",   lambda s: sum(x * tc ** 2 for _, tc, tl, y, x in s)),
    ("SU(2)^2 x X",   lambda s: sum(x * tl ** 2 for _, tc, tl, y, x in s)),
    ("X x grav^2",    lambda s: sum(x for _, tc, tl, y, x in s)),
    ("X^3",           lambda s: sum(x ** 3 for _, tc, tl, y, x in s)),
    ("X^2 x Y",       lambda s: sum(x ** 2 * y for _, tc, tl, y, x in s)),
    ("X x Y^2",       lambda s: sum(x * y ** 2 for _, tc, tl, y, x in s)),
]

P("")
P("=" * 78)
P("STEP 2 -- CONTROL that must NOT fire: the pure Standard-Model channels")
P("=" * 78)
sm_ok = True
for nm, f in CH[:4]:
    vals = {t: f(content(t)) for t in (F(0), F(1), F(-1, 6))}
    ok = set(vals.values()) == {F(0)}
    sm_ok &= ok
    P("  %-14s = %-6s at every X_Q      %s"
      % (nm, vals[F(0)], "OK" if ok else "*** FAILED ***"))
assert sm_ok
P("  >> one full SM generation with no nu_R cancels all four.  The content and")
P("     the state-by-state Cartan bookkeeping are therefore sound.")

P("")
P("=" * 78)
P("STEP 3 -- EVERY channel involving U(1)', as a polynomial in a = X_Q")
P("=" * 78)


def poly(f, deg=3):
    """Exact interpolation of a polynomial of degree <= deg in a."""
    xs = [F(k) for k in range(deg + 1)]
    ys = [f(content(t)) for t in xs]
    # finite differences -> Newton form -> expand into monomial coefficients
    co = [F(0)] * (deg + 1)
    basis = [F(1)] + [F(0)] * deg                       # product (a - x_j)
    div = list(ys)
    for k in range(1, deg + 1):
        for i in range(deg, k - 1, -1):
            div[i] = (div[i] - div[i - 1]) / (xs[i] - xs[i - k])
    for k in range(deg + 1):
        for i in range(deg + 1):
            co[i] += div[k] * basis[i]
        nb = [F(0)] * (deg + 1)                          # basis *= (a - x_k)
        for i in range(deg):
            nb[i + 1] += basis[i]
        for i in range(deg + 1):
            nb[i] -= xs[k] * basis[i]
        basis = nb
    return co


def show(co):
    t = []
    for i in range(len(co) - 1, -1, -1):
        if co[i]:
            t.append("%s%s" % (co[i], {0: "", 1: " a", 2: " a^2", 3: " a^3"}[i]))
    return " + ".join(t) if t else "0"


P("  %-14s %-28s %s" % ("channel", "coefficient", "vanishes at X_Q ="))
roots = {}
for nm, f in CH[4:]:
    co = poly(f)
    # roots over Q of a cubic with rational coefficients: rational root test
    rs = set()
    if any(co[1:]):
        cand = set()
        top = co[max(i for i in range(4) if co[i])]
        for p in range(1, 61):
            for q in range(1, 61):
                cand |= {F(p, q), F(-p, q)}
        cand.add(F(0))
        rs = {r for r in cand if sum(c * r ** i for i, c in enumerate(co)) == 0}
    ans = ("every X_Q" if not any(co) else
           ("NEVER" if not rs else ", ".join(str(r) for r in sorted(rs))))
    roots[nm] = ans
    P("  %-14s %-28s %s" % (nm, show(co), ans))

P("")
P("  CONTROL 3a  SU(2)^2 x X must reproduce su7_anomaly_XQ.py's  A/2 = (3a+1/2)/2")
for t in (F(0), F(1), F(-1, 6), F(-2, 3)):
    assert CH[5][1](content(t)) == (3 * t + H) / 2
P("              verified at four values of a                          : OK")
P("  CONTROL 3b  SU(3)^2 x X identically zero, because X_uc + X_dc = -2X_Q IS")
P("              their eq. (47)")
for t in (F(0), F(1), F(-1, 6)):
    assert CH[4][1](content(t)) == 0
P("              verified at three values of a                         : OK")

P("")
P("=" * 78)
P("STEP 4 -- is there ONE X_Q that cancels them all?")
P("=" * 78)
allch = [nm for nm, _ in CH[4:]]
common = None
for nm, f in CH[4:]:
    co = poly(f)
    if not any(co):
        continue
    s = {t for t in (F(k, 6) for k in range(-30, 31)) if
         sum(c * t ** i for i, c in enumerate(co)) == 0}
    common = s if common is None else (common & s)
P("  X_Q values (searched over k/6, |k| <= 30) killing every U(1)' channel: %s"
  % (sorted(common) if common else "NONE"))
P("")
P("  channel by channel, the obstruction:")
gv = poly(CH[6][1])
P("     X x grav^2   = %s   -- INDEPENDENT of X_Q, and NOT zero" % show(gv))
src = {}
for lab in ("Q", "uc", "dc", "L", "ec"):
    src[lab] = sum(x for l, tc, tl, y, x in content(F(0)) if l == lab)
P("     where does that 1 come from?  sum of U(1)' over each species at a = 0:")
for lab in ("Q", "uc", "dc", "L", "ec"):
    P("        %-4s %s" % (lab, src[lab]))
P("     and at general a the quark part cancels: 6a - 3(a+1/2) - 3(a-1/2) = 0.")
for t in (F(0), F(2), F(-1, 6)):
    assert sum(x for l, tc, tl, y, x in content(t) if l in ("Q", "uc", "dc")) == 0
P("     verified at three values of a                                   : OK")
P("")
P("  >> the U(1)'-gravitational anomaly of the surviving content is 1, sourced")
P("     ENTIRELY by the lepton doublet, and no assignment of the brane quarks'")
P("     charge can touch it.")

P("")
P("=" * 78)
P("STEP 5 -- what state would cancel it, and does their spectrum contain one?")
P("=" * 78)
P("  A left-handed SM-singlet (a right-handed neutrino) of U(1)' charge x adds")
P("  x to X x grav^2, x^3 to X^3, 0 to the nonabelian channels and 0 to X x Y^2.")
P("  Cancelling BOTH the grav and the cubic channel at once needs:")
x3 = poly(CH[7][1])
P("     X x grav^2 : 1 + x = 0            =>  x = -1")
P("     X^3        : (%s) + x^3 = 0" % show(x3))
for t in (F(0), F(-1, 6), F(1, 2)):
    v = sum(x ** 3 for _, tc, tl, y, x in content(t))
    P("        at X_Q = %-6s  X^3 = %-8s  =>  needs x^3 = %-8s  i.e. x = %s"
      % (t, v, -v, "-1  CONSISTENT" if -v == F(-1) else "not -1"))
P("")
P("  So a single SM singlet of U(1)' charge -1 fixes the gravitational channel;")
P("  whether it also fixes the cubic one is a CONDITION ON X_Q, printed above.")

P("")
P("=" * 78)
P("STEP 6 -- the bookkeeping su7_signed_census.py used, side by side")
P("=" * 78)
P("  It summed the SU(2)^2 channel over EVERY bulk zero mode, with no conjugate")
P("  partners.  Below: the same sum, for every channel, on the 21 at their")
P("  eq. (43) parities -- so the size of what the pairing removes is visible.")


def bulk(wts, eta, etap):
    out = []
    for n, g in zero_modes(wts, eta, etap):
        out.append((-g, t3c(n), t3l(n), dot(YIND, n), dot(UP1, n)))
    return out


P("")
P("  %-6s %-12s %-7s %-11s %-11s %-11s %-11s %s"
  % ("rep", "(eta,eta')", "#zero", "SU3^2 X", "SU2^2 X", "X grav", "X^3", "SU2^2 Y"))
for r in ("7", "21", "28", "48", "84"):
    for eta, etap in ((-1, 1),) if r == "21" else ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        st = bulk(REPS[r], eta, etap)
        f = lambda k: sum(w * v for w, *rest in st for v in [k(rest)])
        P("  %-6s %-12s %-7d %-11s %-11s %-11s %-11s %s"
          % (r, "(%+d,%+d)" % (eta, etap), len(st),
             f(lambda z: z[3] * z[0] ** 2), f(lambda z: z[3] * z[1] ** 2),
             f(lambda z: z[3]), f(lambda z: z[3] ** 3),
             f(lambda z: z[2] * z[1] ** 2)))
P("")
P("  CONTROL 6a  the 21 at their eq. (43) parities must give SU2^2 X = 1/4,")
P("              the value su7_signed_census.py reports")
st21 = bulk(REPS["21"], -1, 1)
v = sum(w * x * tl ** 2 for w, tc, tl, y, x in st21)
P("              computed: %s                                          : %s"
  % (v, "OK" if v == F(1, 4) else "*** FAILED ***"))
assert v == F(1, 4)
P("")
P("  >> and the other channels of that SAME multiplet are NOT zero, so the")
P("     census's single-channel condition was never the whole condition even on")
P("     its own bookkeeping.")

P("")
P("=" * 78)
P("STEP 7 -- WHAT U(1)' THEN IS, and why that answers their question structurally")
P("=" * 78)
P("  Three channels -- SU(2)^2 X, X^2 Y, X Y^2 -- independently force X_Q = -1/6,")
P("  and the remaining two force one SM singlet of U(1)' charge -1.  Write out the")
P("  whole assignment and compare it with T3L, Y and B-L:")
P("")
BL = {"Q": F(1, 3), "uc": F(-1, 3), "dc": F(-1, 3), "L": F(-1), "ec": F(1),
      "nc": F(1)}
full = content(F(-1, 6)) + [("nc", F(0), F(0), F(0), F(-1))]
P("  %-5s %-9s %-9s %-9s %-9s %s"
  % ("field", "T3L", "Y", "B-L", "U(1)'", "T3L + Y - (B-L)"))
ok7 = True
seen = set()
for lab, tc, tl, y, x in full:
    if (lab, tl) in seen:
        continue
    seen.add((lab, tl))
    pred = tl + y - BL[lab]
    ok7 &= (pred == x)
    P("  %-5s %-9s %-9s %-9s %-9s %-9s %s"
      % (lab, tl, y, BL[lab], x, pred, "" if pred == x else "*** MISMATCH ***"))
assert ok7
P("")
P("  >> U(1)'  =  T3L + Y - (B-L)   on every field, exactly, with no freedom left.")
P("     The anomaly conditions do not merely PICK a number; they identify the")
P("     leftover U(1) as a combination of T3L, Y and B-L and nothing else.")
P("")
P("  And the four dimension-6 |dB|=1 operators are invariant under all three.")
P("  Computed, not asserted -- their B-L and Y and T3L, from the charges above:")
OPS = [("QQQL", [("Q", 1), ("Q", 1), ("Q", 1), ("L", 1)]),
       ("Q uc* dc* L", [("Q", 1), ("uc", -1), ("dc", -1), ("L", 1)]),
       ("uc uc dc ec", [("uc", 1), ("uc", 1), ("dc", 1), ("ec", 1)]),
       ("QQ uc* ec*", [("Q", 1), ("Q", 1), ("uc", -1), ("ec", -1)])]
YY = {"Q": F(1, 6), "uc": F(-2, 3), "dc": F(1, 3), "L": -H, "ec": F(1)}
P("  %-14s %-8s %-8s %-8s %s" % ("operator", "Y", "B-L", "U(1)'", "verdict"))
for nm, fs in OPS:
    y = sum(s * YY[f] for f, s in fs)
    bl = sum(s * BL[f] for f, s in fs)
    x = y - bl                      # T3L cancels in any SU(2)_L invariant
    P("  %-14s %-8s %-8s %-8s %s"
      % (nm, y, bl, x, "ALLOWED" if x == 0 else "forbidden"))
    assert bl == 0 and x == 0
P("")
P("  >> all four conserve B-L (textbook: the SM dimension-6 baryon-number")
P("     operators are B-L neutral), and all four are Y-neutral, so ANY U(1) in")
P("     span{T3L, Y, B-L} is blind to them.")
P("")
P("  >> THE ANSWER TO KOMORI-MARU'S OPEN QUESTION IS THEREFORE STRUCTURAL, not")
P("     a coincidence at one point of a free parameter: their leftover U(1)'")
P("     cannot forbid proton decay, because anomaly cancellation forces it into")
P("     the span of T3L, Y and B-L, and the dangerous operators are neutral")
P("     under every member of that span.")

P("")
P("=" * 78)
P("STEP 8 -- does their own spectrum contain the SM singlet of charge -1?")
P("=" * 78)
P("  Requirement: a 4D left-handed Weyl, colour singlet, SU(2)_L singlet, Y = 0,")
P("  U(1)' = -1.  That is a RIGHT-HANDED NEUTRINO.  The paper mentions none.")
P("  Search every rep it uses, at every parity choice, over its zero modes:")
P("")
P("  The condition is on the KEPT set, not on one state: with signed charges")
P("  {x_i} of the SM-singlet zero modes one chooses NOT to pair away,")
P("     sum x_i = -1   (X x grav^2)     and     sum x_i^3 = -1   (X^3).")
P("  Two singlets of charge -1/2 satisfy the first and NOT the second, so the")
P("  two channels are independent conditions and the search must be over subsets.")
P("")


def sm_singlets(r, eta, etap):
    """(signed U(1)' charge) of every colour/weak-singlet, Y=0 zero mode."""
    out = []
    for n, g in zero_modes(REPS[r], eta, etap):
        if any(n[i] for i in range(5)) or dot(YIND, n) != 0:
            continue
        out.append((n, -g * dot(UP1, n), g))   # -g = +1 for a left-handed mode
    return out


P("  %-6s %-12s %-22s %-6s %s"
  % ("rep", "(eta,eta')", "component", "chir", "signed U(1)'"))
for r in ("7", "21", "28", "48", "84"):
    for eta in (1, -1):
        for etap in (1, -1):
            for n, sx, g in sm_singlets(r, eta, etap):
                if sx == 0:
                    continue
                P("  %-6s %-12s %-22s %-6s %s"
                  % (r, "(%+d,%+d)" % (eta, etap), str(n),
                     {-1: "L", 1: "R"}[g], sx))
P("")
P("  >> only the pure index-7 components are SM singlets with Y = 0, so the")
P("     available charges are exactly  7 -> -+1/2,  28 -> -+1,  84 -> -+3/2,")
P("     and the 48's is 0.  Nothing else in their reps can do this job.")

P("")
P("  THEIR TABLE 1, with both conditions imposed on the kept singlets")
P("  (parities as printed in their table caption, (eta_R, eta'_R)):")
P("")
T1F = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)], "3.8"),
       ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)], "2.0"),
       ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], "7.5"),
       ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], "6.1"),
       ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], "3.8")]
P("  %-5s %-9s %-30s %s" % ("case", "1/R5", "singlet charges available", "solvable?"))
for tag, cont, iR in T1F:
    pool = []
    for r, a, b, m in cont:
        pool += [sx for _ in range(m) for _, sx, _g in sm_singlets(r, a, b) if sx]
    sol = None
    for k in range(1, len(pool) + 1):
        for c in combinations(range(len(pool)), k):
            v = [pool[i] for i in c]
            if sum(v) == F(-1) and sum(x ** 3 for x in v) == F(-1):
                sol = v
                break
        if sol:
            break
    P("  %-5s %-9s %-30s %s"
      % (tag, iR, str(sorted(pool)),
         "YES  keep %s" % sol if sol else "NO -- no subset works"))
P("")
P("  >> the two rows that cannot do it are exactly the two whose only singlets")
P("     come from a 7 (charge -+1/2) or an 84 (-+3/2).  A charge -1 singlet needs")
P("     a 28, and the 28 appears in cases (1), (2) and (3).")
P("")
P("  >> the state is in the spectrum as an UNWANTED zero mode -- exactly the")
P("     class the paper's sentence pairs away with a conjugate brane fermion.")
P("     Pairing it away is what leaves the anomaly.  KEEPING one such state is")
P("     not optional: it is what makes the model consistent, and it is a")
P("     RIGHT-HANDED NEUTRINO, which the paper does not discuss.")
P("")
P("  SCOPE: the paper writes one generation explicitly.  For N generations the")
P("  lepton doublets give X grav^2 = N and X^3 = N, so N singlets of charge -1")
P("  are needed -- the same statement, N times, and the same reps supply it.")

P("")
P("=" * 78)
P("WHAT IS AND IS NOT ESTABLISHED")
P("=" * 78)
P("  ESTABLISHED, and it is arithmetic on their own equations:")
P("    - with the paper's own pairing prescription the surviving 4D chiral")
P("      content is the SM one, and then NO value of X_Q cancels every U(1)'")
P("      anomaly channel; the gravitational one is X_Q-independent and equal to 1.")
P("    - so su7_signed_census.py's chain (Table 1 -> X_Q -> proton exponent)")
P("      does not survive the pairing sentence: S(added) is cancelled by the")
P("      conjugate brane fermions the paper introduces.")
P("    - the three channels that CAN cancel all force the same X_Q = -1/6, which")
P("      is exactly A = 0: the point at which the U(1)' stops protecting.")
P("    - at that point U(1)' = T3L + Y - (B-L) identically, and the four")
P("      dimension-6 operators are neutral under it for a structural reason.")
P("    - the two channels that cannot cancel demand one SM singlet of charge -1")
P("      per generation, i.e. a right-handed neutrino, which their spectrum")
P("      contains and their prescription pairs away.")
P("  NOT ESTABLISHED:")
P("    - that the localised anomaly must cancel on localised matter: the gate")
P("      (von Gersdorff-Quiros) says mixed U(1) anomalies are the evadable class,")
P("      and a Green-Schwarz two-form remains available for all of them.")
P("    - WHERE each conjugate fermion sits: the paper says 'on the fixed points'")
P("      without saying which, and a localised (as opposed to integrated) anomaly")
P("      is fixed-point by fixed-point.  What is computed here is the INTEGRATED")
P("      4D anomaly, which is the one that cannot depend on that choice.")
P("    - anything about coefficients, <phi>, or a lifetime.")
P("")
P("DONE")
