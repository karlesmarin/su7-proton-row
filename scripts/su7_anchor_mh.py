#!/usr/bin/env python3
"""The anchor, attacked through the column nobody tested: m_h.

Carles Marin + Claude (AI assistant).  2026-08-04.
Source: Komori & Maru, arXiv:2503.04090 (text in ../_papers/SU7_GGHU_2503.04090.txt).

Three parts.

A. THE '7' OF EQ. (68), REBUILT FROM THEIR SECTION 3.1.
   SU7_VACUUM.md recorded it as "the one number that cannot be rebuilt without
   their 3.1 in full".  3.1 is now read.  Eq. (62) lists eight (N, r, P6, P5, P5')
   entries -- four for A_mu (N=4) and four for A_5,A_6 (N=1), the second quartet
   being the first with ALL THREE parities flipped -- and eq. (67) annihilates
   every P6 = -1 entry.  Assembling them must give {2, 4, 7}.  No ghosts appear
   anywhere in their 3.1; N=4 for A_mu and N=1 for A_5,6 (their footnote 3).

B. THE EXACT BRACE, IN CLOSED FORM.
   Their eq. (66) summed and combined with the +2/(pi n)^6 term collapses to
       4*Sum_m + 2/(pi n)^6 + 3 k P6/(4 (pi n)^5)
           = 3k/(4 (pi n)^5) * [ P6 + g(pi k n) ],
       g(x) = coth x + x csch^2 x + (2/3) x^2 coth x csch^2 x,   g(inf) = 1,
   which is eq. (67).  Verified against direct summation.  This is what the
   "exact eqs. (62)-(65)" repair really is, and it bounds its size.

C. THE m_h TEST -- new, and parameter-free.
   Their eq. (80)  m_h^2 = 4 pi^2 g4^2 R5^3 R6 * d2V/da2|min
   with V = C*F, C = 3k/(64 pi^8 R5^6), k = R5/R6, gives
       m_h^2 = 3 g4^2 F''(a_min) / (16 pi^6 R5^2),
   and R6 and k cancel identically.  Their eq. (82) 1/R5 = 2*80.4/a_min then gives

       K = m_h * a_min / sqrt(F''(a_min)) = 160.8 * g4 * sqrt(3/(16 pi^6))

   -- THE SAME NUMBER FOR ALL FIVE ROWS OF THEIR TABLE 1.

   K is invariant under F -> lambda F for any content-INDEPENDENT lambda, so it
   does not test the normalisation (every normalisation repair was already
   excluded); it tests the SHAPE of F content by content, using a column of their
   table that this programme never read.  Constant K would say their potential is
   ours and only the location of their minimum disagrees; non-constant K says the
   content-dependence itself is wrong.
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


_BCACHE = {}


def basis(alpha, s, c, d=0):
    a = np.atleast_1d(alpha)
    key = (a.tobytes(), s, c, d)
    hit = _BCACHE.get(key)
    if hit is not None:
        return hit
    ph = np.outer(a, c * math.pi * _nf)
    f = (np.cos, lambda t: -np.sin(t), lambda t: -np.cos(t),
         lambda t: np.sin(t))[d % 4](ph)
    out = (f * (_w * _sgn[s] * (c * math.pi * _nf) ** d)).sum(axis=1)
    if a.size > 1000:                       # cache only the big fixed grids
        _BCACHE[key] = out
    return out


GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]        # eq. (68) x (1/C)


def terms(rep, eta, etap, r4=True):
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


def V(content, alpha, d=0, r4=True, lam=1.0):
    out = sum(m * basis(alpha, s, c, d) for m, s, c in GAUGE)
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap, r4):
            out = out + lam * m * mult * basis(alpha, s, c, d)
    return out


GRID = np.linspace(0.0, 1.0, 40001)


def minimise(content, r4=True, lam=1.0):
    v = V(content, GRID, r4=r4, lam=lam)
    i = int(np.argmin(v))
    lo, hi = GRID[max(i - 1, 0)], GRID[min(i + 1, len(GRID) - 1)]
    for _ in range(25):
        xs = np.linspace(lo, hi, 15)
        j = int(np.argmin(V(content, xs, r4=r4, lam=lam)))
        lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 14)]
    return float(0.5 * (lo + hi))


T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.8, 3.8),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)], 0.081, 125.5, 2.0),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], 0.021, 125.1, 7.5),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], 0.026, 126.4, 6.1),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.2, 3.8)]

REC = {"source": "Komori & Maru, arXiv:2503.04090", "NMAX": NMAX, "steps": {}}

# ----------------------------------------------------------------------------
P("=" * 78)
P("A -- THE '7' OF EQ. (68), REBUILT FROM SECTION 3.1")
P("=" * 78)
P("  eq. (62), verbatim, as (N, r, P6, P5, P5'):")
# (label, N, r, P6, P5, P5', multiplicity)  -- eq. (62), first bracket then second
EQ62 = [("A_mu   (1,3)", 4, 3, +1, +1, +1, 1),
        ("A_mu   (1,2)", 4, 2, +1, -1, -1, 2),
        ("A_mu   (1,2)", 4, 2, +1, -1, +1, 2),
        ("A_mu   (3,2)", 4, 2, -1, +1, -1, 6),
        ("A_5,A_6 (1,3)", 1, 3, -1, -1, -1, 1),
        ("A_5,A_6 (1,2)", 1, 2, -1, +1, +1, 2),
        ("A_5,A_6 (1,2)", 1, 2, -1, +1, -1, 2),
        ("A_5,A_6 (3,2)", 1, 2, +1, -1, +1, 6)]
P("  %-14s %3s %3s %3s %3s %4s %5s  %-14s %s"
  % ("field/multiplet", "N", "r", "P6", "P5", "P5'", "mult", "eq.(67) 3(1+P6)/4",
     "-> contributes"))
acc = {}                                    # (cos argument, sign) -> weight
for lab, N, r, p6, p5, p5p, mult in EQ62:
    w67 = 3 * (1 + p6) / 4.0                # eq. (67), in units of k/(pi n)^5
    s = p5 * p5p
    contrib = N * mult * w67
    key = (r - 1, s)
    acc[key] = acc.get(key, 0.0) + contrib
    P("  %-14s %3d %3d %+3d %+3d %+4d %5d  %-17s  %s"
      % (lab, N, r, p6, p5, p5p, mult, "%.1f" % w67,
         "-- killed --" if contrib == 0 else
         "%.0f x %s^n cos(%d pi n a)" % (contrib, "(+1)" if s > 0 else "(-1)",
                                         r - 1)))
P("")
# eq. (62)'s overall prefactor is -1/(128 pi^3); eq. (68) is normalised to
# -3k/(128 pi^8 R5^6), i.e. the weights above divided by 3/2 (= 6/4, the surviving
# eq.(67) factor) and then by 2.
got = {k: v / (1.5 * 2.0) for k, v in acc.items() if v != 0}
P("  divide by the surviving eq.(67) factor 6/4 and by 2 (the 1/128 -> 3k/128 pi^8):")
for (c, s), v in sorted(got.items()):
    P("     %s^n cos(%d pi n a) : %.1f" % ("(+1)" if s > 0 else "(-1)", c, v))
WANT = {(2, 1): 2.0, (1, 1): 4.0, (1, -1): 7.0}
ok = got == WANT
P("")
P("  their eq. (68):  {2 cos(2 pi n a) + 4 cos(pi n a) + 7 (-1)^n cos(pi n a)}")
P("  rebuilt       :  %s                                     %s"
  % ({("(-1)^n" if s < 0 else "(+1)^n", c): v for (c, s), v in sorted(got.items())},
     "PASS" if ok else "*** FAIL ***"))
assert ok, got
# where the 7 comes from
mu7 = 4 * 2 * 1.5 / 3.0                      # A_mu  (1,2)(+,-,+), N=4, mult 2
a567 = 1 * 6 * 1.5 / 3.0                     # A_5,6 (3,2)(+,-,+), N=1, mult 6
P("")
P("  >> THE 7 SPLITS AS 7 = %.0f + %.0f, AND NO GHOST ENTERS:" % (mu7, a567))
P("     %.0f  from A_mu (N=4) on the two colour-SINGLET doublets (1,2)+(1,2bar)"
  % mu7)
P("        carrying (P6,P5,P5') = (+,-,+) in their eq. (57);")
P("     %.0f  from A_5,A_6 (N=1, their footnote 3) on the SIX COLOURED states"
  % a567)
P("        (3,2bar)+(3bar,2), whose P6 = -1 in (57) flips to +1 for the")
P("        extra-dimensional components -- which is why eq. (67) keeps them")
P("        for A_5,6 and kills them for A_mu, and vice versa for the doublets.")
P("     So the P6 split is a RELABELLING of the same eight multiplets, not an")
P("     extra ingredient.  Their eq. (68) is confirmed exactly.")
assert abs(mu7 - 4) < 1e-12 and abs(a567 - 3) < 1e-12

P("")
P("  CONTROL THAT COULD HAVE FAILED -- do NOT flip P6 between A_mu and A_5,6")
P("  (i.e. give the A_5,6 quartet the same parities as the A_mu quartet):")
acc_b = {}
for lab, N, r, p6, p5, p5p, mult in EQ62[:4]:
    for NN in (4, 1):
        w67 = 3 * (1 + p6) / 4.0
        acc_b[(r - 1, p5 * p5p)] = acc_b.get((r - 1, p5 * p5p), 0.0) + NN * mult * w67
bad = {k: v / 3.0 for k, v in acc_b.items() if v != 0}
P("     gives %s  instead of {2, 4, 7}                  %s"
  % (sorted(bad.values()), "control fires" if bad != WANT else "*** INERT ***"))
assert bad != WANT
REC["steps"]["seven_rebuilt"] = {
    "coefficients": {"cos2": 2, "cos1": 4, "alt": 7}, "split": {"A_mu": 4, "A_56": 3},
    "ghosts": False, "verdict": "PASS"}

# ----------------------------------------------------------------------------
P("")
P("=" * 78)
P("B -- THE EXACT BRACE IN CLOSED FORM, AND THE SIZE OF THE k >> 1 ERROR")
P("=" * 78)


def brace_exact_sum(n, k, p6, MMAX=200000):
    m = np.arange(1, MMAX + 1, dtype=float)
    s = (4.0 / ((math.pi * n) ** 2 + (math.pi * m / k) ** 2) ** 3).sum()
    return s + 2.0 / (math.pi * n) ** 6 + 3.0 * k * p6 / (4.0 * (math.pi * n) ** 5)


def g(x):
    return (1.0 / math.tanh(x) + x / math.sinh(x) ** 2
            + (2.0 / 3.0) * x * x * (1.0 / math.tanh(x)) / math.sinh(x) ** 2)


def brace_closed(n, k, p6):
    return 3.0 * k / (4.0 * (math.pi * n) ** 5) * (p6 + g(math.pi * k * n))


P("  brace(n,k,P6) = 3k/(4 (pi n)^5) * [P6 + g(pi k n)],")
P("     g(x) = coth x + x csch^2 x + (2/3) x^2 coth x csch^2 x")
P("")
P("  the P6 = -1 brace is a cancellation of two O(1) terms down to O(e^-2 pi k n),")
P("  so it is compared against the P6 = +1 SCALE at the same (n,k), not against")
P("  itself -- below ~1e-16 of that scale the direct sum is pure double-precision")
P("  noise (it even goes negative) and only the closed form is meaningful.")
P("")
P("  %-4s %-6s %-4s %-14s %-14s %s"
  % ("n", "k", "P6", "direct sum", "closed form", "|diff|/scale(n,k)"))
worst = 0.0
for k in (1.2, 3.0, 10.0):
    for n in (1, 2, 5):
        scale = brace_exact_sum(n, k, +1)
        for p6 in (+1, -1):
            a = brace_exact_sum(n, k, p6)
            b = brace_closed(n, k, p6)
            rel = abs(a - b) / abs(scale)
            worst = max(worst, rel)
            P("  %-4d %-6.1f %+4d %-14.6e %-14.6e %.1e" % (n, k, p6, a, b, rel))
P("  worst disagreement, scaled: %.1e                                    %s"
  % (worst, "PASS" if worst < 1e-13 else "*** FAIL ***"))
assert worst < 1e-13
P("")
P("  size of the k >> 1 truncation, g(pi k n) - 1, at n = 1:")
for k in (0.3, 0.5, 0.8, 1.2, 2.0, 3.0, 5.0, 10.0):
    P("     k = %-5.1f  g(pi k) - 1 = %+.4e   (%7.3f %% of the P6=+1 brace)"
      % (k, g(math.pi * k) - 1, 100 * (g(math.pi * k) - 1) / 2))
P("")
P("  >> RETRACTION.  SU7_VACUUM.md (line 50 and the repair table), EVIDENCE.md")
P("     (two rows) and STATE_2026-08-02.md (two places) all record 'the exact")
P("     eqs. (62)-(65) move everything the right way but need k ~ 1.2'.  Five")
P("     copies; the memory never carried it. [[no-claim-lives-in-one-place]]")
P("     THAT NUMBER WAS NEVER COMPUTED.  su7_anchor_hunt.py line 167 prints the")
P("     literal string '~1.2 for all -- against their own k >> 1' on every row;")
P("     its loop fits only lambda and c, and `res` has no 'k' key at all.  A")
P("     conclusion line written without the data. [[the-conclusion-line-is-a-prediction]]")
P("")
P("     The computed answer, above, is different AND it excludes the repair more")
P("     strongly.  The correction is exponentially small in k: 0.02 % at k = 2,")
P("     1.5 % at k = 1.2.  It only becomes O(1) for k < ~0.5, i.e. 1/R6 < 1/R5 --")
P("     the OPPOSITE regime to theirs, not merely a milder version of it.  So no")
P("     value of k consistent with a 6D model of this shape can close a 30-100 %")
P("     gap in a_min, and the honest verdict is not 'needs an unphysical k' but")
P("     'the k >> 1 truncation is not where the residual lives'.")
REC["steps"]["exact_brace"] = {
    "closed_form": "3k/(4(pi n)^5) * [P6 + g(pi k n)]",
    "worst_rel": worst,
    "g_minus_1": {str(k): g(math.pi * k) - 1 for k in (1.2, 2.0, 3.0, 5.0, 10.0)}}

# ----------------------------------------------------------------------------
P("")
P("=" * 78)
P("C -- THE m_h TEST.  K = m_h * a / sqrt(F''(a)) MUST BE ROW-INDEPENDENT")
P("=" * 78)
P("  from their eqs. (80) + (82), with V = C F and C = 3k/(64 pi^8 R5^6):")
P("     m_h^2 = 3 g4^2 F''(a_min) / (16 pi^6 R5^2),   1/R5 = 2*80.4/a_min")
P("  so K = m_h a_min / sqrt(F''(a_min)) = 2*80.4*g4*sqrt(3/(16 pi^6)) for EVERY row.")
P("  K is invariant under F -> lambda F, so it does not test the normalisation.")
P("")
CONST = 2 * MW * math.sqrt(3.0 / (16.0 * math.pi ** 6))
P("  2*80.4*sqrt(3/(16 pi^6)) = %.6f, so K = %.6f * g4" % (CONST, CONST))
P("")
P("  --- evaluated at THEIR published a_min ---")
P("  %-5s %-9s %-9s %-13s %-11s %-11s %s"
  % ("case", "a theirs", "m_h", "F''(a)", "K", "implied g4", "F'(a)/F''(a)"))
Ks, rows = [], []
for tag, cont, a_t, mh_t, iR_t in T1:
    f2 = float(V(cont, np.array([a_t]), d=2)[0])
    f1 = float(V(cont, np.array([a_t]), d=1)[0])
    K = mh_t * a_t / math.sqrt(f2) if f2 > 0 else float("nan")
    Ks.append(K)
    rows.append(dict(case=tag, a_theirs=a_t, mh=mh_t, F2=f2, F1=f1, K=K,
                     g4=K / CONST))
    P("  %-5s %-9.4f %-9.1f %-13.5e %-11.5f %-11.5f %+.4f"
      % (tag, a_t, mh_t, f2, K, K / CONST, f1 / f2))
spread = (max(Ks) - min(Ks)) / (sum(Ks) / len(Ks))
P("")
P("  K spread (max-min)/mean = %.4f  (%.2f %%)" % (spread, 100 * spread))
P("  mean implied g4 = %.4f   (SM SU(2)_L coupling at the weak scale ~ 0.65)"
  % (sum(Ks) / len(Ks) / CONST))
P("")
P("  --- the same quantity at OUR a_min, for contrast ---")
P("  %-5s %-9s %-13s %-11s %s" % ("case", "a ours", "F''(a)", "K", "implied g4"))
Ko = []
for tag, cont, a_t, mh_t, iR_t in T1:
    a_o = minimise(cont)
    f2 = float(V(cont, np.array([a_o]), d=2)[0])
    K = mh_t * a_o / math.sqrt(f2) if f2 > 0 else float("nan")
    Ko.append(K)
    P("  %-5s %-9.4f %-13.5e %-11.5f %.5f" % (tag, a_o, f2, K, K / CONST))
spread_o = (max(Ko) - min(Ko)) / (sum(Ko) / len(Ko))
P("  K spread at our minima = %.4f  (%.2f %%)" % (spread_o, 100 * spread_o))
REC["steps"]["mh_test"] = {"const": CONST, "at_theirs": rows,
                          "spread_theirs": spread, "spread_ours": spread_o,
                          "g4_mean": sum(Ks) / len(Ks) / CONST}

P("")
P("  --- and the third, independent read: is a_theirs a stationary point of F? ---")
P("  a true minimum has F'(a) = 0.  Scale F' by the curvature to get a length:")
P("  d = -F'(a)/F''(a) is the Newton step from a_theirs to the nearby extremum.")
for tag, cont, a_t, mh_t, iR_t in T1:
    f1 = float(V(cont, np.array([a_t]), d=1)[0])
    f2 = float(V(cont, np.array([a_t]), d=2)[0])
    P("  %-5s  a = %.4f   Newton step %+.5f   ->  %.5f   (our min %.5f)"
      % (tag, a_t, -f1 / f2, a_t - f1 / f2, minimise(cont)))

P("")
P("=" * 78)
P("D -- FALSIFYING THE TEST ITSELF, THEN INVERTING IT")
P("=" * 78)
P("  K = m_h * R(a) with R(a) = a/sqrt(F''(a)).  If R were nearly constant over")
P("  the relevant window, K would be ~m_h for ANY a and the agreement of rows")
P("  (1),(2),(5) would be vacuous.  So: how much does R move?")
P("")
P("  %-5s %-11s %-11s %-11s %-11s %-8s %s"
  % ("case", "R(a_theirs)", "R(0.02)", "R(0.10)", "R(0.30)", "[.02,.30]", "[.04,.12]"))
SCAN = np.linspace(0.02, 0.30, 400)
for tag, cont, a_t, mh_t, iR_t in T1:
    f2s = V(cont, SCAN, d=2)
    ok_m = f2s > 0
    Rs = np.where(ok_m, SCAN / np.sqrt(np.abs(f2s)), np.nan)
    f2t = float(V(cont, np.array([a_t]), d=2)[0])
    Rt = a_t / math.sqrt(f2t) if f2t > 0 else float("nan")

    def Rat(x):
        v = float(V(cont, np.array([x]), d=2)[0])
        return x / math.sqrt(v) if v > 0 else float("nan")
    fin = Rs[np.isfinite(Rs)]
    sub = SCAN[(SCAN >= 0.04) & (SCAN <= 0.12)]
    f2b = V(cont, sub, d=2)
    Rb = sub[f2b > 0] / np.sqrt(f2b[f2b > 0])
    P("  %-5s %-11.5f %-11.5f %-11.5f %-11.5f x%-7.1f x%.2f"
      % (tag, Rt, Rat(0.02), Rat(0.10), Rat(0.30),
         (fin.max() / fin.min()) if len(fin) else float("nan"),
         (Rb.max() / Rb.min()) if len(Rb) else float("nan")))
P("")
P("  >> READ THE LAST COLUMN, NOT THE ONE BEFORE IT.  Over the FULL window R does")
P("     move by x4-x130, but that is driven by the ends, where F'' -> 0 and R blows")
P("     up.  Over the physically relevant sub-window a in [0.04, 0.12] R is FLAT to")
P("     a few per cent.  So for a row whose a_min lands in that sub-window, K ~ m_h")
P("     x 0.011 almost automatically, and the agreement of (1),(2),(5) is WEAK")
P("     evidence on its own: it says only that their a sits in the flat region.")
P("     What is NOT weak is the failure mode.  Rows (3) and (4) do not sit there:")
P("     row (4) has R = 0.0333, three times the plateau, and row (3) has F'' < 0,")
P("     so their eq. (80) returns NO REAL m_h AT ALL at their own published a_min.")
P("     That statement is independent of how flat R is.")
P("")
P("  INVERSION -- the m_h column determines a_min WITHOUT minimising anything.")
P("  Solve  m_h * a/sqrt(F''(a)) = 2*80.4*g4*sqrt(3/(16 pi^6))  for a.")
P("  This is a SECOND, independent determination of a_min from their own table.")
P("")
inv = []
for tag, cont, a_t, mh_t, iR_t in T1:
    row = {"case": tag, "a_theirs": a_t, "a_min_ours": minimise(cont), "roots": {}}
    P("  %s   a theirs = %.4f   a from minimising F = %.4f"
      % (tag, a_t, row["a_min_ours"]))
    for g4 in (0.60, 0.63, 0.66):
        tgt = CONST * g4 / mh_t
        xs = np.linspace(1e-3, 1.0, 20001)
        f2s = V(cont, xs, d=2)
        Rs = np.where(f2s > 0, xs / np.sqrt(np.abs(f2s)), np.nan)
        d = Rs - tgt
        lo, hi = d[:-1], d[1:]
        j = np.nonzero(np.isfinite(lo) & np.isfinite(hi) & (lo * hi < 0))[0]
        roots = list(xs[j] - lo[j] * (xs[j + 1] - xs[j]) / (hi[j] - lo[j]))
        row["roots"]["g4=%.2f" % g4] = roots
        P("      g4 = %.2f -> a = %s" % (g4, ", ".join("%.4f" % r for r in roots)
                                         if roots else "NO SOLUTION"))
    inv.append(row)
REC["steps"]["inversion"] = inv
P("")
P("  >> read the three 48-free rows (1),(2),(5) against the two 48-bearing rows")
P("     (3),(4).  The only structural difference between the two groups is the 48,")
P("     i.e. their eq. (75).")

P("")
P("=" * 78)
P("E -- THE DECIDING TEST: does one g4 survive the per-row fermion rescaling?")
P("=" * 78)
P("  SU7_VACUUM.md rejected the fermion rescaling lam because the five fitted")
P("  values (0.85, 0.96, 0.81, 0.80, 0.85) are not one number.  That used only")
P("  the a_min column.  With m_h there is a SECOND equation per row and only ONE")
P("  new unknown in total (g4), so five rows give 10 equations for 6 unknowns.")
P("")
P("  For each row: solve argmin F(lam) = a_theirs for lam, then read g4 off m_h.")
P("  If the five g4 agree -- and agree with the SM value -- then lam is not a fudge")
P("  and the residual is a real, content-dependent factor multiplying the fermion")
P("  sector.  If they scatter, the SHAPE is wrong and no rescaling can fix it.")
P("")


def lam_for(cont, a_t):
    f = lambda L: minimise(cont, lam=L) - a_t
    lo, hi = 0.05, 5.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


P("  %-5s %-9s %-9s %-9s %-13s %-11s %s"
  % ("case", "a theirs", "lam", "check a", "F''(a;lam)", "K", "implied g4"))
g4s, erows = [], []
for tag, cont, a_t, mh_t, iR_t in T1:
    lam = lam_for(cont, a_t)
    if lam is None:
        P("  %-5s %-9.4f  no lam in [0.05,5] reproduces a_theirs" % (tag, a_t))
        erows.append({"case": tag, "lam": None})
        continue
    chk = minimise(cont, lam=lam)
    f2 = float(V(cont, np.array([a_t]), d=2, lam=lam)[0])
    K = mh_t * a_t / math.sqrt(f2) if f2 > 0 else float("nan")
    g4 = K / CONST
    g4s.append(g4)
    erows.append({"case": tag, "lam": lam, "a_check": chk, "F2": f2, "K": K, "g4": g4})
    P("  %-5s %-9.4f %-9.5f %-9.5f %-13.5e %-11.5f %.5f"
      % (tag, a_t, lam, chk, f2, K, g4))
P("")
if g4s:
    mean = sum(g4s) / len(g4s)
    sprd = (max(g4s) - min(g4s)) / mean
    P("  five implied g4: %s" % ", ".join("%.4f" % x for x in g4s))
    P("  mean = %.4f, spread = %.2f %%" % (mean, 100 * sprd))
    P("  SM SU(2)_L coupling: g2(M_Z) = 0.652, and it runs DOWN with scale;")
    P("  at 2-8 TeV g2 ~ 0.63.  1/R5 for these rows is 2.0-7.5 TeV.")
    REC["steps"]["lam_plus_mh"] = {"rows": erows, "g4_mean": mean, "g4_spread": sprd}
P("")
P("  CONTROL -- the same fit with the m_h column SCRAMBLED (rows reversed).")
P("  A test that cannot fail would give the same spread either way.")
g4b = []
for (tag, cont, a_t, _mh, _iR), (_t, _c, _a, mh_s, _i) in zip(T1, T1[::-1]):
    lam = lam_for(cont, a_t)
    if lam is None:
        continue
    f2 = float(V(cont, np.array([a_t]), d=2, lam=lam)[0])
    if f2 > 0:
        g4b.append(mh_s * a_t / math.sqrt(f2) / CONST)
if g4b:
    mb = sum(g4b) / len(g4b)
    P("  scrambled g4: %s" % ", ".join("%.4f" % x for x in g4b))
    P("  scrambled spread = %.2f %%" % (100 * (max(g4b) - min(g4b)) / mb))
    P("  (m_h varies only 125.1-126.8 GeV, so this control is WEAK by construction")
    P("   -- it bounds how much of any agreement is carried by m_h at all.)")

os.makedirs("outputs", exist_ok=True)
os.makedirs("paper_data", exist_ok=True)
with open("paper_data/su7_anchor_mh.json", "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("wrote paper_data/su7_anchor_mh.json")
