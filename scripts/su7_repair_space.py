#!/usr/bin/env python3
"""The anchor residual as a SPACE, and what it can and cannot reach.

Carles Marin + Claude (AI assistant).  2026-08-04.
Source: Komori & Maru, arXiv:2503.04090.
Follows ANCHOR_SECTION_31.md s6 and SU7_VACUUM.md.

STATE COMING IN.
  Their Table 1 a_min is not reproduced from their eqs. (68),(72)-(76).  The
  gauge sector is exact (s1), the k >> 1 truncation is irrelevant (s2), the
  fermion sector is transcribed correctly three independent ways (s6), the
  structure constants, the KK spectrum and all four branchings are derived
  (s7, s8).  A global or per-row rescaling is excluded; per-representation and
  per-channel reweightings are rejected by the two columns together.

THE QUESTION THIS SCRIPT ASKS, AND IT IS NOT THE SAME QUESTION.
  Not "what is the missing ingredient" -- s6 bounded that and could not name it.
  It is: **the headline of this whole line rests on D, computed from the same
  formulas that fail the anchor.  Does it survive the repair?**  Instead of one
  repair we take the WHOLE FAMILY of repairs and ask over which of them the
  verdict "case (2) is the unique row" is invariant.  A conclusion that holds on
  a region is a different object from one that holds at a point.

  PART A  -- the headline as an exact region in repair space.
  PART B  -- and one more negative for the anchor itself: no single extra
             channel s^n cos(c pi n a)/n^p repairs both columns either, and the
             c it would need is measured rather than argued.
"""
import json
import math
import os
from fractions import Fraction as Fr

import numpy as np

P = lambda *a: print(*a, flush=True)
NMAX = 600
MW = 80.4
_n = np.arange(1, NMAX + 1)
_nf = _n.astype(float)
_sgn = {1: np.ones(NMAX), -1: (-1.0) ** _n}


_BC = {}


def basis(alpha, s, c, d=0, p=5):
    """d-th derivative of  sum_n s^n cos(c pi n a) / n^p   at alpha."""
    key = (round(float(alpha), 12), s, c, d, p)
    if key in _BC:
        return _BC[key]
    a = np.atleast_1d(float(alpha))
    ph = np.outer(a, c * math.pi * _nf)
    f = (np.cos, lambda t: -np.sin(t), lambda t: -np.cos(t),
         lambda t: np.sin(t))[d % 4](ph)
    out = (f * ((_nf ** -p) * _sgn[s] * (c * math.pi * _nf) ** d)).sum(axis=1)
    _BC[key] = out
    return out


GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]
REPS = ("7", "28", "48", "84")


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


T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.8),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)], 0.081, 125.5),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], 0.021, 125.1),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], 0.026, 126.4),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.2)]

REC = {"source": "Komori & Maru, arXiv:2503.04090", "steps": {}}

# --------------------------------------------------------------- entry control
P("=" * 78)
P("ENTRY CONTROL -- terms() against the actual Wilson line (eq. 77)")
P("=" * 78)
WILSON = {"7":  {(1, -1): 2},
          "28": {(2, +1): 2, (1, +1): 2, (1, -1): 8},
          "48": {(2, +1): 2, (1, +1): 4, (1, -1): 16},
          "84": {(3, -1): 2, (2, +1): 8, (2, -1): 2, (1, +1): 8, (1, -1): 24}}
ok = True
for rep in REPS:
    want = {k: v // 2 for k, v in WILSON[rep].items()}
    got = {}
    for m, s, c in terms(rep, 1, 1):
        got[(c, s)] = got.get((c, s), 0) + m
    ok &= (got == want)
P("  all four representations reproduce the Wilson-line count:  %s"
  % ("PASS" if ok else "*** FAIL ***"))
assert ok
P("  (the same instrument as su7_content_dependence.py, re-asserted here so this")
P("   script cannot silently drift from it.)")

# ============================================================== PART A
P("")
P("=" * 78)
P("PART A -- THE HEADLINE AS A REGION, NOT AS A POINT")
P("=" * 78)
P("  The headline is: case (2) is the unique row of their Table 1 that both")
P("  supplies the nu_R its own anomaly cancellation demands AND survives")
P("  donating an 84(+,+) to host a third lepton generation.")
P("")
P("  It has two halves and they are NOT equally exposed:")
P("    (i)  which rows can supply a U(1)' = -1 SM singlet: rows (2) and (3).")
P("         Pure group theory on their eqs. (69),(70) + parities -- it is a")
P("         statement about WHICH COMPONENTS EXIST, and no repair of the")
P("         POTENTIAL can touch it.  (su7_anomaly_channels.py; outputs/.)")
P("    (ii) of those two, (3) loses EWSB on donation and (2) does not.")
P("         This comes from D, computed from the very eqs. (68),(72)-(76) that")
P("         fail the anchor.  THIS is the exposed half, and it is what follows.")
P("")

D_TERM = {}          # per (rep, s) exact rational D contribution


def Dof_w(content, w):
    """D with a per-representation weight w[R].  Exact rational."""
    tot = Fr(0)
    for m, s, c in GAUGE:
        tot += Fr(m).limit_denominator() * c ** 2 * (1 if s == 1 else Fr(-3, 4))
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            tot += w[rep] * m * mult * c ** 2 * (1 if s == 1 else Fr(-3, 4))
    return tot


ONE = {R: Fr(1) for R in REPS}
DG = Dof_w([], ONE)
P("  gauge sector          D = %s" % DG)
for lab, cont in (("7(+,+)", [("7", 1, 1, 1)]), ("7(+,-)", [("7", 1, -1, 1)]),
                  ("28(+,+)", [("28", 1, 1, 1)]), ("28(+,-)", [("28", 1, -1, 1)]),
                  ("48(+,+)", [("48", 1, 1, 1)]), ("84(+,+)", [("84", 1, 1, 1)])):
    D_TERM[lab] = Dof_w(cont, ONE) - DG
    P("  %-21s D = %-6s   (per multiplet)" % (lab, D_TERM[lab]))

P("")
P("  D IS EXACTLY LINEAR IN THE WEIGHTS, and the 48's coefficient is EXACTLY 0.")
P("  Proof by the arithmetic itself, not by evaluation: the 48(+,+) contributes")
P("  m c^2 = 1*4 + 2*1 = 6 at s = +1 and 8*1 = 8 at s = -1, and 6 - (3/4)*8 = 0.")
P("  So w(48) multiplies zero -- for EVERY value of w(48).")
chk = [Dof_w([("48", 1, 1, 1)], {**ONE, "48": Fr(q)}) - DG for q in (0, 1, 5, 100)]
P("     control: D(48) at w(48) = 0, 1, 5, 100  ->  %s   %s"
  % (chk, "PASS" if set(chk) == {Fr(0)} else "*** FAIL ***"))
assert set(chk) == {Fr(0)}
P("")
P("  >> THE SINGLE LARGEST REPAIR THE a_min COLUMN ASKS FOR IS w(48) = 5.59")
P("     (ANCHOR_SECTION_31.md s6, system (I), residual 0.0064).  IT IS INVISIBLE")
P("     TO THE HEADLINE.  Not approximately -- identically.")
P("")

P("  The two exposed inequalities, in exact eighths.  After donating one 84:")
P("     row (3):  D = -27/8 + 2 w(28) + (5/4) w(84)          must be < 0")
P("     row (2):  D = -27/8 + 2 w(28) + (15/4) w(84)         must be > 0")
P("  w(7) and w(48) do not appear in either.  The headline lives in a PLANE.")
P("")


def wedge(w28, w84):
    w = {**ONE, "28": Fr(w28).limit_denominator(10 ** 6),
         "84": Fr(w84).limit_denominator(10 ** 6)}
    out = {}
    for tag, cont, *_ in T1:
        red = [(r, e, ep, m - 1 if (r, e, ep) == ("84", 1, 1) else m)
               for r, e, ep, m in cont]
        red = [x for x in red if x[3]]
        out[tag] = (Dof_w(cont, w), Dof_w(red, w))
    return out


P("  The wedge, sliced along w(28) = w(84) = w -- which must reproduce the")
P("  published lam_crit of SU7_VACUUM.md, and that is a control that can fail:")
lo = Fr(27, 8) / (Fr(2) + Fr(15, 4))    # (2)  alive:  (2 + 15/4) w > 27/8
hi = Fr(27, 8) / (Fr(2) + Fr(5, 4))     # (3)  dead :  (2 +  5/4) w < 27/8
P("     (2) survives for w > %s = %.4f      [SU7_VACUUM.md lam_crit(2) = 0.587]"
  % (lo, float(lo)))
P("     (3) dies     for w < %s = %.4f      [SU7_VACUUM.md lam_crit(3) = 1.038]"
  % (hi, float(hi)))
okc = abs(float(lo) - 0.587) < 5e-4 and abs(float(hi) - 1.038) < 5e-4
P("     control against the published values: %s"
  % ("PASS" if okc else "*** FAIL ***"))
assert okc
P("     >> the headline holds on the whole interval w in (%.3f, %.3f) -- a factor"
  % (float(lo), float(hi)))
P("        1.77 wide, and it CONTAINS w = 1.")
P("")

P("  And the same wedge in the full (w28, w84) plane, evaluated at every repair")
P("  this programme has actually fitted:")
FITTED = [("their formulas, w = 1", 1.0, 1.0),
          ("(I) a_min only        s6", 1.081, 0.827),
          ("(II) curvature only   s6", 1.037, 1.030),
          ("(I)+(II) both columns s6", 0.825, 0.923),
          ("per-row lam, lowest   s4", 0.800, 0.800),
          ("per-row lam, highest  s4", 0.962, 0.962),
          ("EDGE: w at which (2) dies", float(lo), float(lo))]
P("  %-26s %-8s %-8s %-11s %-11s %s"
  % ("repair", "w(28)", "w(84)", "D(2) don.", "D(3) don.", "headline"))
AREC = []
for lab, a, b in FITTED:
    d = wedge(a, b)
    d2, d3 = float(d["(2)"][1]), float(d["(3)"][1])
    keep = d2 > 0 and d3 < 0
    P("  %-26s %-8.3f %-8.3f %-11.4f %-11.4f %s"
      % (lab, a, b, d2, d3, "HOLDS" if keep else "*** BROKEN ***"))
    AREC.append({"repair": lab, "w28": a, "w84": b, "D2_donated": d2,
                 "D3_donated": d3, "headline": keep})
NUR = {"(1)": False, "(2)": True, "(3)": True, "(4)": False, "(5)": False}
PERROW = []
for tag, cont, a0, mh0 in T1:
    d = wedge(1.0, 1.0)[tag]
    PERROW.append({"case": tag, "a_theirs": a0, "mh_theirs": mh0,
                   "D_theirs": float(d[0]), "D_donated": float(d[1]),
                   "nu_R": NUR[tag],
                   "headline": bool(NUR[tag] and d[1] > 0)})
REC["steps"]["headline_region"] = {"w_diagonal_interval": [float(lo), float(hi)],
                                   "fitted": AREC, "per_row": PERROW,
                                   "coef": {"D2_donated": [-27 / 8, 2.0, 3.75],
                                            "D3_donated": [-27 / 8, 2.0, 1.25]}}
P("")
P("  ORDERING FACT, and it is why there are three regions and not four:")
P("     D(2) - D(3) = (15/4 - 5/4) w(84) = (5/2) w(84) > 0 for any w(84) > 0.")
P("     So row (3) can NEVER survive the donation while row (2) dies.  The plane")
P("     has exactly three regions: both alive, the headline wedge, both dead.")
P("")
P("  >> every repair the anchor data has ever asked for keeps the headline.")
P("     The one that comes closest to breaking it is the curvature-only fit")
P("     (w(84) = 1.030 against the ceiling 1.038 at w(28) = 1.037) -- and that")
P("     system is exactly determined, so it says nothing.  Recorded anyway.")

# ---- the per-channel repair, killed by an independent instrument
P("")
P("-" * 78)
P("AND THE PER-CHANNEL REPAIR IS KILLED A THIRD TIME, BY THE VACUUM")
P("-" * 78)
P("  s6 rejected it twice: its control fires (it fits scrambled data equally")
P("  well, ratio 0.97) and it needs g4 = 2.39.  Both are arguments about the")
P("  FIT.  D is an independent instrument -- it never sees a_min or m_h at all.")
P("")
CH = [(1, 1), (1, -1), (2, 1), (2, -1), (3, -1)]
CHW = dict(zip(CH, (0.13, 0.28, 0.20, 0.20, 0.20)))
P("  The channel weights s6 measured are all in 0.13-0.28 where eq. (72) says 1.")
P("  Take the whole band [0.13, 0.28] and ask what D their OWN UNMODIFIED rows")
P("  would have -- rows that, by their own Table 1, DO break electroweak")
P("  symmetry.  D <= 0 means a = 0 is a minimum and there is no EWSB at all.")
P("")


def D_channel(content, x):
    tot = Fr(0)
    for m, s, c in GAUGE:
        tot += Fr(m).limit_denominator() * c ** 2 * (1 if s == 1 else Fr(-3, 4))
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            tot += (Fr(x).limit_denominator(10 ** 6) * m * mult * c ** 2
                    * (1 if s == 1 else Fr(-3, 4)))
    return tot


P("  %-6s %-14s %-14s %-14s %s"
  % ("case", "D at w=1", "D at x=0.28", "D at x=0.13", "verdict at x <= 0.28"))
chan_rows = []
for tag, cont, *_ in T1:
    d1, dhi, dlo = (float(Dof_w(cont, ONE)), float(D_channel(cont, 0.28)),
                    float(D_channel(cont, 0.13)))
    chan_rows.append({"case": tag, "D_w1": d1, "D_028": dhi, "D_013": dlo})
    P("  %-6s %-14.4f %-14.4f %-14.4f %s"
      % (tag, d1, dhi, dlo, "EWSB" if dhi > 0 else "*** NO EWSB ***"))
dead = all(r["D_028"] <= 0 for r in chan_rows)
P("")
P("  >> at the largest weight the channel fit allows, %s of the five rows lose"
  % ("ALL" if dead else "some"))
P("     electroweak symmetry breaking ALTOGETHER -- including the rows whose")
P("     published a_min and m_h the fit was constructed to reproduce.  A repair")
P("     that reproduces two numbers of a row by destroying the phenomenon that")
P("     row exists to exhibit is not a repair.")
P("     The instrument is independent: D never sees a_min or m_h.")
REC["steps"]["per_channel_vacuum"] = {"rows": chan_rows, "all_dead": bool(dead)}

# ============================================================== PART B
P("")
P("=" * 78)
P("PART B -- ONE MORE NEGATIVE FOR THE ANCHOR: NO SINGLE EXTRA CHANNEL EITHER")
P("=" * 78)
P("  s6 ended: whatever is missing changes the STRUCTURE -- a cosine argument")
P("  outside {1,2,3} x pi n a, an n-dependence that is not 1/n^5, or multiplets")
P("  not in (69),(70).  All three are the same object: ONE extra term")
P("      delta(a) = v * sum_n s^n cos(c pi n a) / n^p")
P("  with (c, p, s) free.  Add it and the system stays LINEAR, so it is a scan")
P("  over (c, p, s) with a least-squares solve at each point -- no optimiser.")
P("")


def rowdata(alphas=None):
    out = []
    for i, (tag, cont, a0, mh0) in enumerate(T1):
        a_t = a0 if alphas is None else alphas[i]
        r = {"tag": tag, "a": a_t, "mh": mh0, "cont": cont}
        for d in (1, 2):
            r["G%d" % d] = float(sum(m * basis(a_t, s, c, d)[0] for m, s, c in GAUGE))
            r["F%d" % d] = r["G%d" % d] + float(sum(
                m * mult * basis(a_t, s, c, d)[0]
                for rep, e_, ep, mult in cont for m, s, c in terms(rep, e_, ep)))
        r["c"] = (mh0 * a_t / (2 * MW)) ** 2 * 16 * math.pi ** 6 / 3
        out.append(r)
    return out


def lstsq(A, b):
    A, b = np.array(A, float), np.array(b, float)
    nrm = np.linalg.norm(np.c_[A, b], axis=1)
    nrm[nrm == 0] = 1.0
    x, *_ = np.linalg.lstsq(A / nrm[:, None], b / nrm, rcond=None)
    return x, float(np.linalg.norm((A / nrm[:, None]) @ x - b / nrm))


# ---- B0: what the correction must DO, measured
P("-" * 78)
P("B0 -- THE DEMAND, MEASURED BEFORE ANYTHING IS FITTED")
P("-" * 78)
P("  At their published a_i the missing piece must supply")
P("      delta'(a_i)  = -F'(a_i)                    (make a_i stationary)")
P("      delta''(a_i) = c_i/g4^2 - F''(a_i)         (and give their m_h)")
P("  with g4 the 4D SU(2)_L coupling.  Taking the SM value g4 = 0.6300:")
P("")
G4 = 0.6300
RD = rowdata()
P("  %-6s %-8s %-12s %-12s %-12s %-12s %s"
  % ("case", "a", "F'(a)", "demand d'", "F''(a)", "demand d''", "|d''| a / |d'|"))
B0 = []
for r in RD:
    d1 = -r["F1"]
    d2 = r["c"] / G4 ** 2 - r["F2"]
    ratio = abs(d2) * r["a"] / abs(d1) if d1 else float("inf")
    B0.append({"case": r["tag"], "a": r["a"], "F1": r["F1"], "F2": r["F2"],
               "demand1": d1, "demand2": d2, "ratio": ratio})
    P("  %-6s %-8.3f %-12.4f %-12.4f %-12.4f %-12.4f %.4f"
      % (r["tag"], r["a"], r["F1"], d1, r["F2"], d2, ratio))
P("")
_rr = [b["ratio"] for b in B0]
P("  READ THE LAST COLUMN, AND READ IT AS AN EXCLUSION AND NOTHING MORE.")
P("  Any term  v sum_n s^n cos(c pi n a)/n^p  whose argument c pi n a is small")
P("  at n = 1 obeys  delta' = a * delta''  IDENTICALLY (expand the cosine: both")
P("  are -v (c pi)^2 sum_n s^n/n^(p-2), one of them times a).  Such a term")
P("  therefore enters that column at 1 on EVERY row, whatever v, c, p, s are.")
P("  The measured demands range over a factor %.1f, from %.2f to %.2f."
  % (max(_rr) / min(_rr), min(_rr), max(_rr)))
P("  >> NO SINGLE SMALL-ARGUMENT TERM SUPPLIES ALL FIVE DEMANDS.  That is all")
P("     this column shows.  It does NOT say what c is -- the scan measures that.")
REC["steps"]["demand"] = B0

# ---- B1 / B2: the scan
P("")
P("-" * 78)
P("B1/B2 -- THE SCAN")
P("-" * 78)
P("  B1: one universal extra term, the same in every row.  Unknowns (v, 1/g4^2):")
P("      2 unknowns, 10 equations.")
P("  B2: the extra term carried per representation, multiplicity v_R and sign")
P("      tracking each multiplet's own eta*eta'.  Unknowns (v_7,v_28,v_48,v_84,")
P("      1/g4^2): 5 unknowns, 10 equations.")
P("")
P("  THE CONTROL, AND IT IS THE WHOLE TEST.  Scanning (c, p, s) is itself a")
P("  search, so a control that fixes (c,p,s) at the real winner and re-solves on")
P("  scrambled data measures nothing.  The control here re-runs THE ENTIRE SCAN")
P("  on each of the 119 non-trivial permutations of their a_min column and takes")
P("  the best residual of each.  Best-over-search against best-over-search is")
P("  the only comparison that can fail. [[a-control-that-cannot-fail]]")
P("")

CS = [round(0.5 * i, 1) for i in range(1, 81)]        # 0.5 .. 40.0
PS = [3, 4, 5, 6, 7]
SS = [1, -1]


def solve_B1(rows, c, p, s):
    A, b = [], []
    for r in rows:
        A.append([basis(r["a"], s, c, 1, p)[0], 0.0]);  b.append(-r["F1"])
        A.append([basis(r["a"], s, c, 2, p)[0], -r["c"]]);  b.append(-r["F2"])
    return lstsq(A, b)


def solve_B2(rows, c, p, sgn):
    A, b = [], []
    for r in rows:
        for d in (1, 2):
            coef = []
            for R in REPS:
                v = 0.0
                for rep, e_, ep, mult in r["cont"]:
                    if rep == R:
                        v += mult * basis(r["a"], sgn * e_ * ep, c, d, p)[0]
                coef.append(v)
            A.append(coef + ([0.0] if d == 1 else [-r["c"]]))
            b.append(-r["F%d" % d])
    return lstsq(A, b)


BASE, _ = None, None
A0, b0 = [], []
for r in RD:
    A0.append([0.0]);  b0.append(-r["F1"])
    A0.append([-r["c"]]);  b0.append(-r["F2"])
BASE = lstsq(A0, b0)[1]
P("  baseline: no extra term at all, only 1/g4^2 free   ||res|| = %.4f" % BASE)
P("")
from itertools import permutations
A_REAL = [r["a"] for r in RD]
PERMS = [list(q) for q in permutations(A_REAL) if list(q) != A_REAL]


def scan(fn, rows):
    best = None
    for c in CS:
        for p in PS:
            for s in SS:
                x, res = fn(rows, c, p, s)
                if best is None or res < best[0]:
                    best = (res, c, p, s, x)
    return best


WIN = {}
for name, fn in (("B1 universal", solve_B1), ("B2 per-rep", solve_B2)):
    res, c, p, s, x = scan(fn, RD)
    WIN[name] = (res, c, p, s, x)
    u = x[-1]
    g4 = 1 / math.sqrt(u) if u > 0 else float("nan")
    P("  %-14s best (c, p, s) = (%.1f, %d, %+d)   ||res|| = %.5f   g4 = %s"
      % (name, c, p, s, res, ("%.3f" % g4) if u > 0 else "IMAGINARY"))
    P("  %-14s weights: %s   (eq. (72) says the fermion terms carry 1)"
      % ("", ["%.4f" % v for v in x[:-1]]))
    scr = sorted(scan(fn, rowdata(alphas=q))[0] for q in PERMS)
    nbet = sum(1 for v in scr if v <= res)
    P("  %-14s CONTROL, full scan re-run on all %d scrambles of their a column:"
      % ("", len(PERMS)))
    P("  %-14s   best %.5f (%.2f x real)   median %.5f   worst %.5f"
      % ("", scr[0], scr[0] / res, scr[len(scr) // 2], scr[-1]))
    P("  %-14s   scrambles reaching the real residual: %d of %d"
      % ("", nbet, len(scr)))
    verdict = ("A REAL DISCRIMINATION" if scr[0] > 3 * res else
               "*** NOT A DISCRIMINATION -- the family fits anything ***")
    P("  %-14s -> %s" % ("", verdict))
    REC["steps"].setdefault("scan", {})[name] = {
        "c": c, "p": p, "s": s, "resnorm": res, "g4": g4 if u > 0 else None,
        "weights": [float(v) for v in x[:-1]], "baseline": BASE,
        "scramble_best": scr[0], "scramble_median": scr[len(scr) // 2],
        "n_scrambles_beating_real": nbet, "n_scrambles": len(scr),
        "scramble_all": [float(v) for v in scr]}
    # THE GUARD THAT LIED, computed and archived rather than described.  This is
    # the control this directory used everywhere else: fix (c,p,s) at the winner
    # and re-solve on the reversed a column.  It does not re-run the search, so
    # it never pays for the freedom the search used, and it returns a comforting
    # ratio on a test that is dead.  [[the-guard-itself-can-be-the-liar]]
    _, naive = fn(rowdata(alphas=A_REAL[::-1]), c, p, s)
    P("  %-14s THE OLD CONTROL, for the record: (c,p,s) FIXED at the winner and"
      % "")
    P("  %-14s   the a column merely reversed -> ||res|| = %.4f, ratio %.2f"
      % ("", naive, naive / res if res else float("inf")))
    P("  %-14s   i.e. it would have reported %s.  It does not re-run the search."
      % ("", "'informative'" if naive > 3 * res else "'not a discrimination'"))
    REC["steps"]["scan"][name]["naive_fixed_winner_control"] = {
        "resnorm": naive, "ratio": naive / res if res else None}
    P("")
    if name == "B1 universal":
        P("  The residual landscape in c, at the winning (p, s) -- a single sharp")
        P("  minimum is a measurement; a forest of comparable ones is a wiggle")
        P("  fitting four points.")
        land = sorted((fn(RD, cc, p, s)[1], cc) for cc in CS)
        P("     %s" % "  ".join("c=%.1f:%.4f" % (cc, v) for v, cc in land[:8]))
        under = [cc for v, cc in land if v < 3 * res]
        P("     values of c within 3x of the winner: %d of %d  ->  %s"
          % (len(under), len(CS), sorted(under)[:12]))
        REC["steps"]["scan"]["landscape_c"] = {"within3x": sorted(under),
                                               "ngrid": len(CS)}
        P("")

# ---- B3: the n-dependence alone
P("-" * 78)
P("B3 -- 'AN n-DEPENDENCE THAT IS NOT 1/n^5', ON ITS OWN")
P("-" * 78)
P("  Keep the gauge sector at 1/n^5 (it is exact, s1) and let the WHOLE fermion")
P("  sector run as 1/n^q with one overall weight.  Two unknowns, ten equations --")
P("  the tightest test in this file.")
P("")
def solve_q(rows, q):
    A, b = [], []
    for r in rows:
        for d in (1, 2):
            fv = float(sum(m * mult * basis(r["a"], s, c, d, q)[0]
                           for rep, e_, ep, mult in r["cont"]
                           for m, s, c in terms(rep, e_, ep)))
            A.append([fv, 0.0 if d == 1 else -r["c"]])
            b.append(-r["G%d" % d])
    return lstsq(A, b)


QS = [q / 10.0 for q in range(30, 81, 5)]
P("  %-8s %-10s %-10s %s" % ("q", "weight", "g4", "||res||"))
B3 = []
for q in QS:
    x, res = solve_q(RD, q)
    g4 = 1 / math.sqrt(x[1]) if x[1] > 0 else float("nan")
    B3.append({"q": q, "w": float(x[0]), "g4": g4, "res": res})
    P("  %-8.1f %-10.4f %-10s %.4f"
      % (q, x[0], ("%.3f" % g4) if x[1] > 0 else "IMAG", res))
bq = min(B3, key=lambda r: r["res"])
scr3 = sorted(min(solve_q(rowdata(alphas=qq), q)[1] for q in QS) for qq in PERMS)
P("")
P("  best q = %.1f at weight %.4f, g4 = %.3f, ||res|| = %.4f"
  % (bq["q"], bq["w"], bq["g4"], bq["res"]))
P("  CONTROL, same scan over q on all %d scrambles: best %.4f (%.2f x real), "
  "median %.4f" % (len(PERMS), scr3[0], scr3[0] / bq["res"], scr3[len(scr3) // 2]))
P("  %s"
  % ("  >> A REAL DISCRIMINATION" if scr3[0] > 3 * bq["res"] else
     "  >> *** NOT A DISCRIMINATION -- two free parameters fit anything ***"))
P("  And the weight it wants is %.4f where their eq. (72) says 1: the fermion"
  % bq["w"])
P("  sector would have to be %.0fx weaker than their own formula." % (1 / bq["w"]))
REC["steps"]["n_power"] = {"scan": B3, "best": bq,
                           "scramble_best": scr3[0],
                           "scramble_median": scr3[len(scr3) // 2],
                           "scramble_all": [float(v) for v in scr3]}

# ---- B4: the winners handed to the independent instrument
P("")
P("-" * 78)
P("B4 -- AND THE WINNERS HANDED TO THE VACUUM, WHICH NEVER SAW a_min OR m_h")
P("-" * 78)
P("  Their Table 1 rows exist because they BREAK electroweak symmetry.  Any")
P("  repair must leave V''(0) < 0 (a = 0 a maximum) on all five.  V''(0) is")
P("  computed numerically here because a 1/n^q term changes the zeta values the")
P("  closed-form D is built from.")
P("")


def V2_at_origin(cont, extra=None):
    """V''(0) for a row, optionally plus one extra channel (v, c, p, s-rule)."""
    t = float(sum(m * basis(1e-9, s, c, 2)[0] for m, s, c in GAUGE))
    t += float(sum(m * mult * basis(1e-9, s, c, 2)[0]
                   for rep, e_, ep, mult in cont for m, s, c in terms(rep, e_, ep)))
    if extra:
        v, c, p, srule = extra
        for rep, e_, ep, mult in cont:
            t += (v[rep] if isinstance(v, dict) else v) * mult * float(
                basis(1e-9, srule * e_ * ep, c, 2, p)[0])
    return t


P("  %-6s %-16s %-16s %s" % ("case", "V''(0) theirs", "V''(0) + B2 winner", "EWSB"))
B4 = []
r2, c2, p2, s2, x2 = WIN["B2 per-rep"]
vd = {R: float(x2[i]) for i, R in enumerate(REPS)}
for tag, cont, *_ in T1:
    a0, a1 = V2_at_origin(cont), V2_at_origin(cont, (vd, c2, p2, s2))
    B4.append({"case": tag, "V2_theirs": a0, "V2_repaired": a1})
    P("  %-6s %-16.4f %-16.4f %s"
      % (tag, a0, a1, "kept" if a1 < 0 else "*** LOST ***"))
P("")
P("  (a_min itself is not recomputed here: ANCHOR_SECTION_31.md s6's instrument")
P("   note applies -- argmin jumps between local minima and its Jacobian is")
P("   meaningless.  V''(0) is a sign, not a fit, which is why it is the one")
P("   thing this programme states about modified content.)")
REC["steps"]["winners_vs_vacuum"] = B4

os.makedirs("paper_data", exist_ok=True)
with open("paper_data/su7_repair_space.json", "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("wrote paper_data/su7_repair_space.json")
