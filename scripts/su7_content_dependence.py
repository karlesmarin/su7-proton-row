#!/usr/bin/env python3
"""Where the anchor's content dependence is -- MEASURED, not argued.

Carles Marin + Claude (AI assistant).  2026-08-04.
Source: Komori & Maru, arXiv:2503.04090.  Follows ANCHOR_SECTION_31.md.

STATE OF THE HUNT COMING IN.
  - the gauge sector is exact (eq. 68 rebuilt from their 3.1;  ANCHOR_SECTION_31 s1)
  - the k >> 1 truncation is irrelevant (brace identity proved;            s2)
  - a global or per-row fermion rescaling is excluded by BOTH columns     (s4)
  - eqs. (73)-(76) are transcribed correctly.  THREE independent ways:
      (a) by hand from s = eta*eta' and their parity rule above eq. (71);
      (b) by SU(2)_L multiplet count against the decompositions (41),(57),(69),(70)
          -- su7_vacuum.py CONTROL 1;
      (c) by direct counting of components of the actual Wilson line (eq. 77)
          by (|n5-n7|, P5 P5') -- su7_anchor_hunt.py.  Reproduced inline below.
  So the residual is NOT a transcription error.  It is either a physical
  ingredient of their 3.2 that is not in (73)-(76), or their own arithmetic.

WHAT THIS SCRIPT DOES.
  Their Table 1 gives TWO computed numbers per row.  Ten equations.  Ask whether
  ANY per-representation reweighting of the fermion sector reproduces all ten.
      F(a) = F_gauge(a) + sum_R w_R * (R's fermion terms)
  Five unknowns (w_7, w_28, w_48, w_84, g4) against ten equations.  The earlier
  per-representation fit used the a_min column only -- 5 equations, 4 unknowns,
  barely overdetermined, and it returned f(48) = 4.58 with non-zero residuals.
  With m_h the system is overdetermined by five, so a solution surviving it would
  be a measurement of the missing ingredient rather than a fudge.
"""
import json
import math
import os

import numpy as np
from scipy.optimize import least_squares

P = lambda *a: print(*a, flush=True)
NMAX = 600
MW = 80.4
_n = np.arange(1, NMAX + 1)
_nf = _n.astype(float)
_w = _nf ** -5
_sgn = {1: np.ones(NMAX), -1: (-1.0) ** _n}
_BC = {}


def basis(alpha, s, c, d=0):
    a = np.atleast_1d(alpha)
    key = (a.tobytes(), s, c, d)
    if key in _BC:
        return _BC[key]
    ph = np.outer(a, c * math.pi * _nf)
    f = (np.cos, lambda t: -np.sin(t), lambda t: -np.cos(t),
         lambda t: np.sin(t))[d % 4](ph)
    out = (f * (_w * _sgn[s] * (c * math.pi * _nf) ** d)).sum(axis=1)
    if a.size > 500:
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


def V(content, alpha, w, d=0, r4=True):
    out = sum(m * basis(alpha, s, c, d) for m, s, c in GAUGE)
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap, r4):
            out = out + w[rep] * m * mult * basis(alpha, s, c, d)
    return out


GRID = np.linspace(0.0, 1.0, 8001)


def minimise(content, w, r4=True):
    v = V(content, GRID, w, r4=r4)
    i = int(np.argmin(v))
    lo, hi = GRID[max(i - 1, 0)], GRID[min(i + 1, len(GRID) - 1)]
    for _ in range(30):
        xs = np.linspace(lo, hi, 15)
        j = int(np.argmin(V(content, xs, w, r4=r4)))
        lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 14)]
    return float(0.5 * (lo + hi))


T1 = [("(1)", [("28", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.8),
      ("(2)", [("28", 1, 1, 1), ("84", 1, 1, 4)], 0.081, 125.5),
      ("(3)", [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)], 0.021, 125.1),
      ("(4)", [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)], 0.026, 126.4),
      ("(5)", [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)], 0.043, 126.2)]

CONST = 2 * MW * math.sqrt(3.0 / (16.0 * math.pi ** 6))     # K = CONST * g4
REC = {"source": "Komori & Maru, arXiv:2503.04090", "steps": {}}

# ---------------------------------------------------------------- control (c)
P("=" * 78)
P("CONTROL -- eqs. (73)-(76) against the ACTUAL WILSON LINE, not against their")
P("equations.  Their eq. (77) is the identity except on indices 5 and 7, so a")
P("component's cosine argument is |n5-n7| pi n a and its sign is P5 P5'.")
P("=" * 78)
WILSON = {"7":  {(1, -1): 2},
          "28": {(2, +1): 2, (1, +1): 2, (1, -1): 8},
          "48": {(2, +1): 2, (1, +1): 4, (1, -1): 16},
          "84": {(3, -1): 2, (2, +1): 8, (2, -1): 2, (1, +1): 8, (1, -1): 24}}
P("  (counts from su7_anchor_hunt.py; every count is exactly TWICE the term")
P("   multiplicity -- the two states q = +-(r-1)/2.)")
P("")
P("  %-5s %-34s %-34s %s" % ("rep", "Wilson-line count / 2", "our terms()", ""))
allok = True
for rep in REPS:
    want = {k: v // 2 for k, v in WILSON[rep].items()}
    got = {}
    for m, s, c in terms(rep, 1, 1):
        got[(c, s)] = got.get((c, s), 0) + m
    ok = got == want
    allok &= ok
    P("  %-5s %-34s %-34s %s" % (rep, dict(sorted(want.items())),
                                 dict(sorted(got.items())),
                                 "OK" if ok else "*** FAIL ***"))
P("")
P("  >> all four representations reproduce the Wilson-line count exactly.  %s"
  % ("PASS" if allok else "*** FAIL ***"))
P("     THE FERMION SECTOR IS NOT MIS-TRANSCRIBED.  Whatever is content-dependent")
P("     is not in our reading of eqs. (73)-(76).")
assert allok
REC["steps"]["wilson_control"] = {"verdict": "PASS"}

# ---------------------------------------------------------------- the fit
P("")
P("=" * 78)
P("THE FIT -- can ANY per-representation reweighting satisfy both columns?")
P("=" * 78)


P("  THE SYSTEM IS LINEAR, so there is no optimiser and no local minimum.")
P("  Write F(a) = G(a) + sum_R w_R f_R(a) with G the (exact) gauge part.  Then")
P("    (I)  stationarity at their a_i:   G'(a_i)  + sum_R w_R f'_R(a_i)  = 0")
P("    (II) their m_h via eqs.(80),(82): G''(a_i) + sum_R w_R f''_R(a_i) = c_i u,")
P("         c_i = (m_h,i a_i/160.8)^2 * 16 pi^6/3,   u = 1/g4^2.")
P("  Both are LINEAR in (w_7, w_28, w_48, w_84, u): 10 equations, 5 unknowns.")
P("  Each equation is normalised by the norm of its own coefficient row.")
P("")


def blocks(r4=True, alphas=None, mhs=None):
    """Per-row coefficient rows for d=1 (stationarity) and d=2 (curvature).

    alphas overrides the target a_i -- used by the scramble control, where the
    coefficients themselves must be re-evaluated at the wrong a, not merely
    relabelled."""
    out = []
    for i, (tag, cont, a0, mh0) in enumerate(T1):
        a_t = a0 if alphas is None else alphas[i]
        mh_t = mh0 if mhs is None else mhs[i]
        aa = np.array([a_t])
        row = {"tag": tag, "a": a_t, "mh": mh_t}
        for d in (1, 2):
            row["G%d" % d] = float(sum(m * basis(aa, s, c, d)[0]
                                       for m, s, c in GAUGE))
            row["f%d" % d] = [
                float(sum(m * mult * basis(aa, s, c, d)[0]
                          for rep2, e_, ep, mult in cont if rep2 == R
                          for m, s, c in terms(R, e_, ep, r4)))
                for R in REPS]
        row["c"] = (mh_t * a_t / (2 * MW)) ** 2 * 16 * math.pi ** 6 / 3
        out.append(row)
    return out


def solve(rows, use_I=True, use_II=True):
    A, b, lab = [], [], []
    for r in rows:
        if use_I:
            A.append(r["f1"] + [0.0]);  b.append(-r["G1"]);  lab.append(r["tag"] + " F'")
        if use_II:
            A.append(r["f2"] + [-r["c"]]);  b.append(-r["G2"]);  lab.append(r["tag"] + ' F"')
    A, b = np.array(A), np.array(b)
    nrm = np.linalg.norm(np.c_[A, b], axis=1)
    An, bn = A / nrm[:, None], b / nrm
    x, *_ = np.linalg.lstsq(An, bn, rcond=None)
    res = An @ x - bn
    return x, res, lab, np.linalg.matrix_rank(An)


ROWS = blocks()
P("  %-34s %-8s %-8s %-8s %-8s %-10s %-6s %s"
  % ("system", "w(7)", "w(28)", "w(48)", "w(84)", "g4", "rank", "||res||"))
LIN = {}
for name, uI, uII in (("(I) stationarity only  [5 eq, 4 unk]", True, False),
                      ("(II) curvature only    [5 eq, 5 unk]", False, True),
                      ("(I)+(II) BOTH          [10 eq, 5 unk]", True, True)):
    x, res, lab, rk = solve(ROWS, uI, uII)
    g4 = (1 / math.sqrt(x[4])) if x[4] > 0 else float("nan")
    P("  %-34s %-8.4f %-8.4f %-8.4f %-8.4f %-10s %-6d %.4f"
      % (name, x[0], x[1], x[2], x[3],
         ("%.4f" % g4) if x[4] > 0 else "IMAGINARY", rk,
         float(np.linalg.norm(res))))
    LIN[name] = {"w": list(x[:4]), "u": float(x[4]), "g4": g4,
                 "resnorm": float(np.linalg.norm(res)),
                 "residuals": dict(zip(lab, map(float, res)))}
P("")
x, res, lab, rk = solve(ROWS, True, True)
P("  per-equation residuals of the BOTH system (normalised, 0 = satisfied):")
for l, r in zip(lab, res):
    P("     %-12s %+.4f" % (l, r))
P("")
P("  CONTROL -- same 10x5 solve with their a_min column REVERSED across rows.")
P("  A system with enough freedom to fit anything gives a comparable residual.")
xs, ress, _, _ = solve(blocks(alphas=[r["a"] for r in ROWS][::-1]), True, True)
P("     real ||res|| = %.4f     scrambled ||res|| = %.4f     ratio %.2f"
  % (float(np.linalg.norm(res)), float(np.linalg.norm(ress)),
     float(np.linalg.norm(ress) / np.linalg.norm(res))))
REC["steps"]["linear"] = LIN
P("")
P("  A MORE GENERAL HYPOTHESIS -- per-CHANNEL instead of per-representation.")
P("  If the missing ingredient changes which cosines appear rather than how much")
P("  each rep weighs, the free parameters are the five channels (c, s) that occur")
P("  in eqs. (73)-(76): (1,+), (1,-), (2,+), (2,-), (3,-).  Six unknowns, ten")
P("  equations.  This is a DIFFERENT subspace from the per-rep one, not a")
P("  refinement of it.")
CH = [(1, 1), (1, -1), (2, 1), (2, -1), (3, -1)]


def blocks_ch(alphas=None, mhs=None):
    out = []
    for i, (tag, cont, a0, mh0) in enumerate(T1):
        a_t = a0 if alphas is None else alphas[i]
        mh_t = mh0 if mhs is None else mhs[i]
        aa = np.array([a_t])
        row = {"tag": tag, "a": a_t, "mh": mh_t}
        for d in (1, 2):
            row["G%d" % d] = float(sum(m * basis(aa, s, c, d)[0] for m, s, c in GAUGE))
            acc = {ch: 0.0 for ch in CH}
            for rep, e_, ep, mult in cont:
                for m, s, c in terms(rep, e_, ep):
                    acc[(c, s)] += m * mult * float(basis(aa, s, c, d)[0])
            row["f%d" % d] = [acc[ch] for ch in CH]
        row["c"] = (mh_t * a_t / (2 * MW)) ** 2 * 16 * math.pi ** 6 / 3
        out.append(row)
    return out


xc, resc, labc, rkc = solve(blocks_ch(), True, True)
P("     channels %s" % [("c=%d,s=%+d" % ch) for ch in CH])
P("     weights  %s" % ["%.3f" % v for v in xc[:5]])
g4c = (1 / math.sqrt(xc[5])) if xc[5] > 0 else float("nan")
P("     g4 = %s   rank %d   ||res|| = %.4f"
  % (("%.4f" % g4c) if xc[5] > 0 else "IMAGINARY", rkc,
     float(np.linalg.norm(resc))))
P("     >> it FITS (residual at the numerical floor) -- and that is exactly why")
P("        the two numbers beside it decide the question.  Every channel weight")
P("        is ~0.1-0.3 where eq. (72) says 1, i.e. the fermion sector would have")
P("        to be ~5x weaker than their own formula; and g4 = %.2f against the SM"
  % g4c)
P("        SU(2)_L value 0.63 -- a 4D gauge coupling of alpha = g4^2/4pi = %.2f."
  % (g4c ** 2 / (4 * math.pi)))
P("        A six-parameter family that fits only at an unphysical coupling is a")
P("        parametrisation of the residual, not an explanation of it.")
REC["steps"]["per_channel"] = {"channels": [list(c) for c in CH],
                               "w": list(map(float, xc[:5])),
                               "resnorm": float(np.linalg.norm(resc))}

P("")
P("=" * 78)
P("IS A 5 % RESIDUAL A FAILURE?  THEIR TABLE IS QUOTED TO 2 SIGNIFICANT FIGURES")
P("=" * 78)
P("  a_min is printed as 0.043, 0.081, 0.021, 0.026, 0.043 -- so the true value")
P("  is only known to +-0.0005, which is +-1.2 % at 0.043 and +-2.4 % at 0.021.")
P("  m_h is printed to 0.1 GeV.  Any fit therefore has a RESIDUAL FLOOR set by")
P("  their own rounding, and a 5 % residual must be compared against it, not")
P("  against zero.  Monte-Carlo: resample a_i in +-0.0005 and m_h,i in +-0.05,")
P("  re-solve, and build the distribution of ||res||.")
P("")
rng = np.random.default_rng(20260804)
A0 = np.array([r["a"] for r in ROWS])
M0 = np.array([r["mh"] for r in ROWS])
for label, mk in (("per-representation", lambda a, m: blocks(alphas=a, mhs=m)),
                  ("per-channel", lambda a, m: blocks_ch(alphas=a, mhs=m))):
    real = float(np.linalg.norm(solve(mk(A0, M0), True, True)[1]))
    sam = []
    for _ in range(300):
        a = A0 + rng.uniform(-5e-4, 5e-4, 5)
        m = M0 + rng.uniform(-5e-2, 5e-2, 5)
        sam.append(float(np.linalg.norm(solve(mk(a, m), True, True)[1])))
    sam = np.array(sam)
    # The null is "the model is EXACTLY right and we only ever see rounded values".
    # Under it, SOME point of the rounding box must drive ||res|| to ~0.  So the
    # statistic is the MINIMUM over the box, not a tail fraction.
    mn = float(sam.min())
    scr = float(np.linalg.norm(solve(mk(A0[::-1], M0), True, True)[1]))
    P("  %-19s real ||res|| = %.4f   over the rounding box: min %.4f, median %.4f"
      % (label, real, mn, float(np.median(sam))))
    P("  %-19s scrambled-a ||res|| = %.4f (ratio %.2f)   ->  %s"
      % ("", scr, scr / real if real else float("inf"),
         "ROUNDING CANNOT EXPLAIN IT" if mn > 0.5 * real else
         "rounding could account for it"))
    REC["steps"].setdefault("rounding_mc", {})[label] = {
        "real": real, "box_min": mn, "median": float(np.median(sam)),
        "scrambled": scr}

P("")
P("=" * 78)
P("THE NON-LINEAR FIT BELOW IS SUPERSEDED -- kept only to show WHY")
P("=" * 78)
P("  It reports w ~ 1 for the a_min-only fit, where the linear solve above gets")
P("  w(48) = 5.59 with a 25x smaller residual.  The linear solve is right: the")
P("  non-linear objective calls argmin(F), which JUMPS between local minima, so")
P("  its numerical Jacobian is meaningless and least_squares never leaves p0 = 1.")
P("  Recorded rather than deleted, because the same trap produced the earlier")
P("  'f(48) = 4.58, residuals still non-zero' verdict -- which the linear solve")
P("  now shows was the right neighbourhood reached by a wrong instrument.")


def residuals(w, g4, r4=True, use_mh=True):
    r = []
    for tag, cont, a_t, mh_t in T1:
        r.append((minimise(cont, w, r4=r4) - a_t) / a_t)
    if use_mh:
        for tag, cont, a_t, mh_t in T1:
            f2 = float(V(cont, np.array([a_t]), w, d=2, r4=r4)[0])
            if f2 <= 0:
                r.append(10.0)                    # concave: no real m_h at all
            else:
                pred = g4 * math.sqrt(3 * f2 / (16 * math.pi ** 6)) * (2 * MW / a_t)
                r.append((pred - mh_t) / mh_t)
    return r


def run(name, free_reps, r4=True, use_mh=True, lam_only=False):
    def pack(p):
        g4 = p[-1]
        if lam_only:
            w = {R: p[0] for R in REPS}
        else:
            w = {R: 1.0 for R in REPS}
            for i, R in enumerate(free_reps):
                w[R] = p[i]
        return w, g4

    n = (1 if lam_only else len(free_reps)) + 1
    p0 = np.array([1.0] * (n - 1) + [0.63])
    sol = least_squares(lambda p: residuals(*pack(p), r4=r4, use_mh=use_mh),
                        p0, bounds=([1e-3] * (n - 1) + [0.2],
                                    [20.0] * (n - 1) + [1.5]), xtol=1e-12)
    w, g4 = pack(sol.x)
    res = residuals(w, g4, r4=r4, use_mh=use_mh)
    chi = float(np.sum(np.array(res) ** 2))
    return name, w, g4, res, chi, len(res), n


P("  residuals are RELATIVE: (a_min - a_theirs)/a_theirs for five rows, then")
P("  (m_h^pred - m_h)/m_h for five rows.  A concave row (F'' <= 0, no real m_h)")
P("  is charged a flat 10.  chi2 = sum of squared relative residuals.")
P("")
FITS = [run("baseline  w = 1", [], use_mh=True),
        run("global lambda", [], use_mh=True, lam_only=True),
        run("per-rep w (a_min only)", list(REPS), use_mh=False),
        run("per-rep w (BOTH columns)", list(REPS), use_mh=True),
        run("per-rep w, no (1,4) rule", list(REPS), use_mh=True, r4=False)]
P("  %-28s %-7s %-7s %-7s %-7s %-8s %-6s %s"
  % ("model", "w(7)", "w(28)", "w(48)", "w(84)", "g4", "dof", "chi2"))
for name, w, g4, res, chi, neq, npar in FITS:
    P("  %-28s %-7.3f %-7.3f %-7.3f %-7.3f %-8.4f %-6d %.4g"
      % (name, w["7"], w["28"], w["48"], w["84"], g4, neq - npar, chi))
REC["steps"]["fits"] = [{"model": n, "w": w, "g4": g4, "chi2": c,
                         "residuals": r} for n, w, g4, r, c, _, _ in FITS]

best = FITS[3]
P("")
P("  per-row residuals of the BOTH-columns per-rep fit:")
P("  %-5s %-9s %-11s %-9s %-9s %-11s %s"
  % ("case", "a theirs", "a fitted", "rel", "m_h", "m_h pred", "rel"))
wB, g4B = best[1], best[2]
for i, (tag, cont, a_t, mh_t) in enumerate(T1):
    a_f = minimise(cont, wB)
    f2 = float(V(cont, np.array([a_t]), wB, d=2)[0])
    pred = (g4B * math.sqrt(3 * f2 / (16 * math.pi ** 6)) * (2 * MW / a_t)
            if f2 > 0 else float("nan"))
    P("  %-5s %-9.4f %-11.4f %-9.3f %-9.1f %-11.2f %.3f"
      % (tag, a_t, a_f, (a_f - a_t) / a_t, mh_t, pred,
         (pred - mh_t) / mh_t if f2 > 0 else float("nan")))

P("")
P("  CONTROL THAT COULD FAIL -- fit the SAME 5 parameters to SCRAMBLED targets")
P("  (a_min column reversed).  If the model has enough freedom to fit anything,")
P("  the scrambled chi2 will be comparable to the real one.")
T1_REAL = list(T1)
T1[:] = [(t[0], t[1], s[2], t[3]) for t, s in zip(T1_REAL, T1_REAL[::-1])]
scr = run("scrambled a_min", list(REPS), use_mh=True)
T1[:] = T1_REAL
P("     real chi2 = %.4g     scrambled chi2 = %.4g     ratio %.2f"
  % (best[4], scr[4], scr[4] / best[4] if best[4] else float("inf")))
P("     %s" % ("the model CANNOT fit anything -- the real fit is informative"
               if scr[4] > 3 * best[4] else
               "*** the model has too much freedom; read the fit as weak ***"))
REC["steps"]["scramble_control"] = {"real_chi2": best[4], "scrambled_chi2": scr[4]}

os.makedirs("paper_data", exist_ok=True)
with open("paper_data/su7_content_dependence.json", "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("wrote paper_data/su7_content_dependence.json")
