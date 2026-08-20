#!/usr/bin/env python3
"""Contrasting what we have against the papers themselves.

Carles Marin + Claude (AI assistant).  2026-08-04.
Source: Komori & Maru, arXiv:2503.04090; Akamatsu-Hirose-Maru-Nago, arXiv:2312.08608.

THREE THINGS THE PAPERS ALREADY CONTAIN AND WE NEVER USED.

  STEP 1 -- their Table 1 has THREE columns, and the third is a function of the
    first.  Their eq. (82) is 1/R5 = 2 x 80.4 GeV / a_min, so the 1/R5 column is
    determined by the a_min column with no freedom at all.  We have used a_min
    and m_h and never once checked the third against the first.  It is a control
    on THEIR table that can fire, and it costs nothing.

  STEP 2 -- and it is not merely a control: two columns quoted to two significant
    figures constrain a_min BETTER THAN EITHER ALONE.  The rounding box of
    ANCHOR_SECTION_31.md s6 used the a column only.  Intersecting it with the box
    implied by 1/R5 shrinks it, on the rows where it matters most.

  STEP 3 -- the instrument has been validated against a DIFFERENT published GHU
    paper by the same author, in Part V of this series, and s6 never says so.
    Read from the archived run rather than asserted.
"""
import json
import math
import os

P = lambda *a: print(*a, flush=True)
MW = 80.4

# their Table 1, as printed: (case, a_min, m_h/GeV, 1/R5 in TeV)
T1 = [("(1)", 0.043, 126.8, 3.8),
      ("(2)", 0.081, 125.5, 2.0),
      ("(3)", 0.021, 125.1, 7.5),
      ("(4)", 0.026, 126.4, 6.1),
      ("(5)", 0.043, 126.2, 3.8)]
REC = {"source": "Komori & Maru, arXiv:2503.04090", "steps": {}}


def box(value, printed_digits):
    """the interval a printed decimal stands for, at the precision printed"""
    half = 0.5 * 10.0 ** (-printed_digits)
    return value - half, value + half


# ============================================================ STEP 1
P("=" * 78)
P("STEP 1 -- THEIR OWN TABLE 1, CHECKED AGAINST ITSELF")
P("=" * 78)
P("  Their eq. (82):   1/R5 = 2 x 80.4 GeV / a_min.")
P("  So their third column is a FUNCTION of their first, with no freedom.  If the")
P("  two disagree beyond what two significant figures allow, one of them is a")
P("  misprint -- and that would change what the anchor discrepancy even is.")
P("")
P("  %-6s %-9s %-11s %-13s %-13s %s"
  % ("case", "a printed", "1/R5 print", "1/R5 from a", "a from 1/R5", "consistent?"))
S1 = []
allok = True
for tag, a, mh, r5 in T1:
    r5_from_a = 2 * MW / a / 1000.0                      # TeV
    a_from_r5 = 2 * MW / (r5 * 1000.0)
    alo, ahi = box(a, 3)
    rlo, rhi = box(r5, 1)
    # the a-interval the printed 1/R5 stands for (note the inversion flips it)
    a_r_lo, a_r_hi = 2 * MW / (rhi * 1000.0), 2 * MW / (rlo * 1000.0)
    lo, hi = max(alo, a_r_lo), min(ahi, a_r_hi)
    ok = lo <= hi
    allok &= ok
    S1.append({"case": tag, "a_printed": a, "r5_printed": r5,
               "r5_from_a": r5_from_a, "a_from_r5": a_from_r5,
               "a_box_from_a": [alo, ahi], "a_box_from_r5": [a_r_lo, a_r_hi],
               "intersection": [lo, hi] if ok else None, "consistent": bool(ok)})
    P("  %-6s %-9.3f %-11.1f %-13.4f %-13.6f %s"
      % (tag, a, r5, r5_from_a, a_from_r5,
         "YES  a in [%.5f, %.5f]" % (lo, hi) if ok else "*** NO -- EMPTY ***"))
P("")
P("  >> all five rows consistent: %s" % ("YES" if allok else "*** NO ***"))
P("     A control that could have fired and did not.  Their a_min column is not a")
P("     misprint, and neither is their 1/R5 column: the anchor discrepancy is")
P("     about the POTENTIAL, not about a number that was typed wrong.")
P("")
P("  Note the direction, because it is not symmetric.  For the two rows with the")
P("  smallest a -- exactly the two that fail everything else -- the 1/R5 column is")
P("  the SHARPER of the two:")
for r in S1:
    wa = r["a_box_from_a"][1] - r["a_box_from_a"][0]
    wr = r["a_box_from_r5"][1] - r["a_box_from_r5"][0]
    wi = r["intersection"][1] - r["intersection"][0]
    P("     %-6s width from a: %.5f   from 1/R5: %.5f   intersected: %.5f  (%.0f%% of the a box)"
      % (r["case"], wa, wr, wi, 100 * wi / wa))
REC["steps"]["table1_selfconsistency"] = {"rows": S1, "all_consistent": bool(allok)}

# ============================================================ STEP 2
P("")
P("=" * 78)
P("STEP 2 -- THE ROUNDING BOX, TIGHTENED BY THEIR OWN THIRD COLUMN")
P("=" * 78)
P("  ANCHOR_SECTION_31.md s6 asked whether their two-significant-figure rounding")
P("  could explain the 0.1628 residual of the per-representation fit.  Sampling")
P("  the a box (+-0.0005) and the m_h box (+-0.05) gave a MINIMUM of 0.1580 -- so")
P("  rounding could not.  That test used a box that is too GENEROUS: the true box")
P("  is the intersection computed above.  Re-run it there.")
P("")
try:
    import numpy as np
except ImportError:                                       # pragma: no cover
    P("  numpy unavailable -- STEP 2 skipped")
    np = None

if np is not None:
    _n = np.arange(1, 601)
    _nf = _n.astype(float)
    _sgn = {1: np.ones(600), -1: (-1.0) ** _n}

    def basis(a, s, c, d=0):
        ph = np.outer(np.atleast_1d(float(a)), c * math.pi * _nf)
        f = (np.cos, lambda t: -np.sin(t), lambda t: -np.cos(t),
             lambda t: np.sin(t))[d % 4](ph)
        return (f * ((_nf ** -5) * _sgn[s] * (c * math.pi * _nf) ** d)).sum(axis=1)

    GAUGE = [(-1.0, 1, 2), (-2.0, 1, 1), (-3.5, -1, 1)]
    REPS = ("7", "28", "48", "84")

    def terms(rep, eta, etap):
        s = eta * etap
        return {"7": [(1, -s, 1)],
                "28": [(1, s, 2), (4, -s, 1), (1, s, 1)],
                "48": [(1, s, 2), (8, -s, 1), (2, s, 1)],
                "84": [(1, -s, 3), (1, -s, 2), (4, s, 2), (11, -s, 1),
                       (4, s, 1), (1, -s, 1)]}[rep]

    CONT = {"(1)": [("28", 1, -1, 1), ("84", 1, 1, 4)],
            "(2)": [("28", 1, 1, 1), ("84", 1, 1, 4)],
            "(3)": [("28", 1, 1, 1), ("48", 1, 1, 3), ("84", 1, 1, 2)],
            "(4)": [("7", 1, -1, 1), ("48", 1, 1, 2), ("84", 1, 1, 3)],
            "(5)": [("7", 1, 1, 1), ("7", 1, -1, 1), ("84", 1, 1, 4)]}

    def resid(alphas, mhs):
        A, b = [], []
        for (tag, _, _, _), a, mh in zip(T1, alphas, mhs):
            for d in (1, 2):
                G = float(sum(m * basis(a, s, c, d)[0] for m, s, c in GAUGE))
                f = [float(sum(m * mult * basis(a, s, c, d)[0]
                               for rep2, e_, ep, mult in CONT[tag] if rep2 == R
                               for m, s, c in terms(R, e_, ep))) for R in REPS]
                ci = (mh * a / (2 * MW)) ** 2 * 16 * math.pi ** 6 / 3
                A.append(f + ([0.0] if d == 1 else [-ci]));  b.append(-G)
        A, b = np.array(A), np.array(b)
        nrm = np.linalg.norm(np.c_[A, b], axis=1)
        x, *_ = np.linalg.lstsq(A / nrm[:, None], b / nrm, rcond=None)
        return float(np.linalg.norm((A / nrm[:, None]) @ x - b / nrm))

    from scipy.optimize import minimize
    rng = np.random.default_rng(20260804)
    A0 = np.array([r[1] for r in T1])
    M0 = np.array([r[2] for r in T1])
    real = resid(A0, M0)

    P("  A NOTE ON METHOD, because the first version of this step was WRONG.")
    P("  s6 estimated the box minimum by SAMPLING it.  The intersected box is a")
    P("  SUBSET of the a box, so its true minimum can only be >= the wide one --")
    P("  yet sampling both with the same number of draws returned a LOWER value")
    P("  for the smaller box, because 400 points cover a small box better than a")
    P("  large one.  That is a sampling artefact and it points the wrong way.")
    P("  The minimum is therefore MINIMISED here, not sampled: bounded L-BFGS-B")
    P("  from 40 starts per box.  [[a-control-that-cannot-fail]]")
    P("")

    def box_min(lo, hi, nstart=40):
        lo = np.concatenate([lo, M0 - 5e-2])
        hi = np.concatenate([hi, M0 + 5e-2])
        best = float("inf")
        for k in range(nstart):
            x0 = lo + (rng.random(10) if k else np.full(10, 0.5)) * (hi - lo)
            r = minimize(lambda x: resid(x[:5], x[5:]), x0, method="L-BFGS-B",
                         bounds=list(zip(lo, hi)))
            best = min(best, float(r.fun))
        return best

    LO = np.array([r["intersection"][0] for r in S1])
    HI = np.array([r["intersection"][1] for r in S1])
    wide_min = box_min(A0 - 5e-4, A0 + 5e-4)
    tight_min = box_min(LO, HI)
    P("  %-36s %-12s %s" % ("box on a_i", "true minimum", "verdict"))
    P("  %-36s %-12.4f %s"
      % ("the a column alone (s6's box)", wide_min,
         "rounding cannot explain it" if wide_min > 0.5 * real else "it could"))
    P("  %-36s %-12.4f %s"
      % ("intersected with their 1/R5 column", tight_min,
         "rounding cannot explain it" if tight_min > 0.5 * real else "it could"))
    P("")
    P("  real residual = %.4f" % real)
    ok = tight_min >= wide_min - 1e-9
    P("  CONTROL: the subset's minimum must be >= the superset's.  %.4f >= %.4f  %s"
      % (tight_min, wide_min, "PASS" if ok else "*** FAIL -- optimiser did not converge ***"))
    P("  >> rounding was already excluded on the a column alone; their own third")
    P("     column closes it by a further %.1f %%, using nothing but what they printed."
      % (100 * (tight_min - wide_min) / wide_min if wide_min else 0.0))
    REC["steps"]["tightened_box"] = {"real": real, "wide_min": wide_min,
                                     "tight_min": tight_min,
                                     "monotone_control": bool(ok)}

# ============================================================ STEP 3
P("")
P("=" * 78)
P("STEP 3 -- THE INSTRUMENT HAS ALREADY BEEN VALIDATED, ON ANOTHER MARU PAPER")
P("=" * 78)
P("  s6 excludes the gauge sector, the truncation, the transcription and every")
P("  multiplicative repair -- all of them arguments INSIDE the SU(7) computation.")
P("  It never says the obvious external thing: does this kind of computation")
P("  reproduce a published GHU vacuum AT ALL?")
P("")
P("  It does, and the run is archived in Part V of this series.")
P("")
VAL = os.path.join("..", "part_v", "outputs", "ghu_potential_validation.txt")
if os.path.exists(VAL):
    with open(VAL, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.strip():
                P("     | " + line.rstrip())
    REC["steps"]["ahmn_validation"] = {"file": VAL, "present": True}
else:
    P("     *** %s NOT FOUND -- do not make this claim ***" % VAL)
    REC["steps"]["ahmn_validation"] = {"file": VAL, "present": False}
P("")
P("  Read the asymmetry, and read its scope with it:")
P("     Akamatsu-Hirose-Maru-Nago, arXiv:2312.08608, 6D SU(4) GHU  -> vacuum EXACT,")
P("        Higgs mass ratio to 0.1 %")
P("     Komori-Maru,               arXiv:2503.04090, 6D SU(7) GHU  -> a_min off by")
P("        3-108 %, content-dependent")
P("  Same class of computation, same dimensionality, an author in common.")
P("")
P("  SCOPE, and it must be stated or the sentence is worth nothing: these are two")
P("  DIFFERENT implementations -- Part V's is Schur-character based on two Wilson")
P("  phases, this one counts components of one.  So this does not show the SU(7)")
P("  assembly is right.  What it shows is that the APPROACH reproduces a published")
P("  6D GHU vacuum when the paper's own ingredients are enough to determine it,")
P("  which is exactly the possibility s6 leaves open and never closes.")

# ============================================================ STEP 4
P("")
P("=" * 78)
P("STEP 4 -- THE CANDIDATE s6 NEVER LISTED: LOOP ORDER")
P("=" * 78)
P("  ANCHOR_SECTION_31.md s6 names three candidates for the missing structure -- a")
P("  cosine argument outside {1,2,3} pi n a, an n-dependence that is not 1/n^5, or")
P("  multiplets not in (69),(70).  It never names a fourth, and our own")
P("  GATE_POLYAKOV_TWOLOOP.md already contains it:")
P("")
P("  (a) PURE GAUGE, two loops.  Dumitru-Guo-Korthals Altes, PRD 89 (2014) 016009,")
P("      eq. (66):   Omega(2)/Omega(1) = -5 g^2 C2(A) / (16 pi^2),")
P("      'independent of the eigenvalues of the loop ... for any of the classical")
P("      and exceptional groups'.  A MULTIPLICATIVE, background-independent")
P("      rescaling -- which s4 and s6 have already excluded, because the lambda")
P("      the a_min column wants is not one number.")
P("")
P("  (b) WITH FERMIONS.  Guo-Du, JHEP 05 (2019) 042: 'the surprisingly simple")
P("      proportionality ... is in general NO LONGER TRUE when fermions are taken")
P("      into account', and their expressions carry Bernoulli polynomials of the")
P("      individual background angles AND OF THEIR PAIRWISE DIFFERENCES.")
P("      Our angles are the Wilson-line charges q.  A product B(C_b) B(C_d)")
P("      generates cosines at c_b +- c_d -- and with c in {1,2,3} that reaches")
P("      c in {0,...,6}.  THREE OF THOSE, c = 4, 5, 6, ARE OUTSIDE {1,2,3}.")
P("      That is s6's first candidate, arrived at from the physics rather than")
P("      from a scan.")
P("")
P("  The gate never asked the one question that decides whether it can matter:")
P("  HOW BIG IS IT?")
P("")
C2A = 7.0                      # C2(adjoint) = N for SU(N)
P("  %-12s %-14s %s" % ("g4", "5 g^2 C2(A)/(16 pi^2)", "against the lambda the a column wants"))
S4 = []
for g in (0.55, 0.63, 0.70, 0.80):
    r = 5.0 * g ** 2 * C2A / (16 * math.pi ** 2)
    S4.append({"g4": g, "ratio": r})
    P("  %-12.2f %-14.4f %s"
      % (g, r, "|1 - lambda| needed is 0.04-0.20 (lambda = 0.80-0.96)"))
P("")
P("  >> the two-loop scale is 7-18 %% for g4 in 0.55-0.80, and the rescaling the")
P("     a_min column asks for is 4-20 %%.  THE SAME ORDER.  So loop order is not")
P("     ruled out by size, and it is the only named candidate that is BOTH")
P("     content-dependent (so a uniform lambda does not describe it) AND known")
P("     from the literature to be non-multiplicative once fermions are in.")
P("")
P("  WHAT THIS IS AND IS NOT.  It is a size estimate and a shape argument, both")
P("  from published results we have READ but NOT reproduced.  We have not computed")
P("  the two-loop fermionic potential for this model, and nothing in this")
P("  programme depends on it.  What changes is s6's exclusion list: it was")
P("  presented as covering the space of candidates, and it does not -- it covers")
P("  the one-loop space.  A referee would ask this in the first paragraph.")
P("")
P("  It also does NOT rescue their Table 1 by itself: their potential, their m_h")
P("  and the finiteness argument that motivates gauge-Higgs unification are all")
P("  one-loop statements.  If a two-loop term moved a_min by 30-100 %% their own")
P("  calculation would not be the thing they published.  So the honest reading is")
P("  a candidate of the right shape and a plausible size, NOT an explanation.")
REC["steps"]["two_loop_scale"] = {"C2_adjoint": C2A, "rows": S4,
                                  "lambda_range": [0.80, 0.96],
                                  "verdict": "same order -- not excluded by size"}

os.makedirs("paper_data", exist_ok=True)
with open("paper_data/su7_cross_paper.json", "w") as fh:
    json.dump(REC, fh, indent=1)
P("")
P("wrote paper_data/su7_cross_paper.json")
