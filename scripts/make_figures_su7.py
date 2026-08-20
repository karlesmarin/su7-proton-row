#!/usr/bin/env python3
"""make_figures_su7.py - the four figures of the SU(7) GGHU line.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Every number drawn here is READ from paper_data/su7_repair_space.json,
su7_realisable.json and su7_wedge_direction.json, each produced by the script of
the same name and archived in outputs/.
Nothing is recomputed and nothing is transcribed.  [[save-the-outputs-not-just-the-scripts]]

  fig_repair_space     the headline as an exact REGION of repair space, and the
                       two planes whose signs cut it out.  The point of the
                       figure is that the conclusion is not a point estimate.
  fig_anchor_controls  the control that killed Part B.  118 scrambles of their
                       own a_min column, the whole search re-run on each; the
                       real residual sits inside the cloud.
  fig_table1_verdicts  their Table 1 with our verdict columns, colour-coded --
                       including BOTH realisable assignments, one donating one
                       84 and one donating two.
  fig_ratio_line       the verdict as a statement about ONE number, and the two
                       exposures sitting on opposite sides of w = 1.

Palette: validated with the dataviz validator (light surface, categorical, 3
slots) -- lightness band PASS, chroma floor PASS, CVD separation dE 20.7
(deutan) / 22.7 (tritan) PASS, normal-vision floor 24.3 PASS.  The amber carries
a contrast WARN, discharged by direct labels on every region it fills.
Diverging use (the sign of D) is the same two hues with a neutral midpoint --
blue positive, red negative, consistently in all three figures.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ES = "--es" in sys.argv


def T(en, es):
    return es if ES else en


SUF = "_es" if ES else ""
OUT = os.path.join(HERE, "paper")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(HERE, "paper_data", "su7_repair_space.json")) as fh:
    DAT = json.load(fh)
with open(os.path.join(HERE, "paper_data", "su7_realisable.json")) as fh:
    REAL = json.load(fh)
with open(os.path.join(HERE, "paper_data", "su7_wedge_direction.json")) as fh:
    WDIR = json.load(fh)
with open(os.path.join(HERE, "paper_data", "twoloop_wedge.json")) as fh:
    TW = json.load(fh)          # the redrawn Fig. 5: the NET transplanted ratio

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
STROKE, INK, MUTED = "#1F4E79", "#1F2933", "#6B7280"
GRID = "#E3E1DC"
DIVERGE = LinearSegmentedColormap.from_list("d", [RED, "#EFEDE8", BLUE])

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED,
                     "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "xtick.labelsize": 8,
                     "ytick.labelsize": 8})

HR = DAT["steps"]["headline_region"]
G0, G28, G84 = HR["coef"]["D2_donated"]          # D2 = G0 + G28 w28 + G84 w84
_, H28, H84 = HR["coef"]["D3_donated"]
WLO, WHI = HR["w_diagonal_interval"]


def bare(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


# ============================================================ FIG 1
def fig_repair_space():
    fig = plt.figure(figsize=(10.4, 4.35))
    lohi = (0.40, 1.60)

    # ---- panel A: the three regions, categorical
    ax = fig.add_subplot(1, 2, 1)
    g = np.linspace(*lohi, 900)
    X, Y = np.meshgrid(g, g)
    D2 = G0 + G28 * X + G84 * Y
    D3 = G0 + H28 * X + H84 * Y
    cat = np.where(D2 <= 0, 0, np.where(D3 < 0, 1, 2))
    ax.imshow(cat, origin="lower", extent=[*lohi, *lohi], aspect="auto",
              cmap=ListedColormap([RED, AMBER, BLUE]), vmin=0, vmax=2,
              interpolation="nearest", alpha=0.90)
    w = np.linspace(*lohi, 200)
    for coef, lab, lx in (((G28, G84), r"$D_{(2)}=0$", 0.82),
                          ((H28, H84), r"$D_{(3)}=0$", 1.28)):
        ax.plot(w, (-G0 - coef[0] * w) / coef[1], color="white", lw=2.6,
                solid_capstyle="round", zorder=3)
        ax.plot(w, (-G0 - coef[0] * w) / coef[1], color=STROKE, lw=1.1, zorder=4)
        ly = (-G0 - coef[0] * lx) / coef[1]
        assert lohi[0] < ly < lohi[1], (lab, ly)
        ang = np.degrees(np.arctan(-coef[0] / coef[1]))
        ax.text(lx, ly, lab, fontsize=8.2, color=STROKE, zorder=6,
                ha="center", va="center", rotation=ang, rotation_mode="anchor",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none",
                          alpha=0.9))
    ax.plot(lohi, lohi, ls=(0, (4, 3)), color=INK, lw=0.9, zorder=5)
    for v in (WLO, WHI):
        ax.plot([v], [v], "o", ms=5.5, mfc="white", mec=INK, mew=1.2, zorder=7)
    ax.annotate(r"$w=%.3f$" % WLO, (WLO, WLO), textcoords="offset points",
                xytext=(-34, 6), fontsize=7.6, color=INK, zorder=9)
    ax.annotate(r"$w=%.3f$" % WHI, (WHI, WHI), textcoords="offset points",
                xytext=(6, 5), fontsize=7.6, color=INK, zorder=9)

    MK = {"their formulas, w = 1": ("*", 15),
          "(I) a_min only        s6": ("s", 6),
          "(II) curvature only   s6": ("D", 5.4),
          "(I)+(II) both columns s6": ("^", 7),
          "per-row lam, lowest   s4": ("v", 7),
          "per-row lam, highest  s4": ("P", 7.5)}
    NICE = {"their formulas, w = 1": T("their formulas, $w=1$", "sus formulas, $w=1$"),
            "(I) a_min only        s6": r"fit $\alpha_{\min}$",
            "(II) curvature only   s6": T("fit curvature", "ajuste curvatura"),
            "(I)+(II) both columns s6": T("fit both columns", "ajuste ambas columnas"),
            "per-row lam, lowest   s4": r"$\lambda$ min",
            "per-row lam, highest  s4": r"$\lambda$ max"}
    for r in HR["fitted"]:
        if r["repair"] not in MK:
            continue
        mk, ms = MK[r["repair"]]
        ax.plot([r["w28"]], [r["w84"]], mk, ms=ms, mfc="white", mec=INK,
                mew=1.15, zorder=8)
    ax.annotate(T("their formulas, $w=1$", "sus formulas, $w=1$"), (1.0, 1.0),
                textcoords="offset points", xytext=(9, -17), fontsize=8,
                color=INK, zorder=9,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none",
                          alpha=0.85))

    # region labels: fixed points, each ASSERTED to lie in the region it names
    def which(a, b):
        return 0 if G0 + G28 * a + G84 * b <= 0 else (
            1 if G0 + H28 * a + H84 * b < 0 else 2)

    for (px, py), k, lab in (((0.505, 0.455), 0,
                              T("both rows\nlose EWSB",
                                "ambas filas\npierden EWSB")),
                             ((0.655, 1.30), 1,
                              T("CASE (2) UNIQUE", "CASO (2) UNICO")),
                             ((1.41, 1.43), 2,
                              T("both rows keep EWSB\n(no unique row)",
                                "ambas conservan EWSB\n(ninguna fila unica)"))):
        assert which(px, py) == k, (lab, which(px, py))
        ax.text(px, py, lab, fontsize=9.4 if k == 1 else 8.2,
                color="#4A3000" if k == 1 else "white", ha="center",
                va="center", weight="bold", zorder=6, linespacing=1.35)
    ax.set_xlabel(r"$w(\mathbf{28})$", fontsize=10)
    ax.set_ylabel(r"$w(\mathbf{84})$", fontsize=10)
    ax.set_title(T("the headline is a region, not a point",
                   "el titular es una region, no un punto"), fontsize=10)
    ax.set_xlim(*lohi);  ax.set_ylim(*lohi)
    bare(ax)

    # ---- panel B: the two planes whose signs cut it out
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    g2 = np.linspace(*lohi, 40)
    X2, Y2 = np.meshgrid(g2, g2)
    Z2 = G0 + G28 * X2 + G84 * Y2
    Z3 = G0 + H28 * X2 + H84 * Y2
    ax2.plot_surface(X2, Y2, Z2, color=BLUE, alpha=0.55, linewidth=0,
                     antialiased=True, shade=True)
    ax2.plot_surface(X2, Y2, Z3, color=RED, alpha=0.55, linewidth=0,
                     antialiased=True, shade=True)
    ax2.plot_surface(X2, Y2, np.zeros_like(X2), color="#B9B5AC", alpha=0.32,
                     linewidth=0)
    ax2.text(1.62, 1.62, G0 + G28 * 1.6 + G84 * 1.6, r"  $D_{(2)}$",
             color=STROKE, fontsize=9.5, weight="bold")
    ax2.text(1.62, 1.62, G0 + H28 * 1.6 + H84 * 1.6, r"  $D_{(3)}$",
             color=RED, fontsize=9.5, weight="bold")
    ax2.text(0.42, 0.42, 0.25, T("  $D=0$", "  $D=0$"), color="#5B564D",
             fontsize=8.4)
    ax2.set_xlabel(r"$w(\mathbf{28})$", fontsize=9, labelpad=-4)
    ax2.set_ylabel(r"$w(\mathbf{84})$", fontsize=9, labelpad=-4)
    ax2.set_zlabel(T("$D$ after donating one $\\mathbf{84}$",
                     "$D$ tras donar un $\\mathbf{84}$"), fontsize=9, labelpad=-6)
    ax2.set_title(T(r"slopes in $w(\mathbf{84})$: $15/4$ against $5/4$",
                    r"pendientes en $w(\mathbf{84})$: $15/4$ frente a $5/4$"),
                  fontsize=10, pad=-2)
    ax2.tick_params(labelsize=7, pad=-2)
    ax2.view_init(elev=23, azim=-132)
    ax2.set_box_aspect((1, 1, 0.62))

    fig.suptitle(T(r"$D>0$ is electroweak symmetry breaking. "
                   r"$w(\mathbf{7})$ and $w(\mathbf{48})$ do not appear: "
                   r"$D(\mathbf{48})\equiv 0$",
                   r"$D>0$ es ruptura electrodebil. "
                   r"$w(\mathbf{7})$ y $w(\mathbf{48})$ no aparecen: "
                   r"$D(\mathbf{48})\equiv 0$"), fontsize=9.6, y=1.005, color=INK)
    fig.tight_layout(w_pad=2.4)
    save(fig, "fig_repair_space")


# ============================================================ FIG 2
def fig_anchor_controls():
    fig, ax = plt.subplots(figsize=(7.8, 3.45))
    tests = [("B1 universal", T("one universal extra channel\n2 unknowns",
                                "un canal extra universal\n2 incognitas"),
              DAT["steps"]["scan"]["B1 universal"]),
             ("B2 per-rep", T("extra channel per representation\n5 unknowns",
                              "canal extra por representacion\n5 incognitas"),
              DAT["steps"]["scan"]["B2 per-rep"]),
             ("n-power", T(r"fermion sector as $1/n^q$" + "\n2 unknowns",
                           r"sector fermionico como $1/n^q$" + "\n2 incognitas"),
              DAT["steps"]["n_power"])]
    rng = np.random.default_rng(20260804)
    XHI, NT = 5.2, len(tests)
    ax.axvspan(3.0, XHI, color=BLUE, alpha=0.09, lw=0, zorder=0)
    ax.axvline(1.0, color=INK, lw=1.3, zorder=3)
    ax.axvline(3.0, color=BLUE, lw=1.2, ls=(0, (4, 3)), zorder=3)
    yt, yl = [], []
    for i, (_, lab, d) in enumerate(tests):
        y = NT - 1 - i
        scr = np.array(d["scramble_all"])
        real = d.get("resnorm", d.get("best", {}).get("res"))
        rat = scr / real
        assert rat.max() < XHI, ("cloud runs past the axis", rat.max())
        ax.plot([rat.min(), rat.max()], [y, y], color=GRID, lw=9, zorder=1,
                solid_capstyle="round")
        ax.scatter(rat, y + rng.uniform(-0.20, 0.20, rat.size), s=11,
                   facecolor=MUTED, alpha=0.50, linewidth=0, zorder=2)
        nb = int((scr <= real).sum())
        # the ratio is the verdict; the count is secondary and must not read as one
        ax.text(rat.max() * 1.09, y + 0.10, T("best scramble $=%.2f\\times$ real",
                                              "mejor scramble $=%.2f\\times$ real")
                % rat.min(), fontsize=8.2, color=RED, ha="left", va="bottom",
                weight="bold")
        ax.text(rat.max() * 1.09, y - 0.10, T("(%d of %d reach the real fit)",
                                              "(%d de %d alcanzan el ajuste real)")
                % (nb, rat.size), fontsize=7.4, color=MUTED, ha="left", va="top")
        yt.append(y);  yl.append(lab)
    ax.set_xscale("log")
    ax.set_xlim(0.30, XHI)
    ax.set_ylim(-0.70, NT - 0.26)
    ax.set_yticks(yt);  ax.set_yticklabels(yl, fontsize=8.6, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([0.5, 1, 2, 3, 5])
    ax.set_xticklabels([r"$0.5\times$", r"$1\times$", r"$2\times$", r"$3\times$",
                        r"$5\times$"])
    ax.text(1.0, NT - 0.30, T(" the real fit", " el ajuste real"), color=INK,
            fontsize=8.6, ha="left", va="top", weight="bold")
    ax.text(3.08, NT - 0.30,
            T(" a discrimination needs the\n whole cloud beyond here",
              " discriminar exige toda la\n nube mas alla de aqui"),
            color=BLUE, fontsize=8.0, ha="left", va="top", linespacing=1.35)
    ax.set_xlabel(T(r"best $\|\mathrm{res}\|$ of the whole search, "
                    r"$\div$ the real fit's",
                    r"mejor $\|\mathrm{res}\|$ de toda la busqueda, "
                    r"$\div$ el del ajuste real"), fontsize=9.5)
    ax.set_title(T("the control that killed Part B: the search re-run on all 118 "
                   "scrambles of their own $\\alpha_{\\min}$ column",
                   "el control que mato la parte B: la busqueda re-corrida sobre "
                   "las 118 permutaciones de su columna $\\alpha_{\\min}$"),
                 fontsize=9.6, color=INK)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    save(fig, "fig_anchor_controls")


# ============================================================ FIG 3
def fig_table1_verdicts():
    rows = DAT["steps"]["headline_region"]["per_row"]
    content = {"(1)": r"$\mathbf{28}^{(+,-)}+4{\times}\mathbf{84}^{(+,+)}$",
               "(2)": r"$\mathbf{28}^{(+,+)}+4{\times}\mathbf{84}^{(+,+)}$",
               "(3)": r"$\mathbf{28}^{(+,+)}+3{\times}\mathbf{48}^{(+,+)}"
                      r"+2{\times}\mathbf{84}^{(+,+)}$",
               "(4)": r"$\mathbf{7}^{(+,-)}+2{\times}\mathbf{48}^{(+,+)}"
                      r"+3{\times}\mathbf{84}^{(+,+)}$",
               "(5)": r"$\mathbf{7}^{(+,+)}+\mathbf{7}^{(+,-)}"
                      r"+4{\times}\mathbf{84}^{(+,+)}$"}
    cols = [("", 0.030), (T("fermion content", "contenido fermionico"), 0.310),
            (r"$\alpha_{\min}$", 0.075), (r"$m_h$ (GeV)", 0.085),
            (T(r"$\nu_R$ at $-1$", r"$\nu_R$ en $-1$"), 0.095),
            (r"$D$", 0.070), (T(r"$D$ after $-\mathbf{84}$",
                                r"$D$ tras $-\mathbf{84}$"), 0.100),
            (T(r"$D$ after $-2{\times}\mathbf{84}$",
               r"$D$ tras $-2{\times}\mathbf{84}$"), 0.105),
            (T("verdict", "veredicto"), 0.155)]
    x, xs = 0.0, []
    for _, wd in cols:
        xs.append(x);  x += wd
    fig, ax = plt.subplots(figsize=(10.4, 2.55))
    ax.set_xlim(0, x);  ax.set_ylim(0, len(rows) + 1.15)
    ax.axis("off")
    nrm = lambda v: 0.5 + 0.5 * float(np.tanh(v / 1.6))

    for j, (lab, wd) in enumerate(cols):
        ax.text(xs[j] + wd / 2, len(rows) + 0.52, lab, fontsize=8.6,
                ha="center", va="center", color=INK, weight="bold")
    ax.plot([0, x], [len(rows) + 0.14] * 2, color=INK, lw=0.9)

    for i, r in enumerate(rows):
        y = len(rows) - 1 - i + 0.10
        if r["headline"]:
            ax.add_patch(Rectangle((0, y), x, 0.94, facecolor=AMBER,
                                   alpha=0.20, lw=0, zorder=0))
        cells = [r["case"], content[r["case"]], "%.3f" % r["a_theirs"],
                 "%.1f" % r["mh_theirs"], None, None, None, None, None]
        for j, val in enumerate(cells):
            cx = xs[j] + cols[j][1] / 2
            if val is not None:
                ax.text(cx, y + 0.47, val, fontsize=8.4, ha="center",
                        va="center", color=INK)
        # nu_R
        ok = r["nu_R"]
        ax.add_patch(Rectangle((xs[4] + 0.012, y + 0.14), cols[4][1] - 0.024,
                               0.66, facecolor=BLUE if ok else "#EFEDE8",
                               alpha=0.80 if ok else 1.0, lw=0, zorder=1))
        ax.text(xs[4] + cols[4][1] / 2, y + 0.47,
                T("yes", "si") if ok else T("no", "no"), fontsize=8.2,
                ha="center", va="center", zorder=2,
                color="white" if ok else MUTED, weight="bold" if ok else "normal")
        # the two D columns, diverging on the sign
        A2 = {q["n84"]: q["rows"] for q in REAL["assignments"]}[2]
        for j, key in ((5, "D_theirs"), (6, "D_donated"), (7, "two84")):
            if key == "two84":
                cell = A2.get(r["case"])
                if cell is None:
                    continue
                n, d = cell.split("/") if "/" in cell else (cell, "1")
                v = float(n) / float(d)
            else:
                v = r[key]
            ax.add_patch(Rectangle((xs[j] + 0.012, y + 0.14), cols[j][1] - 0.024,
                                   0.66, facecolor=DIVERGE(nrm(v)), lw=0, zorder=1))
            ax.add_patch(Rectangle((xs[j] + 0.012, y + 0.14), 0.0075, 0.66,
                                   facecolor=BLUE if v > 0 else RED, lw=0,
                                   zorder=3))
            ax.text(xs[j] + cols[j][1] / 2 + 0.004, y + 0.47, "%+.3f" % v,
                    fontsize=8.2, ha="center", va="center", zorder=2,
                    color="white" if abs(v) > 1.6 else (INK if v > 0 else RED))
        # verdict
        if r["headline"]:
            txt, col = T("UNIQUE ROW", "FILA UNICA"), "#4A3000"
        elif r["nu_R"]:
            txt, col = T("EWSB lost", "pierde EWSB"), RED
        else:
            txt, col = T(r"no $\nu_R$", r"sin $\nu_R$"), MUTED
        ax.text(xs[8] + cols[8][1] / 2, y + 0.47, txt, fontsize=8.4,
                ha="center", va="center", color=col,
                weight="bold" if r["headline"] else "normal")
        ax.plot([0, x], [y] * 2, color=GRID, lw=0.5, zorder=0)
    ax.set_title(T("Komori-Maru Table 1, with the two conditions imposed. "
                   "$D$ is exact: every entry is a rational in eighths",
                   "Tabla 1 de Komori-Maru con las dos condiciones impuestas. "
                   "$D$ es exacto: cada entrada es un racional en octavos"),
                 fontsize=9.4, color=INK, pad=12)
    fig.tight_layout()
    save(fig, "fig_table1_verdicts")


def fig_ratio_line():
    """The verdict is a statement about ONE number, the fermion-to-gauge weight ratio, and
    section 8 is a statement about where that number may sit: the wedge.  What the wedge does
    not absorb is loop order, so the figure draws the transplanted two-loop ratio against the
    wedge's own ceiling as a function of the gauge coupling.

    NOT drawn any more: the five row-by-row fits to the alpha_min column.  They are a fit and
    not a measurement of w -- section 8 says so and withdraws the argument that used them, and
    a figure may not keep making an argument the text has retracted."""
    ceil = TW["ceiling"]
    xs = np.array([c["g4"] for c in TW["curve"]])
    ws = np.array([c["w"] for c in TW["curve"]])
    gx, gn, wn = TW["crossing_g4"], TW["nominal_g4"], TW["nominal_w"]
    x0, x1 = 0.52, 0.84
    m = (xs >= x0) & (xs <= x1)
    xs, ws = xs[m], ws[m]

    fig, ax = plt.subplots(figsize=(9.7, 3.5))
    ax.set_xlim(x0, x1)
    ax.set_ylim(0.995, max(ws.max(), ceil) + 0.012)
    bare(ax)
    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    # the two horizontals that define the question: one loop, and the wedge's ceiling
    ax.axhline(1.0, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=3)
    ax.text(x0 + 0.004, 1.0 + 0.0015, T("their formulas, one loop: $w=1$",
                                        "sus formulas, un lazo: $w=1$"),
            fontsize=8.2, color=INK, va="bottom")
    ax.axhline(ceil, color=RED, lw=1.2, zorder=3)
    ax.text(x1 - 0.004, ceil - 0.0022,
            T("the wedge's ceiling, $27/26$ --- above this, case (2) is no longer unique",
              "el techo de la cuna, $27/26$ --- por encima, el caso (2) deja de ser unico"),
            fontsize=8.2, color=RED, va="top", ha="right")

    # safe below, lost above
    ax.fill_between(xs, 0.995, np.minimum(ws, ceil), where=ws >= 0, color=AMBER,
                    alpha=0.16, lw=0, zorder=1)
    ax.fill_between(xs, ceil, ws, where=ws > ceil, color=RED, alpha=0.22, lw=0, zorder=2)

    # the transplant itself
    ax.plot(xs, ws, color=STROKE, lw=2.0, zorder=5)
    ax.text(x0 + 0.012, np.interp(x0 + 0.012, xs, ws) - 0.0055,
            T("the transplanted ratio $(1+\\delta_f)/(1+\\delta_b)$",
              "el cociente trasplantado $(1+\\delta_f)/(1+\\delta_b)$"),
            fontsize=8.6, color=STROKE, va="top")

    # the crossing, and the coupling this series actually uses
    ax.plot([gx], [ceil], "o", ms=7, mfc="white", mec=RED, mew=1.6, zorder=7)
    ax.annotate(T("crosses at $g_4=%.4f$" % gx, "cruza en $g_4=%.4f$" % gx),
                xy=(gx, ceil), xytext=(gx - 0.045, ceil + 0.010),
                fontsize=8.4, color=RED,
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.9))
    ax.plot([gn, gn], [0.995, wn], color=MUTED, lw=0.9, ls=(0, (2, 2)), zorder=4)
    ax.plot([gn], [wn], "o", ms=7, mfc=RED, mec="white", mew=1.2, zorder=7)
    ax.annotate(T("$g_4=0.63$, the value this series uses: $w=%.4f$\n"
                  "over the ceiling by $0.13\\,\\%%$" % wn,
                  "$g_4=0.63$, el valor que usa esta serie: $w=%.4f$\n"
                  "por encima del techo un $0.13\\,\\%%$" % wn),
                xy=(gn, wn), xytext=(gn + 0.030, wn + 0.014), fontsize=8.4, color=RED,
                linespacing=1.35,
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.9))

    ax.set_xlabel(T("four-dimensional gauge coupling $g_4$ of the transplanted estimate",
                    "acoplamiento gauge cuatridimensional $g_4$ de la estimacion trasplantada"),
                  fontsize=9.5, color=INK)
    ax.set_ylabel(T("fermion-to-gauge weight ratio",
                    "cociente de pesos fermion/gauge"), fontsize=9.5, color=INK)
    ax.set_title(T("Two-loop safety is open, and this is where its boundary sits",
                   "La seguridad a dos lazos esta abierta, y su frontera esta aqui"),
                 fontsize=9.8, color=INK, pad=12)
    fig.tight_layout()
    save(fig, "fig_ratio_line")


def save(fig, stem):
    p = os.path.join(OUT, stem + SUF + ".pdf")
    fig.savefig(p, bbox_inches="tight")
    fig.savefig(p.replace(".pdf", ".png"), dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s (+ .png)" % p)


if __name__ == "__main__":
    fig_repair_space()
    fig_anchor_controls()
    fig_table1_verdicts()
    fig_ratio_line()
