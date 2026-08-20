# ⚛️ Proton Decay in SU(7) Grand Gauge-Higgs Unification

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22033302-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.22033302)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Gates](https://img.shields.io/badge/gates-10_green-1B6F8C)](scripts/)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](.)

**An obstruction, its minimal escapes, and the one row of their Table 1 that can pay for them.**

**📄 Paper (EN + ES) and every verification script on Zenodo → https://doi.org/10.5281/zenodo.22033302**

> ### 📚 Part **VI** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU***
>   → [github.com/karlesmarin/ghu-su4-completion](https://github.com/karlesmarin/ghu-su4-completion) · [Zenodo 10.5281/zenodo.21432625](https://doi.org/10.5281/zenodo.21432625)
> - **Part II — *Three Gates to a Quark Generation***
>   → [github.com/karlesmarin/su4-sm-cell-criterion](https://github.com/karlesmarin/su4-sm-cell-criterion) · [Zenodo 10.5281/zenodo.21432627](https://doi.org/10.5281/zenodo.21432627)
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential***
>   → [github.com/karlesmarin/centre-parity-selection](https://github.com/karlesmarin/centre-parity-selection) · [Zenodo 10.5281/zenodo.21438226](https://doi.org/10.5281/zenodo.21438226)
> - **Part IV — *Schur Functions at (1,−1,t,t⁻¹)***
>   → [github.com/karlesmarin/schur-nonidentity-o4](https://github.com/karlesmarin/schur-nonidentity-o4) · [Zenodo 10.5281/zenodo.21463000](https://doi.org/10.5281/zenodo.21463000)
> - **Part V — *What the Higgs Potential Cannot See***
>   → [github.com/karlesmarin/higgs-blind-class](https://github.com/karlesmarin/higgs-blind-class) · [Zenodo 10.5281/zenodo.21727094](https://doi.org/10.5281/zenodo.21727094)
> - **Part VI — *Proton Decay in SU(7) Grand Gauge-Higgs Unification*** (this repo)

Komori and Maru (arXiv:2503.04090) leave proton decay to future work in one sentence. This paper takes
that sentence as its question and answers the half that can be answered from their own equations: at
one generation the extra U(1) their model carries **cannot** protect the proton, and the minimal
repairs that could are enumerated, priced, and reduced to a single published row.

Nothing here is a defect in their paper. Their §3 is explicit enough that a reader can rebuild it from
its own equations — which is what §2 of this paper does, and **everything of theirs we checked came out
right**. The one column we cannot reproduce is written up as a question of ours, not as a finding about
them.

## 🙏 With thanks

This work exists because Komori and Maru wrote their model out fully, and because they said plainly what
they had not done — naming proton decay as future work. Stating an open question in public is a generous
thing to do.

Almost nothing here is a new mechanism, and the ones used cost other people a great deal of time:
Gogoladze–Mimura–Nandi (the vector-like adjoint and its `Z₃` cure), Lee–Ma (the family-dependent charge),
Costa–Dobrescu–Fox (the general `U(1)` anomaly solution), Adachi–Maru (reading the curvature at the origin
as a selection rule), Dumitru–Guo–Korthals Altes and Guo–Du (the two-loop halves), and
Hosotani–Maru–Takenaga–Yamashita (the only two-loop Wilson-line potential in this framework). §9 of the
paper says which piece is whose, one by one, and that table is the part we would most like a reader to
check.

## 🧭 What is in it

| | |
|---|---|
| **The obstruction** | Their eqs. (43)–(47) and (76)–(79) force `U(1)' = T₃L + Y − (B−L)` identically, and all four dimension-6 `\|ΔB\|=1` operators have `Y = B−L = 0`. At one generation the symmetry that survives cannot forbid a single one of them. |
| **What the anomalies demand** | Six channels; three force `X_Q = −1/6`, which is exactly `A = 0`. Two do not cancel on the published spectrum and, absent Green–Schwarz, require a Standard-Model singlet of `U(1)'` charge `−1` — **necessary, not merely minimal** — and only the `28` supplies it, so only rows (2) and (3) of their Table 1 can. |
| **The bill, in eighths** | `V''(0) = −π²ζ(3)D` with `D` an exact rational per multiplet: a `48⁽⁺'⁺⁾` is worth exactly `0` and an `84⁽⁺'⁺⁾` exactly `5/4`. Hosting a third lepton generation costs one `84`, which takes case (3) from `9/8` to `−1/8`. **Case (2) is the only published row that satisfies both conditions at once.** |
| **The selection rule** | Protection fails iff `q_φ = \|A_j\|/n`: a countable set with a maximum, so protection is a **half-line and not a tuning**. And `q_φ = 3/2` — the charge of the `84`'s pure index-7 component, present in every row of their table — forbids all four operators to all orders in `⟨φ⟩/M`, for every half-integer brane-quark charge. |
| **The instrument, and where it fails** | The invariant `K = m_h·α_min/√(F″)` must be one number on every row of any such table. We do not reproduce their `α_min` column and say so as an open question rather than a defect: the truncation, the transcription and every multiplicative repair are excluded **at one loop**, and a sixth row is **pre-registered with both candidate values** so the question does not depend on anyone answering it. |
| **Where the result stops** | The verdict survives every multiplicative one-loop deformation: a wedge `(27/46, 27/26)` containing `w = 1` and every repair ever fitted, a factor 1.77 wide, and `D(48) ≡ 0` makes the largest of those repairs invisible to it *identically*. What the wedge does not absorb is loop order. The closest published analogue — 4D and thermal, both halves — lands at `1.0398` against the ceiling `1.0385`, crossing at `g₄ = 0.6205`. **Two-loop safety is left open rather than defused**, and an earlier draft's argument that it fell the safe way is withdrawn inside the paper. |

## 🔍 Reproducing it

```bash
python scripts/su7_anomaly_channels.py     # the six channels, and X_Q = -1/6
python scripts/su7_vacuum.py               # D per multiplet, exact rationals in eighths
python scripts/su7_qphi.py                 # the closed-form selection rule
python scripts/su7_qphi_socratic.py        # and the adversarial audit of it
python scripts/su7_anchor_mh.py            # the m_h invariant, the column nobody tested
python scripts/su7_repair_space.py         # the region result: how much an error could eat
```

Every number printed in the paper is greppable in `outputs/`, and `scripts/check_numbers.py` is the
gate that enforces it — in **both** editions.

## 🔒 The gates

Ten, and they are the reason the paper says what it says:

- `check_numbers.py` — every displayed number backed by an archived run (0 not found, in both editions);
- `check_parity.py` — the two editions carry the same structural elements;
- `check_refs.py` — labels, references, bibitems;
- `check_figures.py` — every figure placed, cited, and with its own images per edition;
- `check_layout.py` — no page with a blank band, and no ink running off the bottom of the sheet. The
  second half was added on publication day: Figure 1's caption was printing past the paper edge in
  both editions, and the filter that stopped the page number counting as content was hiding it;
- `check_letter_map.py` — the page map of the letter to the authors, against the PDF **they were sent**;
- `check_scripts.py` — every script the paper names exists and has an archived output;
- `check_formulas.py` — the display formulas of the two editions compared **body by body**, delimiters
  balanced, and the compiler's `.log` read for dropped glyphs. It found one on its first run;
- `check_reproduces.py` — every cited script is **re-run** and its output diffed against the archive.
  The other eight read files and believe them; this one is the only thing standing behind the sentence
  *every displayed number regenerates from the ancillary scripts*. It found two scripts that had stopped
  running at all, and the numbers they had produced were right;
- `check_prose.py` — a copula running straight into a capitalised word with no full stop between,
  which is what a lead-in left behind by a moved paragraph looks like. The Spanish edition had one:
  *«Los artefactos principales son Un instrumento interactivo…»*. Eight gates and a full read missed it.

## ⚖️ Honesty ledger

**No defect in anyone's paper is claimed anywhere in this work.** The `α_min` discrepancy of §7 is
stated as an open question of ours, with three transcription controls passing, so a residual we cannot
place is more likely ours than theirs. What §2 reports about their `28` and about `π₅ = π₇` are
**unstated data with computed consequences**, not errors.

The mechanisms are other people's and §9 says whose, one by one — including the one this paper
proves and does **not** own: that a real representation under diagonal `±1` parities furnishes no
chiral generation is stated as well known in Gogoladze–Mimura–Nandi, PRD **69** (2004) 075006, which
also names the `Z₃` projection that escapes it. What is ours there is the specialisation to their
parities and the counting.

**What is not claimed**: that case (2) *works* — it is the only candidate left, not a construction;
any `α_min`, `1/R₅` or `m_h` for modified content; any vacuum expectation value, scale or lifetime;
and operators above dimension 6, whose charges were not enumerated.

## 📖 Citation

```bibtex
@misc{marin2026su7proton,
  author = {Mar\'in Mu\~noz, Carles},
  title  = {Proton Decay in {SU(7)} Grand Gauge-Higgs Unification: An Obstruction,
            Its Minimal Escapes, and the One Row of Their Table~1 That Can Pay for Them},
  year   = {2026},
  doi    = {10.5281/zenodo.22033302},
  note   = {Part VI of a series on higher-dimensional gauge-Higgs unification}
}
```

Carles Marín Muñoz · independent researcher · with Claude (Anthropic) as research assistant.
Apache-2.0.
