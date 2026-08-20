#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

THEIR WHOLE FERMION SECTOR, REBUILT FROM THE GROUP.

su7_gauge_from_group.py did this for the gauge sector: their eq. (68) came out of
the adjoint's own (q, P5P5', P6) table instead of being trusted from their eq.
(62). Their FERMION blocks -- eqs. (73)-(76), the 7, the 28, the 48 and the 84 --
were still only checked against their own upstream equations. That is the last
place in the potential where we take their word for a multiplicity.

It matters for one specific reason. The two rows of their Table 1 we fail on are
exactly the two containing a 48, and 48 is the ADJOINT of SU(7). If the term list
for the 48 were wrong -- ours or theirs -- that is precisely where it would show.
So derive it, do not read it.

THE RULE, and it is one rule for all four multiplets:

  the Wilson line is sym(5,7), so in the fundamental only the combinations
  a_+- = (e5 +- e7)/sqrt2 carry charge, q = +-1/2, and both have P5P5' = -1;
  the other five basis vectors are neutral.  A component of Sym^d(7) is a
  multiset of d letters: its charge is the sum and its P5P5' is the product.
  Then a component enters the potential with

      c = 2|q|        and        sign = eta*eta' * P5P5'

Twenty blocks come out of that one rule, so it is overdetermined many times over
and can fail loudly. Nineteen agree on the first run. The twentieth did not, and
it is the useful one: at c = 1 in the 84 the group gives 16 charged pairs where a
naive reading of eq. (76) gives 15. The sixteenth is a_+a_+a_-, and it is not a
missing term of theirs. Sym^3 of the charged doublet is an SU(2) QUADRUPLET, so
one multiplet feeds c = 3 AND c = 1; their own sect. 3.2 list (1 quadruplet, 5
triplets, 15 doublets) already implies 16, and their note under eq. (71) says so.
What was naive is reading "c = 1" as "a doublet" -- our shortcut, not their error.
Section C folds their four published lists and gets the group's numbers exactly.

What this settles, and what it does not: the term lists of eqs. (73)-(76) are not
where our alpha_min discrepancy lives -- not for the 48, not for anything. The
discrepancy has to be in what is done WITH the terms.
"""
import re
from itertools import combinations_with_replacement as multiset

import numpy as np

P = lambda *a: print(*a, flush=True)

# The term lists we are testing must be THE ONES THE PAPER'S NUMBERS COME FROM,
# not a copy of them -- a copy would drift and the test would go quietly vacuous.
# su7_vacuum.py has no __main__ guard (importing it runs the whole minimisation),
# so lift its terms() out of the source and exec exactly that.
_src = open("su7_vacuum.py", encoding="utf-8").read()
_fn = re.search(r"\ndef terms\(.*?\n(?=\n\ndef |\n\nGRID|\Z)", _src, re.S)
assert _fn, "terms() not found in su7_vacuum.py -- the file moved under us"
_ns = {}
exec(_fn.group(0), _ns)
terms = _ns["terms"]
N = 7

# --------------------------------------------------------------- the fundamental
P5P5 = np.diag([1, 1, 1, 1, 1, -1, -1]) @ np.diag([1, 1, 1, -1, -1, -1, 1])
H = np.zeros((N, N))                          # the Wilson line: T^43 = sym(5,7)
H[4, 6] = H[6, 4] = 0.5

P("=" * 78)
P("A -- THE FUNDAMENTAL: WHICH LETTERS CARRY CHARGE")
P("=" * 78)
w, vecs = np.linalg.eigh(H)                   # Hermitian: eigh, never eig
letters = []                                  # (q, p) per basis vector of the 7
for i in range(N):
    v = vecs[:, i]
    p = float(v @ P5P5 @ v)
    assert abs(abs(p) - 1) < 1e-9, "not a P5P5' eigenvector: %s" % p
    letters.append((round(float(w[i]), 9), int(round(p))))
for q, p in sorted(letters, reverse=True):
    P("    q = %+.1f   P5P5' = %+d" % (q, p))
chg = [l for l in letters if l[0] != 0]
P("  charged letters %d, neutral %d" % (len(chg), N - len(chg)))
P("  >> and the two charged ones share P5P5' = -1, which is the precondition")
P("     their sect. 3.2 needs and never states.")
assert len(chg) == 2 and all(p == -1 for _, p in chg)

# ------------------------------------------------------- Sym^d and the term list
DEG = {"7": 1, "28": 2, "48": None, "84": 3}


def blocks_sym(d):
    """(c, sign_relative_to_s) -> multiplicity, for the charged components of
    Sym^d(7).  Pairs (q,-q) are counted once."""
    out = {}
    for comb in multiset(range(N), d):
        q = sum(letters[i][0] for i in comb)
        if q <= 0:                            # count each (q,-q) pair once
            continue
        p = 1
        for i in comb:
            p *= letters[i][1]
        k = (int(round(2 * q)), p)
        out[k] = out.get(k, 0) + 1
    return out


def blocks_adjoint():
    """the adjoint is not a Sym^d: 7 x 7bar minus a singlet.  Charge is the
    DIFFERENCE of the two letters, P5P5' still the product."""
    out = {}
    for i in range(N):
        for j in range(N):
            if i == j:
                continue                      # the Cartan directions: q = 0 anyway
            q = letters[i][0] - letters[j][0]
            if q <= 0:
                continue
            p = letters[i][1] * letters[j][1]
            k = (int(round(2 * q)), p)
            out[k] = out.get(k, 0) + 1
    return out


P("")
P("=" * 78)
P("B -- EQS. (73)-(76), DERIVED AND COMPARED TERM BY TERM")
P("=" * 78)
P("  sign is written relative to s = eta*eta'.  '+s' means the term carries s,")
P("  '-s' means it carries -s.  Our reading of their equations is su7_vacuum.py.")
P("")
def compare(r4):
    ok_all = True
    for rep in ("7", "28", "48", "84"):
        d = DEG[rep]
        got = blocks_adjoint() if d is None else blocks_sym(d)
        # their list, as su7_vacuum reads it, at eta = eta' = +1 so that s = +1
        theirs = {}
        for m, sgn, c in terms(rep, 1, 1, r4=r4):
            theirs[(c, int(sgn))] = theirs.get((c, int(sgn)), 0) + m
        P("  --- %s%s" % (rep, "  (the adjoint)" if d is None else "  = Sym^%d(7)" % d))
        for c, p in sorted(set(got) | set(theirs), reverse=True):
            a, b = got.get((c, p), 0), theirs.get((c, p), 0)
            ok_all &= a == b
            P("      c = %d  sign = %2ss   group %-4d ours %-4d %s"
              % (c, "+" if p > 0 else "-", a, b, "" if a == b else "*** MISMATCH ***"))
        P("      charged states in the multiplet: %d" % (2 * sum(got.values())))
        P("")
    return ok_all


P("  With the eq.(71) quadruplet term ON (r4=True), which is how every number in")
P("  the paper is computed:")
P("")
on = compare(True)
P("  all twenty blocks reproduced: %s" % on)
assert on, "the fermion sector does NOT come out of the group"

P("")
P("  CONTROL -- the same comparison with that term OFF. If it also passed, the")
P("  test would be insensitive to the one term we did not derive, i.e. vacuous.")
P("")
off = compare(False)
P("  reproduced with the term off: %s   (it must be False, or this test is blind)"
  % off)
assert not off, "the control did not fire: the r4 term makes no difference"

P("")
P("  >> So the term is FORCED. Sym^3 of the charged doublet (a_+, a_-) is an")
P("     SU(2) QUADRUPLET, charges 3/2, 1/2, -1/2, -3/2. Its q = 3/2 member is")
P("     a_+a_+a_+ and its q = 1/2 member is a_+a_+a_-, so ONE multiplet feeds")
P("     TWO values of c. That is the whole content of their note under eq. (71),")
P("     and we had been applying it on trust.")

P("")
P("=" * 78)
P("C -- AND THEIR OWN DECOMPOSITION LISTS ALREADY SAID SO")
P("=" * 78)
P("  Their sect. 3.2 lists each multiplet's SU(2)_L content. An SU(2) multiplet")
P("  of dimension d has members q = (d-1)/2, ..., -(d-1)/2, so it contributes a")
P("  charged pair at EVERY c = d-1, d-3, ... down to 1 or 0 -- a quadruplet at")
P("  c = 3 and again at c = 1. Fold their lists that way and compare to B.")
P("")
LISTED = {"7": {2: 1}, "28": {3: 1, 2: 5}, "48": {3: 1, 2: 10},
          "84": {4: 1, 3: 5, 2: 15}}          # their sect. 3.2, verbatim
P("  %-5s %-26s %-18s %-18s" % ("rep", "their SU(2) list", "folded to c", "from the group"))
okD = True
for rep in ("7", "28", "48", "84"):
    fold = {}
    for dim, mult in LISTED[rep].items():
        for c in range(dim - 1, 0, -2):        # every charged pair of the multiplet
            fold[c] = fold.get(c, 0) + mult
    grp = {}
    for (c, _p), m in (blocks_adjoint() if DEG[rep] is None
                       else blocks_sym(DEG[rep])).items():
        grp[c] = grp.get(c, 0) + m
    okD &= fold == grp
    P("  %-5s %-26s %-18s %-18s %s"
      % (rep, LISTED[rep], dict(sorted(fold.items(), reverse=True)),
         dict(sorted(grp.items(), reverse=True)), "" if fold == grp else "*** FAIL ***"))
P("")
P("  their published lists and SU(7) agree on all four multiplets: %s" % okD)
assert okD
P("")
P("  >> So the extra term is not a repair of anything of theirs. Their eq. (76)")
P("     and their sect. 3.2 list are consistent with each other and with the")
P("     group; what is naive is reading c = 1 as 'a doublet'. The 84 has 15")
P("     doublets AND a quadruplet, and the quadruplet's inner pair is the")
P("     sixteenth term at c = 1. Their note under eq. (71) says exactly this.")
P("")
P("  Consequence for our own control: su7_vacuum's CONTROL 1 maps c -> SU(2)")
P("  dimension as d = c+1 and needs r4 OFF to pass. That mapping is only valid")
P("  when no multiplet reaches past its top member, which the 84's quadruplet")
P("  does. The control's CONCLUSION survives -- the last two terms of eq. (76)")
P("  must be 84 terms, since as 48 terms the doublet count cannot be made to")
P("  work either way -- but it passes for a reason narrower than it claims, and")
P("  section C above is the version that does not depend on the c+1 shortcut.")

P("")
P("=" * 78)
P("D -- WHAT THIS EXCLUDES")
P("=" * 78)
adj = blocks_adjoint()
P("  the 48 is the adjoint, and of its 48 generators only %d carry Wilson-line"
  % (2 * sum(adj.values())))
P("  charge -- the other %d are invisible to V_eff at any order in alpha."
  % (48 - 2 * sum(adj.values())))
P("")
P("  Their eq. (76)'s 48 block counts exactly those %d, split %s by P5P5'."
  % (2 * sum(adj.values()),
     " and ".join("%d" % (2 * v) for _, v in sorted(adj.items(), reverse=True))))
P("  So the term list for the 48 is right, and it is right for the same reason")
P("  and by the same rule as the 7, the 28 and the 84.")
P("")
P("  >> The two rows of their Table 1 we fail on are the two containing a 48.")
P("     After this, that failure cannot be a miscount of the 48: the block is")
P("     derived from SU(7) itself. Whatever differs between their alpha_min and")
P("     ours is downstream of the term lists, not inside them.")
P("")
P("  And one thing this changed about OUR OWN reading, not theirs. su7_vacuum's")
P("  r4 flag is justified in its docstring as a second, smaller mass EIGENVALUE")
P("  for the 84's SU(2) quadruplet -- a second eigenvalue of the same states.")
P("  Section B shows that is not what it is: it is the quadruplet's twelfth")
P("  STATE, a_+a_+a_-, which a state count does check and did. The flag was")
P("  right; the reason we gave for it was not.")
