#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

THE ONE SECTOR NEVER CHECKED AGAINST THE GROUP THEORY.

su7_content_dependence.py verifies their eqs. (73)-(76), the FERMION potential,
against the actual Wilson line: it counts components, reads off |n5-n7| and
P5 P5', and reproduces every term.  The GAUGE potential, their eq. (68), has
only ever been rebuilt from their own eq. (62) -- which checks their algebra,
not their content.

That asymmetry matters more than it looks.  su7_loop_order_margin.py shows the
verdict of this paper depends on the RATIO of the fermion weight to the gauge
weight, with 1/26 of margin above; and the wedge of REPAIR_SPACE.md holds the
gauge coefficient FIXED at -27/8 while letting every fermion weight float.  So
the one sector the robustness argument cannot move is the one sector never
verified independently.

This closes it.  The three quantum numbers are good simultaneously -- the
Wilson line of their eq. (77) acts only on indices 5 and 7, and both P5 P5' and
P6 take the SAME value on those two, so all three commute (that is the
precondition their sect. 3.2 needs and never states).  So every one of the 48
adjoint generators carries a definite (q, P5 P5', P6), and their eq. (68) has to
be a sum over that table.
"""
import math
import numpy as np

P = lambda *a: print(*a, flush=True)
N, DIM = 7, 48

P5V = np.diag([1, 1, 1, 1, 1, -1, -1]).astype(float)
P5P = np.diag([1, 1, 1, -1, -1, -1, 1]).astype(float)
P6V = np.diag([1, 1, 1, -1, -1, -1, -1]).astype(float)


def basis():
    T, lab = [], []
    for k in range(2, N + 1):
        for j in range(1, k):
            s = np.zeros((N, N), complex); s[j-1, k-1] = s[k-1, j-1] = 0.5
            T.append(s); lab.append("sym(%d,%d)" % (j, k))
            a = np.zeros((N, N), complex); a[j-1, k-1] = -0.5j; a[k-1, j-1] = 0.5j
            T.append(a); lab.append("antisym(%d,%d)" % (j, k))
        T.append(None); lab.append("diag %d" % k)

    def dg(v, c):
        M = np.zeros((N, N), complex)
        for i, x in enumerate(v):
            M[i, i] = x * c
        return M
    T[2] = dg([1, -1, 0, 0, 0, 0, 0], 0.5)
    T[7] = dg([1, 1, -2, 0, 0, 0, 0], 1/(2*math.sqrt(3)))
    T[14] = dg([0, 0, 0, 1, -1, 0, 0], 0.5)
    T[23] = dg([0, 0, 0, 1, 1, -2, 0], 1/(2*math.sqrt(3)))
    T[34] = dg([0, 0, 0, 1, 1, 1, -3], 1/(2*math.sqrt(6)))
    T[47] = dg([4, 4, 4, -3, -3, -3, -3], 1/math.sqrt(168))
    return T, lab


T, LAB = basis()
f = lambda a, b, c: (-2j * np.trace((T[a] @ T[b] - T[b] @ T[a]) @ T[c])).real
H44 = 43

# ---------------------------------------------------------------- the operators
M = np.array([[f(a, H44, c) for c in range(DIM)] for a in range(DIM)])   # eq. (53)


def conj_matrix(Pm):
    """the matrix of  T -> Pm T Pm  in the T^a basis (2 Tr(T^a Pm T^b Pm))."""
    return np.array([[2 * np.trace(T[a] @ Pm @ T[b] @ Pm).real
                      for b in range(DIM)] for a in range(DIM)])


C55 = conj_matrix(P5V @ P5P)          # P5 P5' acts as one matrix
C6 = conj_matrix(P6V)

P("=" * 78)
P("A -- ARE THE THREE QUANTUM NUMBERS SIMULTANEOUS?  (they must be, or eq. (68)")
P("    cannot be a sum over them at all)")
P("=" * 78)
for nm, X, Y in (("[M, P5P5']", M, C55), ("[M, P6]", M, C6),
                 ("[P5P5', P6]", C55, C6)):
    P("  %-14s max|.| = %.2e" % (nm, np.max(np.abs(X @ Y - Y @ X))))
    assert np.max(np.abs(X @ Y - Y @ X)) < 1e-10
P("  all three commute.  PASS -- and this IS the pi_5 = pi_7 precondition:")
P("  P5P5' on indices 5,7 = (%+d, %+d), P6 on indices 5,7 = (%+d, %+d)"
  % ((P5V @ P5P)[4, 4], (P5V @ P5P)[6, 6], P6V[4, 4], P6V[6, 6]))
P("  had either pair differed, the charge eigenbasis would not be a parity")
P("  eigenbasis and their eq. (68) could not be written this way.")

# ------------------------------------------------- simultaneous classification
P("")
P("=" * 78)
P("B -- ALL 48 ADJOINT GENERATORS, CLASSIFIED")
P("=" * 78)
# 1j*M is hermitian, so eigh -- eig would return a non-orthonormal basis and the
# projections below would silently lose states (it lost six the first time)
w, v = np.linalg.eigh(1j * M)
w = np.real(w)
# refine within each charge eigenspace to diagonalise the two parities as well
table = {}
for q0 in sorted(set(np.round(w, 9))):
    idx = [k for k in range(DIM) if abs(w[k] - q0) < 1e-9]
    B = v[:, idx]
    # simultaneous diagonalisation: P5P5' first, then P6 inside each block
    for Op, other in ((C55, C6),):
        Bre = np.real(B) if np.max(np.abs(np.imag(B))) < 1e-9 else B
        A1 = Bre.conj().T @ Op @ Bre
        e1, U1 = np.linalg.eigh((A1 + A1.conj().T) / 2)
        B = Bre @ U1
        for p1 in (-1, 1):
            jj = [k for k in range(B.shape[1]) if abs(e1[k] - p1) < 1e-6]
            if not jj:
                continue
            B2 = B[:, jj]
            A2 = B2.conj().T @ other @ B2
            e2, U2 = np.linalg.eigh((A2 + A2.conj().T) / 2)
            for p2 in (-1, 1):
                n = sum(1 for k in range(len(e2)) if abs(e2[k] - p2) < 1e-6)
                if n:
                    table[(round(q0, 4), p1, p2)] = table.get((round(q0, 4), p1, p2), 0) + n
P("  %-8s %-10s %-8s %s" % ("q", "P5P5'", "P6", "states"))
for k in sorted(table):
    P("  %-8.1f %+-10d %+-8d %d" % (k[0], k[1], k[2], table[k]))
P("")
P("  total %d = dim adjoint : %s" % (sum(table.values()), sum(table.values()) == DIM))
assert sum(table.values()) == DIM
charged = {k: n for k, n in table.items() if abs(k[0]) > 1e-9}
P("  charged states %d, neutral %d  -- against su7_kk_spectrum.py's 22 and 26"
  % (sum(charged.values()), DIM - sum(charged.values())))
assert sum(charged.values()) == 22

# ------------------------------------------------------------- build eq. (68)
P("")
P("=" * 78)
P("C -- WHAT THEIR EQ. (68) SHOULD BE, FROM THAT TABLE ALONE")
P("=" * 78)
P("  The table in B is group theory and cannot move.  What IS a convention is")
P("  the weight each component carries into the potential: their eqs. (63)-(67)")
P("  close to (3k/4(pi n)^5)[P6 + g(pi k n)], so P6 = +1 and P6 = -1 components")
P("  enter with DIFFERENT constants, and a 6D gauge field splits into A_mu, A_5")
P("  and A_6 whose own P6 is not the generator's.  So do not assume the weights:")
P("  read them off, and let the redundancy be the test.")
P("")
P("  Each term of eq. (68) is a PAIR (q, -q), and c = 2|q|.  Write a for the")
P("  weight of a P6 = +1 pair and b for a P6 = -1 pair.  Three coefficients,")
P("  two unknowns -- OVERDETERMINED, so this can fail.")
P("")
P("  %-6s %-9s %-11s %-11s %s" % ("c", "s=P5P5'", "pairs P6=+", "pairs P6=-",
                                  "theirs"))
THEIRS = {(2, +1): 2, (1, +1): 4, (1, -1): 7}
rowsC = {}
for c in (1, 2):
    for s in (+1, -1):
        pp = table.get((c / 2.0, s, +1), 0) + table.get((-c / 2.0, s, +1), 0)
        pm = table.get((c / 2.0, s, -1), 0) + table.get((-c / 2.0, s, -1), 0)
        pp, pm = pp // 2, pm // 2
        if pp or pm:
            rowsC[(c, s)] = (pp, pm)
            P("  %-6d %+-9d %-11d %-11d %s" % (c, s, pp, pm, THEIRS.get((c, s), "--")))
P("")
from fractions import Fraction as Fr
# a is fixed by the two rows with no P6 = -1 pairs; b by the remaining one
solo = [(k, v) for k, v in rowsC.items() if v[1] == 0]
avals = {Fr(THEIRS[k], v[0]) for k, v in solo}
P("  rows with NO P6 = -1 pairs fix a on their own : %s"
  % sorted(str(x) for x in avals))
P("  CONSISTENCY, and it could have failed: those rows agree on a  : %s"
  % (len(avals) == 1))
assert len(avals) == 1
a = avals.pop()
mixed = [(k, v) for k, v in rowsC.items() if v[1]]
k, (pp, pm) = mixed[0]
b = (Fr(THEIRS[k]) - a * pp) / pm
P("  the remaining row then fixes b                 : %s" % b)
P("")
built = {k: a * v[0] + b * v[1] for k, v in rowsC.items()}
P("  their eq. (68)                    : %s" % {kk: THEIRS[kk] for kk in sorted(THEIRS)})
P("  rebuilt from the table with (a,b) : %s"
  % {kk: str(built[kk]) for kk in sorted(built)})
ok = all(built[kk] == THEIRS[kk] for kk in THEIRS) and set(built) == set(THEIRS)
P("")
if ok:
    P("  >> THEIR EQ. (68) IS REPRODUCED, three coefficients from two weights.")
    P("     a = %s per P6 = +1 pair, b = %s per P6 = -1 pair -- a ratio of %s."
      % (a, b, a / b))
    P("     And it is their own sentence: the %d at (c=1, s=-1) is %d x %s from"
      % (THEIRS[(1, -1)], rowsC[(1, -1)][0], a))
    P("     the P6 = +1 pairs plus %d x %s from the P6 = -1 ones, i.e. 4 + 3 --"
      % (rowsC[(1, -1)][1], b))
    P("     the split their text states and never derives.")
    P("")
    P("  >> So the GAUGE sector is now verified against the Wilson line and the")
    P("     orbifold parities, on the same footing as eqs. (73)-(76).  The one")
    P("     coefficient the repair wedge holds fixed is no longer the one thing")
    P("     in the potential that was only checked against its own upstream")
    P("     equation.")
else:
    P("  >> THEY DIFFER, and every entry is a coefficient the wedge holds FIXED.")
    for kk in sorted(set(THEIRS) | set(built)):
        P("       c=%d s=%+d :  theirs %s   built %s"
          % (kk[0], kk[1], THEIRS.get(kk, 0), built.get(kk, 0)))
