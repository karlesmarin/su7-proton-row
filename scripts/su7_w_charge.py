#!/usr/bin/env python
"""
Authors: Carles Marin + Claude (AI assistant).

IS THEIR EQ. (81) INHERITED, OR IS IT FORCED BY THEIR OWN EQ. (54)?

Komori-Maru introduce the W mass with the words "Noting that the W boson mass
in 5D SU(4) GHU is given by" -- i.e. eq. (81), m_W^2 = alpha^2 / (4 R5^2), is
quoted from a DIFFERENT model.  Their eq. (82), the whole third column of their
Table 1, and the constant in our own invariant K all rest on it.

An inherited relation carries a precondition: that the W's Wilson-line charge is
1/2 in THEIR SU(7) too.  Their eq. (55) already implies the charges on the
adjoint are 1/2 and 1 -- su7_kk_spectrum.py derives {0, +-1/2, +-1} with
multiplicities 26, 20, 2 -- but nothing there says WHICH of them the W is, and
the 2 modes at charge 1 would give m_W = alpha/R5, off by a factor 2.

This decides it from their eq. (54) alone.

  SU(2)_L acts on indices 4,5 (their eq. (41): the (1,2)_{1/2} sits there, and
  T^15 = (1/2)diag(0,0,0,1,-1,0,0) is its T3).  So W^+- are the raising and
  lowering combinations of sym(4,5) and antisym(4,5).  Their eq. (53) makes the
  Wilson-line shift operator M^{ac} = f^{a,44,c} on the adjoint, and a mode of
  charge q gets KK mass (n + q alpha)/R5.  Diagonalise M and read q off the W.
"""
import math
import numpy as np

P = lambda *a: print(*a, flush=True)
N, DIM = 7, 48


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
    lab[2] = lab[7] = "diag SU(3)_C"
    lab[14] = lab[23] = lab[34] = "diag SU(4)"
    lab[47] = "diag U(1)"
    return T, lab


T, LAB = basis()
f = lambda a, b, c: (-2j * np.trace((T[a] @ T[b] - T[b] @ T[a]) @ T[c])).real
H44 = 43                                            # T^44 = sym(5,7), 0-based

iS = LAB.index("sym(4,5)")
iA = LAB.index("antisym(4,5)")
P("=" * 78)
P("A -- SU(2)_L, located from their eq. (41)")
P("=" * 78)
P("  their eq. (41) puts the SM doublet on indices 4,5, and its T3 is")
P("  T^15 = (1/2) diag(0,0,0,1,-1,0,0):")
P("     T[14] diagonal = %s" % np.real(np.diag(T[14])))
assert np.allclose(np.diag(T[14]), [0, 0, 0, .5, -.5, 0, 0])
P("  so W^+- are the raising/lowering pair built from")
P("     %s  (index %d)  and  %s  (index %d)" % (LAB[iS], iS, LAB[iA], iA))
P("")
P("  CONTROL, and it can fail -- the three of them must close on su(2):")
com = lambda a, b: T[a] @ T[b] - T[b] @ T[a]
c1 = np.max(np.abs(com(iS, iA) - 1j * T[14]))
c2 = np.max(np.abs(com(14, iS) - 1j * T[iA]))
P("     [sym, antisym] - i T^15      max|.| = %.2e" % c1)
P("     [T^15, sym]    - i antisym   max|.| = %.2e" % c2)
assert c1 < 1e-12 and c2 < 1e-12
P("     PASS -- they are an su(2), so this really is SU(2)_L.")

P("")
P("=" * 78)
P("B -- THE WILSON-LINE CHARGE OF THE W, FROM THEIR EQ. (54)")
P("=" * 78)
M = np.array([[f(a, H44, c) for c in range(DIM)] for a in range(DIM)])
P("  M^{ac} = f^{a,44,c} antisymmetric : %s" % np.allclose(M, -M.T))
w, v = np.linalg.eig(1j * M)                        # i M hermitian -> real w
w = np.real(w)
P("")
P("  eigenvectors with support on the W directions:")
found = []
for k in range(DIM):
    amp = abs(v[iS, k]) ** 2 + abs(v[iA, k]) ** 2
    if amp > 1e-9:
        found.append((w[k], amp))
found.sort()
for q, amp in found:
    P("     q = %+.6f     weight on (sym(4,5), antisym(4,5)) = %.4f" % (q, amp))
qs = sorted({round(abs(q), 9) for q, _ in found})
P("")
P("  >> |q| carried by the W : %s" % qs)
assert qs == [0.5], qs
P("     The W sits at Wilson-line charge 1/2, NOT 1.  Its lightest KK mass is")
P("     therefore (0 + alpha/2)/R5 = alpha/(2 R5), i.e.")
P("")
P("        m_W = alpha / (2 R5)      -- THEIR EQ. (81), derived in THEIR model")
P("")
P("  Controls that had to fire, or the statement is empty:")
P("     the 2 modes at |q| = 1 exist   : %d modes"
  % sum(1 for x in w if abs(abs(x) - 1) < 1e-9))
P("     had the W been one of them, m_W would be alpha/R5 and their eq. (82)")
P("     would carry a factor 2.  The charge decides it, and it is 1/2.")
P("     the photon must be neutral: T^15 + sqrt3 T^24 is their eq. (78)")
u1em = np.zeros(DIM); u1em[14] = 1.0; u1em[23] = math.sqrt(3)
P("        |M . (T^15 + sqrt3 T^24)| = %.2e   -> uncharged, so it stays massless"
  % np.max(np.abs(M @ u1em)))
assert np.max(np.abs(M @ u1em)) < 1e-12

P("")
P("=" * 78)
P("WHAT THIS CHANGES")
P("=" * 78)
P("  Their eq. (81) is introduced as the 5D SU(4) result.  It is also a theorem")
P("  of their own SU(7): the Wilson line of their eq. (77) gives the W charge")
P("  1/2 and the photon charge 0, both forced by their eq. (54).  So the third")
P("  column of their Table 1, their eq. (82), and the constant 2.245624 in the")
P("  invariant K rest on a relation that is theirs and not merely inherited.")
