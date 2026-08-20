#!/usr/bin/env python3
"""The layer under the anchor: Komori-Maru's KK spectrum, derived instead of assumed.

  Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
  2026-08-04.  Source: Komori & Maru, arXiv:2503.04090v1.

WHY.  Everything checked so far about their section 3 was INTERNAL: eqs. (56),(58),(62)
mutually consistent, eq. (67) proved, eq. (68) rebuilt, Appendix C in agreement.  All of
it sits on eqs. (53)-(55) -- the mass matrix built from the SU(7) structure constants
f^{a,44,b} and its eigenvalues -- which we had taken on trust.  An error there would
propagate everywhere and none of those checks could see it: they all live inside the
building.  So: derive it.

THEIR BASIS, and the false start that found it.  The obvious guess is the nested chain
SU(2) < SU(3) < ... < SU(7), whose diagonal generators sit at 3, 8, 15, 24, 35, 48 --
which is where theirs sit.  It reproduces ELEVEN of the thirteen constants of their
eq. (54) exactly, and fails on the three diagonal ones.  It is the guess that is wrong,
not their paper: it does not satisfy their own eqs. (78) and (79), and it puts T^48 --
which they state commutes with the Wilson line -- on a generator that does not.

Their Cartan is the physical one, SU(3)_C x SU(4) x U(1) with the SU(4) on indices
4,5,6,7, and eqs. (78),(79) fix it completely:

    T^15 = (1/2)      diag(0,0,0, 1,-1, 0, 0)
    T^24 = (1/2sqrt3) diag(0,0,0, 1, 1,-2, 0)
    T^35 = (1/2sqrt6) diag(0,0,0, 1, 1, 1,-3)
    T^48 = (1/sqrt168)diag(4,4,4,-3,-3,-3,-3)      <- the U(1)

The off-diagonal generators coincide with the nested ordering: block k holds, for
j = 1..k-1, the pair sym(j,k), antisym(j,k).  Hence T^44 = sym(5,7) and T^45 =
antisym(5,7), which is exactly the direction their Wilson line (eq. 77) rotates.

Normalisation Tr(T^a T^b) = delta^{ab}/2,  f^{abc} = -2i Tr([T^a,T^b] T^c).
"""
import math
import sys

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
        T.append(None); lab.append("diag %d" % k)          # filled below
    def dg(v, c):
        M = np.zeros((N, N), complex)
        for i, x in enumerate(v):
            M[i, i] = x * c
        return M
    T[2]  = dg([1, -1, 0, 0, 0, 0, 0], 0.5)                        # T^3  SU(3)_C
    T[7]  = dg([1, 1, -2, 0, 0, 0, 0], 1/(2*math.sqrt(3)))         # T^8  SU(3)_C
    T[14] = dg([0, 0, 0, 1, -1, 0, 0], 0.5)                        # T^15 SU(4)
    T[23] = dg([0, 0, 0, 1, 1, -2, 0], 1/(2*math.sqrt(3)))         # T^24 SU(4)
    T[34] = dg([0, 0, 0, 1, 1, 1, -3], 1/(2*math.sqrt(6)))         # T^35 SU(4)
    T[47] = dg([4, 4, 4, -3, -3, -3, -3], 1/math.sqrt(168))        # T^48 U(1)
    lab[2], lab[7] = "diag SU(3)_C", "diag SU(3)_C"
    lab[14] = lab[23] = lab[34] = "diag SU(4)"
    lab[47] = "diag U(1)"
    return T, lab


T, LAB = basis()
f = lambda a, b, c: (-2j * np.trace((T[a] @ T[b] - T[b] @ T[a]) @ T[c])).real
H = 43                                                      # T^44, zero-based

P("=" * 78)
P("A -- THE BASIS, FIXED BY THEIR OWN EQS. (78) AND (79)")
P("=" * 78)
orth = max(abs(np.trace(T[a] @ T[b]) - (0.5 if a == b else 0.0))
           for a in range(DIM) for b in range(DIM))
P("  orthonormality  max |Tr(T^aT^b) - delta/2| = %.2e            %s"
  % (orth, "OK" if orth < 1e-12 else "*** FAIL ***"))
assert orth < 1e-12
e78 = np.real(np.diag(T[14] + math.sqrt(3) * T[23]))
e79 = np.real(np.diag(T[14] - math.sqrt(3)/3 * T[23] + math.sqrt(6)/3 * T[34]))
P("  their eq. (78)  T^15 + sqrt3 T^24              = %s" % np.round(e78, 6))
P("                  they print                      = [0 0 0  1  0 -1  0]   %s"
  % ("MATCH" if np.allclose(e78, [0,0,0,1,0,-1,0]) else "*** FAIL ***"))
P("  their eq. (79)  T^15 -(r3/3)T^24 +(r6/3)T^35   = %s" % np.round(e79, 6))
P("                  they print                      = [0 0 0 .5 -.5 .5 -.5] %s"
  % ("MATCH" if np.allclose(e79, [0,0,0,.5,-.5,.5,-.5]) else "*** FAIL ***"))
assert np.allclose(e78, [0,0,0,1,0,-1,0]) and np.allclose(e79, [0,0,0,.5,-.5,.5,-.5])
P("  T^44 = %s, T^45 = %s  -- the plane their eq. (77) rotates." % (LAB[43], LAB[44]))
P("  T^48 entries on 5 and 7: %+.4f, %+.4f  -> equal, so it commutes with the"
  % (np.real(T[47][4,4]), np.real(T[47][6,6])))
P("  Wilson line, exactly as they state below eq. (77).")

P("")
P("=" * 78)
P("B -- THEIR EQ. (54), EVERY CONSTANT")
P("=" * 78)
THEIRS = {(15,45): -0.5, (16,37): 0.5, (17,36): -0.5, (18,39): 0.5, (19,38): -0.5,
          (20,41): 0.5, (21,40): -0.5, (22,43): 0.5, (23,42): -0.5,
          (24,45): 1/(2*math.sqrt(3)), (33,47): 0.5, (34,46): 0.5,
          (35,45): math.sqrt(2/3)}
mine = {}
for a in range(DIM):
    for c in range(DIM):
        v = f(a, H, c)
        if abs(v) > 1e-10:
            mine[(a+1, c+1)] = v
P("  %-12s %-16s %-16s %-12s %-12s %s"
  % ("f^(a,44,b)", "T^a", "T^b", "computed", "theirs", ""))
bad = []
for (a, c), th in sorted(THEIRS.items()):
    v = mine.get((a, c))
    good = v is not None and abs(v - th) < 1e-10
    if not good:
        bad.append((a, c))
    P("  %-12s %-16s %-16s %-12.6f %-12.6f %s"
      % ("(%d,44,%d)" % (a, c), LAB[a-1], LAB[c-1],
         v if v is not None else float("nan"), th, "match" if good else "*** MISMATCH ***"))
extra = [k for k in mine if k not in THEIRS and (k[1], k[0]) not in THEIRS]
P("")
P("  their 13 constants reproduced      : %d/13" % (13 - len(bad)))
P("  constants they omit that are non-zero: %d %s"
  % (len(extra), extra if extra else "(none)"))
P("  f^{48,44,45} = %.1e  -> zero, which is why the U(1) is absent from their list."
  % f(47, H, 44))
P("  verdict: their eq. (54) is EXACT AND COMPLETE.                        %s"
  % ("PASS" if not bad and not extra else "*** PROBLEM ***"))
assert not bad and not extra

P("")
P("=" * 78)
P("C -- THE KK SPECTRUM OF EQ. (55)")
P("=" * 78)
P("  eq. (53) makes the shift operator on the adjoint  M^{ac} = f^{a,44,c}, and a mode")
P("  of charge q under it acquires KK mass (n + q alpha)/R5.  So diagonalise M.")
M = np.array([[f(a, H, c) for c in range(DIM)] for a in range(DIM)])
P("")
P("  M antisymmetric: %s   (so i M is hermitian and q is real)"
  % np.allclose(M, -M.T))
ev = np.linalg.eigvals(1j * M)
assert max(abs(ev.imag)) < 1e-9
cnt = {}
for x in np.round(ev.real, 6):
    cnt[abs(x) + 0.0] = cnt.get(abs(x) + 0.0, 0) + 1
P("")
P("  %-14s %s" % ("|q|", "modes"))
for x in sorted(cnt):
    P("  %-14.4f %d" % (x, cnt[x]))
P("  total %d = dim adjoint                                                 %s"
  % (sum(cnt.values()), "OK" if sum(cnt.values()) == DIM else "*** FAIL ***"))
P("")
P("  their eq. (55) contains exactly the shifts  alpha  and  alpha/2  and nothing")
P("  else, on towers n and ntilde+1/2.  Charges found: %s"
  % sorted(x for x in cnt if x > 1e-9))
ok = sorted(round(x, 6) for x in cnt if x > 1e-9) == [0.5, 1.0]
P("  >> the only non-zero shift charges are 1/2 and 1: %s" % ("PASS" if ok else "*** FAIL ***"))
assert ok
P("")
P("  So the input to their whole section 3 -- the set of Wilson-line charges carried")
P("  by the 48 gauge modes -- follows from their eq. (54) and nothing else, and it is")
P("  the set their eq. (55) uses.  The layer under the anchor is sound.")
