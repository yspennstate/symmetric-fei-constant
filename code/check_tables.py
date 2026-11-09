"""Check the level, within-level and total bounds on every nonconstant symmetric Boolean function with n <= NMAX.

A symmetric f on n bits is a table (f_0, ..., f_n) of signs, f_w being the value on inputs with w coordinates equal
to -1. Its level-k coefficient is 2^-n sum_w f_w K_k(w), with K_k the Krawtchouk polynomial, computed here in exact
integers. C_N is computed by mpmath quadrature. Prints the largest ratios found and the function attaining H/I.
"""
import itertools
import math
import sys

import mpmath as mp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
mp.mp.dps = 30


def C(N):
    F = lambda x: mp.besseli(0, x) * mp.e ** (-x)
    return mp.quad(lambda t: (F(N * t * (2 - t) / 4) - (1 - t) ** N) / t, [0, mp.mpf(1) / N, 1])


def kraw(n):
    # K[k][w] = sum over inputs x of weight w of chi_S(x), for any fixed S with |S| = k
    return [[sum((-1) ** j * math.comb(k, j) * math.comb(n - k, w - j) for j in range(min(k, w) + 1))
             for w in range(n + 1)] for k in range(n + 1)]


LEVEL = math.log2(4 + math.e)
TOTAL = LEVEL + 2 + 1 / math.log(2)
worst = {"L": 0.0, "B": 0.0, "H": 0.0}
for n in range(1, NMAX + 1):
    K = kraw(n)
    bmax = (1 + float(C(n - 1))) / math.log(2) if n >= 2 else 0.0
    best = (0.0, None)
    for tab in itertools.product((1, -1), repeat=n + 1):
        if len(set(tab)) == 1:
            continue
        c = [sum(tab[w] * K[k][w] for w in range(n + 1)) for k in range(n + 1)]  # 2^n * coefficient
        W = [math.comb(n, k) * c[k] ** 2 / 4 ** n for k in range(n + 1)]
        I = sum(k * W[k] for k in range(n + 1))
        L = sum(w * math.log2(1 / w) for w in W if w > 0)
        B = sum(W[k] * math.log2(math.comb(n, k)) for k in range(n + 1))
        assert abs(sum(W) - 1) < 1e-12
        assert L <= LEVEL * I + 1e-12, (n, tab)
        assert B <= bmax * I + 1e-12, (n, tab)
        assert L + B <= TOTAL * I + 1e-12, (n, tab)
        worst["L"] = max(worst["L"], L / I)
        worst["B"] = max(worst["B"], B / I)
        worst["H"] = max(worst["H"], (L + B) / I)
        if (L + B) / I > best[0] + 1e-12:
            best = ((L + B) / I, tab)
    and_tab = tuple([1] * n + [-1])
    print(f"n={n:2d}  max H/I={best[0]:.6f}  attained by AND/OR: {best[1] in (and_tab, tuple(-x for x in and_tab), and_tab[::-1], tuple(-x for x in and_tab[::-1]))}"
          f"  B bound (1+C_N)/ln2={bmax:.4f}", flush=True)
print(f"largest L/I={worst['L']:.6f} (bound {LEVEL:.6f}); largest B/I={worst['B']:.6f}; largest H/I={worst['H']:.6f} "
      f"(bound {TOTAL:.6f})")
print("all bounds PASS")
