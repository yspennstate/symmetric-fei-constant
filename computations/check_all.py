"""Checks for every number stated in the paper. Run: python3 check_all.py [--nmax 14]

Prints one line per check, ending in PASS or FAIL, and exits nonzero if any check fails.
Needs numpy, scipy and mpmath. Runtime about a minute at the default nmax of 16.
"""
import itertools
import math
import sys
from math import comb, log, log2

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import i0e

mp.mp.dps = 40
NMAX = int(sys.argv[sys.argv.index("--nmax") + 1]) if "--nmax" in sys.argv else 16
fails = []


def report(name, ok, detail=""):
    print("%-58s %s  %s" % (name, "PASS" if ok else "FAIL", detail))
    if not ok:
        fails.append(name)


# 1. The constants.
C_L = mp.log(4 + mp.e, 2)  # Lemma 5.1
C_B = 2 + 1 / mp.log(2)
C = C_L + C_B
K = 12 - mp.mpf(9) / 2 * mp.log(3, 2)
report("log2(4+e) = 2.748092...", abs(C_L - mp.mpf("2.7480923172")) < 1e-9, mp.nstr(C_L, 20))
report("2 + 1/ln 2 = 3.442695...", abs(C_B - mp.mpf("3.4426950409")) < 1e-9, mp.nstr(C_B, 20))
report("Theorem E constant log2(4+e) + 2 + 1/ln 2 < 6.191", C < mp.mpf("6.191"), mp.nstr(C, 20))
report("K = 12 - (9/2) log2 3 = 4.8676687467...", abs(K - mp.mpf("4.86766874675479718345817")) < 1e-20, mp.nstr(K, 25))


def h_pair(t):
    def pl(y):
        return mp.mpf(0) if y == 0 else y * y * mp.log(y * y)
    s_ = mp.sqrt(t)
    return (pl(s_ + 1) + pl(s_ - 1)) / 2 - pl(s_)


k_half = mp.diff(h_pair, 1 + mp.mpf("1e-25"))
report("kappa(1/2) = h'(1) = 2 ln 2", abs(k_half - 2 * mp.log(2)) < 1e-15, mp.nstr(k_half, 20))
gam_up = mp.fsum(mp.mpf(1) / k for k in range(1, 11)) - mp.log(mp.mpf(21) / 2)
report("gamma < H_10 - ln(21/2) < 0.578", mp.euler < gam_up < mp.mpf("0.578"), mp.nstr(gam_up, 12))


# 2. C_N < ln 4 for the envelope integral, with F(x) = exp(-x) I_0(x).
def C_env(N):
    g = lambda t: (i0e(N * t * (2 - t) / 4) - (1 - t) ** N) / t
    pts = sorted({min(0.5, x / N) for x in (0.1, 1, 10, 100)})
    return quad(g, 0, 1, points=pts, limit=400, epsabs=1e-13, epsrel=1e-12)[0]


Ns = list(range(1, 201)) + [500, 1000, 10 ** 4, 10 ** 5, 10 ** 6]
vals = {N: C_env(N) for N in Ns}
worst = max(vals.values())
report("C_N < ln 4 for N = 1..200 and up to 10^6 (quadrature)", worst < log(4),
       "max %.6f at N=%d, ln 4 = %.6f" % (worst, max(vals, key=vals.get), log(4)))
report("C_1, C_2, C_3 below the analytic bounds 1, 59/64, ln 4",
       vals[1] <= 1 and vals[2] <= 59 / 64 and vals[3] < log(4),
       "%.4f %.4f %.4f" % (vals[1], vals[2], vals[3]))
bnd = {N: log(4) - (2 - math.sqrt(2)) / math.sqrt(math.pi * N) + 1 / (2 * N) for N in Ns if N >= 3}
report("C_N <= ln4 - (2-sqrt2)/sqrt(pi N) + 1/(2N), N >= 3", all(vals[N] <= bnd[N] + 1e-10 for N in bnd))


# 3. Bessel facts: sqrt(w) F(w) > 1/sqrt(2 pi) on w >= 1/2, and G(X) < 1/sqrt(2 pi X).
ws = np.concatenate([np.linspace(0.5, 10, 20000), np.geomspace(10, 1e9, 20000)])
h = np.sqrt(ws) * i0e(ws)
report("sqrt(w) e^-w I_0(w) > 1/sqrt(2pi) on [1/2, 1e9]", bool(np.all(h > 1 / math.sqrt(2 * math.pi))),
       "min %.12f" % h.min())
from scipy.special import i1e
Xs = np.geomspace(1e-6, 1e9, 40000)
report("e^-X I_1(X) < 1/sqrt(2 pi X) on [1e-6, 1e9]", bool(np.all(i1e(Xs) < 1 / np.sqrt(2 * math.pi * Xs))))
Xc = mp.mpf(7)
J = mp.quad(lambda u: u * mp.e ** (-u) * (mp.besseli(0, u) - mp.besseli(1, u)), [0, Xc])
report("J(X) = X e^-X I_1(X) at X = 7", abs(J - Xc * mp.e ** (-Xc) * mp.besseli(1, Xc)) < 1e-30)


# 4. Exhaustive symmetric functions: H, L, B, I and the restriction step.
def kraw(n):
    return [[sum((-1) ** j * comb(w, j) * comb(n - w, k - j) for j in range(k + 1)) for w in range(n + 1)]
            for k in range(n + 1)]


def spectrum(prof, K):
    n = len(prof) - 1
    c = [sum(K[w][k] * prof[w] for w in range(n + 1)) for k in range(n + 1)]  # slice sum of chi_S is K_w(|S|)
    return [ck * ck / 4.0 ** n for ck in c]  # squared coefficient at each level


def ent(prof, K):
    n = len(prof) - 1
    sq = spectrum(prof, K)
    W = [comb(n, k) * sq[k] for k in range(n + 1)]
    H = -sum(comb(n, k) * s * log2(s) for k, s in enumerate(sq) if s > 0)
    L = -sum(x * log2(x) for x in W if x > 0)
    B = sum(W[k] * log2(comb(n, k)) for k in range(n + 1))
    I = sum(k * W[k] for k in range(n + 1))
    return H, L, B, I


worst_ratio, worst_L, worst_B, worst_step = (0, None), 0, 0, (0, None)
best_n = {}
for n in range(1, NMAX + 1):
    K = kraw(n)
    Km = kraw(n - 1) if n > 1 else None
    for bits in itertools.product((1, -1), repeat=n):
        prof = (1,) + bits  # the output sign is fixed by f(all ones) = 1
        if len(set(prof)) == 1:
            continue
        H, L, B, I = ent(prof, K)
        if H / I > worst_ratio[0]:
            worst_ratio = (H / I, (n, prof))
        best_n[n] = max(best_n.get(n, 0), H / I)
        worst_L = max(worst_L, L / I)
        worst_B = max(worst_B, B / I)
        if n >= 2:
            def Hs(p):
                if len(set(p)) == 1:
                    return 0.0
                return ent(p, Km)[0]
            I1 = sum(comb(n - 1, w) for w in range(n) if prof[w] != prof[w + 1]) / 2.0 ** (n - 1)
            d = H - (Hs(prof[:-1]) + Hs(prof[1:])) / 2
            if I1 > 0 and d / I1 > worst_step[0]:
                worst_step = (d / I1, (n, prof))
n0, p0 = worst_ratio[1]
report("each n = 2..%d: max H/I equals the AND/OR value" % NMAX,
       all(abs(best_n[n] - ent((1,) * n + (-1,), kraw(n))[0] / ent((1,) * n + (-1,), kraw(n))[3]) < 1e-12
           for n in range(2, NMAX + 1)))
report("all symmetric f, n <= %d: max H/I < 4, attained by AND/OR" % NMAX,
       worst_ratio[0] < 4 and (p0.count(-1) in (1, n0)), "max %.12f at n=%d" % (worst_ratio[0], n0))
report("all symmetric f, n <= %d: L <= log2(4+e) I" % NMAX, worst_L <= float(C_L), "max L/I %.6f" % worst_L)
report("all symmetric f, n <= %d: B <= (2+1/ln2) I" % NMAX, worst_B <= float(C_B), "max B/I %.6f" % worst_B)
or3 = 12 - 4.5 * log2(3)
report("restriction step: max Delta_1/I_1 = 12 - (9/2) log2 3 at n=3",
       abs(worst_step[0] - or3) < 1e-9 and worst_step[1][0] == 3,
       "max %.10f at n=%d, OR_3 %.10f" % (worst_step[0], worst_step[1][0], or3))


# 5. AND/OR: exact closed form and the expansion 4 - (4 - 2/ln 2)/n.
def andor(n):
    n = mp.mpf(n)
    W0 = (1 - 2 ** (1 - n)) ** 2
    H = -W0 * 2 * mp.log1p(-2 ** (1 - n)) / mp.log(2) + (2 ** n - 1) * 4 ** (1 - n) * 2 * (n - 1)
    I = n * 2 ** (1 - n)
    return H / I


ok = True
for n in range(2, NMAX + 1):
    prof = (1,) * n + (-1,)
    H, L, B, I = ent(prof, kraw(n))
    ok &= abs(H / I - float(andor(n))) < 1e-10
report("AND/OR closed form matches direct computation, n <= %d" % NMAX, ok)
ok = all(abs(andor(n) - (4 - (4 - 2 / mp.log(2)) / n)) < 2 ** (-n / 2) for n in range(20, 121))
report("AND/OR ratio = 4 - (4 - 2/ln2)/n + O(2^-n/2), n = 20..120", ok, "n=1024: %s" % mp.nstr(andor(1024), 15))


# 6. Stability lemma on random symmetric functions: sum_j lambda^j q_j <= e^-v I_0(v).
rng = np.random.default_rng(1)
ok, cnt = True, 0
for n in (20, 40, 80):
    N = n - 1
    K = kraw(N)
    for _ in range(40):
        D = [0] * (N + 1)
        sel = sorted(rng.choice(N + 1, size=int(rng.integers(1, N + 1)), replace=False))
        for i, w in enumerate(sel):
            D[w] = (-1) ** i
        c = [sum(K[w][k] * D[w] for w in range(N + 1)) for k in range(N + 1)]
        W = np.array([comb(N, k) * float(c[k]) ** 2 for k in range(N + 1)])
        q = W / W.sum()
        for lam in np.linspace(0, 0.999, 60):
            v = N * (1 - lam * lam) / 4
            ok &= float(np.dot(q, lam ** np.arange(N + 1))) <= i0e(v) + 1e-12
            cnt += 1
report("noise stability <= e^-v I_0(v), %d random checks at n = 20, 40, 80" % cnt, ok)


# 7. The asymptotic restriction constant C*.
B10 = mp.e ** (-10) * mp.besseli(0, 10)
z = 40 * (mp.mpf(23) / 48) ** 2
Cstar = (mp.mpf(14) / 5 + mp.mpf(89) / 50 * B10 + 48 / (23 * mp.sqrt(40 * mp.pi)) + mp.e ** (-z) / z) / mp.log(2)
report("C* = 4.636423666... < 4.637", Cstar < mp.mpf("4.637") and abs(Cstar - mp.mpf("4.636423666")) < 1e-8,
       mp.nstr(Cstar, 15))
report("e^-10 I_0(10) < 0.127833337164", B10 < mp.mpf("0.127833337164"), mp.nstr(B10, 15))

print("\n%d checks failed" % len(fails) if fails else "\nall checks passed")
sys.exit(1 if fails else 0)
