"""Numerical checks of the scalar facts behind Proposition 5.1 (C_N < ln 4), at 30 digits.

F(x) = e^-x I_0(x), G(x) = e^-x I_1(x). Checks:
  (a) sqrt(w) F(w) > 1/sqrt(2 pi) on a log grid of w in [1/2, 5e9];
  (b) J(X) = int_0^X u (F - G) du equals X G(X);
  (c) G(X) sqrt(2 pi X) < 1;
  the identity T3 = H_N - ln N - gamma - E1(N);
  C_N by quadrature against the bound ln 4 - (2 - sqrt 2)/sqrt(pi N) + (H_N - ln N - gamma), for N <= 40 and larger N.
"""
import mpmath as mp

mp.mp.dps = 30
F = lambda x: mp.besseli(0, x) * mp.e ** (-x)
G = lambda x: mp.besseli(1, x) * mp.e ** (-x)
c = 1 / mp.sqrt(2 * mp.pi)

ws = [mp.mpf(1) / 2 * mp.mpf(10) ** (k / mp.mpf(20)) for k in range(0, 200)]
gap = min(mp.sqrt(w) * F(w) - c for w in ws)
assert gap > 0
print("(a) min of sqrt(w)F(w) - 1/sqrt(2pi) on the grid:", mp.nstr(gap, 6))

for X in [mp.mpf(1) / 4, 1, 3, 10, 100]:
    J = mp.quad(lambda u: u * (F(u) - G(u)), [0, X])
    assert abs(J - X * G(X)) < mp.mpf(10) ** -25
    assert G(X) * mp.sqrt(2 * mp.pi * X) < 1
print("(b) J(X) = X G(X) and (c) G(X) < 1/sqrt(2 pi X) at X = 1/4, 1, 3, 10, 100: PASS")


def C(N):
    return mp.quad(lambda t: (F(N * t * (2 - t) / 4) - (1 - t) ** N) / t, [0, mp.mpf(1) / N, 1])


for N in list(range(1, 41)) + [100, 1000, 10 ** 5]:
    H = mp.harmonic(N)
    T3 = mp.quad(lambda t: (mp.e ** (-N * t) - (1 - t) ** N) / t, [0, mp.mpf(1) / N, 1])
    assert abs(T3 - (H - mp.log(N) - mp.euler - mp.e1(N))) < mp.mpf(10) ** -20
    CN = C(N)
    bound = mp.log(4) - (2 - mp.sqrt(2)) / mp.sqrt(mp.pi * N) + (H - mp.log(N) - mp.euler)
    assert CN < mp.log(4)
    if N >= 3:
        assert CN < bound < mp.log(4)
    if N in (1, 2, 3, 10, 40, 100, 1000, 10 ** 5):
        print(f"N={N:6d}  C_N={mp.nstr(CN, 6)}  bound={mp.nstr(bound, 6)}  ln4={mp.nstr(mp.log(4), 6)}")
assert C(2) < mp.mpf(59) / 64
print("T3 identity, C_N < bound < ln 4 (N >= 3) and C_N < ln 4 (all N checked): PASS")
print("constant log2(4+e) + 2 + 1/ln 2 =", mp.nstr(mp.log(4 + mp.e, 2) + 2 + 1 / mp.log(2), 20))
