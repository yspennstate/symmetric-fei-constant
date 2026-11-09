"""Exhaustive symmetric restriction-step certificate, using Python integers only.

No third-party imports. Every comparison in the proof uses integer arithmetic.
Table mask bit w=1 means F(w)=-1; strings are printed weight 0 first.
Run: python3 small_dimensions.py --max-n 15
"""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path
from time import perf_counter

BITS = 144
Q = 1 << BITS
TERMS = 56
# R <= 9/[4(2L+1)3^(2L+1)] for 2 atanh(z), 0<=z<=1/3.
assert 9 * Q < 4 * (2 * TERMS + 1) * 3 ** (2 * TERMS + 1)


def ln_ratio(a, b):
    """Enclose ln((b+a)/(b-a)) in units 1/Q; 0<=a/b<=1/3."""
    assert 0 <= 3 * a <= b
    if a == 0:
        return 0, 0
    an, bn = a, b
    a2, b2 = a * a, b * b
    total = 0
    for j in range(TERMS):
        total += (Q * an) // (bn * (2 * j + 1))
        an *= a2
        bn *= b2
    # Each downward-rounded term loses <1/Q; the tail loses <1/Q.
    return 2 * total, 2 * total + 2 * TERMS + 1


LN2_LO, LN2_HI = ln_ratio(1, 3)


@lru_cache(None)
def log2_bounds(k):
    assert k >= 1
    b = k.bit_length() - 1
    p = 1 << b
    if k == p:
        return b * Q, b * Q
    lo, hi = ln_ratio(k - p, k + p)
    return b * Q + Q * lo // LN2_HI, b * Q + (Q * hi + LN2_LO - 1) // LN2_LO


@lru_cache(None)
def entropy_atom(m, a):
    """Bounds for Q*4^m*phi(a/2^m)."""
    a = abs(a)
    if not a:
        return 0, 0
    lo, hi = log2_bounds(a)
    return a * a * (2 * m * Q - 2 * hi), a * a * (2 * m * Q - 2 * lo)


def kraw(m):
    return [[sum((-1) ** j * comb(s, j) * comb(m - s, w - j)
                 for j in range(max(0, w - (m - s)), min(s, w) + 1))
             for w in range(m + 1)] for s in range(m + 1)]


def transform(table, matrix):
    return [sum(f * k for f, k in zip(table, row)) for row in matrix]


def ent_bounds(m, spectrum, multiplicities):
    lo = hi = 0
    for a, c in zip(spectrum, multiplicities):
        l, u = entropy_atom(m, a)
        lo += c * l
        hi += c * u
    return lo, hi


def ratio_bounds(n, A, B, C, cn, cm):
    j = sum(comb(n - 1, s - 1) * A[s] ** 2 for s in range(1, n + 1))
    if j == 0:
        return None
    lo, hi = ent_bounds(n, A, cn)
    l0, h0 = ent_bounds(n - 1, B, cm)
    l1, h1 = ent_bounds(n - 1, C, cm)
    return lo - 2 * (h0 + h1), hi - 2 * (l0 + l1), j


@lru_cache(None)
def factors(k):
    out = []
    p = 2
    while p * p <= k:
        v = 0
        while k % p == 0:
            k //= p
            v += 1
        if v:
            out.append((p, v))
        p += 1
    if k > 1:
        out.append((k, 1))
    return tuple(out)


def signature(n, A, B, C, cn, cm, j):
    # log2(2)=1 is absorbed into coefficient 2.
    coeff = {2: 2 * 4 ** n}
    for spectrum, mult, scale in ((A, cn, -2), (B, cm, 4), (C, cm, 4)):
        for a, c in zip(spectrum, mult):
            for p, v in factors(abs(a)) if a else ():
                coeff[p] = coeff.get(p, 0) + scale * c * a * a * v
    return tuple((p, Fraction(v, j)) for p, v in sorted(coeff.items()) if v)


REFERENCE = {
    1: '10', 2: '100', 3: '1000', 4: '10000', 5: '100000',
    6: '1000000', 7: '10000110', 8: '100000000', 9: '1000000000',
    10: '10000000000', 11: '011111111000', 12: '1000011111100',
    13: '10000000000000', 14: '110000111111000',
    15: '1000000111100000', 16: '11100001111110000',
}


def table_string(mask, n):
    return ''.join(str((mask >> w) & 1) for w in range(n + 1))


def decimal_bound(num, den, digits=18, upper=False):
    z = num * 10 ** digits
    v = -((-z) // den) if upper else z // den
    sign = '-' if v < 0 else ''
    v = abs(v)
    return f'{sign}{v // 10 ** digits}.{v % 10 ** digits:0{digits}d}'


def certify_n(n):
    start = perf_counter()
    mn, mm = kraw(n), kraw(n - 1)
    cn, cm = [comb(n, s) for s in range(n + 1)], [comb(n - 1, s) for s in range(n)]
    F = [1 - 2 * int(c) for c in REFERENCE[n]]
    ra, rb, rc = transform(F, mn), transform(F[:-1], mm), transform(F[1:], mm)
    rl, rh, rj = ratio_bounds(n, ra, rb, rc, cn, cm)
    rsig = signature(n, ra, rb, rc, cn, cm, rj)
    ties = []
    excluded = constants = 0
    runner = None
    max_width = Fraction(0)
    # Start at all +1; Gray transitions change exactly one weight-table entry.
    F = [1] * (n + 1)
    A, B, C = transform(F, mn), transform(F[:-1], mm), transform(F[1:], mm)
    previous = 0
    for counter in range(1 << (n + 1)):
        mask = counter ^ (counter >> 1)
        if counter:
            w = (mask ^ previous).bit_length() - 1
            change = -2 * F[w]
            F[w] = -F[w]
            for s in range(n + 1):
                A[s] += change * mn[s][w]
            for s in range(n):
                if w < n:
                    B[s] += change * mm[s][w]
                if w > 0:
                    C[s] += change * mm[s][w - 1]
        previous = mask
        bounds = ratio_bounds(n, A, B, C, cn, cm)
        if bounds is None:
            constants += 1
            assert all(x == F[0] for x in F)
            continue
        lo, hi, j = bounds
        # Independent physical influence: probability of adjacent-table sign change.
        edge_mass = sum(cm[w] for w in range(n) if F[w] != F[w + 1])
        assert j == edge_mass * (1 << (n + 1))
        assert sum(c * a * a for c, a in zip(cn, A)) == 4 ** n
        assert sum(c * b * b for c, b in zip(cm, B)) == 4 ** (n - 1)
        assert sum(c * c0 * c0 for c, c0 in zip(cm, C)) == 4 ** (n - 1)
        width = Fraction(hi - lo, j)
        if width > max_width:
            max_width = width
        if hi * rj < rl * j:
            excluded += 1
            if runner is None or hi * runner[1] > runner[0] * j:
                runner = hi, j, mask
        else:
            sig = signature(n, A, B, C, cn, cm, j)
            assert sig == rsig, ('uncertified competitor', n, table_string(mask, n))
            ties.append(table_string(mask, n))
    assert excluded + len(ties) + constants == 1 << (n + 1)
    result = {
        'n': n, 'tables': 1 << (n + 1), 'constant_tables': constants,
        'strictly_excluded_tables': excluded, 'maximizers': sorted(ties),
        'reference': REFERENCE[n], 'influence': str(Fraction(rj, 4 ** n)),
        'ratio_lower_numerator': str(rl), 'ratio_upper_numerator': str(rh),
        'ratio_denominator': str(Q * rj),
        'ratio_decimal_enclosure': [decimal_bound(rl, Q * rj), decimal_bound(rh, Q * rj, upper=True)],
        'exact_ratio_prime_log_coefficients': {str(p): str(v) for p, v in rsig},
        'maximum_interval_width_upper': decimal_bound(max_width.numerator, Q * max_width.denominator, 36, True),
        'seconds': round(perf_counter() - start, 3),
    }
    if runner:
        u, j, mask = runner
        gap = Fraction(rl, rj) - Fraction(u, j)
        assert gap > 0
        result['all_nonmaximizers_ratio_upper'] = decimal_bound(u, Q * j, upper=True)
        result['strict_gap_lower'] = decimal_bound(gap.numerator, Q * gap.denominator)
        result['largest_upper_competitor'] = table_string(mask, n)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=16)
    parser.add_argument('--output', default='certificate.json')
    args = parser.parse_args()
    assert 1 <= args.max_n <= 16
    rows = []
    for n in range(1, args.max_n + 1):
        row = certify_n(n)
        rows.append(row)
        print(json.dumps(row), flush=True)
        Path(args.output).write_text(json.dumps({
            'status': 'proved by exhaustive integer certificate',
            'bits': BITS, 'atanh_terms': TERMS,
            'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'rows': rows,
        }, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
