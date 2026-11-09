"""Outward-rounded certificate for the bound R(n) of Proposition 4.5 of the paper.

No numerical maximization over Boolean functions enters the proof. A positive
Bessel series with an explicit remainder supplies its two special-function calls.
"""
import hashlib
import json
import resource
import time
from pathlib import Path

import mpmath as mp

iv = mp.iv
iv.dps = 60


def rat(p, q=1):
    return iv.mpf(p) / q


def bessel_enclosure(v, terms=100):
    x = v * v / 4
    term = rat(1)
    total = term
    for j in range(1, terms + 1):
        term *= x / (j * j)
        total += term
    ratio = x / ((terms + 2) ** 2)
    assert ratio.b < rat(1).a
    remainder = term * x / ((terms + 1) ** 2) / (1 - ratio)
    low = iv.exp(-v) * total
    high = iv.exp(-v) * (total + remainder)
    return iv.mpf([low.a, high.b])


def bound(n):
    n = rat(n)
    N = n - 1
    a = rat(1, 5)
    r = rat(23, 48)
    c = (N / n) * (1 - iv.exp(-2 * a)) / (2 * a)
    vT = N * (1 - iv.exp(-2 * a)) / 4
    v20 = N * (1 - iv.exp(-40 / n)) / 4
    b = iv.sqrt(2 / (iv.pi * N))
    z = 40 * r * r
    A = -rat(1, 2) - iv.ln(a) + a
    K = rat(24, 23) / iv.sqrt(2 * iv.pi * vT) + iv.exp(-2 * r * r * vT)
    terms = {
        "level_zero_nats": b * (3 + iv.ln(1 + 1 / (b * b)) + 2 * iv.ln(2)),
        "tangent_shift_nats": (6 + 4 * iv.ln(n)) / n,
        "affine_base_nats": rat(14, 5),
        "finite_mean_correction_nats": rat(273, 100) / n,
        "J_upper_nats": rat(48, 23) / iv.sqrt(40 * iv.pi) + iv.exp(-z) / z,
        "laplace_rescaling_upper_nats": bessel_enclosure(20 * c) * iv.ln(1 / c),
        "laplace_tail_upper_nats": iv.exp(a) * A * K,
        "exponential_term_upper_nats": rat(89, 50) * bessel_enclosure(v20),
    }
    total = sum(terms.values(), rat(0)) / iv.ln(2)
    return total, terms, {"c": c, "vT": vT, "v20": v20, "E1_a_upper": A}


def main():
    start = time.monotonic()
    # H_6-ln 7 < gamma, so this proves the elementary gamma lower bound.
    gamma_lower = sum((rat(1, j) for j in range(1, 7)), rat(0)) - iv.ln(7)
    assert gamma_lower.a > rat(1, 2).b
    cutoff, terms, parameters = bound(10000)
    safe_decimal = rat(4866169, 1000000)
    target = rat(48723, 10000)
    assert cutoff.b < safe_decimal.a
    assert safe_decimal.b < target.a
    # A second, coarser advertised cutoff is optional corroboration only.
    check_16384, _, _ = bound(16384)
    assert check_16384.b < target.a
    out = {
        "status": "certified with outward-rounded intervals and analytic monotonicity",
        "n0_variables": 10000,
        "bound_bits_per_unit_coordinate_influence": str(cutoff),
        "safe_decimal_upper_bits": "4.866169",
        "target_bits": "4.8723",
        "certified_margin_bits": str(target - cutoff),
        "terms_nats": {k: str(v) for k, v in terms.items()},
        "parameters": {k: str(v) for k, v in parameters.items()},
        "gamma_lower": str(gamma_lower),
        "n16384_bound_bits": str(check_16384),
        "interval_decimal_digits": iv.dps,
        "bessel_series_nonconstant_terms": 100,
        "seconds": time.monotonic() - start,
        "peak_resident_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path("cutoff_certificate.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
