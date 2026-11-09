# The Fourier entropy-influence inequality for symmetric Boolean functions

Yitzchak Shmalo

Every symmetric Boolean function f satisfies H[f^2] <= C I[f] with C = 673/250 + 2 + 1/ln 2 = 6.1347, where H is
the Fourier entropy in bits and I is the total influence. The previous constant for symmetric functions is 12.04
(O'Donnell, Wright and Zhou, 2011). The AND function has H/I = 4 - (4 - 2/ln 2)/n + O(2^-n), so the best constant for
symmetric functions lies between 4 and 6.1347.

The entropy splits into the entropy of the level weights, at most 2.692 I for every Boolean function, and a
within-level term, at most (2 + 1/ln 2) I for symmetric functions. The paper also shows that the entropy cost of
restricting one coordinate is asymptotically at most 4.637 times its influence.

Version 1 (8 November 2026) had the level bound log2(4+e) I and the constant 6.1908.

- `paper/symmetric_fei.pdf`, `paper/symmetric_fei.tex`: the paper.
- `code/check_all.py`: the constants and finite computations in the paper, including every symmetric function with
  n <= 16 (output in `check_all.out`; about a minute).
- `code/level_entropy_certificate.py`: the interval certificate for the level bound 2.692 (output
  `level_entropy_certificate.json`).
- `code/restriction_majorant_certificate.py`: the interval certificate for the majorant in the restriction estimate
  (output `majorant_certificate.json`).
- `code/check_scalar.py`, `code/check_tables.py`: the checks from version 1 (Bessel inequalities at 30 digits, and
  every symmetric function with n <= 14 against the level, within-level and total bounds).

The scripts need Python 3 with mpmath; `check_all.py` also needs numpy and scipy.

The paper is part of the Math Brain project: https://github.com/yspennstate/math-brain-public
