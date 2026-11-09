# The Fourier entropy-influence inequality for symmetric Boolean functions with constant 6.2

Yitzchak Shmalo

Every symmetric Boolean function f satisfies H[f^2] <= C I[f] with C = log2(4+e) + 2 + 1/ln 2 = 6.1908..., where H is
the Fourier entropy in bits and I is the total influence. The previous constant for symmetric functions is 12.04
(O'Donnell, Wright and Zhou, 2011). The AND function has H/I = 4 - (4 - 2/ln 2)/n + O(2^-n), so the best constant for
symmetric functions lies between 4 and 6.1908.

- `paper/symmetric_fei.pdf`, `paper/symmetric_fei.tex`: the paper.
- `code/check_scalar.py`: the Bessel inequalities and C_N < ln 4 at 30 digits (output in `check_scalar.out`).
- `code/check_tables.py`: every nonconstant symmetric function with n <= 14, checked against the level, within-level
  and total bounds (output in `check_tables.out`).

Both scripts need Python 3 and mpmath. `python code/check_tables.py 14` takes a few minutes.

The paper is part of the Math Brain project: https://github.com/yspennstate/math-brain-public
