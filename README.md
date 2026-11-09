# The Fourier entropy-influence inequality for symmetric Boolean functions

Yitzchak Shmalo

Every symmetric Boolean function f satisfies H[f^2] <= K I[f] with K = 12 - (9/2) log2 3 = 4.8677, where H is the
Fourier entropy in bits and I is the total influence. The previous constant for symmetric functions is 12.04
(O'Donnell, Wright and Zhou, 2011). The AND function has H/I = 4 - (4 - 2/ln 2)/n + O(2^-n), so the best constant for
symmetric functions lies between 4 and 4.8677.

The proof is an induction on one coordinate. Restricting a coordinate of a symmetric function costs at most K times
the influence of that coordinate, in every dimension, with equality for OR on three variables, so K is the best
constant this induction can give. The main tool is a noise-stability bound for the derivative of a symmetric function,
Stab_lambda[g] <= e^-v I_0(v) ||g||^2. For large n the restriction step costs at most 4.637 times the influence, and
the paper also gives a bound with constant 6.191 that uses no computation.

The restriction step is computer-assisted: every symmetric function with n <= 15 is checked in exact arithmetic,
dimensions 16 to 9999 are reduced to finitely many one-variable inequalities with rational multipliers checked in
outward-rounded interval arithmetic, and an explicit bound covers n >= 10^4.

Version 1 (8 November 2026) had the constant 6.1908 and version 2 had 6.1347.

- `paper/symmetric_fei.pdf`, `paper/symmetric_fei.tex`: the paper.
- `code/restriction/run_all.sh`: reruns the certificates for the restriction step (about a minute):
  `small_dimensions.py` (n <= 15, exact), `certificate.py` (16 <= n <= 9999, multipliers in `duals.json` and
  `block_*.json`), `large_cutoff_certificate.py` (n >= 10^4) and `audit.py` (an independent recomputation at 75 digits).
- `code/restriction_majorant_certificate.py`: the interval certificate for the majorant used for large n (output
  `majorant_certificate.json`).
- `code/check_all.py`: the other constants and finite computations in the paper, including every symmetric function
  with n <= 16 (output in `check_all.out`).

The scripts need Python 3 with mpmath; `check_all.py` also needs numpy and scipy.

The paper is part of the Math Brain project: https://github.com/yspennstate/math-brain-public
