# The Fourier entropy–influence inequality for symmetric Boolean functions

Yitzchak Shmalo

[Paper (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) · [Declarations](DECLARATIONS.md)

The paper proves the Fourier entropy–influence inequality for every symmetric Boolean function with constant $931487/200000=4.657435$, improving the constant 12.04 of O'Donnell, Wright and Zhou (ICALP 2011). All logarithms in the statements are to base 2, $H(f)$ is the Fourier entropy and $I(f)$ the total influence.

## Main results

- **Theorem A.** Every symmetric Boolean function satisfies

  $$H(f)\le c\,I(f),\qquad c=\frac{931487}{200000}=4.657435.$$

- **Theorem B.** (i) Restricting one coordinate of a symmetric function costs at most $K=12-\tfrac92\log_2 3=4.8677$ times the influence of that coordinate: $\Delta_1(f)\le K\,I_1(f)$ in every dimension, with equality for $\mathrm{OR}_3$. (ii) In high dimension the step is cheaper: $\limsup_n\sup_f\Delta_1/I_1<4.637$.

- **Averaging.** If $U_d$ bounds $\Delta_1/I_1$ for symmetric functions in dimension $d$, then $H(f)\le\frac1n(U_1+\dots+U_n)\,I(f)$. A single uniform bound cannot beat $K$ because of $\mathrm{OR}_3$, but the average can. Theorem A comes from this average, seeded by a block-restriction census through dimension 20 and the dimension-specific bounds above it.

- **Theorem C.** If $f$ is symmetric on $n\ge2$ variables and $g$ is its derivative in one coordinate, then for $0\le\lambda\le1$

  $$\mathrm{Stab}_\lambda[g]\le e^{-v}I_0(v)\,\|g\|_2^2,\qquad v=\frac{(n-1)(1-\lambda^2)}{4}.$$

- **Theorem D.** Without any computation, $H(f)\le\bigl(\log_2(4+e)+2+1/\ln 2\bigr)I(f)<6.191\,I(f)$ for every symmetric $f$.

- **Theorem E.** For OR and AND, $H/I=4-(4-2/\ln 2)/n+O(2^{-n})$. The optimal constant for symmetric functions lies in $[4,\,4.657435]$.

The finite parts are computer-assisted. Every symmetric function is checked in exact arithmetic for $n\le15$. Dimensions $16$ to $255$ reduce to one-variable inequalities with rational multipliers, and dimensions $256$ to $8{,}388{,}608$ to $65{,}139$ rational cells in 23 blocks, all verified in outward-rounded interval arithmetic. An analytic bound decreasing in $n$ covers every larger dimension. The averaged ledger is summed in exact rational arithmetic; its largest value, $4.65743447\ldots$, is at $n=449$.

Version 1 (8 November 2026) had the constant 6.1908, version 2 had 6.1347 and version 3 had 4.8677.

## Computations

Python needs `mpmath`; the census needs a C++17 compiler. Each command uses one thread.

The one-step bounds, from `computations/restriction_v4`:

```
python3 certify.py --replay
python3 audit.py
```

`certify.py` checks the 32,520 coefficient inequalities in dimensions 16 to 255 and all 65,139 cells in the 23 blocks, with the sharper shift bound integrated over each cell. It also evaluates the analytic tail at $D=8{,}388{,}608$, giving $U_d\le4.654016341712388$ beyond it. Inputs are `source_inputs/finite_duals.json`, `proposals_improved_shift.json` and the compressed `cells_*.json.gz`; outputs are `finite_certificate.json`, `blocks_replay.json`, `tail_certificate.json` and `certification_metrics.json`. Without `--replay` it rebuilds the coverings and writes `blocks_certificate.json`. `audit.py` is a separately written check at 90-digit interval precision with different entropy formulas and a different primitive for the shift integral; it writes `audit.json`.

The block seed, from `computations/block_census`:

```
c++ -O3 -std=c++17 enumerate.cpp -o enumerate
./enumerate 20
python3 audit.py
```

`enumerate.cpp` runs through every symmetric truth table up to dimension 20 (up to complementing the output), checks Parseval and the influence identity in integers, encloses logarithms in outward fixed-point arithmetic, and excludes every nonmaximizer in a second pass. `audit.py` proves the surviving ties as exact prime-log expressions and refines them with 144-bit rational logarithms. Together they give $B_{k,n}$ for $2\le k\le5$, $k\le n\le20$, and the seed $C_{20}<3.944265467369$.

The averaged ledger, from `computations/restriction_v4`:

```
python3 rebuild.py
```

`rebuild.py` seeds at dimension 20, adds each of the 8,388,588 later bounds in exact rational arithmetic, checks the sums again by segment, checks that each constant segment peaks at an endpoint, and repeats the binding value with 85-digit decimals. It writes `ledger_certificate.json` and `head_ledger.csv.gz`. `rebuild.py --verify` replays the interval certificates first.

The two non-symmetric examples in the Questions section, from `computations`:

```
python3 check_examples_v4.py
```

It computes integer Walsh spectra on the full 6- and 16-dimensional cubes and certifies $\Delta_i/I_i=5.19487\ldots$ for the majority of three AND pairs and $6.00375\ldots$ for the OR of eight AND pairs, at every coordinate. Output `examples_v4.json`.

The sharp one-step theorem for $n\le15$ and the pointwise majorant of the analytic range are certified as before: `sh restriction/run_all.sh` and `python3 restriction_majorant_certificate.py` from `computations`. `check_all.py` runs further numerical checks of the constants and of every symmetric function with $n\le16$; those are checks, not part of the proof. Its output is in `check_all.out`.

The paper is part of the Math Brain project, whose public repository is https://github.com/yspennstate/math-brain-public.
