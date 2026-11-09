# The Fourier entropy–influence inequality for symmetric Boolean functions

Yitzchak Shmalo

[Paper (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) · [Declarations](DECLARATIONS.md)

The paper proves the Fourier entropy–influence inequality for every symmetric Boolean function with constant $4650431/1000000=4.650431$, improving the constant 12.04 of O'Donnell, Wright and Zhou (ICALP 2011). All logarithms in the statements are to base 2, $H(f)$ is the Fourier entropy and $I(f)$ the total influence.

## Main results

- **Theorem A.** Every symmetric Boolean function satisfies

  $$H(f)\le c\,I(f),\qquad c=\frac{4650431}{1000000}=4.650431.$$

- **Theorem B.** (i) Restricting one coordinate of a symmetric function costs at most $K=12-\tfrac92\log_2 3=4.8677$ times the influence of that coordinate: $\Delta_1(f)\le K\,I_1(f)$ in every dimension, with equality for $\mathrm{OR}_3$. (ii) In high dimension the step is cheaper: $\limsup_n\sup_f\Delta_1/I_1<4.637$.

- **Averaging.** If $U_d$ bounds $\Delta_1/I_1$ for symmetric functions in dimension $d$, then $H(f)\le\frac1n(U_1+\dots+U_n)\,I(f)$. A single uniform bound cannot beat $K$ because of $\mathrm{OR}_3$, but the average can. Theorem A comes from this average, seeded by a block-restriction census through dimension 24 and the dimension-specific bounds above it.

- **Theorem C.** If $f$ is symmetric on $n\ge2$ variables and $g$ is its derivative in one coordinate, then for $0\le\lambda\le1$

  $$\mathrm{Stab}_\lambda[g]\le e^{-v}I_0(v)\,\|g\|_2^2,\qquad v=\frac{(n-1)(1-\lambda^2)}{4}.$$

- **Theorem D.** Without any computation, $H(f)\le\bigl(\log_2(4+e)+2+1/\ln 2\bigr)I(f)<6.191\,I(f)$ for every symmetric $f$.

- **Theorem E.** For OR and AND, $H/I=4-(4-2/\ln 2)/n+O(2^{-n})$. The optimal constant for symmetric functions lies in $[4,\,4.650431]$.

The finite parts are computer-assisted. Every symmetric function is checked in exact arithmetic for $n\le15$. Dimensions $16$ to $255$ reduce to one-variable inequalities with rational multipliers, and dimensions $256$ to $33{,}554{,}432$ to $75{,}638$ rational cells in 24 blocks, all verified in outward-rounded interval arithmetic. An analytic bound decreasing in $n$ covers every larger dimension. The averaged ledger is summed in exact rational arithmetic; its largest value, $4.650430788786238\ldots$, is at $n=517$.

Version 1 (8 November 2026) had the constant 6.1908, version 2 had 6.1347, version 3 had 4.8677 and version 4 had 4.657435.

## Computations

Python needs `mpmath`; the census needs a C++17 compiler. Each command uses one thread. Run the commands from the repository root.

Rebuild and audit the dimension-24 census, the additional continuous block and analytic tail, and the exact rational ledger:

```
python3 computations/reproduce.py
```

The runner compiles the census programs in a temporary directory, removes the executables afterwards, and runs the following certificates in order. Individual steps are:

```
(cd computations/block_census_24 && c++ -O3 -std=c++17 enumerate.cpp -o enumerate && ./enumerate 24 && rm enumerate)
(cd computations/block_census_24 && python3 audit.py)
(cd computations/block_census_24 && python3 prepare_check.py && c++ -O3 -std=c++17 check_census.cpp -o check_census && ./check_census && rm check_census)
python3 computations/certify_extension.py
python3 computations/check_extension.py
python3 computations/rebuild_ledger.py
python3 computations/audit_ledger.py
```

`block_census_24/enumerate.cpp` checks every symmetric truth table up to dimension 24, up to complementing the output. It checks Parseval and the influence identity in integers, encloses logarithms in outward fixed-point arithmetic, and excludes every nonmaximizer in a second pass. Entropy comparisons use exact signed 192-bit products. `audit.py` proves the surviving ties as exact prime-log expressions and refines them with 144-bit rational logarithms. `check_census.cpp` independently repeats every new table comparison in dimensions 21 through 24, with binary traversal, directly multiplied Krawtchouk polynomials, 48-bit logarithm enclosures and overflow-checked 128-bit products. The certificates cover all 86 pairs $(k,n)$ with $2\le k\le5$ and $k\le n\le24$, and give $C_{24}<3.953557670092$. Outputs include `exhaustive.jsonl`, `certified_maxima.json`, `audit.json` and `exclusion_recheck.jsonl`.

`certify_extension.py` uses the final dual in `restriction_v4/proposals_improved_shift.json` with one additional positive multiplier to certify all 10,499 cells of the block $8{,}388{,}609\le n\le33{,}554{,}432$, giving $U_n\le4.588910281347744$. It evaluates Proposition 4.5 at $n_0=33{,}554{,}432$ with $a=1/45$, giving $U_n\le4.647117187337907$ for every $n\ge n_0$. `check_extension.py` rechecks the cells with different entropy formulas and a different primitive for the shift integral, at 90-digit interval precision, then checks the tail by termwise outward intervals and 85-digit integral quadrature. Inputs and outputs are in `computations/restriction_extension`.

`rebuild_ledger.py` seeds at dimension 24 and adds each of the 33,554,408 later bounds in exact rational arithmetic. It checks the sums again by segment and verifies that each constant segment peaks at an endpoint. The largest bound is $4.650430788786238\ldots$ at $n=517$, below $4650431/1000000$. `audit_ledger.py` independently reconstructs the seed as shortest paths in $nC_n$, reads the one-step bounds from the individual certificates, checks the segment endpoints and all 16,360 saved head rows, and repeats the binding value and cutoff interface at 100 decimal digits. Outputs are `restriction_extension/ledger_certificate.json`, `head_ledger.csv.gz` and `ledger_audit.json`.

The earlier one-step certificates remain in `computations/restriction_v4`:

```
(cd computations/restriction_v4 && python3 certify.py --replay && python3 audit.py)
```

These check the 32,520 coefficient inequalities in dimensions 16 through 255 and the 65,139 cells in the first 23 blocks. Their earlier tail and dimension-20 ledger remain available for comparison. The earlier census remains in `computations/block_census`:

```
(cd computations/block_census && c++ -O3 -std=c++17 enumerate.cpp -o enumerate && ./enumerate 20 && rm enumerate && python3 audit.py)
```

The two non-symmetric examples in the Questions section, from the repository root:

```
(cd computations && python3 check_examples_v4.py)
```

It computes integer Walsh spectra on the full 6- and 16-dimensional cubes and certifies $\Delta_i/I_i=5.19487\ldots$ for the majority of three AND pairs and $6.00375\ldots$ for the OR of eight AND pairs, at every coordinate. Output `examples_v4.json`.

The sharp one-step theorem for $n\le15$ and the pointwise majorant of the analytic range are certified as before: `(cd computations && sh restriction/run_all.sh)` and `python3 computations/restriction_majorant_certificate.py`. `check_all.py` runs further numerical checks of the constants and of every symmetric function with $n\le16$; those are checks, not part of the proof. Run it with `python3 computations/check_all.py`.

The paper is part of the Math Brain project, whose public repository is https://github.com/yspennstate/math-brain-public.

Verify the supplied files before rerunning certificates or compiling the paper. Rebuilding updates timing fields, compressed outputs and PDF metadata.

```
shasum -a 256 -c SHA256SUMS
tectonic paper/main.tex
```
