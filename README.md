# Cyclic (r,δ)-LRCs from structured defining sets: search code and data

Code and data for the computer search in

> [Authors], *Cyclic (r,δ) Locally Recoverable Codes from Structured Defining Sets*, [journal, year].

## Contents

| Path | Description |
|---|---|
| `magma/search.m` | Exhaustive parameter search (Magma) |
| `magma/verify_table.m` | Rebuilds and checks every code of the table of best codes per length |
| `python/analyse.py` | Produces the counts, the table of best codes per length and the figure of the computational section from `data/` |
| `python/compare_literature.py` | Compares every optimal parameter set of global type with the known constructions listed in the literature table of the paper |
| `data/search_q<q>.csv` | All parameter sets found over GF(q), one line each |
| `results/summary.txt` | Counts quoted in the paper |
| `results/table_best_per_length.{csv,tex}` | Optimal code of largest dimension for every length N ∤ q ± 1 |
| `results/figure_rdall.tex` | Standalone source of the figure (pgfplots; compile with pdflatex) |
| `results/literature_comparison.csv` | Status of each optimal parameter set of global type: given by a cyclic construction, only by a non-cyclic one, or new |

## Construction

For N = sn with N | q^m − 1 and β a primitive N-th root of unity in GF(q^m), the code C is the cyclic code of length N over GF(q) with defining set Z = Z_loc ∪ Z_glob, where Z_loc is the union of the q-cyclotomic cosets of ln + h (0 ≤ l ≤ s−1, 0 ≤ h ≤ δ−2) and Z_glob the union of the q-cyclotomic cosets of a + i (0 ≤ i ≤ b−1). In Magma, β = γ^((q^m−1)/N) with γ = `PrimitiveElement(GF(q^m))`. The repair sets are S_j = {j, j+s, …, j+(n−1)s}.

## Search range

q ∈ {3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19}; every length N ≤ 128 at its minimal extension degree m = ord_N(q) ≤ 3; every factorization N = sn with s, n ≥ 2; every 2 ≤ δ ≤ n, 0 ≤ a ≤ N−1, b ≥ 1 satisfying the conditions of the paper, with K > r and K > r_H. There is no restriction on the dimension K.

## Data files

One file per field, `data/search_q<q>.csv`, sorted by (m, N, s, δ, b, a). Columns:

| Column | Meaning |
|---|---|
| `q, m, N, s, n, delta, r, a, b` | parameters of the construction, r = n − δ + 1 |
| `t` | \|Z_loc ∩ Z_glob\| |
| `K` | dimension |
| `d_BCH` | BCH bound of Z |
| `d` | minimum distance (empty if not determined) |
| `d_lower, d_upper` | proven bounds on d (equal to d when d is known) |
| `d_source` | `MinimumWeight`, `certificate` (d_BCH = B_H, hence d = B_H), `corollary` (d = 4 by the corollary of the paper), `time limit` (MinimumWeight did not finish), `not computed` |
| `r_H, delta_H` | locality pair induced by the auxiliary code B_H |
| `delta_H_source` | `exact`, or `lower bound` when d_min(B_H) was not found within the time limit |
| `B, B_H` | Singleton-type bound at (r, δ) and at (r_H, δ_H) |
| `G_H` | B_H − d |
| `type` | `D` if Z_glob ⊆ Z_loc (direct sum of local codes), `G` otherwise |
| `Z_loc, Z_glob` | defining sets, space-separated |

## How the data were produced

All runs used `magma/search.m` with Magma V2.16-11. The minimum distance was computed with `MinimumWeight` under a time limit: in a first pass for every code (600 seconds per code), and in a second pass, which covered the lengths not reached in the first, only for N ≤ 60 (60 seconds per code). In the second pass the computation of d_min(B_H) was also limited to 60 seconds; for 3528 codes (q = 13, N = 122 and q = 19, N = 120) δ_H is therefore a lower bound, and the induced pair is still a valid locality. Long lengths were split into pieces with the options `StartN` and `SplitN`; the pieces were concatenated and duplicate parameter sets removed. For codes whose d was not determined, `d_lower` is d_BCH.

## Reproducing the results

Search, for one field (about minutes for small q, days for q = 19):

```
magma -b q:=11 magma/search.m
```

Magma stays open when the script has finished; type `quit;` to leave.

Parameters: `MaxN` (default 128), `MaxLenD` (default 60), `MaxTimeD` in seconds (default 60), `OutDir` (default `.`).

Counts, table and figure (Python 3, standard library only):

```
python3 python/analyse.py
python3 python/compare_literature.py
```

Check of the table (prints and writes `verify_table_output.txt`; the last line must read `49 codes checked, 0 mismatches.`):

```
magma magma/verify_table.m
```

The literature comparison gives 140 parameter sets of global type attained by a cyclic construction, 71 attained only by non-cyclic constructions and 9 attained by none of the constructions compared. The other constructions listed in the paper (locality r = 2, repeated-root lengths with gcd(N,q) > 1, lengths dividing q +- 1, and curve-based codes, which need q to be a square or give shorter lengths for prime q) cannot give any of these parameter sets; the script checks these conditions as well.

## License

MIT, see `LICENSE`.
