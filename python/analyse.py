#!/usr/bin/env python3
"""Summary, table and figure of the paper from data/search_q*.csv.

Usage:  python3 python/analyse.py            (run from the repository root)

Writes
  results/summary.txt                  counts quoted in Section 4
  results/table_best_per_length.csv    optimal code of largest dimension
  results/table_best_per_length.tex      for every length N not dividing q-1, q+1
  results/figure_rdall.tex             figure of the paper (standalone LaTeX source)

Only the Python standard library is used.
"""
import csv
import glob
import os
from collections import Counter, defaultdict
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
FIGURE_FIELDS = (3, 5, 7, 9, 11, 13, 16, 17, 19)
STYLE = {1: ("blue", "square*", "1.1pt"), 2: ("red", "*", "1.0pt"), 3: ("teal", "triangle*", "1.2pt")}
INT_COLS = ("q", "m", "N", "s", "n", "delta", "r", "a", "b", "t", "K", "d_BCH",
            "d_lower", "d_upper", "r_H", "delta_H", "B", "B_H")


def load():
    rows = []
    for path in sorted(glob.glob(os.path.join(DATA, "search_q*.csv"))):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                for c in INT_COLS:
                    r[c] = int(r[c])
                r["d"] = int(r["d"]) if r["d"] else None
                rows.append(r)
    return rows


def is_optimal(r):
    """Singleton-type optimal at (r_H, delta_H), with the restrictions of the paper."""
    q, N, d = r["q"], r["N"], r["d"]
    return (d is not None and d == r["B_H"] and d >= 3 and r["r_H"] >= 2
            and r["K"] > r["r_H"] and (q - 1) % N != 0 and (q + 1) % N != 0)


def optimal_sets(rows):
    sets = defaultdict(list)
    for r in rows:
        if is_optimal(r):
            sets[(r["q"], r["N"], r["K"], r["d"], r["r_H"], r["delta_H"])].append(r)
    return sets


def summary(rows, sets):
    src = Counter(r["d_source"] for r in rows)
    known = sum(1 for r in rows if r["d"] is not None)
    types = Counter(v[0]["type"] for v in sets.values())
    lines = [
        f"parameter sets                         {len(rows)}",
        f"minimum distance known                 {known}",
        f"  by MinimumWeight                     {src['MinimumWeight']}",
        f"  by the certificate d_BCH = B_H       {src['certificate']}",
        f"  by the corollary                     {src['corollary']}",
        f"MinimumWeight did not finish           {src['time limit']}",
        f"  of which N <= 60                     {sum(1 for r in rows if r['d_source'] == 'time limit' and r['N'] <= 60)}",
        f"not computed (N > 60)                  {src['not computed']}",
        f"delta_H replaced by a lower bound      {sum(1 for r in rows if r['delta_H_source'] == 'lower bound')}",
        f"optimal parameter sets (paper)         {len(sets)}",
        f"  type G / type D                      {types['G']} / {types['D']}",
        f"  with N > 60                          {sum(1 for k in sets if k[1] > 60)}",
        f"  with m = 3                           {sum(1 for v in sets.values() if v[0]['m'] == 3)}",
    ]
    return "\n".join(lines) + "\n"


def best_per_length(sets):
    by_length = defaultdict(list)
    for key, codes in sets.items():
        by_length[key[:2]].append(codes)
    table = []
    for (q, N), groups in sorted(by_length.items()):
        codes = [c for g in groups for c in g]
        codes.sort(key=lambda r: (-r["K"], -r["d"], r["r_H"], r["type"] != "G",
                                  r["s"], r["delta"], r["b"], r["a"]))
        table.append(dict(codes[0], n_sets=len(groups)))
    return table


def write_table(table):
    cols = ["q", "m", "N", "K", "d", "d_source", "s", "n", "delta", "a", "b",
            "r_H", "delta_H", "type", "n_sets", "Z_loc", "Z_glob"]
    with open(os.path.join(RESULTS, "table_best_per_length.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(table)
    out = [r"\begin{tabular}{@{}rrrrrcccr@{}}", r"\toprule",
           r"$q$ & $m$ & $N$ & $K$ & $d$ & $(s,n,\delta,a,b)$ & $(r_H,\delta_H)$ & Type & $\#$ \\",
           r"\midrule"]
    prev = None
    for r in table:
        if prev is not None and r["q"] != prev:
            out.append(r"\midrule")
        prev = r["q"]
        mark = r"$^{\dagger}$" if r["d_source"] == "corollary" else ""
        out.append(f"{r['q']} & {r['m']} & {r['N']} & {r['K']} & {r['d']}{mark} & "
                   f"$({r['s']},{r['n']},{r['delta']},{r['a']},{r['b']})$ & "
                   f"$({r['r_H']},{r['delta_H']})$ & {r['type']} & {r['n_sets']} \\\\")
    out += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(RESULTS, "table_best_per_length.tex"), "w") as f:
        f.write("\n".join(out) + "\n")


def upper_hull(points):
    """Upper convex hull, left to right, of points (x, y) with distinct x."""
    hull = []
    for p in sorted(points):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) >= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return hull


def fmt(points):
    return " ".join(f"({float(x):.4f},{float(y):.4f})" for x, y in points)


def write_figure(rows):
    """Standalone LaTeX source of the figure of the paper (compile with pdflatex).

    For each field and each m, the vertices of the upper convex hull are joined by
    a solid line as long as d/N increases (codes not dominated in r_H/N and d/N);
    the remaining hull vertices, worse in both respects, are shaded and joined by a
    dashed line. All codes with known minimum distance are drawn. Over each field, for every value
    of r_H/N, the code with the largest d/N over all three regimes is drawn in full
    colour and the other codes are shaded. The lines are geometric envelopes."""
    style = {1: ("blue", "square*", "1.3pt"), 2: ("red", "*", "1.2pt"),
             3: ("teal!80!black", "triangle*", "1.5pt")}
    shade = "!30"
    out = [r"""\documentclass[10pt]{standalone}
\usepackage{mathptmx}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{groupplots}
\usetikzlibrary{calc}
\begin{document}
\begin{tikzpicture}
\begin{groupplot}[
  group style={group size=3 by 3, horizontal sep=1.1cm, vertical sep=1.0cm,
               xlabels at=edge bottom, ylabels at=edge left, xticklabels at=edge bottom},
  width=4.9cm, height=3.1cm, scale only axis,
  xmin=0, xmax=0.5, ymin=0, ymax=1,
  xtick={0,0.1,...,0.5}, ytick={0,0.2,...,1},
  xticklabels={,0.1,0.2,0.3,0.4,0.5},
  xlabel={$r_{H}/N$}, ylabel={$d/N$},
  grid=major, grid style={densely dotted,gray!40},
  tick label style={font=\scriptsize}, label style={font=\footnotesize},
  title style={font=\footnotesize, yshift=-1.5pt},
]"""]
    OPEN = {"square*": "square", "*": "o", "triangle*": "triangle"}
    for i, q in enumerate(FIGURE_FIELDS):
        opt = r"title={$q=%d$}" % q
        if i == 0:
            opt += r", legend to name=rdlegend, legend columns=6, legend style={font=\footnotesize, draw=none, column sep=8pt}"
        out.append(r"\nextgroupplot[%s]" % opt)
        if i == 0:
            for m in (1, 2, 3):
                c, mk, ms = style[m]
                out.append(r"\addlegendimage{color=%s, thick, mark=%s, mark size=%s}\addlegendentry{$m=%d$}" % (c, mk, ms, m))
            out.append(r"\addlegendimage{only marks, mark=*, mark size=1.2pt, color=gray}\addlegendentry{best at fixed $r_H/N$}")
            out.append(r"\addlegendimage{only marks, mark=*, mark size=1.2pt, color=gray%s}\addlegendentry{other codes}" % shade)
            out.append(r"\addlegendimage{color=gray, thick, dashed}\addlegendentry{convex hull}")
        layers = {"shaded": [], "dashed": [], "solid": [], "good": []}
        for m in (1, 2, 3):
            best = {}
            for r in rows:
                if r["q"] == q and r["m"] == m and r["d"] is not None:
                    x, y = Fraction(r["r_H"], r["N"]), Fraction(r["d"], r["N"])
                    best[x] = max(best.get(x, y), y)
            if not best:
                continue
            c, mk, ms = style[m]
            hull = upper_hull(best.items())
            k = 0
            while k + 1 < len(hull) and hull[k + 1][1] > hull[k][1]:
                k += 1
            good, rest = hull[:k + 1], hull[k:]
            pts = sorted({(Fraction(r2["r_H"], r2["N"]), Fraction(r2["d"], r2["N"]))
                          for r2 in rows if r2["q"] == q and r2["m"] == m and r2["d"] is not None})
            # over the same field, all three regimes together: for each value of r_H/N
            # the code with the largest d/N is highlighted, the other codes are shaded
            top = {}
            for r2 in rows:
                if r2["q"] == q and r2["d"] is not None:
                    x2, y2 = Fraction(r2["r_H"], r2["N"]), Fraction(r2["d"], r2["N"])
                    top[x2] = max(top.get(x2, y2), y2)
            front = sorted(p for p in pts if p[1] == top[p[0]])
            shaded = sorted(p for p in pts if p[1] < top[p[0]])
            if shaded:
                layers["shaded"].append(r"\addplot[only marks, mark=%s, mark size=%s, color=%s%s, forget plot] coordinates {%s};" % (mk, ms, c, shade, fmt(shaded)))
            if len(rest) >= 2:
                layers["dashed"].append(r"\addplot[color=%s, thick, dashed, no marks, forget plot] coordinates {%s};" % (c, fmt(rest)))
            if len(good) >= 2:
                layers["solid"].append(r"\addplot[color=%s, thick, no marks, forget plot] coordinates {%s};" % (c, fmt(good)))
            layers["good"].append(r"\addplot[only marks, mark=%s, mark size=%s, color=%s, forget plot] coordinates {%s};" % (mk, ms, c, fmt(front)))
        for key in ("shaded", "dashed", "solid", "good"):
            out += layers[key]
    out.append(r"\end{groupplot}")
    out.append(r"\node[anchor=south] at ($(group c2r1.north)+(0,0.45cm)$) {\pgfplotslegendfromname{rdlegend}};")
    out.append(r"\end{tikzpicture}")
    out.append(r"\end{document}")
    with open(os.path.join(RESULTS, "figure_rdall.tex"), "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    rows = load()
    sets = optimal_sets(rows)
    text = summary(rows, sets)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as f:
        f.write(text)
    table = best_per_length(sets)
    write_table(table)
    write_figure(rows)
    print(text, end="")
    print(f"table rows                             {len(table)}")


if __name__ == "__main__":
    main()
