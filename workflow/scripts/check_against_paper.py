#!/usr/bin/env python3
"""Compare the regenerated outputs with the manuscript's copies (read only).

  in:  results/tables/*.tex, results/figures/disorder_vs_stability_points.csv,
       the manuscript folder (PAPER_DIR, default: the Overleaf folder) tables/*.tex
  out: results/check_report.txt (also printed)

Tables: unified diff of each regenerated table against tables/<same name>.tex.
Fig. 2: the plotted points are checked against Table 1 as printed in the
manuscript (OQMD column and P_disorder, at the printed precision), since the PDF
itself cannot be compared byte for byte.  Differences are reported, never fixed.
The manuscript's revision markup is resolved first: \\new{x} counts as x, and
\\old{x} is dropped (a row whose number is \\old{...} is a deleted row).
"""
import difflib
import os
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, 'workflow/src')
from synthesizability.formula import ordered

PAPER_DIR = Path(os.environ.get('PAPER_DIR', Path.home() / 'ufdb/Apps/Overleaf/Synthesizability Project'))
TABLES = ['sample_table.tex', 'database_table.tex', 'sample_phases_table.tex']
OUT = Path('results/check_report.txt')


def _group_end(s, i):
    """Index just past the brace group that opens at s[i] == '{'."""
    depth = 0
    for j in range(i, len(s)):
        depth += {'{': 1, '}': -1}.get(s[j], 0)
        if depth == 0:
            return j + 1
    raise ValueError(f'unbalanced braces: {s[i:i + 60]}')


def strip_markup(s):
    r"""Resolve revision markup: \new{x} -> x, \old{x} -> nothing."""
    out, i = [], 0
    while i < len(s):
        if s.startswith(r'\new{', i) or s.startswith(r'\old{', i):
            end = _group_end(s, i + 4)
            if s.startswith(r'\new{', i):
                out.append(strip_markup(s[i + 5:end - 1]))
            i = end
        else:
            out.append(s[i])
            i += 1
    return ''.join(out)


def manuscript_table(tex):
    r"""Manuscript table with deleted rows (number cell \old{...}) removed and markup resolved."""
    return strip_markup('\n'.join(l for l in tex.splitlines() if not l.lstrip().startswith(r'\old{')))


def table_rows(tex):
    """{row number: cells} of the numbered rows of a LaTeX table; empty cells are kept."""
    rows = {}
    for line in tex.splitlines():
        cells = [c.strip() for c in re.split(r'(?<!\\)&', line.strip().removesuffix('\\\\'))]
        if cells[0].isdigit():
            rows[int(cells[0])] = cells
    return rows


def main():
    report, n_diff = [f'Manuscript folder: {PAPER_DIR}', ''], 0
    for name in TABLES:
        ours, theirs = Path('results/tables') / name, PAPER_DIR / 'tables' / name
        if not theirs.exists():
            report += [f'== {name}: not in the manuscript folder', '']
            continue
        diff = list(difflib.unified_diff(manuscript_table(theirs.read_text()).splitlines(), ours.read_text().splitlines(),
                                         f'manuscript/{name}', f'regenerated/{name}', n=0, lineterm=''))
        n_diff += bool(diff)
        report += [f'== {name}: ' + ('identical' if not diff else f'{sum(l[0] in "+-" for l in diff[2:])} lines differ')]
        report += diff + ['']

    pts = pd.read_csv('results/figures/disorder_vs_stability_points.csv')
    printed = table_rows(manuscript_table((PAPER_DIR / 'tables/sample_table.tex').read_text()))
    formula_of = {n: re.sub(r'\\ch\{(\w+)\}.*', r'\1', c[1]) for n, c in printed.items()}
    by_formula = {f: printed[n] for n, f in formula_of.items()}
    bad = []
    for _, p in pts.iterrows():
        cells = by_formula.get(ordered(p.formula))
        if cells is None:
            bad.append(f'{p.formula}: not in the manuscript Table 1')
            continue
        if f'{p.hull_meV:.1f}' != cells[3] or f'{p.disorder_probability:.2f}' != cells[6]:
            bad.append(f'{p.formula}: figure ({p.disorder_probability:.2f}, {p.hull_meV:.1f}) '
                       f'vs Table 1 ({cells[6]}, {cells[3]})')
    n_diff += bool(bad)
    report += ['== Fig. 2 points vs manuscript Table 1 (P_disorder, OQMD): '
               + ('all 34 agree' if not bad else f'{len(bad)} differ')] + bad
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('\n'.join(report) + '\n')
    print('\n'.join(report))
    print(f'\n{n_diff} of {len(TABLES) + 1} checks show differences -> {OUT}')


if __name__ == '__main__':
    main()
