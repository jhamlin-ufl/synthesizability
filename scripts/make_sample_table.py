#!/usr/bin/env python3
"""Sample table for the paper: formula, hull distance, disorder parameter, XRD outcome.

  results/publication_ready/sample_table.tex

One row per classified target, in the order the figures use (solid solutions,
multiphase, predicted structure found; sample number within each group).  The
sample set, outcome calls, caveats and footnotes all come from
synthesizability.paper_samples, so this table, the XRD figures and the
disorder-vs-stability plot describe the same samples the same way.

The output is a complete RevTeX ``table*`` environment meant to be ``\\input``
into the manuscript.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, 'src')

from synthesizability.formula import sort_by_electronegativity
from synthesizability.paper_samples import (
    OUTCOME_LABEL, OUTCOMES, TABLE_COMMENT, load_classified_samples,
)

OUT_TEX = Path('results/publication_ready/sample_table.tex')
HULL_CSV = Path('data/processed/hull_distances.csv')
WRAP_AT = 56          # characters; longer outcome texts are broken into two lines
LABEL = 'tab:xrd_outcomes'

CAPTION = (
    'Arc-melted targets with an XRD outcome, grouped by outcome class; solid solutions are '
    'ordered by decreasing disorder parameter, the other groups by decreasing energy relative to the hull. '
    'List is the candidate list the target came from (B, C and J as defined in the Methods). '
    r'$\Delta E_{\mathrm{hull}}$ is the energy of the lowest-energy OQMD entry at the target composition relative to the '
    'OQMD convex hull of all other compositions, in meV/atom, so that negative values mean the '
    'database predicts the compound to be stable; the Materials Project and Alexandria views of the same '
    'targets are tabulated in the Supplement. '
    r'$P_{\mathrm{disorder}}$ is the predicted probability that the composition forms a '
    'chemically disordered structure. '
    'The outcome column names the phase or mixture identified by whole-pattern fitting; '
    'numbers in parentheses after a site are refined site fractions, and percentages are '
    'weight fractions.'
)


def list_code(prediction_list):
    """'B (two active binaries, ...)' -> 'B'; the Gao et al. list -> 'J'."""
    head = str(prediction_list).split()[0]
    return {'B': 'B', 'C': 'C', "James's": 'J'}.get(head, head)


def wrap(text, width=WRAP_AT):
    """Break a long outcome text into two lines at the space nearest the middle,
    using makecell (RevTeX's tabular cannot wrap p{} columns on current TeX Live)."""
    if len(text) <= width:
        return text
    cut = max((i for i, ch in enumerate(text) if ch == ' ' and i <= len(text) // 2 + 8), default=None)
    if cut is None:
        return text
    return r'\makecell[l]{' + text[:cut] + r' \\ ' + text[cut + 1:] + '}'


def fmt_hull(value):
    return '' if value is None or value != value else f'{value:.1f}'


def fmt_formula(formula):
    return r'\ch{' + sort_by_electronegativity(formula) + '}'


def main():
    df = load_classified_samples()
    unknown = set(TABLE_COMMENT) - set(df['sample_number'])
    if unknown:
        raise SystemExit(f'TABLE_COMMENT names samples not in the table: {sorted(unknown)}')
    missing = set(df['sample_number']) - set(TABLE_COMMENT)
    if missing:
        raise SystemExit(f'TABLE_COMMENT lacks entries for: {sorted(missing)}')

    # Plain tabular with booktabs rules and no paragraph column: RevTeX's tabular
    # patches choke on p{} columns on current TeX Live (Overleaf: "Extra \\or").
    # The outcome texts are kept short enough for a single line instead.
    hull = pd.read_csv(HULL_CSV).set_index('formula')
    lines = [
        r'\begin{table*}[!t]',
        r'\centering',
        r'\footnotesize',
        r'\caption{' + CAPTION + '}',
        r'\label{' + LABEL + '}',
        r'\vspace{4pt}',
        r'\setlength{\tabcolsep}{6pt}',
        r'\renewcommand{\arraystretch}{1.0}',
        r'\begin{tabular}{rllrrl}',
        r'\toprule',


        r'\# & Formula & List & \makecell{$\Delta E_{\mathrm{hull}}$ \\ (meV/atom)} & $P_{\mathrm{disorder}}$ & XRD outcome \\',
        r'\midrule',
    ]
    # Row order within each group: solid solutions by decreasing disorder
    # parameter, the other two groups by decreasing OQMD energy relative to the hull.
    order = {'SS': ('disorder_probability', False), 'MP': ('oqmd_stability', False), 'P': ('oqmd_stability', False)}
    first = True
    row = 0
    for outcome in OUTCOMES:
        sub = df[df['XRD Result'] == outcome]
        if sub.empty:
            continue
        col, asc = order[outcome]
        sub = sub.sort_values([col, 'sample_number'], ascending=[asc, True])
        if not first:
            lines.append(r'\midrule')
        first = False
        lines.append(r'\multicolumn{6}{l}{\textit{' + OUTCOME_LABEL[outcome] + '}} \\\\')
        for _, r in sub.iterrows():
            num = int(r['sample_number'])
            h = hull.loc[r['formula']]
            row += 1
            lines.append(
                f'{row} & {fmt_formula(r["formula"])} & {list_code(r["prediction_list"])} & '
                f'{fmt_hull(h["oqmd_hull_meV"])} & '
                f'{r["disorder_probability"]:.2f} & {wrap(TABLE_COMMENT[num])} \\\\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']

    OUT_TEX.parent.mkdir(parents=True, exist_ok=True)
    OUT_TEX.write_text('\n'.join(lines) + '\n')
    counts = '  '.join(f'{o}: {(df["XRD Result"] == o).sum()}' for o in OUTCOMES)
    print(f'{len(df)} rows ({counts})\nSaved: {OUT_TEX}')


if __name__ == '__main__':
    main()
