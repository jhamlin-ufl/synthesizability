#!/usr/bin/env python3
"""Cross-database hull distances for the paper's targets, structure-matched.

  results/publication_ready/hull_cross_database.pdf

x: OQMD hull distance of the predicted entry (meV/atom, hull of all other
   compositions, see compute_hull_distances.py).
y: the same quantity in the Materials Project and in Alexandria (PBE), but a
   point is drawn only where that database's ground state at the composition
   is the same structure as the OQMD prediction (``mp_same`` / ``alex_same`` in
   database_agreement.csv, from make_database_table.py).  Targets whose
   structure differs, or for which the database has no entry or no hull
   value, are counted in the log so the omissions stay visible.

Reads results/publication_ready/database_agreement.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, 'src')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from synthesizability.formula import sort_by_electronegativity

IN_CSV = Path('results/publication_ready/database_agreement.csv')
OUT = Path('results/publication_ready/hull_cross_database.pdf')

# Outcome class -> colour and marker, the same as plot_disorder_vs_stability.py.
STYLE = {
    'SS': ('#1565C0', 'o', 'Solid Solution'),
    'MP': ('#E65100', 's', 'Multi-Phase'),
    'P':  ('#2E7D32', '^', 'Predicted'),
}
# Database -> filled or hollow marker.
DBS = [
    ('mp', 'Materials Project', True),
    ('alex', 'Alexandria (PBE)', False),
]
MARKER_SIZE = 30


def main():
    df = pd.read_csv(IN_CSV)
    fig, ax = plt.subplots(figsize=(3.4, 3.4), constrained_layout=True)
    notes = []
    xs, ys = [], []
    counts = {}
    for key, label, filled in DBS:
        same = df[f'{key}_same'].map(lambda v: v is True or str(v) == 'True')
        y = df[f'{key}_hull_meV']
        x = df['oqmd_hull_meV']
        shown = same & y.notna() & x.notna()
        xs += list(x[shown]); ys += list(y[shown])
        counts[key] = int(shown.sum())
        for k, (_, r) in enumerate(df[shown].iterrows()):
            colour, marker, _ = STYLE[r['outcome']]
            ax.scatter([r['oqmd_hull_meV']], [r[f'{key}_hull_meV']], marker=marker, s=MARKER_SIZE,
                       facecolors=colour if filled else 'white', edgecolors=colour, linewidths=0.9, zorder=3)
            dy = 3 if k % 2 == 0 else -6
            ax.annotate(sort_by_electronegativity(r['formula']), (r['oqmd_hull_meV'], r[f'{key}_hull_meV']),
                        fontsize=4, xytext=(2.5, dy), textcoords='offset points', color=colour)
        n_diff = int(((df[f'{key}_same'].notna()) & (~same)).sum())
        n_nohull = int((same & y.isna()).sum())
        n_noentry = int(df[f'{key}_same'].isna().sum())
        notes.append(f'{label}: {int(shown.sum())} shown; omitted {n_diff} with a different structure, '
                     f'{n_noentry} with no entry to compare, {n_nohull} matched but without a hull value')

    lo = 10 * np.floor((min(xs + ys) - 4) / 10)
    hi = 10 * np.ceil((max(xs + ys) + 4) / 10)
    ax.plot([lo, hi], [lo, hi], color='0.6', lw=0.8, ls='--', zorder=0)
    ax.axhline(0, color='0.85', lw=0.6, zorder=0)
    ax.axvline(0, color='0.85', lw=0.6, zorder=0)
    ax.set_xlabel(r'OQMD $\Delta E_{\mathrm{hull}}$ (meV/atom)')
    ax.set_ylabel(r'$\Delta E_{\mathrm{hull}}$ in other database (meV/atom)')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal')
    from matplotlib.lines import Line2D
    outcome_handles = [Line2D([0], [0], linestyle='none', marker=m, markerfacecolor=c, markeredgecolor=c,
                              markersize=5, label=lab) for c, m, lab in STYLE.values()]
    db_handles = [Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='black' if filled else 'white',
                         markeredgecolor='black', markersize=5, label=f'{label} (n={counts[key]})')
                  for key, label, filled in DBS]
    leg1 = ax.legend(handles=outcome_handles, fontsize=6, loc='upper left', frameon=False, title='XRD outcome',
                     title_fontsize=6)
    ax.add_artist(leg1)
    ax.legend(handles=db_handles, fontsize=6, loc='lower right', frameon=False, title='Database', title_fontsize=6)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print('\n'.join(notes))
    print(f'Saved {OUT}')


if __name__ == '__main__':
    main()
