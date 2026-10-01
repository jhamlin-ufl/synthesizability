#!/usr/bin/env python3
"""
Publication figure: OQMD hull distance vs. disorder parameter for the arc-melted
targets, coloured by XRD outcome.

  results/publication_ready/disorder_vs_stability.pdf

Points are the 33 non-Diffusion-Model targets that have an XRD classification
(SS / MP / P).  An "X" is drawn over any point whose most thermodynamically
stable OQMD polymorph (polymorph_rank == 1) is dynamically unstable, i.e. has
imaginary phonon modes.

Each point is labelled with its target formula.  Samples that formed the
predicted structure only with a qualification -- a rediscovery, an off-hull
polymorph, a disordered variant -- carry that caveat under the label.

The sample set, the XRD outcome calls and the caveats all come from
synthesizability.paper_samples, so this figure and the XRD comparison grid
always describe the same 33 samples.
"""
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, 'src')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from synthesizability.formula import formula_to_mathtext
from synthesizability.paper_samples import CAVEAT, PARTIAL, load_classified_samples

try:
    from adjustText import adjust_text
except ImportError:  # pragma: no cover
    raise SystemExit('adjustText is required for the formula labels: '
                     'pip install adjustText')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PHONON_CSV  = Path('data/external/phonon_stability/MANIFEST_w_stability.csv')
OUT_PDF     = Path('results/publication_ready/disorder_vs_stability.pdf')

XLIM = (-0.11, 1.09)   # margin so edge labels are not clipped
YLIM = (-0.08, 0.02)
AXES_RECT = [0.1530, 0.1219, 0.8080, 0.8433]   # matches the original figure
FIGSIZE = (7.2, 5.0)

# XRD outcome → (colour, marker, legend label)
STYLE = {
    'SS': ('#1565C0', 'o', 'Solid Solution'),
    'MP': ('#E65100', 's', 'Multi-Phase'),
    'P':  ('#2E7D32', '^', 'Predicted'),
}
MARKER_SIZE = 55        # pt^2
EDGE_WIDTH  = 0.5
ALPHA       = 0.85

LINESPACING = 1.15      # tightens the caveat under its formula

# Hand corrections to the automatic label placement, keyed by sample number and
# given as (dx, dy) in axis units -- disorder parameter and eV/atom, the same
# numbers you read off the axes.  adjust_text optimises collisions but not
# readability; these are the cases where a person overrules it.  Any leader line
# follows its label.  Set VERBOSE_LABELS=1 to print where every label landed.
LABEL_NUDGE = {
    554: (-0.015, +0.004),   # Gd2Fe3Ge5: off the YCoSi3 marker to its right
}

# Absolute label positions in axis units, for the labels adjust_text cannot
# place well.  Both of these sit in the crowded low-disorder corner, where the
# optimiser flings them across the plot and the leader line then crosses
# unrelated markers.
LABEL_PLACE = {
    507: (+0.1502, -0.0757),   # YSnCu: directly under its own marker
    539: (+0.0600, -0.0350),   # YCoSi3: below, the nearest gap that fits it
}

X_SIZE  = 34            # pt^2 for the dynamic-instability overlay
X_WIDTH = 1.1

LABEL_SIZE  = 5.5       # pt — formula labels
LABEL_COLOR = '#333333'
LEADER_WIDTH = 0.35     # pt — leader lines from label to point
MIN_LEADER_PX = 34      # below this label-to-point distance, drop the leader


def caveat_to_mathtext(caveat):
    """'disordered variant' -> italic mathtext, parenthesised."""
    return '$\\mathit{(' + caveat.replace(' ', '\\ ') + ')}$'


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dynamic_stability():
    """Phonon verdict for the most thermodynamically stable polymorph."""
    phonon = pd.read_csv(PHONON_CSV)
    ground = phonon[phonon['polymorph_rank'] == 1]

    dup = ground['target_formula'].duplicated()
    if dup.any():
        raise ValueError(f"multiple rank-1 rows for "
                         f"{sorted(ground.loc[dup, 'target_formula'])}")

    verdict = ground['phonon_status'].str.strip().str.lower()
    unknown = set(verdict) - {'stable', 'unstable'}
    if unknown:
        raise ValueError(f"unrecognised phonon_status values: {sorted(unknown)}")

    return pd.Series(verdict.values == 'unstable',
                     index=ground['target_formula'].values,
                     name='dynamically_unstable')


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_figure(df):
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.7, alpha=0.4,
               zorder=1)

    partial = df['sample_number'].isin(PARTIAL)
    for outcome, (colour, marker, _) in STYLE.items():
        g = df[(df['XRD Result'] == outcome) & ~partial]
        ax.scatter(g['disorder_probability'], g['oqmd_stability'],
                   s=MARKER_SIZE, c=colour, marker=marker, alpha=ALPHA,
                   edgecolors='black', linewidths=EDGE_WIDTH, zorder=2)

    # Partial successes: the marker of their outcome class, half filled.
    for _, r in df[partial].iterrows():
        colour, marker, _ = STYLE[r['XRD Result']]
        ax.plot(r['disorder_probability'], r['oqmd_stability'], linestyle='none',
                marker=marker, markersize=np.sqrt(MARKER_SIZE),
                markerfacecolor=colour, markerfacecoloralt='white',
                fillstyle='bottom', markeredgecolor='black',
                markeredgewidth=EDGE_WIDTH, alpha=ALPHA, zorder=2)

    unstable = df[df['dynamically_unstable']]
    ax.scatter(unstable['disorder_probability'], unstable['oqmd_stability'],
               s=X_SIZE, c='black', marker='x', linewidths=X_WIDTH, zorder=3)

    # The caveat is a second line of the formula label rather than a text of
    # its own, so adjust_text reserves room for it and any leader line lands
    # on the whole block.
    def label(row):
        text = formula_to_mathtext(row['formula'])
        if isinstance(row['caveat'], str):
            text += '\n' + caveat_to_mathtext(row['caveat'])
        return text

    texts = [
        ax.text(r['disorder_probability'], r['oqmd_stability'], label(r),
                fontsize=LABEL_SIZE, color=LABEL_COLOR, zorder=4,
                ha='center', va='center', linespacing=LINESPACING,
                # invisible on the white ground; masks the dashed zero line
                # where a label happens to sit on it
                bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none',
                          alpha=0.75))
        for _, r in df.iterrows()
    ]

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(np.arange(0.0, 1.01, 0.2))
    ax.set_xlabel('Disorder Parameter', fontsize=11)
    ax.set_ylabel('Stability (eV/atom)', fontsize=11)
    ax.tick_params(labelsize=10)

    counts = df.loc[~partial, 'XRD Result'].value_counts()
    handles = [
        Line2D([0], [0], linestyle='none', marker=marker, color=colour,
               markeredgecolor='black', markeredgewidth=EDGE_WIDTH,
               markersize=np.sqrt(MARKER_SIZE), alpha=ALPHA,
               label=f'{label} (n={counts.get(outcome, 0)})')
        for outcome, (colour, marker, label) in STYLE.items()
    ]
    if partial.any():
        colour, marker, _ = STYLE['P']
        handles.append(
            Line2D([0], [0], linestyle='none', marker=marker,
                   markerfacecolor=colour, markerfacecoloralt='white',
                   fillstyle='bottom', markeredgecolor='black',
                   markeredgewidth=EDGE_WIDTH, markersize=np.sqrt(MARKER_SIZE),
                   alpha=ALPHA, label=f'Partial success (n={int(partial.sum())})'))
    handles.append(
        Line2D([0], [0], linestyle='none', marker='x', color='black',
               markeredgewidth=X_WIDTH, markersize=np.sqrt(X_SIZE),
               label=f'Dynamically unstable (n={int(df["dynamically_unstable"].sum())})')
    )
    legend = ax.legend(handles=handles, loc='lower right', fontsize=9,
                       framealpha=0.8)

    # Nudge labels off their markers and off each other.  The legend is passed
    # as an obstacle so labels do not end up underneath it.
    fig.canvas.draw()
    texts, arrows = adjust_text(texts, ax=ax,
                x=df['disorder_probability'].to_numpy(),
                y=df['oqmd_stability'].to_numpy(),
                objects=[legend.get_frame()],
                expand=(1.6, 1.75),
                force_text=(0.55, 0.7),
                force_static=(1.1, 1.2),
                force_pull=(0.035, 0.035),
                max_move=28,
                iter_lim=3000,
                arrowprops=dict(arrowstyle='-', color=LABEL_COLOR,
                                lw=LEADER_WIDTH, shrinkA=1, shrinkB=3))

    apply_label_nudges(fig, ax, df, texts, arrows)

    return fig


def apply_label_nudges(fig, ax, df, texts, arrows):
    """Shift hand-listed labels vertically, dragging their leader lines along."""
    for name, table in (('LABEL_NUDGE', LABEL_NUDGE), ('LABEL_PLACE', LABEL_PLACE)):
        unknown = set(table) - set(df['sample_number'])
        if unknown:
            raise ValueError(f'{name} names samples not in the figure set: '
                             f'{sorted(unknown)}')
    both = set(LABEL_NUDGE) & set(LABEL_PLACE)
    if both:
        raise ValueError(f'samples in both LABEL_NUDGE and LABEL_PLACE: '
                         f'{sorted(both)}')

    arrow_for = {id(a.patchA): a for a in arrows
                 if getattr(a, 'patchA', None) is not None}

    for text, (_, row) in zip(texts, df.iterrows()):
        num = row['sample_number']
        if num in LABEL_PLACE:
            text.set_position(LABEL_PLACE[num])
        elif num in LABEL_NUDGE:
            dx, dy = LABEL_NUDGE[num]
            x, y = text.get_position()
            text.set_position((x + dx, y + dy))
        else:
            continue
        arrow = arrow_for.get(id(text))
        if arrow is None:
            continue
        x, y = text.get_position()
        point = (row['disorder_probability'], row['oqmd_stability'])
        px = ax.transData.transform((x, y)) - ax.transData.transform(point)
        if (px[0] ** 2 + px[1] ** 2) ** 0.5 < MIN_LEADER_PX:
            arrow.remove()
        else:
            # patchA is the text itself, so the tail re-clips to the new box.
            arrow.set_positions((x, y), point)

    if os.environ.get('VERBOSE_LABELS'):
        print('Label positions (x, y) after nudging:')
        for text, (_, row) in zip(texts, df.iterrows()):
            x, y = text.get_position()
            print(f"  {int(row['sample_number']):04d} {row['formula']:<12} "
                  f"label=({x:+.4f}, {y:+.5f})  "
                  f"point=({row['disorder_probability']:+.4f}, "
                  f"{row['oqmd_stability']:+.5f})")


def main():
    df = (load_classified_samples()
          .sort_values('sample_number')
          .reset_index(drop=True))
    unstable = load_dynamic_stability()

    missing = sorted(set(df['formula']) - set(unstable.index))
    if missing:
        raise ValueError(f"no phonon data for: {missing}")

    df['dynamically_unstable'] = df['formula'].map(unstable)

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig = make_figure(df)
    fig.savefig(OUT_PDF)
    plt.close(fig)

    counts = df['XRD Result'].value_counts()
    print(f"Plotted {len(df)} samples "
          f"(SS: {counts.get('SS', 0)}  MP: {counts.get('MP', 0)}  "
          f"P: {counts.get('P', 0)})")
    print(f"Dynamically unstable ground states: {int(df['dynamically_unstable'].sum())}")
    for _, r in df[df['dynamically_unstable']].sort_values('formula').iterrows():
        print(f"  {int(r['sample_number']):04d}  {r['formula']:<12} {r['XRD Result']}")
    print('Caveats flagged on the plot:')
    for _, r in df[df['caveat'].notna()].sort_values('sample_number').iterrows():
        print(f"  {int(r['sample_number']):04d}  {r['formula']:<12} {r['caveat']}")
    print(f"Saved: {OUT_PDF}")


if __name__ == '__main__':
    main()
