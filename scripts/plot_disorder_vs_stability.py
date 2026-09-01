#!/usr/bin/env python3
"""
Publication figure: OQMD hull distance vs. disorder parameter for the arc-melted
targets, coloured by XRD outcome.

  results/publication_ready/disorder_vs_stability.pdf

Points are the 33 non-Diffusion-Model targets that have an XRD classification
(SS / MP / P).  An "X" is drawn over any point whose most thermodynamically
stable OQMD polymorph (polymorph_rank == 1) is dynamically unstable, i.e. has
imaginary phonon modes.

Each point is labelled with its target formula.

Sample classification follows scripts/plot_xrd_grids.py so the two figures
always describe the same set of samples.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from adjustText import adjust_text
except ImportError:  # pragma: no cover
    raise SystemExit('adjustText is required for the formula labels: '
                     'pip install adjustText')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SYNTH_CSV   = Path('data/processed/synthesis_data.csv')
FINAL_CSV   = Path('data/temp/Arc Melted Samples - Final Analysis.csv')
SHEET1_CSV  = Path('data/temp/Arc Melted Samples - Sheet1.csv')
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

X_SIZE  = 34            # pt^2 for the dynamic-instability overlay
X_WIDTH = 1.1

LABEL_SIZE  = 5.5       # pt — formula labels
LABEL_COLOR = '#333333'
LEADER_WIDTH = 0.35     # pt — leader lines from label to point


# Allen electronegativity (configuration energy).  Chosen over Pauling because
# it is the scale under which Si, Ge and Sn sit above the late transition
# metals, giving the ordering conventional for these intermetallics
# (ZrPt5Si, not ZrSiPt5).  La and Gd are estimates: Allen's tabulation omits
# the f block, but both are the most electropositive element in every formula
# here, so the ordering does not depend on their exact values.
ALLEN_ELECTRONEGATIVITY = {
    'Al': 1.613, 'Au': 1.92,  'B':  2.051, 'C':  2.544, 'Co': 1.84,
    'Cr': 1.65,  'Cu': 1.85,  'Fe': 1.80,  'Gd': 1.20,  'Ge': 1.994,
    'Hf': 1.16,  'Ir': 1.68,  'La': 1.10,  'Mn': 1.75,  'Mo': 1.47,
    'Nb': 1.41,  'Ni': 1.88,  'Pd': 1.58,  'Pt': 1.72,  'Ru': 1.54,
    'Sc': 1.19,  'Si': 1.916, 'Sn': 1.824, 'Ti': 1.38,  'Y':  1.12,
    'Zr': 1.32,
}


def parse_formula(formula):
    return [(el, n) for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula) if el]


def sort_by_electronegativity(formula):
    """'ZrSiPt5' -> 'ZrPt5Si': elements in order of increasing electronegativity."""
    parts = parse_formula(formula)
    missing = [el for el, _ in parts if el not in ALLEN_ELECTRONEGATIVITY]
    if missing:
        raise KeyError(f'no Allen electronegativity for {missing} '
                       f'(formula {formula}); add it to ALLEN_ELECTRONEGATIVITY')
    parts.sort(key=lambda p: (ALLEN_ELECTRONEGATIVITY[p[0]], p[0]))
    return ''.join(el + n for el, n in parts)


def formula_to_mathtext(formula):
    """'ZrSiPt5' -> '$\\mathrm{ZrPt_{5}Si}$': reordered, with subscripts."""
    return ('$\\mathrm{'
            + ''.join(f'{el}_{{{n}}}' if n else el
                      for el, n in parse_formula(sort_by_electronegativity(formula)))
            + '}$')


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_classified_samples():
    """The 33 targets with an XRD outcome, as used by plot_xrd_grids.py."""
    final = pd.read_csv(FINAL_CSV)
    sheet1 = pd.read_csv(SHEET1_CSV)
    synth = pd.read_csv(SYNTH_CSV)

    final['Sample Number'] = pd.to_numeric(final['Sample Number'], errors='coerce')
    sheet1['Sample Number'] = pd.to_numeric(sheet1['Sample Number'], errors='coerce')
    sheet1 = sheet1.dropna(subset=['Sample Number'])
    sheet1['Sample Number'] = sheet1['Sample Number'].astype(int)

    nums = final['Sample Number'].dropna().astype(int).tolist()
    sub = synth[synth['sample_number'].isin(nums) &
                (synth['prediction_list'] != 'Diffusion Model')].copy()
    sub = sub.join(sheet1.set_index('Sample Number')[['XRD Result']],
                   on='sample_number')
    return sub[sub['XRD Result'].isin(STYLE)].copy()


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

    for outcome, (colour, marker, _) in STYLE.items():
        g = df[df['XRD Result'] == outcome]
        ax.scatter(g['disorder_probability'], g['oqmd_stability'],
                   s=MARKER_SIZE, c=colour, marker=marker, alpha=ALPHA,
                   edgecolors='black', linewidths=EDGE_WIDTH, zorder=2)

    unstable = df[df['dynamically_unstable']]
    ax.scatter(unstable['disorder_probability'], unstable['oqmd_stability'],
               s=X_SIZE, c='black', marker='x', linewidths=X_WIDTH, zorder=3)

    texts = [
        ax.text(r['disorder_probability'], r['oqmd_stability'],
                formula_to_mathtext(r['formula']),
                fontsize=LABEL_SIZE, color=LABEL_COLOR, zorder=4,
                ha='center', va='center',
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

    counts = df['XRD Result'].value_counts()
    handles = [
        Line2D([0], [0], linestyle='none', marker=marker, color=colour,
               markeredgecolor='black', markeredgewidth=EDGE_WIDTH,
               markersize=np.sqrt(MARKER_SIZE), alpha=ALPHA,
               label=f'{label} (n={counts.get(outcome, 0)})')
        for outcome, (colour, marker, label) in STYLE.items()
    ]
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
    adjust_text(texts, ax=ax,
                x=df['disorder_probability'].to_numpy(),
                y=df['oqmd_stability'].to_numpy(),
                objects=[legend.get_frame()],
                expand=(1.45, 1.6),
                force_text=(0.5, 0.65),
                force_static=(0.95, 1.05),
                force_pull=(0.04, 0.04),
                max_move=26,
                iter_lim=3000,
                arrowprops=dict(arrowstyle='-', color=LABEL_COLOR,
                                lw=LEADER_WIDTH, shrinkA=1, shrinkB=3))

    return fig


def main():
    df = load_classified_samples()
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
    print(f"Saved: {OUT_PDF}")


if __name__ == '__main__':
    main()
