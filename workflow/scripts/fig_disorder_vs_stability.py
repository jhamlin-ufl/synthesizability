#!/usr/bin/env python3
"""Fig. 2: OQMD Delta H_f - C(x) against P_disorder, by XRD outcome category.

  in:  tables/{paper_samples,outcomes}.csv, computed/p_disorder.csv,
       results/hull_energies.csv (recomputed OQMD value, as in Table 1),
       results/phonon_flags.csv
  out: results/figures/disorder_vs_stability.{pdf,png}, disorder_vs_stability_points.csv

Symbols: one colour family per outcome family (green = a predicted structure formed,
blue = solid solution formed, dark grey = on-hull structure not formed); within a
family the first category is filled and the second open; a split marker shows the
two categories that the data cannot distinguish; * / dagger / double dagger after a
label mark previously reported structures, compositions and series; a red X with a
thin white outline marks a target whose every structure below the hull of the other
compositions (OQMD, MP or Alexandria PBE) is dynamically unstable.  The category
legend sits above the axes in three family columns.  Formula labels are placed by
adjustText, with hand placements in LABEL_PLACE, then drawn under the symbols and
pushed clear of every symbol.
"""
import sys
from pathlib import Path

sys.path.insert(0, 'workflow/src')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.colors import to_rgb
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.transforms import Affine2D, Bbox, ScaledTranslation
import matplotlib.patheffects as pe

from synthesizability.formula import mathtext
from synthesizability.samples import CATEGORY_NAME, load_targets

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PHONON_CSV = Path('results/phonon_flags.csv')
HULL_CSV = Path('results/hull_energies.csv')
OUT_PDF = Path('results/figures/disorder_vs_stability.pdf')
OUT_POINTS = Path('results/figures/disorder_vs_stability_points.csv')

# colour families checked for deuteranopia, protanopia and tritanopia (Machado 2009)
GREEN, BLUE, DARK = '#1A9850', '#3B6FB6', '#404040'
# category -> (colour, marker, filled)
STYLE = {
    1: (GREEN, 'o', True),
    2: (GREEN, 'o', False),
    3: (BLUE,  's', True),
    4: (BLUE,  's', False),
    5: (DARK,  '^', True),
    6: (DARK,  '^', False),
}
FAMILY_HEAD = ['A predicted structure formed', 'Solid solution formed', 'On-hull structure not formed']
SHORT = {1: 'on-hull structure', 2: 'above-hull polymorph',
         3: 'disordered parent of on-hull', 4: 'other solid solution',
         5: 'multiphase', 6: 'unresolved'}
SIZE = 95            # pt^2, as scatter's s
EDGE_OPEN = 1.6      # pt, thick edge of the open markers
EDGE_FILLED = 0.9    # pt, darker-shade edge of the filled markers
DARKEN = 0.55        # fill * DARKEN = edge colour of a filled marker
REPORTED_MARK = {'yes': '*', 'composition': '†', 'series': '‡', 'no': ''}

X_SIZE = SIZE        # pt^2 for the dynamic-instability overlay (red X, white outline)
X_WIDTH = 1.4
X_COLOR = '#D7191C'
X_HALO = [pe.withStroke(linewidth=2.4, foreground='white')]

XLIM = (-0.11, 1.09)
YLIM = (-80, 25)                # meV/atom
FIGSIZE = (7.2, 5.75)
AXES_RECT = [0.1530, 0.106, 0.8080, 0.735]
LEGEND_FONTSIZE = 6.8

LABEL_SIZE = 5.5
LABEL_COLOR = '#333333'
LEADER_WIDTH = 0.35   # pt
MIN_LEADER_PX = 16    # label centre closer than this to its point: no leader

# Hand placements (data coordinates: disorder parameter, eV/atom; drawn in
# meV/atom) for the
# labels adjust_text cannot place unambiguously; every other label is placed by
# adjust_text.  Every label whose centre is farther than MIN_LEADER_PX from its
# point gets a leader line, so a label never has to be guessed to a symbol.
LABEL_PLACE = {
    513: (-0.050, +0.0038),    # LaRuC2
    555: (-0.050, -0.0018),    # GdFeGe2
    554: (-0.050, -0.0095),    # Gd2Fe3Ge5
    539: (-0.050, -0.0163),    # YCoSi3
    553: (-0.045, -0.0232),    # GdFeGe3
    510: (+0.031, +0.0158),    # NbNi2B2
    537: (+0.100, -0.0022),    # GdSn6Co4
    512: (+0.105, -0.0094),    # ScCrC2
    514: (+0.072, -0.0292),    # LaNiSi3
    493: (+0.112, -0.0140),    # TiPt5Si
    491: (+0.128, -0.0222),    # TiPt5Ge
    540: (+0.242, -0.0108),    # YAl3Cu2
    480: (+0.238, -0.0340),    # GdPt2Ge
    561: (+0.252, -0.0440),    # ZrPt5Si
    538: (+0.085, -0.0452),    # YCoSi
    507: (+0.205, -0.0691),    # YSnCu
    494: (+0.165, +0.0012),    # Ti2Pt8Si
    492: (+0.318, -0.0045),    # Ti2Pt8Ge
    472: (+0.312, -0.0192),    # ScPdCu4
    556: (+0.312, -0.0288),    # Gd3Sn6Co
    483: (+0.695, -0.0048),    # La5Al3Co2
    469: (+0.749, +0.0215),    # ScPd2Cu
    477: (+0.772, +0.0068),    # ScPd5Au2
    473: (+0.800, -0.0046),    # Sc3Pd2Cu
    475: (+0.842, +0.0118),    # La2PdAu2
    474: (+0.955, +0.0062),    # Sc2PdCu
    545: (+0.962, -0.0084),    # Al4FeCo3
    535: (+0.968, -0.0016),    # Al8FeCo7
    466: (+0.830, -0.0175),    # Hf2MoIr (point at the recomputed -12.0, not the stored -18.6)
    465: (+0.985, -0.0228),    # Pd6Cu5Ge
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_points():
    """The targets with x = P_disorder, y = recomputed OQMD Delta H_f - C(x) (meV/atom)
    and the dynamical-stability cross."""
    df = load_targets()
    hull = pd.read_csv(HULL_CSV)[['sample_number', 'oqmd_hull_meV']]
    flags = pd.read_csv(PHONON_CSV)[['sample_number', 'all_below_hull_unstable']]
    df = df.merge(hull, on='sample_number', how='left').merge(flags, on='sample_number', how='left')
    if df[['oqmd_hull_meV', 'all_below_hull_unstable']].isna().any().any():
        raise ValueError('missing hull energy or phonon flag')
    return df.rename(columns={'oqmd_hull_meV': 'hull_meV', 'all_below_hull_unstable': 'dynamically_unstable'})


# ---------------------------------------------------------------------------
# Markers
# ---------------------------------------------------------------------------

def darker(colour):
    return tuple(c * DARKEN for c in to_rgb(colour))


def face_edge(cat):
    colour, _, filled = STYLE[cat]
    if filled:
        return colour, darker(colour), EDGE_FILLED
    return 'white', colour, EDGE_OPEN


def marker_path(cat):
    """Unit-sized path of the category's marker, as scatter would draw it."""
    ms = MarkerStyle(STYLE[cat][1])
    return ms.get_path().transformed(ms.get_transform())


def point_frame(ax, x, y):
    """Transform from marker units (extent 1 = sqrt(SIZE) pt) to the point."""
    side_in = np.sqrt(SIZE) / 72.0
    return (Affine2D().scale(side_in) + ax.figure.dpi_scale_trans
            + ScaledTranslation(x, y, ax.transData))


def draw_marker(ax, x, y, cat, alt=None, zorder=3):
    """One sample.  With `alt`, the marker is split into `cat` and `alt`, the two
    categories the data cannot distinguish, both in the shape of `cat`."""
    frame = point_frame(ax, x, y)
    path = marker_path(cat)
    if alt is None:
        fc, ec, lw = face_edge(cat)
        p = PathPatch(path, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
        p.set_transform(frame)
        p.set_clip_on(False)
        ax.add_patch(p)
        return
    # geometric halves (clip paths proved unreliable here): the marker path as
    # a polygon, cut at x = 0; each half keeps its own face and edge colours
    # and the edge along the cut is the divider
    poly = MPath(path.to_polygons()[0])
    halves = ((cat, poly.clip_to_bbox(Bbox([[-1, -1], [0, 1]]))),
              (alt, poly.clip_to_bbox(Bbox([[0, -1], [1, 1]]))))
    for c, half in halves:
        fc, ec, lw = face_edge(c)
        p = PathPatch(half, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
        p.set_transform(frame)
        p.set_clip_on(False)
        ax.add_patch(p)


def legend_marker(cat):
    fc, ec, lw = face_edge(cat)
    return Line2D([0], [0], linestyle='none', marker=STYLE[cat][1],
                  markerfacecolor=fc, markeredgecolor=ec, markeredgewidth=lw,
                  markersize=np.sqrt(SIZE))


def legend_split(cat, alt):
    """Two Line2D halves plus an outline, for HandlerTuple(ndivide=1)."""
    fc1, ec1, lw1 = face_edge(cat)
    fc2, ec2, lw2 = face_edge(alt)
    m = STYLE[cat][1]
    left = Line2D([0], [0], linestyle='none', marker=m, fillstyle='left',
                  markerfacecolor=fc1, markerfacecoloralt='none',
                  markeredgecolor='none', markeredgewidth=0, markersize=np.sqrt(SIZE))
    right = Line2D([0], [0], linestyle='none', marker=m, fillstyle='right',
                   markerfacecolor=fc2, markerfacecoloralt='none',
                   markeredgecolor='none', markeredgewidth=0, markersize=np.sqrt(SIZE))
    edge = Line2D([0], [0], linestyle='none', marker=m, markerfacecolor='none',
                  markeredgecolor=ec2 if STYLE[alt][2] is False else ec1,
                  markeredgewidth=max(lw1, lw2) * 0.6, markersize=np.sqrt(SIZE))
    return (left, right, edge)


def push_labels_clear(fig, ax, texts, df, steps=60):
    """Move any label whose box comes within a symbol's radius of a symbol away from it."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pts = ax.transData.transform(df[['disorder_probability', 'hull_meV']].to_numpy())
    rad = 0.75 * np.sqrt(SIZE) / 72 * fig.dpi
    inv = ax.transData.inverted()
    for _ in range(steps):
        moved = False
        for text in texts:
            bb = text.get_window_extent(renderer=renderer).expanded(1.08, 1.15)
            for p in pts:
                cx, cy = min(max(p[0], bb.x0), bb.x1), min(max(p[1], bb.y0), bb.y1)
                if np.hypot(cx - p[0], cy - p[1]) < rad:
                    v = np.array([(bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2]) - p
                    v = v / (np.hypot(*v) or 1) * 2.0
                    text.set_position(inv.transform(ax.transData.transform(text.get_position()) + v))
                    moved = True
        if not moved:
            break
        fig.canvas.draw()


def blank():
    return Line2D([0], [0], linestyle='none', marker='none')


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_figure(df):
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.7, alpha=0.4, zorder=1)

    for _, r in df.iterrows():
        alt = int(r['alt']) if r['split'] else None
        draw_marker(ax, r['disorder_probability'], r['hull_meV'], int(r['category']), alt)

    unstable = df[df['dynamically_unstable']]
    ax.scatter(unstable['disorder_probability'], unstable['hull_meV'],
               s=X_SIZE, c=X_COLOR, marker='x', linewidths=X_WIDTH, zorder=4, path_effects=X_HALO)

    texts = [
        ax.text(r['disorder_probability'], r['hull_meV'],
                mathtext(r['formula']) + REPORTED_MARK[r['reported']],
                fontsize=LABEL_SIZE, color=LABEL_COLOR, zorder=2.6,   # under the symbols
                ha='center', va='center',
                # invisible on the white ground; masks the dashed zero line
                # where a label happens to sit on it
                bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none', alpha=0.75))
        for _, r in df.iterrows()
    ]

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(np.arange(0.0, 1.01, 0.2))
    ax.set_xlabel(r'$P_{\mathrm{disorder}}$', fontsize=11)
    ax.set_ylabel(r'$\Delta H_{\mathrm{f}} - C(x)$ in OQMD (meV/atom)', fontsize=11)
    ax.tick_params(labelsize=10)

    # category legend: three family columns above the axes, bold heads
    handles, labels = [], []
    for i, head in enumerate(FAMILY_HEAD):
        handles += [blank(), legend_marker(2 * i + 1), legend_marker(2 * i + 2)]
        labels += [head, SHORT[2 * i + 1], SHORT[2 * i + 2]]
    top = fig.legend(handles=handles, labels=labels, ncol=3, loc='upper center',
                     bbox_to_anchor=(0.557, 0.995), fontsize=LEGEND_FONTSIZE + 0.6, frameon=False,
                     handleheight=1.7, columnspacing=2.2, handletextpad=0.5)
    for t, lab in zip(top.get_texts(), labels):
        if lab in FAMILY_HEAD:
            t.set_fontweight('bold')
            t.set_x(-20)
    # conventions in the empty lower right
    x_handle = Line2D([0], [0], linestyle='none', marker='x', color=X_COLOR,
                      markeredgewidth=X_WIDTH, markersize=np.sqrt(X_SIZE), path_effects=X_HALO)
    h2 = [legend_split(3, 1), x_handle, blank(), blank(), blank()]
    l2 = ['split: the data cannot distinguish the two categories',
          'X: every below-hull structure dynamically unstable',
          '* structure previously reported at the target composition',
          '† composition previously reported without a structure',
          '‡ ordered member of a previously reported, partly disordered series']
    legend = ax.legend(handles=h2, labels=l2, loc='lower right', fontsize=LEGEND_FONTSIZE,
                       framealpha=0.85, handler_map={tuple: HandlerTuple(ndivide=1, pad=0)},
                       handleheight=1.7)

    # Automatic placement first (the legend is an obstacle), then the hand
    # placements, then a leader line for every displaced label.
    fig.canvas.draw()
    unknown = set(LABEL_PLACE) - set(df['sample_number'])
    if unknown:
        raise ValueError(f'LABEL_PLACE names samples not in the figure: {sorted(unknown)}')
    adjust_text(texts, ax=ax,
                x=df['disorder_probability'].to_numpy(),
                y=df['hull_meV'].to_numpy(),
                objects=[legend.get_frame()],
                expand=(1.6, 1.75), force_text=(0.55, 0.7), force_static=(1.1, 1.2),
                force_pull=(0.035, 0.035), max_move=28, iter_lim=3000)
    for text, (_, r) in zip(texts, df.iterrows()):
        num = int(r['sample_number'])
        if num in LABEL_PLACE:
            px, py = LABEL_PLACE[num]
            text.set_position((px, 1000 * py))
    push_labels_clear(fig, ax, texts, df)
    for text, (_, r) in zip(texts, df.iterrows()):
        pos = text.get_position()
        point = (r['disorder_probability'], r['hull_meV'])
        px = ax.transData.transform(pos) - ax.transData.transform(point)
        if np.hypot(*px) < MIN_LEADER_PX:
            continue
        leader = FancyArrowPatch(posA=pos, posB=point, patchA=text, arrowstyle='-',
                                 shrinkA=1, shrinkB=3, color=LABEL_COLOR,
                                 lw=LEADER_WIDTH, zorder=5, transform=ax.transData)
        ax.add_patch(leader)

    return fig


def main():
    df = load_points()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig = make_figure(df)
    fig.savefig(OUT_PDF, metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(OUT_PDF.with_suffix('.png'), dpi=200)
    plt.close(fig)
    df[['sample_number', 'formula', 'disorder_probability', 'hull_meV', 'category', 'alt', 'reported',
        'dynamically_unstable']].to_csv(OUT_POINTS, index=False)
    print(f'{len(df)} points -> {OUT_PDF}; crossed: ' + ', '.join(df.loc[df.dynamically_unstable, 'formula']))


if __name__ == '__main__':
    main()
