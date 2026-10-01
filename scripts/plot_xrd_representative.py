#!/usr/bin/env python3
"""Body-text figure: one representative XRD comparison per kind of outcome.

  results/publication_ready/xrd_representative.pdf

Six full-width panels, two across by three down on an A4 page, each showing
what was predicted (red, bottom), what was measured (blue) and, where we have
one, the structure we believe formed (green, top).  The six are chosen in
``REPRESENTATIVE`` below to span the ways a prediction can play out, from a
solid solution with no trace of the predicted ordering to a probable new
compound.  Every other target appears in the appendix grids drawn by
plot_xrd_comparison_grid.py; the two scripts share their panel machinery
through synthesizability.xrd_figure so the same sample looks the same in both.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, 'src')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from synthesizability.paper_samples import load_classified_samples
from synthesizability.xrd_figure import (
    CM, FIG_W_CM, BG_WINDOW, FWHM, TWO_THETA_RANGE,
    build_record, draw_panel, label_outer_axes, page_height_cm, panel_letter,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUT_PDF = Path('results/publication_ready/xrd_representative.pdf')

# (sample number, panel title) in reading order: left to right, top to bottom.
# The title is the kind of outcome the panel stands for.
REPRESENTATIVE = [
    # Al4FeCo3: predicted as an Im-3m ordering of Al, Fe and Co on bcc.  What
    # formed is B2 with Al on one site (0.92) and Fe/Co sharing the other; the
    # Fe/Co ordering of the prediction is undetectable (no scattering contrast).
    (545, 'B2 solid solution: Al ordered, Fe and Co mixed', None),
    # Sc2CuPd: predicted full Heusler.  B2 with a 0.89 Sc site fits; the L2_1
    # (111) line at 23.6 deg is absent.  The inset shows the missing line.
    (474, r'Predicted L2$_1$ ordering absent: B2 solid solution', (21.5, 26)),
    # LaNiSi3: low disorder probability, and a clean multiphase mixture of
    # known La-Ni-Si phases (Jade fit shown).
    (514, 'Multiphase mixture of known phases', None),
    # CoSiY: the target composition formed, but as the known Pnma polymorph
    # ~74 meV above the predicted P4/nmm ground state.
    (538, 'Known polymorph of the target formed', None),
    # CuSnY: the predicted P6_3mc structure formed, and is PDF 04-014-1199.
    (507, 'Predicted structure formed; already known', None),
    # TiSiPt5: the predicted P4/mmm structure appears to have formed, with a
    # Pt-rich fcc impurity; no PDF entry.
    (493, 'Predicted structure formed; probably novel', None),
]

NCOLS = 2
CAPTION_LINES = 9             # this caption has six panels to describe

# Margins around the axes block, in cm.
LEFT_CM = 1.0
RIGHT_CM = 0.25
TOP_CM = 0.45                 # the panel titles sit in the top gap
BOTTOM_CM = 0.8
HSPACE = 0.20
WSPACE = 0.06

FONTSIZE = 7.5


def make_figure(records):
    n = len(records)
    nrows = -(-n // NCOLS)
    fig_h_cm = page_height_cm(CAPTION_LINES)
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(FIG_W_CM * CM, fig_h_cm * CM),
                             sharex=True, squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)

    for i, (rec, title, inset) in enumerate(records, start=1):
        ax = axes[(i - 1) // NCOLS, (i - 1) % NCOLS]
        ax.set_visible(True)
        draw_panel(ax, rec, panel_letter(i), fontsize=FONTSIZE,
                   label_traces=True, title=title, inset=inset, minor=5)

    label_outer_axes(axes, fontsize=FONTSIZE)
    fig.subplots_adjust(left=LEFT_CM / FIG_W_CM, right=1 - RIGHT_CM / FIG_W_CM,
                        top=1 - TOP_CM / fig_h_cm, bottom=BOTTOM_CM / fig_h_cm,
                        hspace=HSPACE, wspace=WSPACE)
    return fig, (FIG_W_CM, fig_h_cm)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT_PDF)
    args = parser.parse_args()

    df = load_classified_samples().set_index('sample_number', drop=False)
    print(f'Range {TWO_THETA_RANGE[0]:g}-{TWO_THETA_RANGE[1]:g} deg, '
          f'simulated FWHM {FWHM} deg, SNIP window {BG_WINDOW} deg\n')

    records = []
    for num, title, inset in REPRESENTATIVE:
        if num not in df.index:
            raise SystemExit(f'{num:04d} is not among the classified samples')
        rec = build_record(df.loc[num])
        if rec is None:
            raise SystemExit(f'{num:04d} has no usable measurement')
        records.append((rec, title, inset))

    fig, (w_cm, h_cm) = make_figure(records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    plt.close(fig)
    print(f'\nSaved: {args.out}  ({len(records)} panels, {w_cm:.2f} x {h_cm:.2f} cm)')


if __name__ == '__main__':
    main()
