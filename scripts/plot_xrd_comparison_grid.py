#!/usr/bin/env python3
"""Appendix figures: measured vs predicted XRD for every arc-melted target.

  results/publication_ready/xrd_appendix_SS.pdf   solid solutions
  results/publication_ready/xrd_appendix_MP.pdf   multiphase mixtures
  results/publication_ready/xrd_appendix_P.pdf    predicted structure found

One figure per XRD outcome -- the three point types of the disorder-vs-stability
plot -- and one panel per sample within it, so that the six samples singled out
in the body of the paper (plot_xrd_representative.py) can be seen in the company
of the rest of their class.  What a panel shows (file choice, background
removal, simulated structures, the green solution trace) is decided in
synthesizability.xrd_figure, shared with the body figure.

Every figure is the width of the A4 text block and every panel row has the same
height, so the three figures read as a set.  The height of each figure follows
from its row count, capped at the page.  Column counts are chosen per outcome
below to keep the panels wide: a diffraction pattern wants width more than
height.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, 'src')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from synthesizability.paper_samples import OUTCOMES, OUTCOME_LABEL, load_classified_samples
from synthesizability.xrd_figure import (
    CM, FIG_W_CM, BG_WINDOW, FWHM, TWO_THETA_RANGE,
    build_records, draw_panel, label_outer_axes, legend_handles, page_height_cm,
    panel_letter,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUT_DIR = Path('results/publication_ready')
OUT_NAME = 'xrd_appendix_{outcome}.pdf'

# Columns per outcome figure.  12 SS panels in 2 columns and 16 MP panels in 3
# fill a page at six rows each; the 5 P panels take half a page in 2 columns,
# with the legend in the spare cell.
NCOLS = {'SS': 2, 'MP': 3, 'P': 2}

CAPTION_LINES = 6             # room left under a full-page figure
FULL_PAGE_ROWS = 6            # rows a full-page figure holds: sets the row height

# Margins around the axes block, in cm.  TOP_LEGEND_CM replaces TOP_CM when the
# grid has no spare cell to hold the legend.
LEFT_CM = 1.0
RIGHT_CM = 0.25
TOP_CM = 0.25
TOP_LEGEND_CM = 0.85
BOTTOM_CM = 0.75
HSPACE = 0.16
WSPACE = 0.07

FONTSIZE = 6.5


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def layout_slots(n, ncols, nrows):
    """Map panel index -> (row, col), leaving the spare slots in the upper right.

    The top row is the short one, so the empty cells and the legend that goes in
    them end up above the grid rather than stranded at the bottom.
    """
    if n > ncols * nrows:
        raise ValueError(f'{n} panels do not fit in a {nrows}x{ncols} grid')
    n_top = n - (nrows - 1) * ncols          # panels in the (partial) top row
    slots = []
    for idx in range(n):
        if idx < n_top:
            slots.append((0, idx))
        else:
            j = idx - n_top
            slots.append((1 + j // ncols, j % ncols))
    return slots, n_top


def figure_geometry(nrows, legend_band):
    """(fig_w_in, fig_h_in, subplots_adjust kwargs) for an nrows-deep grid."""
    page_h = page_height_cm(CAPTION_LINES)
    top = TOP_LEGEND_CM if legend_band else TOP_CM
    # Row height of a full page with the plain top margin; the legend band, if
    # any, comes out of the page height rather than growing the figure.
    row_h = (page_h - TOP_CM - BOTTOM_CM) / FULL_PAGE_ROWS
    fig_h = min(page_h, nrows * row_h + top + BOTTOM_CM)
    adjust = dict(left=LEFT_CM / FIG_W_CM, right=1 - RIGHT_CM / FIG_W_CM,
                  top=1 - top / fig_h, bottom=BOTTOM_CM / fig_h,
                  hspace=HSPACE, wspace=WSPACE)
    return FIG_W_CM * CM, fig_h * CM, adjust


def make_figure(records, ncols):
    n = len(records)
    nrows = math.ceil(n / ncols)
    slots, n_top = layout_slots(n, ncols, nrows)
    legend_band = n_top == ncols            # no spare cell for the legend

    fig_w, fig_h, adjust = figure_geometry(nrows, legend_band)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharex=True,
                             squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)

    for i, (rec, (r, c)) in enumerate(zip(records, slots), start=1):
        ax = axes[r, c]
        ax.set_visible(True)
        draw_panel(ax, rec, panel_letter(i), fontsize=FONTSIZE)

    label_outer_axes(axes, fontsize=FONTSIZE)
    fig.subplots_adjust(**adjust)

    handles = legend_handles(records)
    if legend_band:
        fig.legend(handles=handles, loc='upper center', ncols=len(handles),
                   bbox_to_anchor=(0.5, 1.0), fontsize=FONTSIZE + 1,
                   frameon=False, handlelength=1.6, columnspacing=1.6)
    else:
        empty = [axes[0, c].get_position() for c in range(n_top, ncols)]
        centre = (0.5 * (empty[0].x0 + empty[-1].x1),
                  0.5 * (empty[0].y0 + empty[0].y1))
        fig.legend(handles=handles, loc='center', bbox_to_anchor=centre,
                   fontsize=FONTSIZE + 1, frameon=False, handlelength=1.6,
                   labelspacing=0.7, borderaxespad=0)
    return fig, (nrows, ncols), (fig_w / CM, fig_h / CM)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outcome', choices=OUTCOMES, action='append',
                        help='build only this outcome (repeatable); default all')
    parser.add_argument('--out-dir', type=Path, default=OUT_DIR)
    args = parser.parse_args()
    outcomes = args.outcome or list(OUTCOMES)

    df = load_classified_samples()
    print(f'{len(df)} classified samples '
          + '  '.join(f'{o}: {(df["XRD Result"] == o).sum()}' for o in OUTCOMES))
    print(f'Range {TWO_THETA_RANGE[0]:g}-{TWO_THETA_RANGE[1]:g} deg, '
          f'simulated FWHM {FWHM} deg, SNIP window {BG_WINDOW} deg\n')

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for outcome in outcomes:
        print(f'== {outcome}: {OUTCOME_LABEL[outcome]}')
        records = build_records(df[df['XRD Result'] == outcome])
        if not records:
            raise SystemExit(f'no records to plot for {outcome}')

        fig, (nrows, ncols), (w_cm, h_cm) = make_figure(records, NCOLS[outcome])
        out = args.out_dir / OUT_NAME.format(outcome=outcome)
        fig.savefig(out)
        plt.close(fig)

        missing = [r['sample_number'] for r in records if r['predicted'] is None]
        if missing:
            print(f'  No OQMD structure for: {missing}')
        print(f'  Saved: {out}  ({len(records)} panels, {nrows}x{ncols} grid, '
              f'{w_cm:.2f} x {h_cm:.2f} cm)\n')


if __name__ == '__main__':
    main()
