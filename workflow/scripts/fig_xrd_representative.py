#!/usr/bin/env python3
"""Body-text figure: representative XRD panels, one per outcome the text discusses.

  results/publication_ready/xrd_representative.pdf

Six panels, two across by three down (PANELS below; the list is the only place
that decides which samples appear, and the layout also works with five, the
legend then taking the empty slot).  Each panel stacks, from the top: the
measured pattern with the whole-pattern main fit (the refinement behind
the phases and weight fractions printed on the panel), the difference curve on
the same scale, and short rows of Bragg positions, one per phase (the target
phase in green).  A small picture of the formed crystal structure
(results/publication_ready/structure_insets/<NNNN>.png, made from the refined
cell) sits in an empty part of each panel.  The disordered-parent panel also
carries the pattern simulated from the predicted ordered structure, with a
zoom on its missing superlattice line.

Which fit is the main fit for each panel is a data file,
workflow/figure_inputs/xrd_representative_fits.csv (one row per sample; the
choices and the R_wp numbers behind them are in its note column).  Fits come
from three kinds of source:
  refinements / verdict_extra  pipeline exports (<fit>.txt + <fit>_ticks.json),
                               statistics from the pipeline tables;
  order_ladder                 GSAS-II order-parameter ladder exports, Bragg
                               positions from the CIF beside the export;
  fit_folder                   a refinement folder (fit.txt, ticks.json,
                               summary.json with record.wR and
                               record.state.wt in tick-file phase order).

Text, the structure pictures and the zoom inset are placed automatically: each
is tried at candidate positions in display coordinates and kept only where it
touches no data point, fit, difference curve, tick row or other item; if an
item fits nowhere the panel's headroom grows and the panel is redrawn.  The
script checks the final layout again and fails if anything overlaps.

The output is vector throughout except the structure pictures.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator

from synthesizability.parsers.xrd import parse_xrd_file
from synthesizability.xrd_figure import (
    CM, COLOR_SOLUTION, FIG_W_CM, FWHM, XLABEL, YLABEL,
    _mathtext_space_group, find_oqmd_cif, panel_letter, space_group_symbol,
)
from synthesizability.xrd_simulator import peak_list, simulate_pattern
from xrd_panel_helpers import (                                   # noqa: E402
    COLOR_FIT, COLOR_TICK_OTHER, REPORTED_MARK, load_fit, lookup_fit_stats, parse_phases,
)

def load_samples():
    """Paper samples keyed by the button numbers the panels use (target and record button)."""
    ps = pd.read_csv('tables/paper_samples.csv', dtype={'sample': str, 'record_sample': str})
    tg = pd.read_csv('databases/snapshot/targets.csv')
    ps = ps.merge(tg[['formula', 'oqmd_entry_id']], on='formula', how='left')
    rows = []
    for _, r in ps.iterrows():
        for num in {int(r['sample']), int(r['record_sample'])}:
            rows.append(dict(r, sample_number=num))
    return pd.DataFrame(rows).drop_duplicates('sample_number').set_index('sample_number', drop=False)


CATEGORY_NAME = dict(pd.read_csv('tables/paper_samples.csv')[['category', 'category_name']]
                     .drop_duplicates().itertuples(index=False, name=None))


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUT_PDF = Path('results/figures/xrd_representative.pdf')
FITS_CSV = Path('workflow/figure_inputs/xrd_representative_fits.csv')
PHASEID = Path('.')
STRUCTURE_DIR = Path('workflow/figure_inputs/structure_insets')
SAMPLES_DIR = Path('samples')   # samples/<formula>/buttons/<NNNN>/<scan>

# Reading order; one entry per panel.  ``sample`` is the sample_number of the
# row in representative_fits.csv (and of the sample table: 0615 is a remake of
# paper sample 0540, so its panel is keyed 540 with fit_sample 615).  ``name``
# is the formula as the text writes it.  ``show_wt=False`` leaves the weight
# fractions off the phase list; ``split_list`` places the lines after the
# numbered phases (unindexed marker, R_wp) as a block of their own; ``unindexed`` marks observed lines of an
# unidentified phase (angles in deg).  ``predicted`` adds the pattern
# simulated from the OQMD prediction below the ticks; ``inset`` = (lo, hi)
# zooms on a predicted line that is absent.  The structure picture is
# STRUCTURE_DIR/<fit_sample>.png unless ``structure`` names another file.
PANELS = [
    dict(sample=480, name=r'GdPt$_2$Ge'),
    dict(sample=493, name=r'TiPt$_5$Si'),
    dict(sample=561, name=r'ZrPt$_5$Si'),
    dict(sample=539, name=r'YCoSi$_3$'),
    # 0615: an unidentified phase is not modelled, so no weight fractions;
    # its clear lines (positive residual > 5 % of the strongest line where the
    # fit calculates almost nothing) are marked.
    # ``diff_scale`` shrinks a large difference curve (labelled on the panel)
    # to leave room for the structure picture.
    dict(sample=540, name=r'YAl$_3$Cu$_2$', show_wt=False, split_list=True, diff_scale=0.2,
         unindexed=[24.28, 29.06, 35.20, 36.54, 41.22, 42.22, 46.28, 60.14]),
    # 0466 annealed (zero-background scan; replaced 0610/0478 Sc2Pd5Au):
    # the zoom covers the absent L2_1 (111) Mo/Ir-order line and
    # the weak (200) Hf-order line (``inset_lines`` labels them; ``inset_fit``
    # also draws the fit there).  The two lines no phase indexes are marked.
    dict(sample=466, name=r'Hf$_2$MoIr', predicted=True, inset=(21.5, 28.5),
         inset_size=(0.20, 0.24), inset_fit=True, inset_wrap=True, split_list=True,
         inset_lines=[(23.6, '(111)\nMo/Ir'), (27.3, '(200)\nHf')],
         unindexed=[33.46, 36.60], unindexed_label='unindexed lines'),
]

TWO_THETA = (5.0, 120.0)      # common window; fits that ran to 150 deg are clipped
NCOLS = 2
NROWS = 3

# The manuscript includes the figure at the text width, capped at 0.86 of the
# text height with the aspect ratio kept (revtex4-2), so the height must stay
# below about 1.15 times the width for the width to govern.
FIG_H_CM = 16.5

LEFT_CM = 0.85
RIGHT_CM = 0.15
TOP_CM = 1.05                 # legend row, then the first row of panel titles
TOP_CM_NO_LEGEND = 0.5        # five panels: the legend takes the empty slot
BOTTOM_CM = 0.95
HSPACE = 0.16
WSPACE = 0.07

FONTSIZE = 8.5
RASTER_DPI = 450
PNG_DPI = 400
YLABEL_NET = r'$I - I_{\mathrm{bg}}$ (arb. units)'   # every panel: fitted background subtracted
SMALL = 7.5

# Colours (validated with the dataviz palette check, light background: all
# pairs pass the normal-vision floor; red fit vs green ticks is the one
# deuteranope-close pair and they differ in mark type).
COLOR_DIFF = '#0072B2'        # difference curve: blue
COLOR_PRED = '#7B2D8E'        # simulated prediction: purple
COLOR_UNINDEXED = '#6D4C41'   # markers on unindexed lines: dark brown
MS_UNINDEXED = 3.6

MS_DATA = 2.6                 # data crosses
MEW_DATA = 0.55
LW_FIT = 1.0
LW_DIFF = 0.75
LW_TICK = 0.7
LW_SIM = 0.9

# Vertical layout in data units (net pattern scaled to a maximum of 1), except
# the tick rows, which are fixed fractions of the panel height.
TICK_FRAC = 0.028             # height of one row of Bragg positions
TICK_GAP_FRAC = 0.010         # gap between rows
GAP = 0.03                    # between bands
SIM_FRAC = 0.14               # simulated prediction band, fraction of panel height
HEADROOM0 = 0.03              # starting headroom above the strongest point
HEADROOM_STEP = 0.05
HEADROOM_MAX = 2.0

# Structure picture: top-right corner of every panel (right of the strongest
# line).  All pictures share one scale (pixels per angstrom), so an atom is the
# same size in every panel: the largest scale at which every panel fits without
# pushing its pattern down by more than STRUCT_HEADROOM_ALLOW over what its other
# items need.  The key goes under the picture or beside it, whichever fits.
STRUCT_K_MIN, STRUCT_K_MAX, STRUCT_K_STEP = 3.0, 14.0, 0.5   # pixels (100 dpi) per angstrom
STRUCT_K_SPREAD = 1.3         # largest / smallest picture scale
STRUCT_ANCHOR = (0.992, 0.985)  # its top-right corner, axes fraction (all panels)
STRUCT_HEADROOM_ALLOW = 0.45
KEY_ROW_PX = 28               # key row under each picture (pixels at 100 dpi)
KEY_FONT = 6.5
KEY_DOT_PX = 7.0              # element sphere diameter in the key
KEY_GAP_PX = 5.0
KEY_AXIS_PX = 10.0            # length of an in-plane axis arrow
KEY_LINE_PX = 11.0            # line spacing of a key set beside the picture
KEY_AXIS_BLOCK_PX = 28.0      # room for the axes marker above that column
DIFF_LABEL_GAP = 0.09         # room above a shrunk difference curve for its label
BAND_ABOVE = 0.15             # shaded zoom window reaches this far above the data
INSET_W, INSET_H = 0.27, 0.27  # zoom inset (axes fraction, plot area only)
PAD_PX = 2.5                  # clearance around every placed item (pixels at 100 dpi)


# ---------------------------------------------------------------------------
# Main fits
# ---------------------------------------------------------------------------

def ladder_stats(row):
    """(wR, {key: 1.0}) for a single-phase fit of the GSAS-II order-parameter
    ladder (workflow/scripts/phaseid_order_ladder.py), from results/phaseid/order_ladder.csv."""
    df = pd.read_csv(PHASEID / 'order_ladder.csv')
    s, eta = int(row['fit_sample']), float(row['source_key'])
    hit = df[(df['sample'] == s) & (np.isclose(df['eta'], eta)) & (df['failed'] != True)]  # noqa: E712
    if hit.empty:
        raise KeyError(f'{s:04d}: no converged ladder row at eta {eta} in order_ladder.csv')
    return float(hit.iloc[0]['wR']), {str(row['target_key']): 1.0}


def load_fit_folder(fit_path):
    """A refinement folder: fit.txt, ticks.json, summary.json.  The weight
    fractions in summary.json follow the phase order of ticks.json (both are
    written from the same GSAS-II phase list)."""
    folder = fit_path.parent
    x, yo, yc, bg, diff = np.loadtxt(fit_path, unpack=True)
    ticks = {k: np.asarray(v, float)
             for k, v in json.loads((folder / 'ticks.json').read_text()).items()}
    record = json.loads((folder / 'summary.json').read_text())['record']
    wt = record.get('state', {}).get('wt') or [1.0] * len(ticks)
    if len(wt) != len(ticks):
        raise ValueError(f'{folder}: {len(wt)} weight fractions for {len(ticks)} phases')
    total = float(sum(wt))
    fractions = {k: float(w) / total for k, w in zip(ticks, wt)}
    return x, yo, yc, bg, diff, ticks, float(record['wR']), fractions


def load_scan(name, fit_sample):
    """(two_theta, counts) of the measured scan ``name`` of button ``fit_sample``."""
    hits = sorted(SAMPLES_DIR.glob(f'*/buttons/{fit_sample:04d}/{name}')) or sorted(SAMPLES_DIR.glob(f'*/buttons/*/{name}'))
    if not hits:
        raise FileNotFoundError(f'scan {name} not found under {SAMPLES_DIR}')
    d = parse_xrd_file(hits[0])
    return np.asarray(d['two_theta'], float), np.asarray(d['intensity'], float)


def structure_path(cfg, row):
    name = cfg.get('structure') or f'{int(row["fit_sample"]):04d}.png'
    return STRUCTURE_DIR / name


def build_record(cfg, row, sample):
    num = cfg['sample']
    fit_spec = str(row['fit'])
    rec = dict(sample_number=num, fit_sample=int(row['fit_sample']), name=cfg['name'],
               category=int(sample['category']), reported=str(sample['reported']),
               state=str(row['state']), inset=cfg.get('inset'), predicted=None,
               structure=structure_path(cfg, row), split_list=cfg.get('split_list', False),
               inset_size=cfg.get('inset_size'), show_wt=cfg.get('show_wt', True),
               unindexed=cfg.get('unindexed') or [], diff_scale=cfg.get('diff_scale', 1.0),
               unindexed_label=cfg.get('unindexed_label', 'unindexed phase'),
               inset_fit=cfg.get('inset_fit', False), inset_wrap=cfg.get('inset_wrap', False),
               inset_lines=cfg.get('inset_lines') or [])
    phases = parse_phases(str(row['phases']))
    fit_path = PHASEID / fit_spec
    target = str(row['target_key']) if pd.notna(row['target_key']) else None

    if row['source_table'] == 'order_ladder':
        # The ladder exports carry no tick file; the Bragg positions come
        # from the refined cell's CIF beside the export.
        x, yo, yc, bg, diff = np.loadtxt(fit_path, unpack=True)
        wR, fractions = ladder_stats(row)
        angles, _ = peak_list(fit_path.with_suffix('.cif'), two_theta_range=TWO_THETA)
        ticks = {target: np.asarray(angles, float)}
    elif row['source_table'] == 'fit_folder':
        x, yo, yc, bg, diff, ticks, wR, fractions = load_fit_folder(fit_path)
    else:
        x, yo, yc, bg, diff, ticks = load_fit(fit_path)
        wR, fractions = lookup_fit_stats(dict(row, _fit_path=fit_path))

    keys = [k for k, _ in phases]
    if set(keys) != set(fractions):
        raise ValueError(f'{num:04d}: phases in {FITS_CSV.name} {keys} differ from the fit '
                         f'{sorted(fractions)}')
    if set(keys) != set(ticks):
        raise ValueError(f'{num:04d}: phases {keys} differ from the tick file {sorted(ticks)}')
    if target is not None and target not in keys:
        raise ValueError(f'{num:04d}: target_key {target!r} is not among the phases')

    # Scale on the displayed part of the refined range so the strongest
    # visible line sets the zoom.
    lo, hi = TWO_THETA
    m = (x >= lo) & (x <= hi)
    x, yo, yc, bg, diff = x[m], yo[m], yc[m], bg[m], diff[m]
    net_obs, net_calc = yo - bg, yc - bg
    scale = max(float(net_obs.max()), float(net_calc.max()))

    # Measured points of the scan outside the refined range (the fit is not
    # extended there): background held at its value at the nearer end of the
    # refined range.
    sx, sy = load_scan(row['scan'], rec['fit_sample'])
    below = (sx >= max(lo, cfg.get('data_from', lo))) & (sx < x[0] - 1e-6)
    above = (sx > x[-1] + 1e-6) & (sx <= hi)
    x_obs = np.concatenate([sx[below], x, sx[above]])
    obs = np.concatenate([sy[below] - bg[0], net_obs, sy[above] - bg[-1]])
    rec.update(x=x, x_obs=x_obs, obs=obs / scale, obs_fit=net_obs / scale,
               calc=net_calc / scale, diff=diff / scale * rec['diff_scale'],
               ticks=ticks, phases=phases, fractions=fractions, target=target, wR=wR,
               fit_range=(float(x[0]), float(x[-1])))

    msg = (f"  {num:04d} ({rec['fit_sample']:04d}) cat {rec['category']} {row['scan']:<34} "
           f"fit {x[0]:.1f}-{x[-1]:.1f}, data from {x_obs[0]:.1f}  R_wp {wR:5.2f}  "
           + ' + '.join(f'{k[:26]} {100 * fractions[k]:.1f}' for k, _ in phases))
    if cfg.get('predicted'):
        # What was predicted: the most stable OQMD structure at the composition.
        cif = find_oqmd_cif(sample['formula'], sample['oqmd_entry_id'])
        if cif is None:
            raise FileNotFoundError(f'{num:04d}: no OQMD CIF for {sample["formula"]}')
        spg = space_group_symbol(cif)
        rec['predicted'] = simulate_pattern(cif, two_theta_range=TWO_THETA, fwhm=FWHM)
        rec['predicted_label'] = ('OQMD prediction'
                                  + (f' ({_mathtext_space_group(spg)})' if spg else ''))
        msg += f'  | predicted {Path(cif).name} {spg}'
    if not rec['structure'].exists():
        msg += f'  | no structure picture {rec["structure"]}'
    print(msg)
    return rec


# ---------------------------------------------------------------------------
# Geometry helpers (display coordinates)
# ---------------------------------------------------------------------------

def _segment_boxes(ax, x, y, pad):
    """Bounding boxes [x0, x1, y0, y1] of the polyline segments, in pixels."""
    xy = ax.transData.transform(np.column_stack([x, y]))
    x0 = np.minimum(xy[:-1, 0], xy[1:, 0]) - pad
    x1 = np.maximum(xy[:-1, 0], xy[1:, 0]) + pad
    y0 = np.minimum(xy[:-1, 1], xy[1:, 1]) - pad
    y1 = np.maximum(xy[:-1, 1], xy[1:, 1]) + pad
    return np.column_stack([x0, x1, y0, y1])


def _hits(boxes, rect):
    """True if rect [x0, x1, y0, y1] overlaps any box."""
    if len(boxes) == 0:
        return False
    return bool(np.any((boxes[:, 0] < rect[1]) & (boxes[:, 1] > rect[0])
                       & (boxes[:, 2] < rect[3]) & (boxes[:, 3] > rect[2])))


class Panel:
    """One panel: draws the bands for a given headroom, then places items."""

    def __init__(self, ax, rec, letter):
        self.ax, self.rec, self.letter = ax, rec, letter
        self.headroom = HEADROOM0
        self.struct_k = STRUCT_K_MIN
        self.key_side = 'below'

    # -- bands ---------------------------------------------------------------
    def layout(self):
        rec = self.rec
        n = len(rec['phases'])
        y = GAP
        self.y_pred = None
        sim = SIM_FRAC if rec['predicted'] is not None else 0.0
        if sim:
            self.y_pred = y
            y += GAP
        # Everything above the ticks, in data units.
        dmin, dmax = float(rec['diff'].min()), float(rec['diff'].max())
        base = min(0.0, float(rec['obs'].min()), float(rec['calc'].min()))
        # A shrunk difference curve carries a label just above it.
        self.label_gap = DIFF_LABEL_GAP if rec['diff_scale'] != 1.0 else 0.0
        above = GAP + (dmax - dmin) + GAP + self.label_gap + (-base) + 1.0 + self.headroom
        frac = n * (TICK_FRAC + TICK_GAP_FRAC) + sim
        total = (y + above) / (1.0 - frac)
        self.th, self.tg = TICK_FRAC * total, TICK_GAP_FRAC * total
        self.sim_h = sim * total
        y += self.sim_h
        y_ticks_top = y + n * (self.th + self.tg) - self.tg
        self.y_res = y_ticks_top + GAP - dmin
        self.y_res_top = self.y_res + dmax
        self.y_fit = self.y_res_top + GAP + self.label_gap - base
        self.y_top = self.y_fit + 1.0 + self.headroom
        self.y_ticks_top = y_ticks_top

    def draw(self):
        ax, rec = self.ax, self.rec
        ax.cla()
        self.layout()
        lo, hi = TWO_THETA
        x = rec['x']
        self.curves = []                          # (x, y, pad_px) obstacles

        if rec['predicted'] is not None:
            pt, pi = rec['predicted']
            yy = pi / pi.max() * self.sim_h + self.y_pred
            ax.plot(pt, yy, color=COLOR_PRED, lw=LW_SIM)
            self.curves.append((pt, yy, 1.0))

        # Bragg positions, one short row per phase, the main phase on top
        # (next to the data), the target phase in green.
        phases = sorted(rec['phases'], key=lambda kd: -rec['fractions'][kd[0]])
        self.phase_order = phases
        self.tick_rows = []
        for i, (key, _) in enumerate(phases):
            y1 = self.y_ticks_top - i * (self.th + self.tg)
            y0 = y1 - self.th
            tt = rec['ticks'][key]
            tt = tt[(tt >= lo) & (tt <= hi)]
            colour = COLOR_SOLUTION if key == rec['target'] else COLOR_TICK_OTHER
            ax.vlines(tt, y0, y1, color=colour, lw=LW_TICK)
            ax.text(-0.006, (y0 + y1) / 2, str(i + 1), transform=ax.get_yaxis_transform(),
                    fontsize=SMALL - 1, color=colour, ha='right', va='center', clip_on=False)
            self.tick_rows.append((y0, y1))

        # Difference curve, then data and fit.
        ax.axhline(self.y_res, color='0.85', lw=0.3, zorder=0)
        ax.plot(x, rec['diff'] + self.y_res, color=COLOR_DIFF, lw=LW_DIFF)
        ax.plot(rec['x_obs'], rec['obs'] + self.y_fit, color='black', lw=0, marker='x', ms=MS_DATA,
                mew=MEW_DATA, alpha=0.9, zorder=2)
        ax.plot(x, rec['calc'] + self.y_fit, color=COLOR_FIT, lw=LW_FIT,
                solid_joinstyle='round', zorder=3)
        # Unindexed lines: an open down-triangle just above the measured
        # pattern (its highest point under the marker's width).
        self.marker_pts = []
        if rec['unindexed']:
            bbw = ax.get_window_extent().width
            half = 0.5 * (TWO_THETA[1] - TWO_THETA[0]) * (MS_UNINDEXED * 100 / 72) / bbw
            for a in rec['unindexed']:
                m = np.abs(rec['x_obs'] - a) <= max(half, 0.1)
                self.marker_pts.append((a, float(rec['obs'][m].max()) + self.y_fit))
            off = matplotlib.transforms.offset_copy(ax.transData, fig=ax.figure,
                                                    y=MS_UNINDEXED / 2 + 1.5, units='points')
            ax.plot([a for a, _ in self.marker_pts], [y for _, y in self.marker_pts], lw=0,
                    marker='v', ms=MS_UNINDEXED, mfc='none', mec=COLOR_UNINDEXED, mew=0.7,
                    transform=off, zorder=4, clip_on=False)
        self.curves += [(x, rec['diff'] + self.y_res, 1.0),
                        (rec['x_obs'], rec['obs'] + self.y_fit, MS_DATA / 2 * 100 / 72 + 0.5),
                        (x, rec['calc'] + self.y_fit, 1.0)]

        if rec['inset']:
            # Shade the zoomed window up to just above the data there (not the
            # full height, which would block the structure picture's corner).
            ilo, ihi = rec['inset']
            mw = (rec['x_obs'] >= ilo) & (rec['x_obs'] <= ihi)
            self.band_top = self.y_fit + float(rec['obs'][mw].max()) + BAND_ABOVE
            ax.fill_between([ilo, ihi], 0, self.band_top, color='0.5', alpha=0.10, zorder=0,
                            lw=0)

        ax.set_xlim(lo, hi)
        ax.set_ylim(0.0, self.y_top)
        ax.set_yticks([])
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}' if lo < v < hi else ''))
        ax.tick_params(axis='x', which='major', direction='in', top=True, length=3.5,
                       labelbottom=False, labelsize=FONTSIZE, pad=2)
        ax.tick_params(axis='x', which='minor', direction='in', top=True, length=2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.6)

        # Titles: panel letter and formula left, outcome category right.
        mark = REPORTED_MARK[rec['reported']]
        ax.set_title(f'({self.letter}) {rec["name"]}{mark}', loc='left', fontsize=FONTSIZE + 0.5,
                     fontweight='bold', pad=3)
        ax.set_title(CATEGORY_NAME[rec['category']], loc='right', fontsize=SMALL,
                     color=COLOR_SOLUTION, pad=3.5)

    # -- structure key -----------------------------------------------------------
    def _key_width(self, key, renderer):
        """Width in pixels of the key row: axes marker, then one sphere and
        label per element."""
        w = 2 * KEY_AXIS_PX + 4
        for item in key:
            t = self.ax.text(0, 0, item['label'], fontsize=KEY_FONT)
            w += KEY_DOT_PX + 2 + t.get_window_extent(renderer).width + KEY_GAP_PX
            t.remove()
        return w

    def _key_column_size(self, key, renderer):
        """Width and height in pixels of the key set beside the picture: axes
        marker on top, then one sphere and label per line."""
        w = 2 * KEY_AXIS_PX + 16
        for item in key:
            t = self.ax.text(0, 0, item['label'], fontsize=KEY_FONT)
            w = max(w, KEY_DOT_PX + 3 + t.get_window_extent(renderer).width)
            t.remove()
        return w, KEY_AXIS_BLOCK_PX + KEY_LINE_PX * len(key)

    def _draw_axes_marker(self, side, org):
        ax = self.ax
        inv = ax.transAxes.inverted()
        for name, a in side['axes'].items():
            u = np.array([a['dx'], a['dy']])
            tip = org + u * KEY_AXIS_PX * max(a['frac'], 0.4)
            ax.annotate('', xy=inv.transform(tip), xytext=inv.transform(org),
                        xycoords='axes fraction', textcoords='axes fraction',
                        arrowprops=dict(arrowstyle='-|>,head_length=0.25,head_width=0.12',
                                        lw=0.5, color='0.2', shrinkA=0, shrinkB=0), zorder=6,
                        gid='key')
            ax.text(*inv.transform(tip + u * 4.0), f'${name}$', transform=ax.transAxes,
                    fontsize=KEY_FONT - 0.5, ha='center', va='center', zorder=6, gid='key')

    def _key_dot(self, item, x, y):
        ax = self.ax
        inv = ax.transAxes.inverted()
        c = item['colour']
        kw = (dict(markerfacecolor=c[0], markerfacecoloralt=c[1], fillstyle='left')
              if isinstance(c, list) else dict(markerfacecolor=c))
        ax.plot(*inv.transform((x, y)), transform=ax.transAxes, marker='o',
                ms=KEY_DOT_PX * 72 / 100, mec='0.25', mew=0.4, lw=0, zorder=6,
                clip_on=False, **kw)

    def _draw_key_column(self, side, x0, y1):
        """Key beside the picture: axes marker at the top, elements below,
        left-aligned at x0 (pixels; y1 = top)."""
        ax = self.ax
        inv = ax.transAxes.inverted()
        self._draw_axes_marker(side, np.array([x0 + KEY_AXIS_PX + 8, y1 - KEY_AXIS_BLOCK_PX + 10]))
        y = y1 - KEY_AXIS_BLOCK_PX - KEY_LINE_PX / 2
        for item in side['key']:
            self._key_dot(item, x0 + KEY_DOT_PX / 2, y)
            ax.text(*inv.transform((x0 + KEY_DOT_PX + 3, y)), item['label'], transform=ax.transAxes,
                    fontsize=KEY_FONT, ha='left', va='center', zorder=6, gid='key')
            y -= KEY_LINE_PX

    def _draw_key(self, side, x0, x1, yc):
        """Axes marker at the left of the row, element key right-aligned
        (pixel coordinates, row centred on yc)."""
        ax = self.ax
        inv = ax.transAxes.inverted()
        self._draw_axes_marker(side, np.array([x0 + KEY_AXIS_PX + 2, yc - 4]))
        # element key, laid out right to left
        x = x1
        for item in reversed(side['key']):
            t = ax.text(*inv.transform((x, yc)), item['label'], transform=ax.transAxes,
                        fontsize=KEY_FONT, ha='right', va='center', zorder=6, gid='key')
            x -= t.get_window_extent(ax.figure.canvas.get_renderer()).width + 2 + KEY_DOT_PX / 2
            self._key_dot(item, x, yc)
            x -= KEY_DOT_PX / 2 + KEY_GAP_PX

    # -- obstacles -------------------------------------------------------------
    def obstacles(self):
        ax = self.ax
        boxes = [_segment_boxes(ax, cx, cy, pad) for cx, cy, pad in self.curves]
        # Tick rows (full width) and, for the zoom panel, the shaded window.
        bb = ax.get_window_extent()
        rows = []
        for y0, y1 in self.tick_rows:
            p0 = ax.transData.transform([(TWO_THETA[0], y0), (TWO_THETA[1], y1)])
            rows.append([bb.x0, bb.x1, p0[0, 1], p0[1, 1]])
        if self.rec['inset']:
            ilo, ihi = self.rec['inset']
            p = ax.transData.transform([(ilo, 0), (ihi, self.band_top)])
            rows.append([p[0, 0], p[1, 0], p[0, 1], p[1, 1]])
        r = MS_UNINDEXED / 2 * 100 / 72
        for a, y in getattr(self, 'marker_pts', []):
            px, py = ax.transData.transform((a, y))
            py += (MS_UNINDEXED / 2 + 1.5) * 100 / 72
            rows.append([px - r - 1, px + r + 1, py - r - 1, py + r + 1])
        if rows:
            boxes.append(np.asarray(rows, float))
        return np.vstack(boxes)

    # -- items -----------------------------------------------------------------
    def items(self, renderer):
        """Items to place: (name, width_px, height_px, region, cost(x0, y1), draw)."""
        rec, ax = self.rec, self.ax
        out = []

        # Phase list with fractions and R_wp.
        lines = []
        for i, (key, display) in enumerate(self.phase_order):
            # Name and fraction on one line, the structure description under
            # it, so the block is narrow enough to sit beside the picture.
            wt = f'  {100 * rec["fractions"][key]:.0f} wt%' if rec['show_wt'] else ''
            colour = COLOR_SOLUTION if key == rec['target'] else 'black'
            name, _, desc = display.partition(', ')
            lines.append((f'{i + 1}  {name}{wt}', colour))
            if desc:
                lines.append((f'     {desc}', colour))
        extra = []
        if rec['unindexed']:
            extra.append((f'\u25bd  {rec["unindexed_label"]}', COLOR_UNINDEXED))
        extra.append((rf'$R_{{wp}}$ = {rec["wR"]:.1f} %', 'black'))
        # ``split_list``: the extra lines form a second block placed on its own.
        blocks = ([('phases', lines), ('phases_extra', extra)] if rec['split_list']
                  else [('phases', lines + extra)])
        for bname, blines in blocks:
            sizes = []
            for text, _ in blines:
                t = ax.text(0, 0, text, fontsize=SMALL, transform=ax.transAxes)
                e = t.get_window_extent(renderer)
                sizes.append((e.width, e.height))
                t.remove()
            lh = max(h for _, h in sizes) * 1.12
            w = max(wd for wd, _ in sizes)
            h = lh * len(blines)

            def draw_text(x0, y1, lines=blines, lh=lh, w=w):
                # x0, y1: left and top of the block in pixels.
                right = x0 + w > ax.get_window_extent().x0 + 0.5 * ax.get_window_extent().width
                inv = ax.transAxes.inverted()
                for j, (text, colour) in enumerate(lines):
                    px = x0 + w if right else x0
                    ax.text(*inv.transform((px, y1 - j * lh)), text, transform=ax.transAxes,
                            fontsize=SMALL, color=colour, ha='right' if right else 'left',
                            va='top', zorder=5)
            out.append((bname, w, h, 'fit', lambda fx, fy: (1 - fy) + 0.6 * (1 - fx), draw_text))

        # Structure picture, with its element key and axes in a row under it
        # (drawn here, so they are the same size in every panel).
        bb = ax.get_window_extent()
        if rec['structure'].exists():
            img = plt.imread(rec['structure'])
            ih, iw = img.shape[:2]
            side = json.loads(rec['structure'].with_suffix('.json').read_text())
            sh = self.struct_k * side['height_A']
            iw_px = sh * iw / ih
            if self.key_side == 'below':
                sw = max(iw_px, self._key_width(side['key'], renderer))
                ph = sh + KEY_ROW_PX
            else:
                kw, kh = self._key_column_size(side['key'], renderer)
                sw = iw_px + 4 + kw
                ph = max(sh, kh)

            def draw_img(x0, y1, img=img, sw=sw, sh=sh, side=side, iw_px=iw_px):
                inv = ax.transAxes.inverted()
                below = self.key_side == 'below'
                xi = x0 + (sw - iw_px) / 2 if below else x0
                (ax0, ay0), (ax1, ay1) = inv.transform([(xi, y1 - sh), (xi + iw_px, y1)])
                iax = ax.inset_axes([ax0, ay0, ax1 - ax0, ay1 - ay0])
                iax.imshow(img, interpolation='lanczos')
                iax.set_axis_off()
                self.image_ax = iax
                if below:
                    self._draw_key(side, x0, x0 + sw, y1 - sh - KEY_ROW_PX / 2)
                else:
                    self._draw_key_column(side, x0 + iw_px + 4, y1)
            out.append(('structure', sw, ph,
                        ('fixed', STRUCT_ANCHOR[0] - sw / bb.width, STRUCT_ANCHOR[1]),
                        lambda fx, fy: 0.0, draw_img))

        # Zoom on the missing predicted line (with room for its title above
        # and its tick labels below).
        if rec['inset']:
            iw_frac, ih_frac = rec.get('inset_size') or (INSET_W, INSET_H)
            iw_px, ih_px = iw_frac * bb.width, ih_frac * bb.height
            title_px, label_px, side_px = 12.0, 11.0, 6.0
            ilo, ihi = rec['inset']
            t = ax.text(0, 0, self._zoom_title(), fontsize=SMALL - 1, transform=ax.transAxes)
            tw = t.get_window_extent(renderer).width
            title_px += t.get_window_extent(renderer).height - 9.0   # a wrapped title is taller
            t.remove()
            side_px = max(side_px, (tw - iw_px) / 2 + 2)

            def draw_inset(x0, y1, iw_px=iw_px, ih_px=ih_px, title_px=title_px,
                           side_px=side_px):
                inv = ax.transAxes.inverted()
                (ax0, ay0), (ax1, ay1) = inv.transform(
                    [(x0 + side_px, y1 - title_px - ih_px), (x0 + side_px + iw_px, y1 - title_px)])
                self._draw_zoom(ax.inset_axes([ax0, ay0, ax1 - ax0, ay1 - ay0]))
            out.append(('zoom', iw_px + 2 * side_px, title_px + ih_px + label_px, 'fit',
                        lambda fx, fy: (1 - fy) + 0.8 * (1 - fx), draw_inset))

        # Label of the simulated prediction, inside its own band.
        if rec['predicted'] is not None:
            t = ax.text(0, 0, rec['predicted_label'], fontsize=SMALL, transform=ax.transAxes)
            e = t.get_window_extent(renderer)
            t.remove()

            def draw_pl(x0, y1, w=e.width):
                inv = ax.transAxes.inverted()
                ax.text(*inv.transform((x0 + w, y1)), rec['predicted_label'],
                        transform=ax.transAxes, fontsize=SMALL, color=COLOR_PRED,
                        ha='right', va='top', zorder=5)
            out.append(('pred_label', e.width, e.height, 'pred',
                        lambda fx, fy: (1 - fy) + 0.5 * (1 - fx), draw_pl))

        # Label of a shrunk difference curve, beside it.
        if rec['diff_scale'] != 1.0:
            label = rf'difference $\times${rec["diff_scale"]:g}'
            t = ax.text(0, 0, label, fontsize=SMALL, transform=ax.transAxes)
            e = t.get_window_extent(renderer)
            t.remove()

            def draw_dl(x0, y1, label=label):
                inv = ax.transAxes.inverted()
                ax.text(*inv.transform((x0, y1)), label, transform=ax.transAxes,
                        fontsize=SMALL, color=COLOR_DIFF, ha='left', va='top', zorder=5)
            out.append(('diff_label', e.width, e.height, 'res',
                        lambda fx, fy: fy + 2.0 * fx, draw_dl))
        return out

    def _zoom_title(self):
        ilo, ihi = self.rec['inset']
        sep = '\n' if self.rec['inset_wrap'] else ' '      # two lines where the slot is narrow
        return f'{ilo:g}–{ihi:g}°,{sep}same scale'

    def _draw_zoom(self, axins):
        rec = self.rec
        ilo, ihi = rec['inset']
        pt, pi = rec['predicted']
        x = rec['x_obs']
        # Same scale as the panel: the simulated pattern normalised as in its
        # band relative to the measured pattern's strongest line.
        m = (pt >= ilo) & (pt <= ihi)
        axins.plot(pt[m], pi[m] / pi.max(), color=COLOR_PRED, lw=LW_SIM)
        mx = (x >= ilo) & (x <= ihi)
        axins.plot(x[mx], rec['obs'][mx], color='black', lw=0, marker='x', ms=MS_DATA,
                   mew=MEW_DATA)
        if rec['inset_fit']:
            mf = (rec['x'] >= ilo) & (rec['x'] <= ihi)
            axins.plot(rec['x'][mf], rec['calc'][mf], color=COLOR_FIT, lw=0.8, zorder=3)
        top = max(float((pi / pi.max())[m].max()), float(rec['obs'][mx].max()))
        # Predicted lines named in the zoom (label just above the simulated peak).
        for t, lab in rec['inset_lines']:
            mm = np.abs(pt - t) < 0.4
            axins.annotate(lab, (t, float((pi / pi.max())[mm].max()) if mm.any() else 0.0),
                           xytext=(0, 2), textcoords='offset points', ha='center', va='bottom',
                           fontsize=SMALL - 2, color=COLOR_PRED, linespacing=0.9)
        axins.set_xlim(ilo, ihi)
        axins.set_ylim(-0.03 * top, (1.8 if rec['inset_lines'] else 1.25) * top)
        axins.set_yticks([])
        axins.xaxis.set_major_locator(MultipleLocator(2))
        axins.tick_params(axis='x', labelsize=SMALL - 1, pad=1, length=2)
        axins.set_title(self._zoom_title(), fontsize=SMALL - 1, pad=2)
        for spine in axins.spines.values():
            spine.set_linewidth(0.5)
        self.zoom_ax = axins

    # -- placement -------------------------------------------------------------
    def place(self, renderer, order=None):
        """Place every item (in ``order``); the plan, or None if one does not fit."""
        ax = self.ax
        bb = ax.get_window_extent()
        obst = self.obstacles()
        y_fit_px = ax.transData.transform((0, self.y_res_top))[1]
        regions = {'fit': (y_fit_px, bb.y1),
                   # beside the difference curve, in the gap kept above it
                   'res': tuple(ax.transData.transform(
                       [(0, self.y_ticks_top),
                        (0, self.y_res_top + GAP / 2 + self.label_gap)])[:, 1])}
        if self.y_pred is not None:
            p = ax.transData.transform([(0, self.y_pred), (0, self.y_pred + self.sim_h + GAP)])
            regions['pred'] = (p[0, 1], p[1, 1])
        # Occupancy bitmap of the plot area (1 px cells) with a summed-area
        # table, so an empty-rectangle test is four lookups.
        W, H = int(np.ceil(bb.width)) + 1, int(np.ceil(bb.height)) + 1
        grid = np.zeros((H + 1, W + 1), np.int32)

        def mark(boxes):
            b = np.asarray(boxes, float)
            x0 = np.clip(np.floor(b[:, 0] - bb.x0), 0, W).astype(int)
            x1 = np.clip(np.ceil(b[:, 1] - bb.x0), 0, W).astype(int)
            y0 = np.clip(np.floor(b[:, 2] - bb.y0), 0, H).astype(int)
            y1 = np.clip(np.ceil(b[:, 3] - bb.y0), 0, H).astype(int)
            ok = (x1 > x0) & (y1 > y0)
            x0, x1, y0, y1 = x0[ok], x1[ok], y0[ok], y1[ok]
            np.add.at(grid, (y0, x0), 1)
            np.add.at(grid, (y0, x1), -1)
            np.add.at(grid, (y1, x0), -1)
            np.add.at(grid, (y1, x1), 1)

        mark(obst)
        plan = []
        items = self.items(renderer)
        if order is not None:
            items = [items[i] for i in order]
        for name, w, h, region, cost, draw in items:
            occ = (np.cumsum(np.cumsum(grid, 0), 1)[:H, :W] > 0).astype(np.int32)
            sat = np.zeros((H + 1, W + 1), np.int64)
            sat[1:, 1:] = occ.cumsum(0).cumsum(1)
            rw = int(np.ceil(w + 2 * PAD_PX))
            rh = int(np.ceil(h + 2 * PAD_PX))
            # Candidate lower-left corners of the padded rectangle (pixels
            # from the plot-area origin), restricted to the region; a fixed
            # item ('fixed', left, top in axes fraction) has one candidate.
            if isinstance(region, tuple):
                _, fl, ft = region
                ys = np.asarray([int(round(ft * bb.height - h - PAD_PX))])
                xs = np.asarray([int(round(fl * bb.width - PAD_PX))])
            else:
                ylo, yhi = regions[region]
                ys = np.arange(max(0, int(ylo - bb.y0)), int(yhi - bb.y0) - rh + 1)
                xs = np.arange(0, W - rw)
            if len(ys) == 0 or len(xs) == 0:
                self.failed = name
                return None
            Y, X = np.meshgrid(ys, xs, indexing='ij')
            count = (sat[Y + rh, X + rw] - sat[Y, X + rw] - sat[Y + rh, X] + sat[Y, X])
            free = count == 0
            if not free.any():
                self.failed = name
                return None
            fx = (X + rw - PAD_PX) / bb.width          # right edge of the item
            fy = (Y + rh - PAD_PX) / bb.height         # top edge of the item
            c = np.where(free, cost(fx, fy), np.inf)
            k = np.unravel_index(np.argmin(c), c.shape)
            px0, py0 = bb.x0 + X[k], bb.y0 + Y[k]
            rect = (px0, px0 + rw, py0, py0 + rh)
            mark([rect])
            plan.append((name, px0 + PAD_PX, py0 + rh - PAD_PX, draw, rect, float(c[k])))
        return plan

    def fit_items(self, renderer, headroom_max):
        """Smallest headroom (from HEADROOM0) at which every item fits; the
        plan, or None.  Leaves the panel drawn at that headroom."""
        self.headroom = HEADROOM0
        while self.headroom <= headroom_max + 1e-9:
            self.draw()
            items = self.items(renderer)
            fixed = [i for i, it in enumerate(items) if isinstance(it[3], tuple)]
            free = [i for i in range(len(items)) if i not in fixed]
            plans = [p for p in (self.place(renderer, fixed + list(order))
                                 for order in itertools.permutations(free))
                     if p is not None]
            if plans:
                return min(plans, key=lambda p: sum(e[5] for e in p))
            self.headroom += HEADROOM_STEP
        return None

    def build(self, renderer, k_cap=STRUCT_K_MAX, draw=True):
        """Largest picture scale up to k_cap (with the key under or beside the
        picture) that costs at most STRUCT_HEADROOM_ALLOW of headroom over the
        smallest picture; returns that scale.  With draw=False nothing is
        placed, so the figure can first agree on one scale for all panels."""
        self.struct_k, self.key_side = STRUCT_K_MIN, 'below'
        plan = self.fit_items(renderer, HEADROOM_MAX)
        if plan is None:
            raise RuntimeError(f'{self.rec["sample_number"]:04d}: {self.failed} does not fit')
        limit = min(self.headroom + STRUCT_HEADROOM_ALLOW, HEADROOM_MAX)
        found = None
        for k in np.arange(k_cap, STRUCT_K_MIN - 1e-9, -STRUCT_K_STEP):
            for side in ('below', 'right'):
                self.struct_k, self.key_side = float(k), side
                bigger = self.fit_items(renderer, limit)
                if bigger is not None:
                    found = (float(k), side, bigger)
                    break
            if found:
                break
        if found is None:
            self.struct_k, self.key_side = STRUCT_K_MIN, 'below'
            plan = self.fit_items(renderer, HEADROOM_MAX)
        else:
            self.struct_k, self.key_side, plan = found
            self.fit_items(renderer, limit)            # leave the panel drawn at this layout
        if draw:
            for _, x0, y1, draw_item, _, _ in plan:
                draw_item(x0, y1)
            self.plan = plan
        return self.struct_k

    def verify(self, renderer):
        """Independent check after drawing: every text artist and inset inside
        the plot area clears every curve, tick row and each other."""
        ax = self.ax
        obst = self.obstacles()
        rects = []
        key_parts = []                          # key row: one group, laid out by _draw_key
        for t in ax.texts:
            if t.get_transform() == ax.get_yaxis_transform():
                continue                        # tick-row numbers sit outside the axes
            e = t.get_window_extent(renderer)
            if t.get_gid() == 'key':
                key_parts.append(e)
                continue
            rects.append((t.get_text()[:30], (e.x0, e.x1, e.y0, e.y1)))
        if key_parts:
            rects.append(('structure key', (min(e.x0 for e in key_parts), max(e.x1 for e in key_parts),
                                            min(e.y0 for e in key_parts), max(e.y1 for e in key_parts))))
        for ch in ax.child_axes:
            e = ch.get_tightbbox(renderer)
            rects.append(('inset', (e.x0, e.x1, e.y0, e.y1)))
        bad = []
        for i, (name, r) in enumerate(rects):
            if _hits(obst, r):
                bad.append(f'{name!r} touches a curve')
            for name2, r2 in rects[i + 1:]:
                if _hits(np.asarray([r2]), r):
                    bad.append(f'{name!r} overlaps {name2!r}')
            bb = ax.get_window_extent()
            if r[0] < bb.x0 or r[1] > bb.x1 or r[2] < bb.y0 or r[3] > bb.y1:
                bad.append(f'{name!r} leaves the plot area')
        return bad


def legend_handles(with_prediction):
    h = [
        Line2D([0], [0], color='black', lw=0, marker='x', ms=5, mew=1.0,
               label='Measured'),
        Line2D([0], [0], color=COLOR_FIT, lw=1.6, label='Fit'),
        Line2D([0], [0], color=COLOR_DIFF, lw=1.3, label='Difference'),
        Line2D([0], [0], color=COLOR_SOLUTION, lw=0, marker='|', ms=7, mew=1.2,
               label='Reflections, target phase'),
        Line2D([0], [0], color=COLOR_TICK_OTHER, lw=0, marker='|', ms=7, mew=1.2,
               label='other phases'),
    ]
    if with_prediction:
        h.append(Line2D([0], [0], color=COLOR_PRED, lw=1.6, label='Prediction (simulated)'))
    return h


def make_figure(records):
    fig, axes = plt.subplots(NROWS, NCOLS, figsize=(FIG_W_CM * CM, FIG_H_CM * CM),
                             sharex=True, squeeze=False)
    # The legend row above the grid is only needed when no slot is spare.
    top_cm = TOP_CM if len(records) == NROWS * NCOLS else TOP_CM_NO_LEGEND
    fig.subplots_adjust(left=LEFT_CM / FIG_W_CM, right=1 - RIGHT_CM / FIG_W_CM,
                        top=1 - top_cm / FIG_H_CM, bottom=BOTTOM_CM / FIG_H_CM,
                        hspace=HSPACE, wspace=WSPACE)
    renderer = fig.canvas.get_renderer()
    panels = [Panel(axes[(i - 1) // NCOLS, (i - 1) % NCOLS], rec, panel_letter(i))
              for i, rec in enumerate(records, start=1)]
    # Nearly one picture scale for all panels: each panel takes its largest
    # scale, but no more than STRUCT_K_SPREAD times the smallest panel's.
    ks = [p.build(renderer, draw=False) for p in panels]
    for p, k in zip(panels, ks):
        p.build(renderer, k_cap=min(k, STRUCT_K_SPREAD * min(ks)))
    print('  structure pictures (px/angstrom at 100 dpi, key position): '
          + ', '.join(f'({p.letter}) {p.struct_k:g} {p.key_side}' for p in panels))

    spare = list(axes.flat[len(records):])
    for ax in spare:
        ax.axis('off')
    handles = legend_handles(any(r['predicted'] is not None for r in records))
    if spare:
        spare[0].legend(handles=handles, loc='center', fontsize=FONTSIZE, frameon=False,
                        handlelength=2.0, labelspacing=0.9, handletextpad=0.8)
    else:
        fig.legend(handles=handles, loc='upper center', ncol=len(handles), fontsize=SMALL + 0.5,
                   frameon=False, handlelength=1.5, handletextpad=0.4, columnspacing=1.1,
                   bbox_to_anchor=(0.5 + (LEFT_CM - RIGHT_CM) / 2 / FIG_W_CM, 1.0),
                   borderaxespad=0.15)

    for c in range(NCOLS):
        for rr in range(NROWS - 1, -1, -1):
            if axes[rr, c].axison:
                axes[rr, c].tick_params(labelbottom=True)
                axes[rr, c].set_xlabel(XLABEL, fontsize=FONTSIZE + 1, labelpad=1)
                break
    for rr in range(NROWS):
        axes[rr, 0].set_ylabel(YLABEL_NET, fontsize=FONTSIZE + 0.5, labelpad=7)

    fig.canvas.draw()
    problems = []
    for p in panels:
        problems += [f'({p.letter}) {b}' for b in p.verify(renderer)]
        print(f'  ({p.letter}) headroom {p.headroom:.2f}: '
              + ', '.join(f'{name} at ({(x0 - p.ax.get_window_extent().x0) / p.ax.get_window_extent().width:.2f}, '
                          f'{(y1 - p.ax.get_window_extent().y0) / p.ax.get_window_extent().height:.2f})'
                          for name, x0, y1, _, _, _ in p.plan))
    # Titles: left and right parts of each panel title must not meet.
    for p in panels:
        left, right = p.ax._left_title, p.ax._right_title
        el, er = left.get_window_extent(renderer), right.get_window_extent(renderer)
        if el.x1 + 1 > er.x0:
            problems.append(f'({p.letter}) title parts overlap')
    if problems:
        raise SystemExit('Layout check failed:\n  ' + '\n  '.join(problems))
    print('  layout check: no text or inset touches a curve, a tick row or another item')
    return fig


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT_PDF)
    args = parser.parse_args()

    samples = load_samples()
    fits = pd.read_csv(FITS_CSV).set_index('sample_number', drop=False)
    missing = [c['sample'] for c in PANELS if c['sample'] not in fits.index]
    if missing:
        raise SystemExit(f'{FITS_CSV} has no main fit for {missing}')
    if len(PANELS) > NROWS * NCOLS:
        raise SystemExit(f'{len(PANELS)} panels do not fit a {NROWS} x {NCOLS} grid')

    print(f'Window {TWO_THETA[0]:g}-{TWO_THETA[1]:g} deg, simulated FWHM {FWHM} deg\n')
    records = [build_record(cfg, fits.loc[cfg['sample']], samples.loc[cfg['sample']])
               for cfg in PANELS]

    fig = make_figure(records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Curves, data points and tick rows are rasterized (thousands of markers made
    # the vector PDF slow to display); text, axes and spines stay vector.
    for ax in fig.axes:
        for a in [*ax.lines, *ax.collections, *ax.images]:
            a.set_rasterized(True)
    fig.savefig(args.out, dpi=RASTER_DPI)
    # The manuscript includes a single PNG (loads fast in the PDF viewer).
    fig.savefig(args.out.with_suffix('.png'), dpi=PNG_DPI)
    plt.close(fig)
    print(f'\nSaved: {args.out}  ({len(records)} panels, {FIG_W_CM:.2f} x {FIG_H_CM:.2f} cm)')


if __name__ == '__main__':
    main()
