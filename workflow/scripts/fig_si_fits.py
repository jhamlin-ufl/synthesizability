#!/usr/bin/env python3
"""Supplementary figures: each target's main fit and the contrast fits behind its category.

  in:  tables/si_fits.csv (one row per fit; columns in tables/README.md),
       tables/paper_samples.csv, samples/<formula>/fits/<button>_<fit>/, CIFs in the `sticks` column
  out: results/supplement/fig_S_<row>_<formula>.pdf, one per Table 1 row with a deposited fit;
       results/supplement/fig_S_ascast_annealed.pdf (annealed main fit above the as-cast fit)

Panels: (a) main fit: observed minus background (crosses), fit (red), difference (blue), one
numbered row of Bragg positions per phase (green = target_phase) and weight fractions; (b) the
difference curves of every fit on one scale, labelled with model and R_wp (recomputed from fit.txt
with w = 1/obs and checked against summary.json); (c) zoom windows on the deciding lines: observed
counts and each model's calculated pattern, or the Cu K-alpha stick pattern of the CIF in `sticks`
(strongest line scaled to the strongest observed line) where the refined pattern is unusable.
"""
import argparse
import json
import re
import sys
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, 'workflow/src')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core.periodic_table import Element

from synthesizability.formula import ALLEN, mathtext
from synthesizability.samples import CATEGORY_NAME

FITS_CSV = Path('tables/si_fits.csv')
SAMPLES_CSV = Path('tables/paper_samples.csv')
FITS_DIR = Path('.')   # fit_dir in si_fits.csv is relative to the repository root
OUT_DIR = Path('results/supplement')

FIGSIZE = (7.0, 4.6)            # inches; two per page
FONTSIZE, SMALL = 8.5, 7.0
ROLES = ['FoR', 'C1', 'C2', 'C3']
# Main fit red, first contrast Okabe-Ito sky blue, second contrast black and dashed in the zoom panels:
# mutually distinct also for colour-blind readers (red vs orange was not), and clear of the green target ticks.
ROLE_COLOR = {'FoR': '#D32F2F', 'C1': '#1E88E5', 'C2': 'black', 'C3': '#8E24AA', 'ascast': '#1565C0'}
# zoom panels: main fit widest and lowest; each alternative model thinner and on top, so both stay visible where they overlap
ZOOM_LW = {'FoR': 2.0, 'C1': 1.5, 'C2': 1.5, 'C3': 1.5}
ZOOM_Z = {'FoR': 3, 'C1': 4, 'C2': 5, 'C3': 5}
# zoom panels: alternative models dashed (blue) and dotted (black) so the red main fit shows through the gaps
ROLE_DASH = {'C1': (0, (5, 2.5)), 'C2': (0, (1.2, 1.6)), 'C3': (0, (1.2, 1.6))}
COLOR_DIFF, COLOR_TARGET, COLOR_OTHER = '#0072B2', '#2E7D32', '0.45'
TICK_FRAC, TICK_GAP_FRAC = 0.035, 0.012   # tick rows in (a), fraction of its height
RWP_TOL = 0.05                  # recomputed vs recorded R_wp (percent)
MAX_ZOOMS = 8
plt.rcParams.update({'font.size': FONTSIZE, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})

# Display names that the key would get wrong (database polymorphs filed under a prediction tag).
DISPLAY = {'YCoSi3_pred_OQMD-1799027': r'$\mathrm{YCoSi_{3}}$ Cmmm (SmNiGe$_3$-type)',
           'Y(CoSi)2_mp_mp-5129': r'$\mathrm{YCo_{2}Si_{2}}$',
           'LaOH3_cod_4031381': r'La(OH)$_3$',
           'ScCuPd2_pred_ALEX-agm078275705': r'$\mathrm{ScPd_{2}Cu}$ ordered P4/mmm (Alexandria, above hull)'}
# Phase keys -> display words (panel a). Keys are <formula>_<source>_<id>[_dis]; IDs are dropped.
WORDS = {'pred': 'predicted', 'dis': 'site-mixed', 'oqmd': 'OQMD', 'alex': 'Alexandria', 'mp': 'MP',
         'cod': 'COD', 'L12': r'L1$_2$', 'L21': r'L2$_1$', 'E93': r'E9$_3$',
         'C11b': r'C11$_b$', 'C2m': 'C2/m', 'A1': 'fcc', 'A2': 'bcc'}
DROP = re.compile(r'\d{3,}|agm\d+|mp-\d+|(OQMD|ALEX)-\w+|oqmd\d+\w*|o1|proto|Krachan|[A-Z][a-z]+\d{4}')
SPECIAL = [(re.compile(r'.+_L12_(\w+)_([A-Z][a-z])([A-Z][a-z])_proposed'),
            lambda m: rf'L1$_2$ {m[1]}({m[2]},{m[3]})$_3$ parent'),
           (re.compile(r'.+_pm_(\w+)-(\w+)_x(\d+)'),
            lambda m: f'{fmt_formula(m[1])}-{fmt_formula(m[2])} solid soln. (x = 0.{m[3]})'),
           (re.compile(r'(.+)_([A-Z][a-z]?)([A-Z][a-z]?)mix'), lambda m: f'{phase_label(m[1])}, {m[2]}/{m[3]} mixed')]


def load_fit(folder):
    """Curves, tick positions, R_wp and weight fractions of one fit folder."""
    x, obs, calc, bg, _ = np.loadtxt(folder / 'fit.txt', unpack=True)
    ticks = {k: np.asarray(v, float) for k, v in json.loads((folder / 'ticks.json').read_text()).items()}
    pos = obs > 0
    rwp = 100 * np.sqrt(np.sum((obs - calc)[pos] ** 2 / obs[pos]) / np.sum(obs[pos]))
    s = json.loads((folder / 'summary.json').read_text())
    if 'record' in s:
        recorded, wt = s['record']['wR'], s['record']['state']['wt']
    else:                                       # search-match pipeline export
        recorded = s['rec']['wR']
        wt = [{p['name']: p['wt'] for p in s['rec']['phases']}[k] for k in ticks]
    if abs(recorded - rwp) > RWP_TOL:
        raise ValueError(f'{folder}: R_wp {rwp:.2f} from the curve, {recorded:.2f} in summary.json')
    if len(wt) != len(ticks):
        raise ValueError(f'{folder}: {len(wt)} weight fractions for {len(ticks)} phases')
    return dict(x=x, obs=obs, calc=calc, bg=bg, ticks=ticks, rwp=rwp, wt=np.asarray(wt, float) / np.sum(wt))


@lru_cache
def sticks_ka12(cif):
    """Cu K-alpha1 and K-alpha2 line positions and intensities of a CIF (K-alpha2 at half weight), 5-150 deg."""
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure
    st = Structure.from_file(cif)
    p1 = XRDCalculator('CuKa1').get_pattern(st, two_theta_range=(5, 150))
    p2 = XRDCalculator('CuKa2').get_pattern(st, two_theta_range=(5, 150))
    return np.concatenate([p1.x, p2.x]), np.concatenate([p1.y, 0.5 * np.asarray(p2.y)])


def observed_fwhm(fit):
    """FWHM (deg) of the strongest line of the best fit's calculated net pattern."""
    x, y = fit['x'], fit['calc'] - fit['bg']
    i = int(np.argmax(y)); half = y[i] / 2
    lo = i
    while lo > 0 and y[lo] > half:
        lo -= 1
    hi = i
    while hi < len(y) - 1 and y[hi] > half:
        hi += 1
    return float(np.clip(x[hi] - x[lo], 0.05, 1.0))


def simulated_pattern(cif, x, fwhm):
    """Pattern of a CIF on grid x: Cu K-alpha1/2 lines broadened by pseudo-Voigt profiles (eta 0.5) of the
    given FWHM; strongest point = 1."""
    pos, inten = sticks_ka12(str(cif))
    w = fwhm / 2
    y = np.zeros_like(x, dtype=float)
    for p0, a in zip(pos, inten):
        if x[0] - 5 * fwhm < p0 < x[-1] + 5 * fwhm:
            d = (x - p0) / w
            y += a * (0.5 / (1 + d ** 2) + 0.5 * np.exp(-np.log(2) * d ** 2))
    return y / np.max(y) if np.max(y) > 0 else y


def predicted_cif(model):
    """CIF of the predicted structure of an alternative model (its starting structure)."""
    c = FITS_DIR / model['dir'] / 'start' / f"{model['target']}.cif"
    if not c.exists():
        raise FileNotFoundError(c)
    return c


@lru_cache
def stick_pattern(cif):
    """Cu K-alpha positions and intensities (strongest = 1) of a CIF, 5-150 deg."""
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure
    p = XRDCalculator('CuKa').get_pattern(Structure.from_file(cif), two_theta_range=(5, 150))
    return np.asarray(p.x), np.asarray(p.y) / np.max(p.y)


def fmt_formula(tok):
    """'ZrSiPt5' -> mathtext in the paper's element order; anything else gets subscripted digits."""
    els = re.findall(r'[A-Z][a-z]?', tok)
    if re.fullmatch(r'([A-Z][a-z]?\d*){2,}', tok) and all(e in ALLEN for e in els):
        return mathtext(tok)
    return re.sub(r'(?<=[A-Za-z)])(\d+\.?\d*)', r'$_{\1}$', tok)


def is_formula(tok):
    """At least two element symbols, e.g. 'CaCu5', 'Y3(Al', '(Ni0.5Pt0.5)Ge'; not 'C2)' or 'Pnma'."""
    els = re.findall(r'[A-Z][a-z]?', tok)
    return (re.fullmatch(r'[A-Z(][A-Za-z\d.(),]*', tok) is not None and len(els) >= 2
            and all(Element.is_valid_symbol(e) for e in els))


def fmt_label(text):
    """Legend text from si_fits.csv: formulas typeset, eta -> greek, X_n -> subscript."""
    # no order-parameter values in legends (PI 2026-11-09): 'X order eta 0.71' -> 'partial X order'
    text = re.sub(r'\(?([A-Za-z/]+) order \(?eta \d*\.?\d+\)?\)?', r'partial \1 order', text)
    text = text.replace('partial partial', 'partial')
    text = re.sub(r'([A-Za-z0-9])_(\d+|[a-z])(?![a-z])', r'\1$_{\2}$', text)
    return re.sub(r'[A-Za-z0-9.(),]+', lambda m: fmt_formula(m[0]) if is_formula(m[0]) else m[0], text)


def phase_label(key):
    """'ZrSiPt5_pred_OQMD-2146113' -> 'ZrPt5Si predicted'."""
    if key in DISPLAY:
        return DISPLAY[key]
    for pattern, render in SPECIAL:
        m = pattern.fullmatch(key)
        if m:
            return render(m)
    words = []
    toks = key.split('_')
    predicted = 'pred' in toks
    for tok in toks:
        # a database name next to a known phase reads like 'the database ground state'; name it only for predictions
        if not predicted and tok in ('oqmd', 'alex', 'mp', 'cod'):
            continue
        if m := re.fullmatch(r'eta(0?\.?)(\d+)', tok):
            words.append('partly ordered')
        elif tok in WORDS:
            words.append(WORDS[tok])
        elif not DROP.fullmatch(tok):
            words.append(fmt_formula(tok) if is_formula(tok) else tok)
    return ' '.join(dict.fromkeys(words))


def parse_zoom(texts):
    """Union of 'lo-hi;lo-hi' windows over the given rows, overlapping windows merged."""
    wins = sorted(tuple(float(v) for v in part.split('-')) for t in texts for part in t.split(';') if part.strip())
    merged = []
    for lo, hi in wins:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    if len(merged) > MAX_ZOOMS:
        print(f'  {len(merged)} zoom windows, drawing the first {MAX_ZOOMS}', file=sys.stderr)
    return merged[:MAX_ZOOMS]


def tick_rows(ax, fit, target, x0, x1, y0, span):
    """Numbered rows of Bragg positions below y0; returns the colour of each phase."""
    colours = []
    for i, (key, tt) in enumerate(fit['ticks'].items()):
        colour = COLOR_TARGET if key == target else COLOR_OTHER
        y1 = y0 - (TICK_GAP_FRAC + i * (TICK_FRAC + TICK_GAP_FRAC)) * span
        ax.vlines(tt[(tt >= x0) & (tt <= x1)], y1 - TICK_FRAC * span, y1, color=colour, lw=0.6)
        ax.text(x0 - 0.006 * (x1 - x0), y1 - TICK_FRAC * span / 2, str(i + 1), color=colour,
                fontsize=SMALL - 1, ha='right', va='center')
        colours.append(colour)
    return colours


# targets whose refined weight fractions are not meaningful (text says so): print phases without fractions
NO_FRACTIONS = {'YAl3Cu2', 'ScPdCu4', 'La5Al3Co2', 'Pd6Cu5Ge'}
# Rows where the point is only that the prediction did not form (PI 2026-11-09): panel (b) overlays the measured
# net pattern, the best fit and the pattern of the predicted structure refined alone; panel (c) shows the two
# difference curves on one scale. No zoom windows.
OVERLAY_FORMULAS = {'PtNiGe2', 'ScPdCu4', 'La5Al3Co2', 'TiPt5Ge', 'Ti2Pt8Si', 'LaNiSi3', 'GdSn6Co4', 'NbNi2B2',
                    'GdFeGe2', 'Gd2Fe3Ge5', 'GdFeGe3', 'LaRuC2', 'Pd6Cu5Ge', 'La2PdAu2', 'Ti2Pt8Ge'}
# Rows without a satisfactory description of the pattern: only the measured pattern and the predicted pattern.
PREDONLY_FORMULAS = {'ScPdCu4', 'Pd6Cu5Ge', 'La2PdAu2', 'Ti2Pt8Ge'}
# zoom insets for prediction-only rows: (lo, hi, deciding line of the predicted structure)
INSET = {}
# prediction-only rows drawn as a computed pattern (fixed narrow lines, no background) stacked above the measured counts
STACKED_FORMULAS = {'ScPdCu4'}
STACKED_FWHM = 0.15   # deg, close to the D500 instrumental width


PRED_COLOR = '#1E88E5'      # simulated predicted pattern (blue); the best-fit residual is drawn in neutral gray
RESID_COLOR = '0.35'


def panel_a(ax, fit, target, zooms, colour=ROLE_COLOR['FoR'], show_wt=True, extra=None, show_diff=True, diff_colour=COLOR_DIFF):
    """Main fit: net pattern, fit, difference, numbered tick rows, phase list."""
    scale = np.max(fit['obs'] - fit['bg'])
    x, net, calc, diff = (fit['x'], (fit['obs'] - fit['bg']) / scale,
                          (fit['calc'] - fit['bg']) / scale, (fit['obs'] - fit['calc']) / scale)
    y_diff = np.min(net) - 0.04 - np.max(diff)
    if not show_diff:
        diff = np.zeros_like(diff)
        y_diff = np.min(net) - 0.02
    ax.plot(x, net, lw=0, marker='x', ms=2.4, mew=0.45, color='black', zorder=2)
    ax.plot(x, calc, lw=0.8, color=colour, zorder=3)
    if show_diff:
        ax.plot(x, diff + y_diff, lw=0.6, color=diff_colour)
    if extra is not None:                       # the predicted pattern (max 1), stacked above the data
        extra = dict(extra, y=extra['y'] + np.max(np.maximum(net, calc)) + 0.12)
        ax.plot(extra['x'], extra['y'], lw=0.9, color=extra['colour'], zorder=4)
    y0 = y_diff + np.min(diff) - 0.04          # top of the tick rows
    n = len(fit['ticks'])
    # The phase list sits in the upper corner over the weaker half of the pattern; the top is
    # raised just enough for it to clear the data, and the tick rows take a fixed fraction.
    top_data = np.maximum(net, calc)
    if extra is not None:
        top_data = np.maximum(top_data, np.interp(x, extra['x'], extra['y']))
    half = 0.45 * (x[-1] - x[0])
    left_max, right_max = np.max(top_data[x <= x[0] + half]), np.max(top_data[x >= x[-1] - half])
    height_in = ax.get_position().height * ax.figure.get_figheight()
    frac = (n + 1 + (extra is not None)) * SMALL * 1.35 / 72 / height_in + 0.03
    rows = n * (TICK_FRAC + TICK_GAP_FRAC) + TICK_GAP_FRAC
    c = frac / (1 - rows)
    top = max(top_data.max() + 0.07, (min(left_max, right_max) + 0.03 - c * y0) / (1 - c))
    span = (top - y0) / (1 - rows)
    colours = tick_rows(ax, fit, target, x[0], x[-1], y0, span)
    for lo, hi in zooms:
        ax.axvspan(lo, hi, color='0.5', alpha=0.12, lw=0, zorder=0)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(top - span, top)
    lines = [(f'{i + 1}  {phase_label(k)}' + (('  ' + (f'{100 * w:.0f}' if w >= 0.005 else '<1') + ' wt%') if show_wt else ''), col)
             for i, (k, w, col) in enumerate(zip(fit['ticks'], fit['wt'], colours))]
    lines.append((rf'$R_{{wp}}$ = {fit["rwp"]:.1f} %', 'black'))
    if extra is not None:
        lines.append((extra['text'], extra['colour']))
    right = right_max <= left_max
    for j, (text, colour) in enumerate(lines):
        ax.text(0.985 if right else 0.015, 0.975 - j * (frac - 0.03) / len(lines), text,
                transform=ax.transAxes, ha='right' if right else 'left', va='top', fontsize=SMALL, color=colour)
    ax.set_yticks([])
    ax.set_ylabel(r'$I - I_{\mathrm{bg}}$ (arb. units)')


def panel_b(ax, models):
    """Difference curves, top to bottom in role order, on one scale; stick-only rows as text."""
    fitted = [m for m in models if m['fit'] is not None]
    curves = [(m, (m['fit']['obs'] - m['fit']['calc']) / np.max(m['fit']['obs'] - m['fit']['bg'])) for m in fitted]
    step = 1.1 * max(np.ptp(d) for _, d in curves)
    for k, (m, d) in enumerate(curves):
        ax.axhline(-k * step, color='0.85', lw=0.3, zorder=0)
        ax.plot(m['fit']['x'], d - k * step, lw=0.8, color=ROLE_COLOR[m['role']])
        extra = '; zoom shows its stick pattern' if m['sticks'] else ''
        t = ax.text(0.012, -k * step + 0.22 * step, rf'{fmt_label(m["label"])}   $R_{{wp}}$ = {m["fit"]["rwp"]:.1f} %{extra}',
                    transform=ax.get_yaxis_transform(), fontsize=SMALL, color=ROLE_COLOR[m['role']],
                    va='bottom', bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))
        fit_inside(t, ax)
    for j, m in enumerate(m for m in models if m['fit'] is None):
        ax.text(0.988, 0.06 + 0.12 * j, f'{fmt_label(m["label"])}: zoom panels', transform=ax.transAxes,
                ha='right', fontsize=SMALL, color=ROLE_COLOR[m['role']])
    ax.set_ylim(np.min(curves[-1][1]) - (len(curves) - 1) * step - 0.05 * step,
                max(np.max(curves[0][1]), 0.5 * step) + 0.05 * step)
    ax.set_yticks([])
    ax.set_ylabel(r'$I_{\mathrm{obs}} - I_{\mathrm{calc}}$', fontsize=SMALL)
    ax.set_xlabel(r'2$\theta$ (deg)')


def panel_b_overlay(ax, best, pred):
    """Top: measured net pattern (crosses), best fit (red), predicted structure refined alone (blue dashed).
    Below, on the same scale: the difference curves of the two fits."""
    f, g = best['fit'], pred['fit']
    scale = np.max(f['obs'] - f['bg'])
    ax.plot(f['x'], (f['obs'] - f['bg']) / scale, lw=0, marker='x', ms=1.6, mew=0.35, color='black', zorder=2)
    ax.plot(f['x'], (f['calc'] - f['bg']) / scale, lw=0.9, color=ROLE_COLOR['FoR'], zorder=3)
    # the refined scale of a structure that does not fit collapses, so show its calculated pattern scaled to the
    # strongest observed line (the difference curve below is that of the refined fit)
    pc = g['calc'] - g['bg']
    ax.plot(g['x'], pc / np.max(pc), lw=0.9, color=ROLE_COLOR['C1'], ls=(0, (4, 2)), zorder=4)
    d1 = (f['obs'] - f['calc']) / scale
    d2 = (g['obs'] - g['calc']) / scale
    gap = 0.12
    o1 = -gap - np.max(d1)
    o2 = o1 + np.min(d1) - gap - np.max(d2)
    for d, o, m in [(d1, o1, best), (d2, o2, pred)]:
        ax.axhline(o, color='0.85', lw=0.3, zorder=0)
        ax.plot(m['fit']['x'], d + o, lw=0.6, color=ROLE_COLOR[m['role']])
    for j, (text, colour) in enumerate([(rf'best fit: {fmt_label(re.sub(r"^best (fit|match): ", "", best["label"]))}   $R_{{wp}}$ = {f["rwp"]:.1f} %', ROLE_COLOR['FoR']),
                                        (rf'{fmt_label(pred["label"])} (pattern scaled to the strongest line); refined alone $R_{{wp}}$ = {g["rwp"]:.1f} %', ROLE_COLOR['C1'])]):
        t = ax.text(0.985, 0.97 - 0.09 * j, text, transform=ax.transAxes, ha='right', va='top', fontsize=SMALL, color=colour,
                    bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))
        fit_inside(t, ax)
    ax.set_ylim(o2 + np.min(d2) - 0.05, 1.12)
    ax.set_yticks([])
    ax.set_ylabel(r'$I - I_{\mathrm{bg}}$; differences', fontsize=SMALL)
    ax.set_xlabel(r'2$\theta$ (deg)')


def fit_inside(text, ax, margin=0.012):
    """Shrink a label until it ends inside the axes (long model names ran past the right edge)."""
    renderer = ax.figure.canvas.get_renderer()
    right = ax.get_window_extent(renderer).x1 - margin * ax.get_window_extent(renderer).width
    while text.get_window_extent(renderer).x1 > right and text.get_fontsize() > 5.0:
        text.set_fontsize(text.get_fontsize() - 0.25)


def draw_sticks(ax, cif, ref, lo, hi, colour):
    """Stick pattern on the reference background, strongest line = strongest observed net line."""
    pos, inten = stick_pattern(cif)
    m = (pos >= lo) & (pos <= hi)
    base = np.interp(pos[m], ref['x'], ref['bg'])
    tops = base + inten[m] * np.max(ref['obs'] - ref['bg'])
    ax.vlines(pos[m], base, tops, color=colour, lw=1.3, zorder=4)
    return np.max(tops, initial=0)


def panel_c(ax, models, ref, lo, hi, target=None):
    """Observed counts and each model's calculated pattern (or stick pattern) in one window."""
    m = (ref['x'] >= lo) & (ref['x'] <= hi)
    ax.plot(ref['x'][m], ref['obs'][m], lw=0, marker='x', ms=3, mew=0.55, color='black', zorder=2)
    top = np.max(ref['obs'][m])
    for mod in models:
        if mod['sticks']:
            top = max(top, draw_sticks(ax, mod['sticks'], ref, lo, hi, ROLE_COLOR[mod['role']]))
        else:
            f = mod['fit']
            mf = (f['x'] >= lo) & (f['x'] <= hi)
            ax.plot(f['x'][mf], f['calc'][mf], lw=ZOOM_LW.get(mod['role'], 1.2), color=ROLE_COLOR[mod['role']], zorder=ZOOM_Z.get(mod['role'], 4),
                    ls=ROLE_DASH.get(mod['role'], '-'))
            top = max(top, np.max(f['calc'][mf], initial=0))
    ax.set_xlim(lo, hi)
    ax.set_ylim(top=1.08 * top)
    # Bragg positions of every phase of the main fit, as in panel (a): green = phase at the nominal composition
    y0, y1 = ax.get_ylim()
    span = y1 - y0
    n = len(ref['ticks'])
    for i, (key, tt) in enumerate(ref['ticks'].items()):
        colour = COLOR_TARGET if key == target else COLOR_OTHER
        yt = y0 - (0.04 + i * 0.07) * span
        ax.vlines(tt[(tt >= lo) & (tt <= hi)], yt - 0.05 * span, yt, color=colour, lw=1.0)
    # phases that only an alternative model adds (e.g. the predicted structure): its line positions in that model's colour
    extra = []
    for mod in models:
        f = mod.get('fit')
        if f is None or mod['role'] == 'FoR':
            continue
        for key, tt in f['ticks'].items():
            if key not in ref['ticks']:
                extra.append((tt, ROLE_COLOR[mod['role']]))
    for j, (tt, colour) in enumerate(extra):
        yt = y0 - (0.04 + (n + j) * 0.07) * span
        ax.vlines(tt[(tt >= lo) & (tt <= hi)], yt - 0.05 * span, yt, color=colour, lw=1.0)
    ax.set_ylim(bottom=y0 - (0.04 + (n + len(extra)) * 0.07) * span)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.set_yticks([])


def load_models(group):
    return [dict(role=r.role, label=r.label, sticks=r.sticks, target=r.target_phase, dir=r.fit_dir, zoom=r.zoom,
                 fit=load_fit(FITS_DIR / r.fit_dir) if r.fit_dir else None) for r in group.itertuples()]


def make_figure(row, group, sample):
    zooms = parse_zoom(group['zoom'])
    models = load_models(group)
    head = next(m for m in models if m['fit'] is not None)      # FoR, or the first contrast if pending
    alt = next((m for m in models if m['role'] != 'FoR' and m['fit'] is not None), None)
    if sample['formula'] in MATCHED_FORMULAS:
        return matched_figure(row, models, sample)
    if sample['formula'] in STACKED_FORMULAS and alt is not None:
        return stacked_predonly_figure(row, head, alt, sample)
    if sample['formula'] in PREDONLY_FORMULAS and alt is not None:
        return predonly_figure(row, head, alt, sample)
    if sample['formula'] in OVERLAY_FORMULAS and alt is not None:
        return overlay_figure(row, head, alt, sample)
    fig = plt.figure(figsize=FIGSIZE)
    ncol = 1 if len(zooms) <= 3 else 2
    outer = fig.add_gridspec(1, 2 if zooms else 1, width_ratios=[3.3, ncol] if zooms else None,
                             left=0.055, right=0.985, top=0.93, bottom=0.10, wspace=0.12)
    left = outer[0].subgridspec(2, 1, height_ratios=[2.1, 1], hspace=0.06)
    ax_a = fig.add_subplot(left[0])
    ax_b = fig.add_subplot(left[1], sharex=ax_a)
    plt.setp(ax_a.get_xticklabels(), visible=False)
    panel_a(ax_a, head['fit'], head['target'], zooms, ROLE_COLOR[head['role']],
            show_wt=sample['formula'] not in NO_FRACTIONS)
    panel_b(ax_b, models)
    cat = f'{sample["category"]}' + (f'/{int(sample["alt_category"])}' if pd.notna(sample['alt_category']) else '')
    ax_a.set_title(mathtext(sample["formula"]), loc='left', fontsize=FONTSIZE + 1)
    if head['role'] != 'FoR':
        ax_b.text(0.988, 0.9, f'main fit pending; (a) shows {head["role"]}', transform=ax_b.transAxes,
                  ha='right', va='top', fontsize=SMALL, color='0.3')
    for ax, letter in [(ax_a, 'a'), (ax_b, 'b')]:
        ax.text(-0.035, 1.0, f'({letter})', transform=ax.transAxes, fontsize=FONTSIZE, va='top', ha='right')
    if zooms:
        nrow = -(-len(zooms) // ncol)
        right = outer[1].subgridspec(nrow, ncol, hspace=0.45, wspace=0.22)
        for k, (lo, hi) in enumerate(zooms):
            ax = fig.add_subplot(right[k // ncol, k % ncol])
            panel_c(ax, models, head['fit'], lo, hi, target=head['target'])
            if k == 0:
                ax.set_title('(c)', loc='left', fontsize=FONTSIZE)
            if k % ncol == 0:
                ax.set_ylabel('counts', fontsize=SMALL)
            if k >= len(zooms) - ncol:
                ax.set_xlabel(r'2$\theta$ (deg)', fontsize=SMALL)
    for ax in fig.axes:
        ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / f'fig_S_{row:02d}_{sample["formula"]}.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


# Rows whose contrast models are all refined with the same adjustable parameters (ZrPt5Si, TiPt5Si): (a) the main fit;
# (b) the difference curves of the matched models; (c) the ordered (C1) and site-mixed (C2) models at the
# ordering lines; (d) the ordered model with (C1) and without (C3) the minor phase at its lines.
MATCHED_FORMULAS = {'ZrPt5Si', 'TiPt5Si'}


def matched_figure(row, models, sample):
    by = {m['role']: m for m in models}
    head, c1, c2, c3 = by['FoR'], by['C1'], by['C2'], by['C3']
    zc, zd = parse_zoom([c1['zoom']]), parse_zoom([c3['zoom']])
    fig = plt.figure(figsize=(FIGSIZE[0], FIGSIZE[1] * 1.45))
    outer = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.0], left=0.055, right=0.985, top=0.955, bottom=0.07, hspace=0.3)
    top = outer[0].subgridspec(1, 2, width_ratios=[3.3, 2], wspace=0.12)
    left = top[0].subgridspec(2, 1, height_ratios=[2.1, 1.3], hspace=0.06)
    ax_a = fig.add_subplot(left[0])
    ax_b = fig.add_subplot(left[1], sharex=ax_a)
    plt.setp(ax_a.get_xticklabels(), visible=False)
    panel_a(ax_a, head['fit'], head['target'], zc + zd, ROLE_COLOR['FoR'])
    panel_b(ax_b, [c1, c2, c3])
    ax_a.set_title(mathtext(sample['formula']), loc='left', fontsize=FONTSIZE + 1)
    for ax, letter in [(ax_a, 'a'), (ax_b, 'b')]:
        ax.text(-0.035, 1.0, f'({letter})', transform=ax.transAxes, fontsize=FONTSIZE, va='top', ha='right')
    right = top[1].subgridspec(2, 2, hspace=0.45, wspace=0.22)
    for k, (lo, hi) in enumerate(zc):
        ax = fig.add_subplot(right[k // 2, k % 2])
        panel_c(ax, [c1, c2], c1['fit'], lo, hi, target=c1['target'])
        if k == 0:
            ax.set_title('(c)', loc='left', fontsize=FONTSIZE)
        if k % 2 == 0:
            ax.set_ylabel('counts', fontsize=SMALL)
        if k >= 2:
            ax.set_xlabel(r'2$\theta$ (deg)', fontsize=SMALL)
    bottom = outer[1].subgridspec(1, len(zd), wspace=0.25)
    for k, (lo, hi) in enumerate(zd):
        ax = fig.add_subplot(bottom[k])
        panel_c(ax, [c1, c3], c1['fit'], lo, hi, target=c1['target'])
        if k == 0:
            ax.set_title('(d)', loc='left', fontsize=FONTSIZE)
            ax.set_ylabel('counts', fontsize=SMALL)
        ax.set_xlabel(r'2$\theta$ (deg)', fontsize=SMALL)
    for ax in fig.axes:
        ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / f'fig_S_{row:02d}_{sample["formula"]}.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


def smooth_background(x, y, width_deg=3.0):
    """Rough smooth baseline: rolling 10th percentile over width_deg, then smoothed over the same width."""
    from scipy.ndimage import percentile_filter, uniform_filter1d
    n = max(5, int(round(width_deg / np.median(np.diff(x)))))
    return uniform_filter1d(percentile_filter(y, 10, size=n, mode='nearest'), n, mode='nearest')


def stacked_predonly_figure(row, head, alt, sample):
    """Measured counts below; the pattern computed from the predicted structure (fixed line width, no background)
    above them. No fit."""
    f = head['fit']
    fig, ax = plt.subplots(figsize=(FIGSIZE[0], FIGSIZE[1] * 0.8))
    fig.subplots_adjust(left=0.05, right=0.985, top=0.92, bottom=0.13)
    obs = f['obs']; lo, hi = float(np.min(obs)), float(np.max(obs)); span = hi - lo
    ax.plot(f['x'], obs, lw=0, marker='x', ms=1.8, mew=0.4, color='black', zorder=2, label='measured counts')
    sim = simulated_pattern(predicted_cif(alt), f['x'], STACKED_FWHM)
    name = re.sub(r"\s*(alone|in place of .*)$", "", fmt_label(alt["label"])).replace("ordered ", "")
    ax.plot(f['x'], hi + 0.08 * span + sim * span, lw=0.8, color=PRED_COLOR, zorder=3, label=f'{name}, computed')
    ax.legend(loc='upper right', fontsize=SMALL, frameon=False)
    ax.set_xlim(f['x'][0], f['x'][-1])
    ax.set_ylim(lo - 0.05 * span, hi + 1.2 * span)
    ax.set_yticks([])
    ax.set_ylabel('intensity (arb. units)')
    ax.set_xlabel(r'2$\theta$ (deg)')
    ax.set_title(mathtext(sample['formula']), loc='left', fontsize=FONTSIZE + 1)
    ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / f'fig_S_{row:02d}_{sample["formula"]}.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


def predonly_figure(row, head, alt, sample):
    """Measured counts and the predicted pattern, scaled to the strongest line, on a smooth background; no fit."""
    f, g = head['fit'], alt['fit']
    fig, ax = plt.subplots(figsize=(FIGSIZE[0], FIGSIZE[1] * 0.62))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.9, bottom=0.17)
    bg = smooth_background(f['x'], f['obs'])
    ax.plot(f['x'], f['obs'], lw=0, marker='x', ms=1.8, mew=0.4, color='black', zorder=2, label='measured counts')
    sim = simulated_pattern(predicted_cif(alt), f['x'], observed_fwhm(f))
    ax.plot(f['x'], bg + sim * np.max(f['obs'] - bg), lw=0.9, color=PRED_COLOR, zorder=3,
            label=f'{re.sub(r"\s*(alone|in place of .*)$", "", fmt_label(alt["label"])).replace("ordered ", "")}, simulated')
    ax.legend(loc='upper right', fontsize=SMALL, frameon=False)
    ax.set_xlim(f['x'][0], f['x'][-1])
    ax.set_ylabel('counts')
    if sample['formula'] in INSET:
        lo, hi, mark = INSET[sample['formula']]
        ins = ax.inset_axes([0.62, 0.22, 0.25, 0.48])
        m = (f['x'] >= lo) & (f['x'] <= hi)
        ins.plot(f['x'][m], f['obs'][m], lw=0, marker='x', ms=2.5, mew=0.5, color='black')
        ins.plot(f['x'][m], (bg + sim * np.max(f['obs'] - bg))[m], lw=1.0, color=PRED_COLOR)
        mark = float(f['x'][m][np.argmax(sim[m])])   # the simulated line itself (DFT cell)
        ins.axvline(mark, color='0.5', lw=0.6, ls=':')
        ins.set_xlim(lo, hi)
        ins.set_yticks([])
        ins.tick_params(length=2, width=0.5, labelsize=SMALL - 1)
        ins.text(mark, 0.97, f' {mark:.1f}' + r'$^\circ$', transform=ins.get_xaxis_transform(), fontsize=SMALL - 1, va='top', color='0.3')
    ax.set_xlabel(r'2$\theta$ (deg)')
    ax.set_title(mathtext(sample['formula']), loc='left', fontsize=FONTSIZE + 1)
    ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / f'fig_S_{row:02d}_{sample["formula"]}.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


def overlay_figure(row, head, alt, sample):
    """One panel: best fit with its difference curve, and the predicted pattern (scaled) overlaid."""
    f, g = head['fit'], alt['fit']
    # simulated from the predicted CIF (a failed refinement of a structure that is not there is no guide to its pattern)
    sim = simulated_pattern(predicted_cif(alt), f['x'], observed_fwhm(f))
    name = re.sub(r'\s*(alone|in place of .*)$', '', fmt_label(alt['label'])).replace('ordered ', '')
    extra = dict(x=f['x'], y=sim, colour=PRED_COLOR, text=f'{name}, simulated')
    fig, ax = plt.subplots(figsize=(FIGSIZE[0], FIGSIZE[1] * 0.95))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.93, bottom=0.10)
    # the predicted compound did not form in these rows, so no phase is drawn in the green 'nominal composition' colour
    panel_a(ax, f, None, [], ROLE_COLOR['FoR'], show_wt=sample['formula'] not in NO_FRACTIONS, extra=extra,
            show_diff=True, diff_colour=RESID_COLOR)
    ax.set_xlabel(r'2$\theta$ (deg)')
    cat = f'{sample["category"]}' + (f'/{int(sample["alt_category"])}' if pd.notna(sample['alt_category']) else '')
    ax.set_title(mathtext(sample['formula']), loc='left', fontsize=FONTSIZE + 1)
    ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / f'fig_S_{row:02d}_{sample["formula"]}.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


def net(fit, x0, x1):
    """Observed and calculated net intensity on [x0, x1], scaled to a maximum of 1."""
    m = (fit['x'] >= x0) & (fit['x'] <= x1)
    s = np.max(fit['obs'][m] - fit['bg'][m])
    return fit['x'][m], (fit['obs'][m] - fit['bg'][m]) / s, (fit['calc'][m] - fit['bg'][m]) / s


def ascast_figure(fits, samples):
    """One row per annealed target: annealed main fit above the as-cast fit, and zooms."""
    ascast = fits[fits['role'] == 'ascast']
    # the two isostructural Pt5 compounds side by side, then GdPt2Ge
    order = {'ZrPt5Si': 0, 'TiPt5Si': 1, 'GdPt2Ge': 2}
    ascast = ascast.sort_values('formula', key=lambda c: c.map(order).fillna(99))
    fig = plt.figure(figsize=(7.0, 2.5 * len(ascast) + 0.3))
    outer = fig.add_gridspec(len(ascast), 2, width_ratios=[2.4, 1.6], left=0.04, right=0.985,
                             top=0.97, bottom=0.06, hspace=0.42, wspace=0.08)
    for i, (_, a) in enumerate(ascast.iterrows()):
        f_row = fits[(fits['table_row'] == a['table_row']) & (fits['role'] == 'FoR')].iloc[0]
        pair = [('annealed', f_row, ROLE_COLOR['FoR']), ('as cast', a, ROLE_COLOR['ascast'])]
        fitted = [(name, r, col, load_fit(FITS_DIR / r['fit_dir'])) for name, r, col in pair]
        x0, x1 = max(f['x'][0] for *_, f in fitted), min(f['x'][-1] for *_, f in fitted)
        ax = fig.add_subplot(outer[i, 0])
        for k, (name, r, col, f) in enumerate(fitted):
            x, o, c = net(f, x0, x1)
            off = 1.35 * (1 - k)
            ax.plot(x, o + off, lw=0, marker='x', ms=1.8, mew=0.35, color='black')
            ax.plot(x, c + off, lw=0.7, color=col)
            tt = f['ticks'][r['target_phase']]
            ax.vlines(tt[(tt >= x0) & (tt <= x1)], off - 0.16, off - 0.06, color=COLOR_TARGET, lw=0.6)
            wt = 100 * f['wt'][list(f['ticks']).index(r['target_phase'])]
            ax.text(0.99, (off + 0.75) / 2.55 + 0.05, f'{name}: {mathtext(r["formula"])} {wt:.0f} wt%, '
                    rf'$R_{{wp}}$ = {f["rwp"]:.1f} %', transform=ax.transAxes, ha='right', fontsize=SMALL, color=col)
        ax.set_xlim(x0, x1)
        ax.set_ylim(-0.2, 2.45)
        ax.set_yticks([])
        ax.set_title(mathtext(a["formula"]), loc='left', fontsize=FONTSIZE + 0.5)
        ax.set_xlabel(r'2$\theta$ (deg)', fontsize=SMALL)
        zooms = parse_zoom([a['zoom']])
        right = outer[i, 1].subgridspec(1, len(zooms), wspace=0.22)
        for k, (lo, hi) in enumerate(zooms):
            az = fig.add_subplot(right[k])
            for j, (name, r, col, f) in enumerate(fitted):
                x, o, c = net(f, lo, hi)
                off = 1.15 * (1 - j)
                az.plot(x, o + off, lw=0, marker='x', ms=2.2, mew=0.4, color='black')
                az.plot(x, c + off, lw=1.4, color=col)
                tt = f['ticks'][r['target_phase']]
                az.vlines(tt[(tt >= lo) & (tt <= hi)], off - 0.13, off - 0.04, color=COLOR_TARGET, lw=0.8)
            az.set_xlim(lo, hi)
            az.set_ylim(-0.18, 2.3)
            az.set_yticks([])
            az.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
            az.set_xlabel(r'2$\theta$ (deg)', fontsize=SMALL)
    for ax in fig.axes:
        ax.tick_params(length=2.5, width=0.5, labelsize=SMALL)
    out = OUT_DIR / 'fig_S_ascast_annealed.pdf'
    fig.savefig(out)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('rows', nargs='*', type=int, help='Table 1 rows to draw (default: all, plus the as-cast figure)')
    args = ap.parse_args()
    fits = pd.read_csv(FITS_CSV, dtype={'sample': str}, keep_default_na=False)
    samples = pd.read_csv(SAMPLES_CSV, dtype={'alt_category': 'Int64'}).set_index('table_row')
    bad = fits[~fits['role'].isin(ROLES + ['ascast', 'pending'])]
    if not bad.empty:
        raise ValueError(f'unknown roles: {sorted(set(bad["role"]))}')
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for row, group in fits.groupby('table_row'):
        if group['formula'].iloc[0] != samples.loc[row, 'formula']:
            raise ValueError(f'row {row}: {group["formula"].iloc[0]} in si_fits.csv, {samples.loc[row, "formula"]} in Table 1')
        group = group[group['role'].isin(ROLES)].sort_values('role', key=lambda r: r.map(ROLES.index), kind='stable')
        if (args.rows and row not in args.rows) or not (group['fit_dir'] != '').any():
            continue
        out = make_figure(int(row), group, samples.loc[row])
        print(f'{out}  ' + '  '.join(f'{r.role} {r.fit_dir or r.sticks}' for r in group.itertuples()))
    if not args.rows:
        print(ascast_figure(fits, samples))


if __name__ == '__main__':
    main()
