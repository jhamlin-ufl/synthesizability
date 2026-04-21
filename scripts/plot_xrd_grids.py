#!/usr/bin/env python3
"""
Three PDF grid plots of measured XRD patterns vs OQMD-predicted patterns.
  plot_xrd_ss.pdf   — solid solution samples
  plot_xrd_mp.pdf   — multi-phase samples
  plot_xrd_p.pdf    — predicted-structure-found samples

Each subplot: red (OQMD predicted, bottom) + blue (experimental, above).
Offset is minimum to prevent overlap; both normalized to peak = 1.
"""
import sys, os, tempfile
from pathlib import Path

sys.path.insert(0, str(Path.home() / 'repos/search_match'))
sys.path.insert(0, 'src')

from synthesizability.parsers.xrd import parse_xrd_file, is_xrd_file
from search_match.xrd_simulator import simulate_xrd_pattern

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pymatgen.io.cif import CifParser, CifWriter
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TWO_THETA_RANGE = (30, 115)   # degrees — common to all patterns
FWHM = 0.5                    # degrees — broadening for simulated patterns
GAP = 0.15                    # normalized units between predicted top and exp bottom
NCOLS = 4                     # columns for SS and MP grids
NCOLS_P = 3                   # columns for P grid (6 samples fits 3×2 exactly)
PREFER_PANALYTICAL = {473}    # Sc3CuPd2: prefer Panalytical over Siemens
RAW_DIR = Path('data/raw')
OQMD_DIR = Path('data/external/oqmd_structures')
OUT_DIR = Path('data/temp')

# Caglioti W parameter → FWHM ≈ sqrt(W) at mid angles
W_PARAM = 0.2 * (FWHM / 0.5) ** 2
INSTRUMENT_PARAMS = {'U': 0.0, 'V': 0.0, 'W': W_PARAM, 'X': 0.0, 'Y': 0.0}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_samples():
    fa = pd.read_csv('data/temp/Arc Melted Samples - Final Analysis.csv')
    sheet1 = pd.read_csv('data/temp/Arc Melted Samples - Sheet1.csv')
    synth = pd.read_csv('data/processed/synthesis_data.csv')

    fa['Sample Number'] = pd.to_numeric(fa['Sample Number'], errors='coerce')
    sheet1['Sample Number'] = pd.to_numeric(sheet1['Sample Number'], errors='coerce')
    sheet1 = sheet1.dropna(subset=['Sample Number'])
    sheet1['Sample Number'] = sheet1['Sample Number'].astype(int)

    fa_nums = [482 if n == 561 else n
               for n in fa['Sample Number'].dropna().astype(int).tolist()]

    sub = synth[synth['sample_number'].isin(fa_nums) &
                (synth['prediction_list'] != 'Diffusion Model')].copy()
    sub = sub.join(sheet1.set_index('Sample Number')[['XRD Result']], on='sample_number')
    sub['XRD Result'] = sub['XRD Result'].fillna('--')
    return sub.sort_values('sample_number').reset_index(drop=True)


def find_xrd_file(num):
    dirs = [d for d in RAW_DIR.iterdir() if d.name.startswith(f'{num:04d}_')]
    if not dirs:
        return None
    files = [f for f in sorted(dirs[0].iterdir())
             if f.suffix in ('.xy', '.txt') and is_xrd_file(f)]
    xy  = [f for f in files if f.suffix == '.xy']
    txt = [f for f in files if f.suffix == '.txt']
    if num in PREFER_PANALYTICAL:
        return xy[0] if xy else (txt[0] if txt else None)
    return txt[0] if txt else (xy[0] if xy else None)


def find_oqmd_cif(formula, entry_id):
    cif_dir = OQMD_DIR / formula
    if not cif_dir.exists():
        return None
    cif = cif_dir / f'{int(entry_id)}.cif'
    if cif.exists():
        return cif
    cifs = sorted(cif_dir.glob('*.cif'))
    return cifs[0] if cifs else None


# ---------------------------------------------------------------------------
# XRD helpers
# ---------------------------------------------------------------------------

def load_exp(num):
    """Load experimental pattern clipped to TWO_THETA_RANGE. Returns (tt, ii) or None."""
    fpath = find_xrd_file(int(num))
    if fpath is None:
        return None
    try:
        d = parse_xrd_file(fpath)
    except Exception:
        return None
    tt, ii = d['two_theta'], d['intensity']
    mask = (tt >= TWO_THETA_RANGE[0]) & (tt <= TWO_THETA_RANGE[1])
    if mask.sum() < 10:
        return None
    return tt[mask], ii[mask].astype(float)


def simulate_oqmd(formula, entry_id):
    """Simulate OQMD predicted pattern. Returns (tt, ii) or None."""
    cif_path = find_oqmd_cif(formula, entry_id)
    if cif_path is None:
        return None
    try:
        buf = TWO_THETA_RANGE[1] - TWO_THETA_RANGE[0]
        sim_range = (TWO_THETA_RANGE[0] - 2, TWO_THETA_RANGE[1] + 2)
        pattern = simulate_xrd_pattern(
            str(cif_path),
            two_theta_range=sim_range,
            step_size=0.02,
            wavelength='CuKa',
            background_type='polynomial',
            noise_level=0.0,
            instrument_params=INSTRUMENT_PARAMS,
        )
        tt = pattern['angles']
        ii = pattern['intensities'] - pattern['background']
        ii = np.maximum(ii, 0)
        mask = (tt >= TWO_THETA_RANGE[0]) & (tt <= TWO_THETA_RANGE[1])
        return tt[mask], ii[mask]
    except Exception as e:
        print(f'  WARNING: simulation failed for {formula}: {e}')
        return None


def norm_peak(ii):
    """Normalize so peak = 1, baseline near 0."""
    mn, mx = ii.min(), ii.max()
    return (ii - mn) / (mx - mn) if mx > mn else np.zeros_like(ii)


# ---------------------------------------------------------------------------
# Grid layout: empty cells in upper-right
# ---------------------------------------------------------------------------

def sample_to_grid(idx, n_samples, ncols):
    """
    Map sample index → (row, col) with empty cells in upper-right.
    Top row has (n_samples mod ncols) samples; empty slots fill its right end.
    If n_samples is a multiple of ncols, no empty cells.
    """
    nrows = (n_samples + ncols - 1) // ncols
    n_top = n_samples - (nrows - 1) * ncols   # samples in the (partial) top row
    if idx < n_top:
        return 0, idx
    j = idx - n_top
    return 1 + j // ncols, j % ncols


# ---------------------------------------------------------------------------
# Core plotting function
# ---------------------------------------------------------------------------

def plot_group(samples_df, title, ncols, out_path):
    """
    samples_df: rows with sample_number, formula, oqmd_entry_id, oqmd_stability
    """
    records = []
    for _, row in samples_df.iterrows():
        num = int(row['sample_number'])
        formula = row['formula']
        exp = load_exp(num)
        if exp is None:
            print(f'  SKIP {num} {formula}: no experimental data')
            continue
        oqmd = simulate_oqmd(formula, row['oqmd_entry_id']) if pd.notna(row['oqmd_entry_id']) else None
        records.append((num, formula, exp, oqmd))

    n = len(records)
    if n == 0:
        print(f'  No data for group, skipping.')
        return

    nrows = (n + ncols - 1) // ncols

    # Figure size: fit on 8.5×11 portrait with caption room (~1.5 in)
    # Target: 7 in wide, leaving ~1.5 in for caption → height ≤ 9.5 in
    subplot_h = min(9.0 / nrows, 2.4)
    subplot_w = 7.0 / ncols
    fig_w = subplot_w * ncols
    fig_h = subplot_h * nrows

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(fig_w, fig_h),
                              sharex=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    # Hide all axes first, then show used ones
    for ax in axes.flat:
        ax.set_visible(False)

    for idx, (num, formula, exp, oqmd) in enumerate(records):
        row_i, col_i = sample_to_grid(idx, n, ncols)
        ax = axes[row_i, col_i]
        ax.set_visible(True)

        tt_exp, ii_exp = exp
        ii_exp_n = norm_peak(ii_exp)

        if oqmd is not None:
            tt_oqmd, ii_oqmd = oqmd
            ii_oqmd_n = norm_peak(ii_oqmd)
            offset = 1.0 + GAP   # exp baseline sits above oqmd peak
            ax.plot(tt_oqmd, ii_oqmd_n,        color='#D32F2F', lw=0.8, rasterized=True)
            ax.plot(tt_exp,  ii_exp_n + offset, color='#1565C0', lw=0.8, rasterized=True)
            ax.set_ylim(-0.1, offset + 1.15)
        else:
            ax.plot(tt_exp, ii_exp_n + 1.0 + GAP, color='#1565C0', lw=0.8, rasterized=True)
            ax.set_ylim(-0.1, 2.25)

        ax.set_xlim(*TWO_THETA_RANGE)
        ax.set_yticks([])
        ax.tick_params(axis='x', direction='in', top=True,
                       labelbottom=False, labelsize=7)

        # Formula label — upper left, inside axes
        ax.text(0.03, 0.97, formula,
                transform=ax.transAxes, fontsize=6.5,
                va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))

    # x-axis tick labels only on bottom visible row
    bottom_row = nrows - 1
    for col_i in range(ncols):
        ax = axes[bottom_row, col_i]
        if ax.get_visible():
            ax.tick_params(axis='x', labelbottom=True, labelsize=7)
            ax.set_xlabel(r'2$\theta$ (°)', fontsize=7, labelpad=2)

    # y-axis label on leftmost column, middle row
    mid_row = nrows // 2
    for r in range(nrows):
        ax = axes[r, 0]
        if ax.get_visible():
            ax.set_ylabel('Intensity (arb. units)', fontsize=7)
            break   # only first visible leftmost

    # Legend — small, placed in one of the bottom-left axis
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='#D32F2F', lw=1.2, label='OQMD predicted'),
        Line2D([0], [0], color='#1565C0', lw=1.2, label='Experimental'),
    ]
    # Put legend in the bottom-left subplot
    ax_legend = axes[bottom_row, 0]
    ax_legend.legend(handles=legend_elements, loc='lower left',
                     fontsize=6, framealpha=0.8)

    fig.tight_layout(pad=0.4, h_pad=0.3, w_pad=0.3)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}  ({n} samples, {nrows}×{ncols} grid)')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print('Loading sample list...')
    sub = load_samples()

    ss = sub[sub['XRD Result'] == 'SS'].reset_index(drop=True)
    mp = sub[sub['XRD Result'] == 'MP'].reset_index(drop=True)
    p  = sub[sub['XRD Result'] == 'P' ].reset_index(drop=True)

    print(f'SS: {len(ss)}  MP: {len(mp)}  P: {len(p)}')
    print(f'2-theta range: {TWO_THETA_RANGE[0]}–{TWO_THETA_RANGE[1]}°,  FWHM={FWHM}°')

    print('\nPlotting SS grid...')
    plot_group(ss, 'Solid Solution', NCOLS,   OUT_DIR / 'plot_xrd_ss.pdf')

    print('\nPlotting MP grid...')
    plot_group(mp, 'Multi-Phase',    NCOLS,   OUT_DIR / 'plot_xrd_mp.pdf')

    print('\nPlotting P grid...')
    plot_group(p,  'Predicted',      NCOLS_P, OUT_DIR / 'plot_xrd_p.pdf')

    print('\nDone.')
