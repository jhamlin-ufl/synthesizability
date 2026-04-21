#!/usr/bin/env python3
"""
Plot XRD patterns for all solid-solution (SS) samples.

Produces a grid plot — one subplot per sample, labeled by formula and
XRD result classification from the Sheet1 CSV.

When a sample has two XRD datasets (e.g. both Siemens .txt and Panalytical .xy),
both are shown: one on the left y-axis and one on the right y-axis, each
labeled by instrument.

Samples are taken from the Final Analysis CSV (preferred) and Sheet1 CSV
(XRD classification). For samples replaced by remakes, only the remake is used.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, 'src')
from synthesizability.parsers.xrd import parse_xrd_file, is_xrd_file

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Sample definitions
# SS samples from Sheet1, using Final Analysis numbers (remakes replace originals):
#   447 (MoTaTi2) → 541 (MoTiTa2)
#   448 (MoTa2Ti) → 542 (MoTi2Ta)
#   486 (Al4FeCo3) → 545 (Al4FeCo3)
#   487 (Al8FeCo7) → 535 (Al8FeCo7)
# "two phase SS" (451, 462) included alongside single-phase SS.
# ---------------------------------------------------------------------------

SS_SAMPLES = [
    (449, 'HfTa4Zr'),
    (450, 'HfMoTa2'),
    (451, 'MoTaZr2'),
    (452, 'Hf2MoTi'),
    (453, 'HfMoTi'),
    (454, 'Nb2TaZr'),
    (455, 'NbTa2Zr'),
    (456, 'Hf2MoNb'),
    (457, 'MoNbZr2'),
    (460, 'HfMo2Zr'),
    (461, 'Hf2MoRe3'),
    (462, 'Hf2Nb3Ru'),
    (465, 'Cu5GePd6'),
    (466, 'Hf2MoIr'),
    (469, 'ScCuPd2'),
    (470, 'MoNbTa2'),
    (471, 'MoNb2Ta'),
    (472, 'ScCu4Pd'),
    (473, 'Sc3CuPd2'),
    (474, 'Sc2CuPd'),
    (477, 'ScPd5Au2'),
    (478, 'Sc2Pd5Au'),
    (479, 'MnAl2Au6'),
    (481, 'NiGePt2'),
    (526, 'Re4TaZr'),
    (535, 'Al8FeCo7'),   # remake of 487
    (541, 'MoTiTa2'),    # remake of 447
    (542, 'MoTi2Ta'),    # remake of 448
    (545, 'Al4FeCo3'),   # remake of 486
]

RAW_DIR = Path('data/raw')
SHEET1_CSV = Path('data/temp/Arc Melted Samples - Sheet1.csv')

# Colors for dual-axis plots
COLOR_LEFT = 'C0'    # Siemens (or first dataset)
COLOR_RIGHT = 'C3'   # Panalytical (or second dataset)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sheet1_labels() -> dict[int, str]:
    """Return {sample_num: xrd_result} from Sheet1 CSV."""
    labels = {}
    with open(SHEET1_CSV, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                num = int(row['Sample Number'])
                result = (row.get('XRD Result') or '').strip()
                if result:
                    labels[num] = result
            except (ValueError, KeyError):
                continue
    return labels


def find_xrd_files(sample_num: int) -> list[Path]:
    """Return all valid XRD data files for a sample, sorted .xy before .txt."""
    dirs = [d for d in RAW_DIR.iterdir() if d.name.startswith(f'{sample_num:04d}_')]
    if not dirs:
        return []
    sample_dir = dirs[0]
    xrd_files = [f for f in sorted(sample_dir.iterdir())
                 if f.suffix in ('.xy', '.txt') and is_xrd_file(f)]
    # Put .xy files first
    return sorted(xrd_files, key=lambda f: (0 if f.suffix == '.xy' else 1, f.name))


def normalize(intensity: np.ndarray) -> np.ndarray:
    mn, mx = intensity.min(), intensity.max()
    return (intensity - mn) / (mx - mn) if mx > mn else np.zeros_like(intensity)


def instrument_name(data: dict) -> str:
    return data.get('instrument', 'Unknown')


# ---------------------------------------------------------------------------
# Grid plot
# ---------------------------------------------------------------------------

def plot_grid(records, sheet1_labels: dict[int, str]):
    """
    records: list of (sample_num, formula, [list of (path, parsed_data)])
    """
    n = len(records)
    ncols = 5
    nrows = (n + ncols - 1) // ncols

    # Taller rows to give room between x-axis label and next subplot title
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4 * ncols, 3.5 * nrows),
                              sharex=False, sharey=False)
    axes_flat = axes.flatten()

    for i, (num, formula, datasets) in enumerate(records):
        ax = axes_flat[i]

        xrd_result = sheet1_labels.get(num, 'SS')
        title = f'{num} · {formula}  [{xrd_result}]'
        ax.set_title(title, fontsize=7.5, pad=4)

        if len(datasets) == 1:
            # Single dataset — simple plot
            fpath, data = datasets[0]
            tt = data['two_theta']
            ii = normalize(data['intensity'])
            ax.plot(tt, ii, lw=0.8, color=COLOR_LEFT, rasterized=True)
            ax.set_xlim(tt[0], tt[-1])
            ax.set_ylim(-0.05, 1.15)
            ax.tick_params(axis='both', labelsize=6)
            ax.set_yticks([])
            ax.set_xlabel(r'2θ (°)', fontsize=7, labelpad=2)

        else:
            # Two datasets — left and right y-axes
            fpath0, data0 = datasets[0]
            fpath1, data1 = datasets[1]

            instr0 = instrument_name(data0)
            instr1 = instrument_name(data1)

            # Extend x-range to the union of both patterns
            tt0, ii0 = data0['two_theta'], normalize(data0['intensity'])
            tt1, ii1 = data1['two_theta'], normalize(data1['intensity'])

            ax.plot(tt0, ii0, lw=0.8, color=COLOR_LEFT, rasterized=True,
                    label=instr0)
            ax.set_xlim(min(tt0[0], tt1[0]), max(tt0[-1], tt1[-1]))
            ax.set_ylim(-0.05, 1.15)
            ax.tick_params(axis='both', labelsize=6)
            ax.set_yticks([])
            ax.set_xlabel(r'2θ (°)', fontsize=7, labelpad=2)

            ax2 = ax.twinx()
            ax2.plot(tt1, ii1, lw=0.8, color=COLOR_RIGHT, rasterized=True,
                     label=instr1)
            ax2.set_ylim(-0.05, 1.15)
            ax2.set_yticks([])

            # Instrument labels as small colored text inside the axes
            ax.text(0.03, 0.88, instr0, transform=ax.transAxes,
                    fontsize=6, color=COLOR_LEFT, va='top')
            ax2.text(0.97, 0.88, instr1, transform=ax2.transAxes,
                     fontsize=6, color=COLOR_RIGHT, va='top', ha='right')

    # Hide unused axes
    for j in range(len(records), len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle('XRD Patterns — Solid Solution Samples', fontsize=13)
    # rect leaves space for suptitle at top and prevents clipping
    fig.tight_layout(rect=[0, 0, 1, 0.97], h_pad=3.0, w_pad=1.5)
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_all():
    """Returns list of (sample_num, formula, [(path, data), ...])."""
    records = []
    for num, formula in SS_SAMPLES:
        fpaths = find_xrd_files(num)
        if not fpaths:
            print(f"  WARNING: no XRD file found for {num} ({formula})")
            continue
        datasets = []
        for fpath in fpaths:
            try:
                data = parse_xrd_file(fpath)
                datasets.append((fpath, data))
            except Exception as e:
                print(f"  WARNING: could not parse {fpath}: {e}")
        if datasets:
            records.append((num, formula, datasets))
    return records


if __name__ == '__main__':
    print("Loading Sheet1 labels...")
    sheet1_labels = load_sheet1_labels()

    print("Loading XRD data...")
    records = load_all()
    print(f"Loaded {len(records)} samples "
          f"({sum(1 for *_, ds in records if len(ds) > 1)} with dual datasets).")

    fig = plot_grid(records, sheet1_labels)
    plt.show()
