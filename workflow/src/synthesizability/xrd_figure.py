"""Shared constants and helpers of the XRD figures."""
from __future__ import annotations

import string
from pathlib import Path

import pandas as pd
from pymatgen.core import Structure

OQMD_DIR = Path('workflow/figure_inputs/predicted')   # predicted structures drawn in the figures
A4_CM = (21.0, 29.7)
MARGIN_CM = 1.4
CM = 1 / 2.54                 # cm -> inch
COLOR_SOLUTION = '#2E7D32'
FIG_W_CM = A4_CM[0] - 2 * MARGIN_CM
FWHM = 0.4              # deg, broadening applied to the simulated patterns
XLABEL = r'2$\theta$ (deg)'
YLABEL = r'$I$ (arb. units)'


def _mathtext_space_group(symbol):
    """'P6_3mc' -> '$P6_3mc$', 'Fm-3m' -> '$Fm\\bar{3}m$'."""
    out = symbol.replace('-3', r'\bar{3}').replace('-4', r'\bar{4}') \
                .replace('-6', r'\bar{6}').replace('-1', r'\bar{1}')
    return f'${out}$'


def find_oqmd_cif(formula, entry_id):
    """CIF for the most stable OQMD structure at this composition."""
    cif_dir = OQMD_DIR / formula
    if not cif_dir.exists():
        return None
    if pd.notna(entry_id):
        cif = cif_dir / f'{int(entry_id)}.cif'
        if cif.exists():
            return cif
    cifs = sorted(cif_dir.glob('*.cif'))
    return cifs[0] if cifs else None


def panel_letter(index):
    """1 -> 'a', 2 -> 'b', ...  Figures here never exceed 26 panels."""
    return string.ascii_lowercase[index - 1]


def space_group_symbol(cif_path):
    """Hermann-Mauguin symbol of the structure in a CIF, or None."""
    if cif_path is None:
        return None
    try:
        return Structure.from_file(cif_path).get_space_group_info()[0]
    except Exception:
        return None
