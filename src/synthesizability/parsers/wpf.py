# src/synthesizability/parsers/wpf.py
"""
Parser for Jade WPF-Rietveld Refinement Report files (*.wpf.txt).

These are plain-text reports exported by MDI Jade after whole-pattern fitting
Rietveld refinement. Each file describes one or more phases fit to an XRD
pattern, with per-phase lattice parameters and an optional quantitative
phase analysis section.
"""
import re
from pathlib import Path


# Regex for the phase header line: "Phase #1 [Bragg-R = 7.50%]: Hf2Ta6"
_RE_PHASE_HEADER = re.compile(
    r'^Phase\s+#(\d+)\s+\[Bragg-R\s*=\s*([\d.]+)%\]:\s*(.+)$'
)

# Crystal system line: "Cubic: Im-3m (229), Z=2, cI2 [PDF#04-003-6604] [CSD#457969]"
# Also handles: "Tetragonal: I41/amd (141) ● origin at 2/m, Z=4, ..."
_RE_CRYSTAL_SYSTEM = re.compile(
    r'^\s*(Cubic|Hexagonal|Tetragonal|Orthorhombic|Monoclinic|Triclinic|Trigonal):\s*'
    r'([^\s(]+)\s+\((\d+)\)'          # space group name + number
    r'(?:.*?\[PDF#([\d-]+)\])?'        # optional PDF#
)

# Lattice parameter line: "[x] a  = 3.29236 (0.00089) <2>"
# Also matches angle lines: "[x] α  = 90.000  (2261.044) <2>"
_RE_LATTICE = re.compile(
    r'^\s*\[(.)\]\s+([abcαβγ])\s*=\s*([-\d.]+)\s+\(([-\d.]+)\)'
)

# Unit cell volume: "Unit Cell Volume = 35.688 (17)(Å³)"
_RE_VOLUME = re.compile(
    r'^\s*Unit Cell Volume\s*=\s*([\d.]+)\s+\((\d+)\)\(Å'
)

# Note field with bracketed CIF filename: "Note: [Hf1Ta3_2082282_stab+87meV.cif]"
_RE_NOTE_CIF = re.compile(r'^\s*Note:\s*\[([^\]]+\.cif)\]')

# Chemical formula line
_RE_FORMULA = re.compile(r'^\s*Chemical Formula\s*=\s*(.+)$')

# Quantitative Analysis data row
# Handles "Phase Name" and "Phase Name ● Long Name" formats
_RE_QA_ROW = re.compile(
    r'^\s*(.+?)\s+'
    r'([\d.]+)\s+\(([\d.]+)\)\s+'   # Wt% (sigma)
    r'([\d.]+)\s+\(([\d.]+)\)\s+'   # Vol% (sigma)
    r'([\d.]+)\s+\(([\d.]+)\)'      # DD% (sigma)
)

# Refinement Halted line
_RE_HALTED = re.compile(
    r'^Refinement Halted\s+\(R/E=([\d.]+)\),\s*Round=(\d+),\s*Iter=(\d+),\s*'
    r'P=(\d+),\s*N=(\d+),\s*R=([\d.]+)%\s+\(E=([\d.]+)%.*χ²=([\d.]+)\)'
)


def _format_lattice_param(value: float, sigma: float) -> str:
    """
    Format a lattice parameter using standard crystallographic notation.
    Sigma is rounded to 2 significant figures; value is formatted to the
    same number of decimal places.  E.g. 5.13734 ± 0.00554 → '5.1373(55)'.
    """
    if sigma <= 0:
        return f'{value:.4f}'

    import math
    # Position of first significant figure in sigma (negative = right of decimal)
    first_sig = -int(math.floor(math.log10(sigma)))
    # Keep 2 significant figures of sigma → n_decimals = first_sig + 1
    n_decimals = max(first_sig + 1, 0)

    value_str = f'{value:.{n_decimals}f}'
    sigma_int = round(sigma * 10 ** n_decimals)

    return f'{value_str}({sigma_int})'


def _parse_lattice_lines(lines: list[str]) -> dict:
    """
    Parse the 3 rows of lattice parameter lines (a/α, b/β, c/γ).
    Returns dict with keys 'a','b','c','alpha','beta','gamma',
    each as (value, sigma) tuple, only for refined ([x]) params.
    """
    param_map = {'a': 'a', 'b': 'b', 'c': 'c', 'α': 'alpha', 'β': 'beta', 'γ': 'gamma'}
    lattice = {}
    for line in lines:
        m = _RE_LATTICE.match(line)
        if m:
            refined, param, value, sigma = m.group(1), m.group(2), m.group(3), m.group(4)
            key = param_map.get(param)
            if key and refined == 'x':
                lattice[key] = (float(value), float(sigma))
    return lattice


def _split_into_phase_blocks(lines: list[str]) -> list[tuple[int, list[str]]]:
    """
    Split file lines into phase blocks. Returns list of (line_index, lines)
    for each phase section.
    """
    phase_starts = []
    for i, line in enumerate(lines):
        if _RE_PHASE_HEADER.match(line):
            phase_starts.append(i)

    if not phase_starts:
        return []

    blocks = []
    for idx, start in enumerate(phase_starts):
        end = phase_starts[idx + 1] if idx + 1 < len(phase_starts) else len(lines)
        blocks.append((start, lines[start:end]))

    return blocks


def _parse_phase_block(block_lines: list[str]) -> dict | None:
    """Parse a single phase block into a structured dict."""
    phase = {
        'number': None,
        'name': None,
        'formula': None,
        'crystal_system': None,
        'space_group': None,
        'space_group_number': None,
        'pdf_number': None,
        'bragg_r': None,
        'lattice': {},
        'cell_volume': None,
        'wt_pct': None,
        'wt_pct_sigma': None,
        'vol_pct': None,
        'vol_pct_sigma': None,
        'cif_note': None,
    }

    # Collect all lattice parameter lines (they come in groups of 3)
    lattice_lines = []

    for i, line in enumerate(block_lines):
        # Phase header
        m = _RE_PHASE_HEADER.match(line)
        if m:
            phase['number'] = int(m.group(1))
            phase['bragg_r'] = float(m.group(2))
            phase['name'] = m.group(3).strip()
            continue

        # Formula
        m = _RE_FORMULA.match(line)
        if m:
            phase['formula'] = m.group(1).strip()
            continue

        # Crystal system
        m = _RE_CRYSTAL_SYSTEM.match(line)
        if m:
            phase['crystal_system'] = m.group(1)
            phase['space_group'] = m.group(2)
            phase['space_group_number'] = int(m.group(3))
            if m.group(4):
                phase['pdf_number'] = m.group(4)
            continue

        # Lattice parameter lines (collect them all, parse after)
        if _RE_LATTICE.match(line):
            lattice_lines.append(line)
            continue

        # Unit cell volume
        m = _RE_VOLUME.match(line)
        if m:
            phase['cell_volume'] = float(m.group(1))
            continue

        # CIF note
        m = _RE_NOTE_CIF.match(line)
        if m:
            phase['cif_note'] = m.group(1)
            continue

    phase['lattice'] = _parse_lattice_lines(lattice_lines)

    # For cubic, b and c are constrained equal to a — drop them if present
    if phase['crystal_system'] == 'Cubic':
        phase['lattice'].pop('b', None)
        phase['lattice'].pop('c', None)
    # For hexagonal/trigonal, b = a — drop b if present
    elif phase['crystal_system'] in ('Hexagonal', 'Trigonal'):
        phase['lattice'].pop('b', None)

    return phase if phase['number'] is not None else None


def _parse_quantitative_analysis(lines: list[str], phases: list[dict]) -> None:
    """
    Find the Quantitative Analysis section and populate wt_pct/vol_pct on phases.
    Modifies phases in-place.
    """
    in_qa = False
    header_seen = False

    # Build a lookup from short phase name → phase dict
    name_to_phase = {}
    for p in phases:
        if p['name']:
            name_to_phase[p['name'].strip().lower()] = p
        if p['formula']:
            name_to_phase[p['formula'].strip().lower()] = p

    for line in lines:
        if line.strip() == 'Quantitative Analysis:':
            in_qa = True
            continue

        if not in_qa:
            continue

        # Skip blank lines and the header row
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith('Phase ID'):
            header_seen = True
            continue

        # Stop at XRF line or Density line
        if stripped.startswith('XRF(') or stripped.startswith('Density of Specimen'):
            break

        if not header_seen:
            continue

        m = _RE_QA_ROW.match(line)
        if not m:
            continue

        raw_name = m.group(1).strip()
        wt_pct = float(m.group(2))
        wt_sigma = float(m.group(3))
        vol_pct = float(m.group(4))
        vol_sigma = float(m.group(5))

        # Name may be "FormulaName" or "FormulaName ● Long Name"
        short_name = raw_name.split('●')[0].strip().lower()

        matched = name_to_phase.get(short_name)
        if matched:
            matched['wt_pct'] = wt_pct
            matched['wt_pct_sigma'] = wt_sigma
            matched['vol_pct'] = vol_pct
            matched['vol_pct_sigma'] = vol_sigma


def parse_wpf_file(filepath: Path) -> dict | None:
    """
    Parse a Jade WPF-Rietveld Refinement Report file.

    Returns a dict with keys:
        phases: list of phase dicts (see module docstring)
        r_factor, r_expected, r_over_e, chi2, n_params, n_points

    Returns None if the file cannot be parsed.
    """
    filepath = Path(filepath)
    try:
        text = filepath.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return None

    lines = text.splitlines()

    # Verify this is a Jade WPF report
    if len(lines) < 2 or 'WPF-Rietveld Refinement Report' not in lines[1]:
        return None

    # Parse phase blocks
    phase_blocks = _split_into_phase_blocks(lines)
    if not phase_blocks:
        return None

    phases = []
    for _, block_lines in phase_blocks:
        phase = _parse_phase_block(block_lines)
        if phase:
            phases.append(phase)

    if not phases:
        return None

    # Parse Quantitative Analysis section (multi-phase only)
    _parse_quantitative_analysis(lines, phases)

    # Parse Refinement Halted line (last non-empty line)
    r_factor = r_expected = r_over_e = chi2 = None
    n_params = n_points = None

    for line in reversed(lines):
        m = _RE_HALTED.match(line.strip())
        if m:
            r_over_e = float(m.group(1))
            n_params = int(m.group(4))
            n_points = int(m.group(5))
            r_factor = float(m.group(6))
            r_expected = float(m.group(7))
            chi2 = float(m.group(8))
            break

    return {
        'phases': phases,
        'r_factor': r_factor,
        'r_expected': r_expected,
        'r_over_e': r_over_e,
        'chi2': chi2,
        'n_params': n_params,
        'n_points': n_points,
    }


def format_lattice_param(value: float, sigma: float) -> str:
    """Public helper: format a (value, sigma) lattice param as '3.2924(9)'."""
    return _format_lattice_param(value, sigma)
