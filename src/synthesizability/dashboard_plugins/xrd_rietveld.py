# src/synthesizability/dashboard_plugins/xrd_rietveld.py
"""
Dashboard plugin: XRD Rietveld fit results.

Displays the Jade XRD fit image, a per-phase table of refined lattice
parameters and phase fractions, overall fit quality statistics, and
CIF download links for the reference structures used in the fit.
"""
import sys
from pathlib import Path

sys.path.insert(0, 'src')
from synthesizability.parsers.wpf import parse_wpf_file, format_lattice_param


def _raw_dir(sample_id: str) -> Path:
    return Path('data/raw') / sample_id


def _rel_base(sample_id: str) -> str:
    """Relative URL prefix from a sample detail page to its raw data dir."""
    return f'../../../data/raw/{sample_id}'


# -------------------------------------------------------------------------
# Plugin interface
# -------------------------------------------------------------------------

def get_summary_cards(df) -> list[dict]:
    n = sum(
        1 for row in df.itertuples()
        if list(_raw_dir(row.sample_id).glob('*_XRD_fit.JPG'))
    )
    return [{'label': 'With XRD Fits', 'value': str(n)}]


def get_table_columns(df) -> list[str]:
    return []


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    sample_id = row['sample_id']
    sample_dir = _raw_dir(sample_id)
    rel_base = _rel_base(sample_id)

    jpg_files = sorted(sample_dir.glob('*_XRD_fit.JPG'))
    wpf_files = sorted(sample_dir.glob('*.wpf.txt'))
    cif_files = sorted(sample_dir.glob('*.cif'))

    if not jpg_files and not wpf_files:
        return None

    html = ''

    # --- Fit image ---------------------------------------------------------
    for jpg in jpg_files:
        html += f'''
<div class="plot-container">
    <img src="{rel_base}/{jpg.name}" alt="XRD Rietveld fit" style="max-width:100%;">
</div>'''

    # --- Parse wpf and render tables ---------------------------------------
    wpf_data = None
    if wpf_files:
        wpf_data = parse_wpf_file(wpf_files[0])

    if wpf_data:
        phases = wpf_data['phases']
        has_fractions = any(p['wt_pct'] is not None for p in phases)

        # Overall fit quality
        r = wpf_data['r_factor']
        e = wpf_data['r_expected']
        roe = wpf_data['r_over_e']
        chi2 = wpf_data['chi2']
        quality_parts = []
        if r is not None:
            quality_parts.append(f'R = {r:.2f}%')
        if e is not None:
            quality_parts.append(f'E = {e:.2f}%')
        if roe is not None:
            quality_parts.append(f'R/E = {roe:.2f}')
        if chi2 is not None:
            quality_parts.append(f'χ² = {chi2:.4f}')
        if quality_parts:
            html += f'''
<div class="field">
    <div class="field-label">Fit Quality</div>
    <div class="field-value">{' &nbsp;|&nbsp; '.join(quality_parts)}</div>
</div>'''

        # Phase table
        html += '''
<div class="field">
    <div class="field-label">Phases</div>
    <table class="fit-table">
        <thead>
            <tr>
                <th>Phase</th>
                <th>Formula</th>
                <th>Crystal System</th>
                <th>Space Group</th>
                <th>Lattice Parameters (Å / °)</th>
                <th>Bragg-R</th>'''
        if has_fractions:
            html += '\n                <th>Wt% (σ)</th>'
        html += '''
            </tr>
        </thead>
        <tbody>'''

        for p in phases:
            lattice_str = _format_lattice_html(p)
            sg_str = f'{p["space_group"]} ({p["space_group_number"]})' if p['space_group'] else '—'
            bragg = f'{p["bragg_r"]:.2f}%' if p['bragg_r'] is not None else '—'
            wt_cell = ''
            if has_fractions:
                if p['wt_pct'] is not None:
                    wt_cell = f'\n                <td>{p["wt_pct"]:.1f} ({p["wt_pct_sigma"]:.1f})</td>'
                else:
                    wt_cell = '\n                <td>—</td>'
            html += f'''
            <tr>
                <td>{p["name"] or "—"}</td>
                <td>{p["formula"] or "—"}</td>
                <td>{p["crystal_system"] or "—"}</td>
                <td>{sg_str}</td>
                <td style="font-family:monospace;font-size:0.9em;">{lattice_str}</td>
                <td>{bragg}</td>{wt_cell}
            </tr>'''

        html += '''
        </tbody>
    </table>
</div>'''

    # --- CIF download links -----------------------------------------------
    if cif_files:
        cif_labels = _build_cif_labels(cif_files, wpf_data)
        links = ' '.join(
            f'<a href="{rel_base}/{c.name}" download class="cif-link">'
            f'{cif_labels[c.name]}</a>'
            for c in cif_files
        )
        html += f'''
<div class="field">
    <div class="field-label">Reference CIFs</div>
    <div class="field-value">{links}</div>
</div>'''

    return {'title': 'XRD Rietveld Fit', 'html': html}


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def _build_cif_labels(cif_files: list, wpf_data: dict | None) -> dict:
    """
    Return a dict mapping each CIF filename → human-readable label.

    Matches CIF files to phases using two strategies:
      1. phase['cif_note']     — bracketed filename in the Note: field
      2. phase['pdf_number']   — matches 'PDF Card - {pdf_number}.cif'

    Unmatched files fall back to their stem.
    """
    labels = {c.name: c.stem for c in cif_files}

    if not wpf_data:
        return labels

    for phase in wpf_data.get('phases', []):
        label = phase.get('formula') or phase.get('name') or ''
        if not label:
            continue

        # Match via explicit cif_note
        cif_note = phase.get('cif_note')
        if cif_note and cif_note in labels:
            labels[cif_note] = label
            continue

        # Match via PDF number → "PDF Card - XX-XXX-XXXX.cif"
        pdf_num = phase.get('pdf_number')
        if pdf_num:
            pdf_name = f'PDF Card - {pdf_num}.cif'
            if pdf_name in labels:
                labels[pdf_name] = label

    return labels


def _format_lattice_html(phase: dict) -> str:
    """Build a compact lattice parameter string for display."""
    lattice = phase.get('lattice', {})
    if not lattice:
        return '—'

    crystal_system = phase.get('crystal_system', '')
    parts = []

    # Lengths
    for key in ('a', 'b', 'c'):
        if key in lattice:
            val, sig = lattice[key]
            parts.append(f'{key} = {format_lattice_param(val, sig)}')

    # Angles — only show non-90° ones (or ones with real uncertainty)
    for key, label in (('alpha', 'α'), ('beta', 'β'), ('gamma', 'γ')):
        if key in lattice:
            val, sig = lattice[key]
            # Skip if very close to 90° with trivial sigma (fixed angle)
            if abs(val - 90.0) < 0.01 and sig < 0.01:
                continue
            # Skip 120° gamma for hexagonal (fixed)
            if key == 'gamma' and crystal_system in ('Hexagonal', 'Trigonal') and abs(val - 120.0) < 0.1:
                continue
            parts.append(f'{label} = {format_lattice_param(val, sig)}°')

    return ' &nbsp; '.join(parts) if parts else '—'
