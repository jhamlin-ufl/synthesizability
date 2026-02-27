"""Dashboard plugin for OQMD thermodynamic data."""
import pandas as pd
from pathlib import Path


def get_summary_cards(df) -> list[dict]:
    """Return summary cards for OQMD data."""
    if 'oqmd_stability' not in df.columns:
        return []

    n_with_oqmd = df['oqmd_stability'].notna().sum()
    n_oqmd_stable = (df['oqmd_stability'] < 0.1).sum()
    avg_stability = df['oqmd_stability'].mean()

    return [
        {'label': 'With OQMD Data', 'value': str(n_with_oqmd)},
        {'label': 'OQMD Stable (<0.1 eV)', 'value': str(n_oqmd_stable)},
        {'label': 'Avg E-above-hull', 'value': f'{avg_stability:.3f} eV'},
    ]


def get_table_columns(df) -> list[str]:
    """Return OQMD-related columns present in df."""
    candidates = ['oqmd_formula', 'oqmd_delta_e', 'oqmd_stability',
                  'oqmd_entry_id', 'oqmd_n_polymorphs']
    return [c for c in candidates if c in df.columns]


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    """No dashboard-specific generation needed for OQMD - data is in the dataframe."""
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    """Generate OQMD detail section HTML."""
    if pd.isna(row.get('oqmd_stability')):
        return {
            'title': 'OQMD Thermodynamic Data',
            'html': '<div style="color: #999; font-style: italic;">No OQMD data available for this composition</div>'
        }

    stability = row['oqmd_stability']
    delta_e = row['oqmd_delta_e']
    entry_id = int(row['oqmd_entry_id'])
    formula = row['formula']

    # Stability classification
    if stability < 0.01:
        box_class = "oqmd-stable"
        status_text = "✓ On convex hull - thermodynamically stable"
    elif stability < 0.1:
        box_class = "oqmd-warning"
        status_text = "⚠ Metastable - likely synthesizable"
    else:
        box_class = "oqmd-unstable"
        status_text = "✗ Above hull - competing phases more favorable"

    # Find CIF files
    oqmd_structures_dir = results_dir.parent / 'data' / 'external' / 'oqmd_structures'
    cif_dir = oqmd_structures_dir / formula
    cif_links_html = ''
    if cif_dir.exists():
        for cif_path in sorted(cif_dir.glob('*.cif')):
            rel_path = f'../../../data/external/oqmd_structures/{formula}/{cif_path.name}'
            cif_links_html += f'<a href="{rel_path}" class="cif-link" download>Download CIF</a>\n'

    html = f'''
<div class="{box_class}">
    <strong>{status_text}</strong>
</div>
<div class="field">
    <div class="field-label">Formation Energy (ΔHf)</div>
    <div class="field-value">{delta_e:.4f} eV/atom</div>
</div>
<div class="field">
    <div class="field-label">Energy Above Hull</div>
    <div class="field-value">{stability:.4f} eV/atom</div>
</div>
<div class="field">
    <div class="field-label">OQMD Entry</div>
    <div class="field-value">
        <a href="https://oqmd.org/materials/entry/{entry_id}" class="external-link" target="_blank">
        Entry {entry_id} ↗</a>
    </div>
</div>
'''
    if cif_links_html:
        html += f'''
<div class="field">
    <div class="field-label">Structure Files</div>
    <div class="field-value">{cif_links_html}</div>
</div>
'''

    return {'title': 'OQMD Thermodynamic Data', 'html': html}