"""Dashboard plugin for susceptibility analysis."""
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, "src")
from synthesizability.susceptibility import (
    load_all_chi_data,
    extract_tc_values,
    plot_single_chi_real,
    plot_single_chi_imaginary,
    plot_single_hc2,
)


def _composition_from_sample_id(sample_id: str) -> str:
    parts = sample_id.split('_')
    return parts[2] if len(parts) > 2 else sample_id


def _fit_params_for_sample(row, results_dir: Path) -> dict | None:
    """Load fit parameters for this sample from the hc2 CSV."""
    fit_params_path = results_dir / 'susceptibility' / 'hc2_fit_parameters.csv'
    if not fit_params_path.exists():
        return None
    fit_df = pd.read_csv(fit_params_path)
    composition = _composition_from_sample_id(row['sample_id'])
    match = fit_df[fit_df['Composition'] == composition]
    if len(match) == 0:
        return None
    r = match.iloc[0]
    return {
        'linear': {'Hc2_0': r['Hc2(0) Linear (T)'], 'Tc': r['Tc Linear (K)']},
        'quadratic': {'Hc2_0': r['Hc2(0) Quadratic (T)'], 'Tc': r['Tc Quadratic (K)']},
    }


def get_summary_cards(df) -> list[dict]:
    """Return summary cards for susceptibility data."""
    if 'chi_n_files' not in df.columns:
        return []
    n_with_chi = (df['chi_n_files'] > 0).sum()
    n_superconducting = df['tc_kelvin'].notna().sum() if 'tc_kelvin' in df.columns else 0
    return [
        {'label': 'With χ Data', 'value': str(n_with_chi)},
        {'label': 'Superconducting', 'value': str(n_superconducting)},
    ]


def get_table_columns(df) -> list[str]:
    """Return susceptibility-related columns present in df."""
    candidates = ['tc_kelvin', 'chi_n_files', 'chi_has_high_field', 'chi_fields']
    return [c for c in candidates if c in df.columns]


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    """Generate per-sample susceptibility plots if chi data exists."""
    if row['chi_n_files'] == 0:
        return

    sample_id = row['sample_id']
    composition = _composition_from_sample_id(sample_id)
    sample_dir = Path('data/raw') / sample_id

    chi_data = load_all_chi_data(sample_dir)
    if len(chi_data) == 0:
        return

    # Real part
    real_path = plots_dir / f'{sample_id}_chi_real.png'
    if not real_path.exists():
        fig = plot_single_chi_real(chi_data, composition, sample_id)
        fig.savefig(real_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Imaginary part
    imag_path = plots_dir / f'{sample_id}_chi_imag.png'
    if not imag_path.exists():
        fig = plot_single_chi_imaginary(chi_data, composition, sample_id)
        fig.savefig(imag_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Hc2
    tc_data = extract_tc_values(chi_data)
    if len(tc_data) > 0:
        hc2_path = plots_dir / f'{sample_id}_hc2.png'
        if not hc2_path.exists():
            fig, _ = plot_single_hc2(tc_data, composition)
            fig.savefig(hc2_path, dpi=150, bbox_inches='tight')
            plt.close(fig)


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    """Generate susceptibility detail section HTML."""
    sample_id = row['sample_id']

    chi_real = plots_dir / f'{sample_id}_chi_real.png'
    chi_imag = plots_dir / f'{sample_id}_chi_imag.png'
    hc2 = plots_dir / f'{sample_id}_hc2.png'

    has_plots = any(p.exists() for p in [chi_real, chi_imag, hc2])
    if not has_plots:
        return None

    html = ''

    if chi_real.exists():
        html += f'''
<div class="plot-container">
    <img src="../plots/{chi_real.name}" alt="Real Susceptibility">
</div>'''

    if chi_imag.exists():
        html += f'''
<div class="plot-container">
    <img src="../plots/{chi_imag.name}" alt="Imaginary Susceptibility">
</div>'''

    if hc2.exists():
        html += f'''
<div class="plot-container">
    <img src="../plots/{hc2.name}" alt="Upper Critical Field">
</div>'''

    fit = _fit_params_for_sample(row, results_dir)
    if fit:
        html += f'''
<div class="field">
    <div class="field-label">Fit Parameters</div>
    <table class="fit-table">
        <thead><tr><th>Model</th><th>Hc2(0) (T)</th><th>Tc (K)</th></tr></thead>
        <tbody>
            <tr><td>Linear</td><td>{fit["linear"]["Hc2_0"]:.3f}</td><td>{fit["linear"]["Tc"]:.3f}</td></tr>
            <tr><td>Quadratic</td><td>{fit["quadratic"]["Hc2_0"]:.3f}</td><td>{fit["quadratic"]["Tc"]:.3f}</td></tr>
        </tbody>
    </table>
</div>'''

    return {'title': 'Susceptibility Analysis', 'html': html}