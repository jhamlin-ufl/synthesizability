# src/synthesizability/dashboard_plugins/xrd_comparison.py
"""
Dashboard plugin: measured XRD pattern vs. predicted structure simulation.

For each sample, simulates an XRD pattern from the predicted crystal structure
(OQMD stable phase or Diffusion Model genai structure) and overlays it on every
measured XRD pattern. One comparison plot is generated per measured file.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, 'src')

_OQMD_DIR = Path('data/external/oqmd_structures')
_GENAI_DIR = Path('data/external/genai_structures')

_FWHM = 0.15  # degrees — approximate lab instrument broadening

_ANNEALED_OFFSET = 1.15  # vertical gap between stacked measured patterns


# ---------------------------------------------------------------------------
# CIF lookup
# ---------------------------------------------------------------------------

def _find_cif(formula: str) -> tuple[Path | None, str | None]:
    """Return (cif_path, source_label) for the best available predicted structure."""
    oqmd_dir = _OQMD_DIR / formula
    if oqmd_dir.exists():
        cifs = sorted(oqmd_dir.glob('*.cif'))
        if cifs:
            return cifs[0], 'OQMD'
    genai_dir = _GENAI_DIR / formula
    if genai_dir.exists():
        cifs = sorted(genai_dir.glob('*.cif'))
        if cifs:
            return cifs[0], 'Diffusion Model'
    return None, None


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def _simulate_pattern(cif_path: Path, two_theta_min: float, two_theta_max: float):
    """
    Simulate a powder XRD pattern from a CIF file.

    Returns (two_theta_array, intensity_array) both normalized to max=1.
    Uses CuKα radiation and Gaussian peak broadening with FWHM = _FWHM°.
    """
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    structure = Structure.from_file(str(cif_path))
    calc = XRDCalculator(wavelength='CuKa')
    pattern = calc.get_pattern(structure, two_theta_range=(two_theta_min, two_theta_max))

    if len(pattern.x) == 0:
        return None, None

    two_theta = np.arange(two_theta_min, two_theta_max, 0.02)
    sigma = _FWHM / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    intensities = np.zeros_like(two_theta)
    for pos, intensity in zip(pattern.x, pattern.y):
        intensities += intensity * np.exp(-0.5 * ((two_theta - pos) / sigma) ** 2)

    if intensities.max() == 0:
        return None, None

    return two_theta, intensities / intensities.max()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _make_comparison_plot(
    meas_two_theta: np.ndarray,
    meas_intensity: np.ndarray,
    sim_two_theta: np.ndarray,
    sim_intensity: np.ndarray,
    sim_label: str,
    measured_filename: str,
    annealed: list | None = None,
    x_range: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Return a matplotlib Figure with measured (blue) and simulated (red) patterns.

    Annealed re-measurements are stacked above the as-cast pattern with a fixed
    vertical offset rather than sharing its baseline — drawn on one baseline the
    two measured patterns obscure each other exactly where they differ.
    """
    def _norm(values: np.ndarray) -> np.ndarray:
        peak = values.max()
        return values / peak if peak > 0 else values

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(meas_two_theta, _norm(meas_intensity),
            color='#1f4e79', linewidth=0.8,
            label='Measured (as-cast)' if annealed else 'Measured')
    ax.plot(sim_two_theta, sim_intensity,
            color='#c0392b', linewidth=1.0, alpha=0.85, label=sim_label)

    top = 1.0
    for i, (tt, intensity, label) in enumerate(annealed or []):
        offset = _ANNEALED_OFFSET * (i + 1)
        ax.plot(tt, _norm(intensity) + offset,
                color='#1a7f5a', linewidth=0.8, label=label)
        top = offset + 1.0

    ax.set_xlabel('2θ (°)', fontsize=11)
    ax.set_ylabel('Intensity (normalized, offset)' if annealed else 'Intensity (normalized)',
                  fontsize=11)
    ax.set_xlim(*(x_range or (meas_two_theta.min(), meas_two_theta.max())))
    if annealed:
        ax.set_ylim(0, top + 0.30)
    else:
        ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Only the stacked layout needs a more opaque legend, and keeping the
    # as-cast-only call untouched leaves those plots byte-identical.
    if annealed:
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    else:
        ax.legend(loc='upper right', fontsize=9)
    ax.set_title(measured_filename, fontsize=9, color='#555555', pad=4)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plugin interface
# ---------------------------------------------------------------------------

def get_summary_cards(df) -> list[dict]:
    n = sum(
        1 for formula in df['formula'].dropna().unique()
        if (_OQMD_DIR / formula).exists() or (_GENAI_DIR / formula).exists()
    )
    return [{'label': 'With XRD Comparison', 'value': str(n)}]


def get_table_columns(df) -> list[str]:
    # Flags which samples carry an annealed scan overlaid on their comparison
    # plot. collect_table_columns() drops it if the column is absent.
    return ['annealed_data_available']


def _pattern_curve(pattern) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (two_theta, intensity) for a pattern, or None if unusable."""
    two_theta = pattern.get('two_theta')
    intensity = pattern.get('intensity')
    if two_theta is None or intensity is None or len(two_theta) < 10:
        return None
    return two_theta, intensity


def _plot_is_current(plot_path: Path, sample_id: str,
                     filenames: list[str], cif_path: Path) -> bool:
    """
    True when plot_path exists and postdates every input that feeds it.

    A plain exists() check would leave a cached plot in place after an annealed
    scan is added to the sample folder, so the new curve would never appear.
    """
    if not plot_path.exists():
        return False

    plot_mtime = plot_path.stat().st_mtime
    raw_dir = Path('data/raw') / sample_id
    sources = [raw_dir / name for name in filenames if name] + [cif_path]

    for source in sources:
        try:
            if source.stat().st_mtime > plot_mtime:
                return False
        except OSError:
            continue

    return True


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    patterns = row.get('xrd_patterns') if hasattr(row, 'get') else None
    if patterns is None:
        try:
            patterns = row['xrd_patterns']
        except (KeyError, TypeError):
            patterns = []

    if not patterns:
        return

    cif_path, source = _find_cif(row['formula'])
    if cif_path is None:
        return

    annealed_patterns = [p for p in patterns if p.get('annealed')]
    as_cast_patterns = [p for p in patterns if not p.get('annealed')]

    # Annealed scans ride along on the as-cast plots instead of getting their
    # own. A sample with nothing but an annealed scan still gets plotted.
    primary = as_cast_patterns or annealed_patterns
    overlay = annealed_patterns if as_cast_patterns else []

    overlay_curves = []
    for pattern in overlay:
        curve = _pattern_curve(pattern)
        if curve is not None:
            overlay_curves.append((curve[0], curve[1], 'Measured (annealed)'))

    overlay_names = [p.get('filename', '') for p in overlay]

    for i, meas in enumerate(primary):
        plot_path = plots_dir / f'{row["sample_id"]}_xrd_comparison_{i}.png'
        if _plot_is_current(plot_path, row['sample_id'],
                            [meas.get('filename', '')] + overlay_names, cif_path):
            continue

        curve = _pattern_curve(meas)
        if curve is None:
            continue
        two_theta, intensity = curve

        # Span every curve drawn, so the simulation also covers the annealed
        # scan — it often starts lower in 2θ than the as-cast pattern.
        t_min = min([float(two_theta.min())] + [float(t.min()) for t, _, _ in overlay_curves])
        t_max = max([float(two_theta.max())] + [float(t.max()) for t, _, _ in overlay_curves])

        try:
            sim_tt, sim_int = _simulate_pattern(cif_path, t_min, t_max)
        except Exception:
            continue

        if sim_tt is None:
            continue

        formula = row['formula']
        sim_label = f'{source}: {formula}'

        fig = _make_comparison_plot(
            two_theta, intensity,
            sim_tt, sim_int,
            sim_label,
            meas.get('filename', f'pattern {i}'),
            annealed=overlay_curves,
            x_range=(t_min, t_max),
        )
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close(fig)


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    sample_id = row['sample_id']

    # Collect plots in index order
    plot_files = sorted(
        plots_dir.glob(f'{sample_id}_xrd_comparison_*.png'),
        key=lambda p: int(p.stem.rsplit('_', 1)[-1]),
    )
    if not plot_files:
        return None

    cif_path, source = _find_cif(row['formula'])

    # Build metadata line
    if source == 'OQMD':
        stability = row.get('oqmd_stability') if hasattr(row, 'get') else None
        try:
            stability = row['oqmd_stability']
        except (KeyError, TypeError):
            stability = None

        if stability is not None and not (stability != stability):  # nan check
            sign = '+' if stability >= 0 else ''
            meta = (f'OQMD stable phase &nbsp;·&nbsp; '
                    f'stability = {sign}{stability:.3f} eV/atom above hull')
        else:
            meta = 'OQMD stable phase'
    else:
        meta = 'Diffusion Model predicted structure'

    try:
        patterns = row['xrd_patterns']
    except (KeyError, TypeError):
        patterns = None
    if isinstance(patterns, list) and any(p.get('annealed') for p in patterns):
        meta += ' &nbsp;·&nbsp; annealed scan overlaid (offset above as-cast)'

    html = f'<div class="field-value" style="margin-bottom:8px;">{meta}</div>\n'

    for plot_file in plot_files:
        html += f'''
<div class="plot-container">
    <img src="../plots/{plot_file.name}" alt="XRD comparison" style="max-width:100%;">
</div>'''

    return {'title': 'XRD: Measured vs. Predicted', 'html': html}
