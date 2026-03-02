# src/synthesizability/dashboard_plugins/composition.py
"""
Dashboard plugin for synthesis composition deviation analysis.

Compares measured element masses against formula-derived expected mole fractions,
flagging samples where the composition deviates significantly from the target.
"""

import pandas as pd
from pathlib import Path


# Deviation thresholds for color coding
_THRESHOLD_WARN = 0.02   # >2%  yellow
_THRESHOLD_BAD = 0.05    # >5%  orange
_THRESHOLD_FAIL = 0.10   # >10% red


def _deviation_style(max_dev: float) -> tuple[str, str]:
    """Return (css_border_color, css_background_color) for a given max deviation."""
    if max_dev > _THRESHOLD_FAIL:
        return "#cc0000", "#fff0f0"
    elif max_dev > _THRESHOLD_BAD:
        return "#cc6600", "#fff5e6"
    elif max_dev > _THRESHOLD_WARN:
        return "#ccaa00", "#fffbe6"
    else:
        return "#007700", "#f0fff0"


def get_summary_cards(df) -> list[dict]:
    """Return count of samples with composition flags."""
    if 'composition_ok' not in df.columns:
        return []
    n_flagged = (df['composition_ok'] == False).sum()
    return [
        {'label': 'Composition Flags', 'value': str(n_flagged)},
    ]


def get_table_columns(df) -> list[str]:
    """Return composition-related columns present in df."""
    candidates = ['composition_ok', 'composition_max_deviation']
    return [c for c in candidates if c in df.columns]


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    """No plot generation needed for this plugin."""
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    """
    Render composition deviation detail section.

    Shows a status banner and a per-element table of expected vs measured
    mole fractions with deviations.
    """
    max_dev = row.get('composition_max_deviation')
    expected = row.get('composition_expected_fractions')
    measured = row.get('composition_measured_fractions')

    if pd.isna(max_dev) or expected is None or measured is None:
        return {
            'title': 'Composition Check',
            'html': '<div style="color:#999; font-style:italic;">No composition data available for this sample.</div>'
        }

    border_color, bg_color = _deviation_style(max_dev)

    # Status banner
    pct = max_dev * 100
    if max_dev <= _THRESHOLD_WARN:
        status_text = f"✓ Composition OK — max deviation {pct:.2f}%"
    elif max_dev <= _THRESHOLD_BAD:
        status_text = f"⚠ Minor deviation — max {pct:.2f}%"
    elif max_dev <= _THRESHOLD_FAIL:
        status_text = f"⚠ Significant deviation — max {pct:.2f}%"
    else:
        status_text = f"✗ Large deviation — max {pct:.2f}% — check formula label or weighing"

    banner_html = f"""
<div style="background:{bg_color}; border-left:4px solid {border_color};
            padding:12px; margin:8px 0; border-radius:4px;">
    <strong>{status_text}</strong>
</div>"""

    # Per-element table
    elements = sorted(expected.keys())
    rows_html = ""
    for el in elements:
        exp = expected[el]
        meas = measured.get(el)
        if meas is None:
            diff_str = "—"
            diff_pct_str = "—"
            row_style = ""
        else:
            diff = meas - exp
            diff_pct = diff * 100
            diff_str = f"{diff:+.4f}"
            diff_pct_str = f"{diff_pct:+.2f}%"
            el_dev = abs(diff)
            el_border, el_bg = _deviation_style(el_dev)
            row_style = f' style="background:{el_bg};"' if el_dev > _THRESHOLD_WARN else ""

        rows_html += f"""<tr{row_style}>
    <td><strong>{el}</strong></td>
    <td>{exp:.4f}</td>
    <td>{f'{meas:.4f}' if meas is not None else '—'}</td>
    <td>{diff_str}</td>
    <td>{diff_pct_str}</td>
</tr>\n"""

    table_html = f"""
<table class="fit-table" style="margin-top:12px;">
    <thead>
        <tr>
            <th>Element</th>
            <th>Expected (mol frac)</th>
            <th>Measured (mol frac)</th>
            <th>Δ (mol frac)</th>
            <th>Δ (%)</th>
        </tr>
    </thead>
    <tbody>
        {rows_html}
    </tbody>
</table>"""

    # Legend
    legend_html = """
<div style="margin-top:10px; font-size:0.85em; color:#555;">
    <span style="background:#fffbe6; padding:2px 6px; border-radius:3px;">■ &gt;2%</span>
    &nbsp;
    <span style="background:#fff5e6; padding:2px 6px; border-radius:3px;">■ &gt;5%</span>
    &nbsp;
    <span style="background:#fff0f0; padding:2px 6px; border-radius:3px;">■ &gt;10%</span>
</div>"""

    html = banner_html + table_html + legend_html

    return {'title': 'Composition Check', 'html': html}
