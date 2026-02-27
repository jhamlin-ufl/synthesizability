"""
SuperCon database plugin for the synthesis dashboard.

Displays known superconductors in the same ternary/binary/unary element system
as each sample, with highlights for closest composition and closest Tc.

Data source: MDR SuperCon Datasheet Ver.240322 (NIMS)
DOI: https://doi.org/10.48505/nims.3837
"""

import json
import math
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition


SUPERCON_DIR = Path("results/supercon")


def _load_hits(formula: str) -> list[dict]:
    """Load cached SuperCon hits for a sample formula."""
    path = SUPERCON_DIR / f"{formula}.json"
    if not path.exists():
        return []
    return json.loads(path.read_text())


def _n_elements(formula: str) -> int:
    """Return number of distinct elements in a formula string."""
    try:
        return len(Composition(formula).elements)
    except Exception:
        return 0


def _composition_distance(formula_a: str, formula_b: str) -> float:
    """
    Compute Euclidean distance between two compositions in fractional coordinate space.
    Only considers elements present in formula_a (the sample).
    Returns infinity if either formula cannot be parsed.
    """
    try:
        a = Composition(formula_a).fractional_composition
        b = Composition(formula_b).fractional_composition
        elements = list(a.elements)
        dist = sum((a[el] - b.get(el, 0)) ** 2 for el in elements)
        dist += sum(b[el] ** 2 for el in b.elements if el not in a)
        return math.sqrt(dist)
    except Exception:
        return float("inf")


def get_summary_cards(df) -> list[dict]:
    """Return summary card: number of samples with at least one SuperCon hit."""
    count = 0
    for formula in df["formula"]:
        if _load_hits(formula):
            count += 1
    return [{"label": "Samples w/ SuperCon Hits", "value": str(count)}]


def get_table_columns(df) -> list[str]:
    """No new dataframe columns — data lives in per-sample JSONs."""
    return []


def generate(row, plots_dir, results_dir) -> None:
    """No plot generation needed for this plugin."""
    pass


def get_detail_section(row, plots_dir, results_dir) -> dict | None:
    """
    Render SuperCon hits table with highlight cards for closest composition
    and closest Tc.
    """
    hits = _load_hits(row["formula"])
    if not hits:
        return None

    sample_formula = row["formula"]
    measured_tc = row.get("tc_kelvin")
    has_measured_tc = pd.notna(measured_tc)

    # --- Find highlight entries ---

    # Closest composition: minimize compositional distance, prefer more elements
    def composition_sort_key(h):
        return (_composition_distance(sample_formula, h["formula"]),
                -_n_elements(h["formula"]))

    hits_with_tc = [h for h in hits if h["tc"] is not None]

    closest_comp = min(hits, key=composition_sort_key)

    closest_tc = None
    if has_measured_tc and hits_with_tc:
        closest_tc = min(hits_with_tc, key=lambda h: abs(h["tc"] - measured_tc))

    # --- Highlight cards HTML ---
    def ref_link_html(hit: dict) -> str:
        if hit.get("url"):
            label = hit["doi"] if hit.get("doi") else "Search"
            return f'<a href="{hit["url"]}" target="_blank" class="external-link">{label}</a>'
        return hit.get("journal") or "—"

    def highlight_card(label: str, hit: dict, extra: str = "") -> str:
        tc_str = f"{hit['tc']:.2f} K" if hit["tc"] is not None else "N/A"
        name = hit["compound"] or hit["formula"]
        return f"""
<div style="background:#f0f7ff; border-left:4px solid #0066cc; padding:12px;
            margin:8px 0; border-radius:4px;">
    <strong>{label}:</strong> {name}
    &nbsp;|&nbsp; <code>{hit['formula']}</code>
    &nbsp;|&nbsp; T<sub>c</sub> = {tc_str}
    &nbsp;|&nbsp; {ref_link_html(hit)}
    {f'<br><small style="color:#666;">{extra}</small>' if extra else ''}
</div>"""

    highlights_html = ""
    comp_dist = _composition_distance(sample_formula, closest_comp["formula"])
    highlights_html += highlight_card(
        "Closest Composition", closest_comp,
        extra=f"Compositional distance: {comp_dist:.3f}"
    )

    if closest_tc is not None and closest_tc is not closest_comp:
        delta = abs(closest_tc["tc"] - measured_tc)
        highlights_html += highlight_card(
            "Closest T<sub>c</sub>", closest_tc,
            extra=f"ΔT<sub>c</sub> = {delta:.2f} K from measured {measured_tc:.2f} K"
        )
    elif closest_tc is closest_comp and has_measured_tc:
        # Same entry — just note it
        highlights_html += f"""
<div style="background:#f8f9fa; padding:8px; border-radius:4px; margin:4px 0;
            font-size:0.9em; color:#555;">
    Closest composition and closest T<sub>c</sub> are the same entry.
</div>"""

    # --- Main table (sorted Tc descending) ---
    sorted_hits = sorted(hits, key=lambda h: h["tc"] if h["tc"] is not None else -999, reverse=True)

    rows_html = ""
    for h in sorted_hits:
        tc_str = f"{h['tc']:.2f}" if h["tc"] is not None else "—"
        name = h["compound"] or h["formula"]
        ref_html = ref_link_html(h)

        # Highlight closest comp and closest tc rows in the table
        style = ""
        if h is closest_comp:
            style = ' style="background:#e8f4e8;"'
        elif closest_tc is not None and h is closest_tc:
            style = ' style="background:#fff3cd;"'

        rows_html += f"""<tr{style}>
    <td>{name}</td>
    <td><code>{h['formula']}</code></td>
    <td>{tc_str}</td>
    <td>{ref_html}</td>
</tr>\n"""

    table_id = f"sc_table_{sample_formula.replace('.', '_')}"
    n_total = len(sorted_hits)
    default_show = 10

    html = f"""
{highlights_html}

<div style="margin-top:20px;">
    <p style="color:#555; font-size:0.9em; margin-bottom:8px;">
        {n_total} known superconductor(s) in the
        <strong>{sample_formula}</strong> element system.
        Sorted by T<sub>c</sub> descending.
        <span style="display:inline-block; margin-left:8px;">
            <span style="background:#e8f4e8; padding:2px 6px; border-radius:3px; font-size:0.85em;">■ closest composition</span>
            &nbsp;
            {'<span style="background:#fff3cd; padding:2px 6px; border-radius:3px; font-size:0.85em;">■ closest T<sub>c</sub></span>' if closest_tc is not None else ''}
        </span>
    </p>

    <table id="{table_id}" class="fit-table">
        <thead>
            <tr>
                <th>Compound</th>
                <th>Formula</th>
                <th>T<sub>c</sub> (K)</th>
                <th>Reference</th>
            </tr>
        </thead>
        <tbody>
            {rows_html}
        </tbody>
    </table>

    <div id="{table_id}_controls" style="margin-top:8px; font-size:0.9em; color:#555;">
        Showing <span id="{table_id}_showing">{min(default_show, n_total)}</span> of {n_total} entries.
        <button onclick="toggleSCTable('{table_id}', {n_total}, {default_show})"
                id="{table_id}_btn"
                style="margin-left:10px; padding:3px 10px; cursor:pointer;">
            Show all
        </button>
    </div>
</div>

<script>
(function() {{
    // Hide rows beyond default on load
    const table = document.getElementById('{table_id}');
    const rows = table.tBodies[0].rows;
    for (let i = {default_show}; i < rows.length; i++) {{
        rows[i].style.display = 'none';
    }}
}})();

function toggleSCTable(tableId, nTotal, defaultShow) {{
    const table = document.getElementById(tableId);
    const rows = table.tBodies[0].rows;
    const btn = document.getElementById(tableId + '_btn');
    const showing = document.getElementById(tableId + '_showing');
    const allVisible = rows[defaultShow] && rows[defaultShow].style.display !== 'none';
    for (let i = defaultShow; i < rows.length; i++) {{
        rows[i].style.display = allVisible ? 'none' : '';
    }}
    btn.textContent = allVisible ? 'Show all' : 'Show less';
    showing.textContent = allVisible ? Math.min(defaultShow, nTotal) : nTotal;
}}
</script>
"""

    return {
        "title": "Known Superconductors in Element System",
        "html": html,
    }