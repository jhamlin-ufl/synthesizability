"""
Dashboard plugin for OQMD ternary phase diagram data.

Displays all phases in the unary/binary/ternary chemical spaces for each
sample, with CIF download links. Data lives in per-space JSONs under
data/external/oqmd_ternary_phases/<space>.json.
"""
import json
from itertools import combinations
from pathlib import Path

import pandas as pd

from synthesizability.oqmd import parse_elements_from_formula


TERNARY_DIR = Path("data/external/oqmd_ternary_phases")
DEFAULT_SHOW = 10


def _load_all_entries(formula: str) -> list[dict]:
    """
    Load and merge all entries across all subspaces for a given formula.

    Each entry is augmented with an 'order' field (1=unary, 2=binary, 3=ternary)
    and a 'space' field. Sorted by stability ascending (most stable first),
    with None stability last.
    """
    elements = parse_elements_from_formula(formula)
    all_entries = []

    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = '-'.join(sorted(combo))
            json_path = TERNARY_DIR / f"{space}.json"
            if not json_path.exists():
                continue
            payload = json.loads(json_path.read_text())
            for entry in payload["entries"]:
                all_entries.append({**entry, "space": space, "order": r})

    # Sort by stability: non-null ascending, then nulls
    all_entries.sort(key=lambda e: (
        e["stability"] is None,
        e["stability"] if e["stability"] is not None else 0
    ))
    return all_entries


def _cif_rel_path(space: str, composition_id: str, entry_id: int,
                  stability: float | None) -> str | None:
    """
    Return relative path from a detail page to the CIF file, or None if missing.

    Detail pages live at results/dashboard/samples/<id>.html so the relative
    path to data/external/... is ../../../data/external/...
    """
    from synthesizability.oqmd import make_cif_filename
    filename = make_cif_filename(composition_id, entry_id, stability)
    cif_path = TERNARY_DIR / space / "cifs" / filename
    if not cif_path.exists():
        return None
    return f"../../../data/external/oqmd_ternary_phases/{space}/cifs/{filename}"


def get_summary_cards(df) -> list[dict]:
    """Return total CIF count across all spaces."""
    total = sum(
        1 for p in TERNARY_DIR.rglob("*.cif")
    ) if TERNARY_DIR.exists() else 0
    return [
        {"label": "Ternary Phase CIFs", "value": str(total)},
    ]


def get_table_columns(df) -> list[str]:
    """No dataframe columns owned by this plugin."""
    return []


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    """No plot generation needed — data is in JSONs and CIFs."""
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    """
    Render table of all OQMD phases in the sample's chemical space,
    sorted by stability with expandable rows.
    """
    formula = row["formula"]
    entries = _load_all_entries(formula)

    if not entries:
        return None

    elements = parse_elements_from_formula(formula)
    space_label = "-".join(elements)
    n_total = len(entries)
    n_stable = sum(1 for e in entries if e["stability"] is not None and e["stability"] <= 0)
    n_icsd = sum(1 for e in entries if e["icsd"])

    # Build table rows
    rows_html = ""
    for e in entries:
        comp = e["composition_id"].replace(" ", "")
        entry_id = e["entry_id"]
        stability = e["stability"]
        delta_e = e["delta_e"]
        icsd = e["icsd"]
        space = e["space"]

        # Stability cell with color coding
        if stability is None:
            stab_str = "—"
            stab_style = ""
        else:
            stab_meV = round(stability * 1000)
            stab_str = f"{stab_meV:+d}"
            if stability <= 0:
                stab_style = ' style="color:#28a745; font-weight:bold;"'
            elif stability <= 0.05:
                stab_style = ' style="color:#856404;"'
            else:
                stab_style = ' style="color:#721c24;"'

        delta_e_str = f"{delta_e*1000:.1f}" if delta_e is not None else "—"
        icsd_str = "✓" if icsd else ""

        # CIF link
        cif_rel = _cif_rel_path(space, e["composition_id"], entry_id, stability)
        if cif_rel:
            cif_html = f'<a href="{cif_rel}" class="cif-link" download>CIF</a>'
        else:
            cif_html = "—"

        rows_html += f"""<tr>
    <td><code>{comp}</code></td>
    <td{stab_style}>{stab_str}</td>
    <td>{delta_e_str}</td>
    <td style="text-align:center; color:#28a745;">{icsd_str}</td>
    <td>{space}</td>
    <td><a href="https://oqmd.org/materials/entry/{entry_id}"
           class="external-link" target="_blank">{entry_id}</a></td>
    <td>{cif_html}</td>
</tr>\n"""

    table_id = f"tp_table_{formula.replace('.', '_')}"

    html = f"""
<p style="color:#555; font-size:0.9em; margin-bottom:12px;">
    {n_total} entries in the <strong>{space_label}</strong> element system
    ({n_stable} on hull, {n_icsd} ICSD-tagged).
    Stability in meV/atom. Sorted most stable first.
</p>

<table id="{table_id}" class="fit-table">
    <thead>
        <tr>
            <th>Composition</th>
            <th>Stability (meV/atom)</th>
            <th>ΔE (meV/atom)</th>
            <th>ICSD</th>
            <th>Space</th>
            <th>Entry</th>
            <th>CIF</th>
        </tr>
    </thead>
    <tbody>
        {rows_html}
    </tbody>
</table>

<div id="{table_id}_controls" style="margin-top:8px; font-size:0.9em; color:#555;">
    Showing <span id="{table_id}_showing">{min(DEFAULT_SHOW, n_total)}</span>
    of {n_total} entries.
    <button onclick="toggleTPTable('{table_id}', {n_total}, {DEFAULT_SHOW})"
            id="{table_id}_btn"
            style="margin-left:10px; padding:3px 10px; cursor:pointer;">
        Show all
    </button>
</div>

<script>
(function() {{
    const table = document.getElementById('{table_id}');
    const rows = table.tBodies[0].rows;
    for (let i = {DEFAULT_SHOW}; i < rows.length; i++) {{
        rows[i].style.display = 'none';
    }}
}})();

function toggleTPTable(tableId, nTotal, defaultShow) {{
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
        "title": "OQMD Phases in Element System",
        "html": html,
    }