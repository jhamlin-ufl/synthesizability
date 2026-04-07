"""
Dashboard plugin for OQMD ternary phase diagram data.

For three-element samples, renders an interactive Plotly ternary phase
diagram followed by a sortable table of all entries with CIF download
links.  For binary/unary samples only the table is shown.

Data lives in per-space JSONs under
data/external/oqmd_ternary_phases/<space>.json.
"""
import json
import re
from itertools import combinations
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

from synthesizability.oqmd import parse_elements_from_formula, parse_formula_to_oqmd


TERNARY_DIR = Path("data/external/oqmd_ternary_phases")
DEFAULT_SHOW = 10

# Color scale: blue = very negative (stable), white = 0, red = unstable
COLORSCALE = "RdBu_r"
CMIN, CMAX = -0.5, 0.3


# ---------------------------------------------------------------------------
# Data helpers (shared by figure and table)
# ---------------------------------------------------------------------------

def _load_all_entries(formula: str) -> list[dict]:
    """
    Load and merge all entries across all subspaces for a given formula.
    Entries are augmented with 'order' (1=unary, 2=binary, 3=ternary)
    and 'space'.  Sorted by stability ascending (most stable first).
    """
    elements = parse_elements_from_formula(formula)
    all_entries = []
    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = '-'.join(sorted(combo))
            json_path = TERNARY_DIR / f"{space}.json"
            if not json_path.exists():
                continue
            for entry in json.loads(json_path.read_text())["entries"]:
                all_entries.append({**entry, "space": space, "order": r})

    all_entries.sort(key=lambda e: (
        e["stability"] is None,
        e["stability"] if e["stability"] is not None else 0,
    ))
    return all_entries


def _lowest_per_composition(entries: list[dict]) -> list[dict]:
    """Keep only the lowest-ΔE entry per unique composition_id."""
    by_comp: dict[str, dict] = {}
    for e in entries:
        cid = e["composition_id"]
        if cid not in by_comp or e["delta_e"] < by_comp[cid]["delta_e"]:
            by_comp[cid] = e
    return list(by_comp.values())


def _cif_rel_path(space: str, composition_id: str, entry_id: int,
                  stability: float | None) -> str | None:
    """
    Return relative path from a detail page to the CIF file, or None if missing.
    Detail pages live at results/dashboard/samples/<id>.html.
    """
    from synthesizability.oqmd import make_cif_filename
    filename = make_cif_filename(composition_id, entry_id, stability)
    cif_path = TERNARY_DIR / space / "cifs" / filename
    if not cif_path.exists():
        return None
    return f"../../../data/external/oqmd_ternary_phases/{space}/cifs/{filename}"


# ---------------------------------------------------------------------------
# Ternary figure helpers
# ---------------------------------------------------------------------------

def _parse_composition(composition_id: str) -> dict[str, float]:
    """'Mo2 Ta1 Ti2' → {'Mo': 2.0, 'Ta': 1.0, 'Ti': 2.0}"""
    out = {}
    for token in composition_id.split():
        m = re.match(r"([A-Z][a-z]*)(\d+(?:\.\d+)?)", token)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def _target_fracs_from_formula(formula: str, elements: list[str]) -> dict[str, float]:
    """Return element → fraction dict for the sample formula."""
    oqmd_str = parse_formula_to_oqmd(formula)   # e.g. 'Mo1 Ta2 Ti1'
    counts = _parse_composition(oqmd_str)
    total = sum(counts.values())
    if total == 0:
        return {el: 1.0 / len(elements) for el in elements}
    return {el: counts.get(el, 0.0) / total for el in elements}


def _to_ternary(counts: dict[str, float], elements: list[str]) -> tuple[float, float, float]:
    vals = [counts.get(el, 0.0) for el in elements]
    total = sum(vals)
    if total == 0:
        return (0.0, 0.0, 0.0)
    return tuple(v / total for v in vals)


def _fracs_match(counts: dict[str, float], target_fracs: dict[str, float],
                 elements: list[str], tol: float = 0.005) -> bool:
    fa, fb, fc = _to_ternary(counts, elements)
    return (abs(fa - target_fracs.get(elements[0], 0.0)) < tol and
            abs(fb - target_fracs.get(elements[1], 0.0)) < tol and
            abs(fc - target_fracs.get(elements[2], 0.0)) < tol)


def _make_tooltip(e: dict, is_target: bool = False) -> str:
    stab = e.get("stability")
    stab_str = f"{stab * 1000:+.0f} meV/atom" if stab is not None else "unknown"
    icsd_str = "  [ICSD]" if e.get("icsd") else ""
    target_str = "★ Target composition<br>" if is_target else ""
    return (
        f"{target_str}"
        f"<b>{e['composition_id']}</b>{icsd_str}<br>"
        f"Formation energy ΔE = {e['delta_e'] * 1000:.0f} meV/atom<br>"
        f"Hull distance = {stab_str}"
    )


def _make_ternary_figure(all_entries: list[dict], formula: str,
                         elements: list[str]) -> go.Figure:
    """
    Build the interactive ternary figure.

    Trace indices used by filter buttons:
      0 – non-ICSD, off-hull
      1 – non-ICSD, on-hull
      2 – ICSD, off-hull
      3 – ICSD, on-hull
      4 – target composition
    """
    entries = _lowest_per_composition(all_entries)
    n_total = len(entries)
    n_hull = sum(1 for e in entries
                 if e.get("stability") is not None and e["stability"] <= 1e-6)
    n_icsd = sum(1 for e in entries if e.get("icsd"))

    target_fracs = _target_fracs_from_formula(formula, elements)

    # Separate target entry from the rest (avoid double-plotting)
    target_entry = None
    remaining = []
    for e in entries:
        counts = _parse_composition(e["composition_id"])
        if target_entry is None and _fracs_match(counts, target_fracs, elements):
            target_entry = e
        else:
            remaining.append(e)

    def _split(subset):
        ni_off, ni_on, ic_off, ic_on = [], [], [], []
        for e in subset:
            on = e.get("stability") is not None and e["stability"] <= 1e-6
            if e.get("icsd"):
                (ic_on if on else ic_off).append(e)
            else:
                (ni_on if on else ni_off).append(e)
        return ni_off, ni_on, ic_off, ic_on

    ni_off, ni_on, ic_off, ic_on = _split(remaining)

    def _build_coords(subset, is_target_group=False):
        a, b, c, de, tips = [], [], [], [], []
        for e in subset:
            counts = _parse_composition(e["composition_id"])
            fa, fb, fc = _to_ternary(counts, elements)
            a.append(fa); b.append(fb); c.append(fc)
            de.append(e["delta_e"])
            tips.append(_make_tooltip(e, is_target=is_target_group))
        return a, b, c, de, tips

    fig = go.Figure()

    def _add_trace(subset, name, symbol, size, opacity, outline_color,
                   is_target_group=False):
        outline_width = 1.5 if outline_color != "rgba(0,0,0,0)" else 0
        if not subset:
            fig.add_trace(go.Scatterternary(
                a=[], b=[], c=[], mode="markers", name=name,
                showlegend=False,
                marker=dict(symbol=symbol, size=size, color=[],
                            colorscale=COLORSCALE, cmin=CMIN, cmax=CMAX,
                            opacity=opacity,
                            line=dict(color=outline_color, width=outline_width)),
                hoverinfo="text",
            ))
            return
        a, b, c, de, tips = _build_coords(subset, is_target_group)
        fig.add_trace(go.Scatterternary(
            a=a, b=b, c=c,
            mode="markers",
            name=name,
            marker=dict(
                symbol=symbol, size=size,
                color=de, colorscale=COLORSCALE, cmin=CMIN, cmax=CMAX,
                opacity=opacity,
                line=dict(color=outline_color, width=outline_width),
            ),
            text=tips,
            hoverinfo="text",
            cliponaxis=False,
        ))

    _add_trace(ni_off, "Non-ICSD, off-hull", "circle",   7, 0.55, "rgba(0,0,0,0)")
    _add_trace(ni_on,  "Non-ICSD, on-hull",  "circle",  10, 0.90, "black")
    _add_trace(ic_off, "ICSD, off-hull",     "diamond",  8, 0.55, "rgba(0,0,0,0)")
    _add_trace(ic_on,  "ICSD, on-hull",      "diamond", 11, 0.90, "black")

    # Attach colorbar to the first trace that has data
    for tr in fig.data:
        if tr.marker.color is not None and len(tr.marker.color) > 0:
            tr.marker.colorbar = dict(
                title=dict(text="Formation Energy<br>ΔE (eV/atom)", side="right"),
                thickness=14, len=0.45,
                x=1.02, xanchor="left",
                y=0.98, yanchor="top",
                tickformat=".2f",
            )
            tr.marker.showscale = True
            break

    # Target composition (trace 4)
    if target_entry is not None:
        a, b, c, de, _ = _build_coords([target_entry], is_target_group=True)
        symbol = "star" if not target_entry.get("icsd") else "star-diamond"
        fig.add_trace(go.Scatterternary(
            a=a, b=b, c=c,
            mode="markers",
            name=f"Target: {formula}",
            marker=dict(
                symbol=symbol, size=18,
                color=de, colorscale=COLORSCALE, cmin=CMIN, cmax=CMAX,
                line=dict(color="black", width=1.5),
            ),
            text=[_make_tooltip(target_entry, is_target=True)],
            hoverinfo="text",
            cliponaxis=False,
        ))
    else:
        ta = target_fracs.get(elements[0], 0.0)
        tb = target_fracs.get(elements[1], 0.0)
        tc = target_fracs.get(elements[2], 0.0)
        fig.add_trace(go.Scatterternary(
            a=[ta], b=[tb], c=[tc],
            mode="markers",
            name=f"Target: {formula}",
            marker=dict(
                symbol="star", size=18,
                color="gold",
                line=dict(color="black", width=1.5),
            ),
            text=[f"★ Target: {formula}<br>(no OQMD entry at this composition)"],
            hoverinfo="text",
            cliponaxis=False,
        ))

    # Filter buttons — visibility array: [T0, T1, T2, T3, T4]
    buttons = [
        dict(label="All compositions",
             method="restyle",
             args=[{"visible": [True, True, True, True, True]}]),
        dict(label="On-hull only",
             method="restyle",
             args=[{"visible": [False, True, False, True, True]}]),
        dict(label="ICSD only",
             method="restyle",
             args=[{"visible": [False, False, True, True, True]}]),
    ]

    system_name = "-".join(elements)
    fig.update_layout(
        title=dict(
            text=(
                f"<b>{system_name} — OQMD phases</b>"
                f"<br><sup>{n_total} unique compositions · "
                f"{n_hull} on hull · "
                f"{n_icsd} ICSD entries · "
                f"target: {formula}</sup>"
            ),
            x=0.5, y=0.97, yanchor="top",
        ),
        ternary=dict(
            sum=1,
            aaxis=dict(title=dict(text=f"<b>{elements[0]}</b>", font=dict(size=14)),
                       linewidth=2, ticks="outside", tickfont=dict(size=10),
                       layer="below traces"),
            baxis=dict(title=dict(text=f"<b>{elements[1]}</b>", font=dict(size=14)),
                       linewidth=2, ticks="outside", tickfont=dict(size=10),
                       layer="below traces"),
            caxis=dict(title=dict(text=f"<b>{elements[2]}</b>", font=dict(size=14)),
                       linewidth=2, ticks="outside", tickfont=dict(size=10),
                       layer="below traces"),
        ),
        updatemenus=[dict(
            type="buttons",
            direction="left",
            buttons=buttons,
            x=0.5, xanchor="center",
            y=-0.06, yanchor="top",
            bgcolor="white",
            bordercolor="#aaa",
            font=dict(size=12),
            pad=dict(r=4, t=4),
        )],
        legend=dict(
            x=1.02, xanchor="left",
            y=0.48, yanchor="top",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#ccc", borderwidth=1,
            font=dict(size=11),
        ),
        width=750,
        height=580,
        margin=dict(l=60, r=175, t=80, b=80),
    )

    return fig


# ---------------------------------------------------------------------------
# Plugin interface
# ---------------------------------------------------------------------------

def get_summary_cards(df) -> list[dict]:
    total = sum(1 for _ in TERNARY_DIR.rglob("*.cif")) if TERNARY_DIR.exists() else 0
    return [{"label": "Ternary Phase CIFs", "value": str(total)}]


def get_table_columns(df) -> list[str]:
    return []


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    formula = row["formula"]
    entries = _load_all_entries(formula)

    if not entries:
        return None

    elements = parse_elements_from_formula(formula)
    space_label = "-".join(elements)
    n_total = len(entries)
    n_stable = sum(1 for e in entries
                   if e["stability"] is not None and e["stability"] <= 0)
    n_icsd = sum(1 for e in entries if e["icsd"])

    html = ""

    # --- Ternary diagram (3-element systems only) --------------------------
    if len(elements) == 3:
        fig = _make_ternary_figure(entries, formula, elements)
        html += fig.to_html(
            full_html=False,
            include_plotlyjs="cdn",
            config={"displayModeBar": False},
        )

    # --- Phase table -------------------------------------------------------
    html += f"""
<p style="color:#555; font-size:0.9em; margin:12px 0 8px 0;">
    {n_total} entries in the <strong>{space_label}</strong> element system
    ({n_stable} on hull, {n_icsd} ICSD-tagged).
    Stability in meV/atom. Sorted most stable first.
</p>
"""

    table_id = f"tp_table_{formula.replace('.', '_')}"
    rows_html = ""
    for e in entries:
        comp = e["composition_id"].replace(" ", "")
        entry_id = e["entry_id"]
        stability = e["stability"]
        delta_e = e["delta_e"]
        space = e["space"]

        if stability is None:
            stab_str, stab_style = "—", ""
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
        icsd_str = "✓" if e["icsd"] else ""

        cif_rel = _cif_rel_path(space, e["composition_id"], entry_id, stability)
        cif_html = (f'<a href="{cif_rel}" class="cif-link" download>CIF</a>'
                    if cif_rel else "—")

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

    html += f"""<table id="{table_id}" class="fit-table">
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

    return {"title": "OQMD Phases in Element System", "html": html}
