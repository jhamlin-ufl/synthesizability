"""
Dashboard plugin for ternary phase diagram data.

Sources: OQMD, Materials Project, Alexandria PBE, Alexandria PBEsol.
For three-element samples, renders interactive Plotly ternary phase diagrams
with a source toggle.  For binary/unary samples only the phase table is shown.

OQMD data:       data/external/oqmd_ternary_phases/<space>.json
MP data:         data/external/mp_ternary_phases/<space>.json
Alexandria PBE:  data/external/alexandria_pbe_ternary_phases/<space>.json
Alexandria PBEsol: data/external/alexandria_pbesol_ternary_phases/<space>.json
"""
import json
import re
from itertools import combinations
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

from synthesizability.oqmd import parse_elements_from_formula, parse_formula_to_oqmd


OQMD_DIR        = Path("data/external/oqmd_ternary_phases")
MP_DIR          = Path("data/external/mp_ternary_phases")
ALEX_PBE_DIR    = Path("data/external/alexandria_pbe_ternary_phases")
ALEX_PBESOL_DIR = Path("data/external/alexandria_pbesol_ternary_phases")

# Keep old name working for backward compatibility in other code
TERNARY_DIR = OQMD_DIR

DEFAULT_SHOW = 10

# Color scale: blue = stable (low energy), white = 0, red = unstable
COLORSCALE = "RdBu_r"
# OQMD stability ranges from negative (below hull) to positive (above hull)
# MP and Alexandria: stability ≥ 0 (0 = on hull)
CMIN_OQMD,  CMAX_OQMD  = -0.5, 0.3   # eV/atom
CMIN_MP,    CMAX_MP    =  0.0, 0.5   # eV/atom
CMIN_ALEX,  CMAX_ALEX  =  0.0, 0.5   # eV/atom (same scale as MP)


# ---------------------------------------------------------------------------
# Data helpers (shared by figure and table)
# ---------------------------------------------------------------------------

def _load_all_entries(formula: str, data_dir: Path) -> list[dict]:
    """
    Load and merge all entries across all subspaces for a given formula
    from *data_dir* (either OQMD_DIR or MP_DIR).

    Entries are augmented with 'order' (1=unary, 2=binary, 3=ternary)
    and 'space'.  Sorted by stability ascending (most stable first).
    """
    elements = parse_elements_from_formula(formula)
    all_entries = []
    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = '-'.join(sorted(combo))
            json_path = data_dir / f"{space}.json"
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
        if cid not in by_comp or (e["delta_e"] is not None and
                (by_comp[cid]["delta_e"] is None or
                 e["delta_e"] < by_comp[cid]["delta_e"])):
            by_comp[cid] = e
    return list(by_comp.values())


def _oqmd_cif_rel_path(space: str, composition_id: str, entry_id: int,
                        stability: float | None) -> str | None:
    """Relative path to an OQMD CIF file from a sample detail page."""
    from synthesizability.oqmd import make_cif_filename
    filename = make_cif_filename(composition_id, entry_id, stability)
    cif_path = OQMD_DIR / space / "cifs" / filename
    if not cif_path.exists():
        return None
    return f"../../../data/external/oqmd_ternary_phases/{space}/cifs/{filename}"


def _mp_cif_rel_path(space: str, composition_id: str, mp_id: str,
                     stability: float | None) -> str | None:
    """Relative path to an MP CIF file from a sample detail page."""
    compact = composition_id.replace(" ", "")
    safe_mp_id = mp_id.replace("/", "_")
    if stability is None:
        stab_str = "stabNone"
    else:
        meV = round(stability * 1000)
        stab_str = f"stab+{meV}meV" if meV >= 0 else f"stab{meV}meV"
    filename = f"{compact}_{safe_mp_id}_{stab_str}.cif"
    cif_path = MP_DIR / space / "cifs" / filename
    if not cif_path.exists():
        return None
    return f"../../../data/external/mp_ternary_phases/{space}/cifs/{filename}"


def _alex_cif_rel_path(source: str, space: str, composition_id: str,
                       entry_id: str, stability: float | None) -> str | None:
    """Relative path to an Alexandria CIF file from a sample detail page."""
    data_dir = ALEX_PBE_DIR if source == "alex_pbe" else ALEX_PBESOL_DIR
    compact = composition_id.replace(" ", "")
    if stability is None:
        stab_str = "stabNone"
    else:
        meV = round(stability * 1000)
        stab_str = f"stab+{meV}meV" if meV >= 0 else f"stab{meV}meV"
    filename = f"{compact}_{entry_id}_{stab_str}.cif"
    cif_path = data_dir / space / "cifs" / filename
    if not cif_path.exists():
        return None
    return f"../../../{data_dir.name}/{space}/cifs/{filename}"


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


def _make_tooltip(e: dict, is_target: bool = False, source: str = "oqmd") -> str:
    stab = e.get("stability")
    if source == "oqmd":
        stab_str = f"{stab * 1000:+.0f} meV/atom" if stab is not None else "unknown"
    else:
        # MP and Alexandria: stability ≥ 0, 0 = on hull
        stab_str = f"{stab * 1000:.0f} meV/atom above hull" if stab is not None else "unknown"
    icsd_str = "  [ICSD]" if e.get("icsd") else ""
    target_str = "★ Target composition<br>" if is_target else ""
    de = e.get("delta_e")
    de_str = f"{de * 1000:.0f} meV/atom" if de is not None else "—"
    return (
        f"{target_str}"
        f"<b>{e['composition_id']}</b>{icsd_str}<br>"
        f"Formation energy ΔE = {de_str}<br>"
        f"Hull distance = {stab_str}"
    )


def _make_ternary_figure(all_entries: list[dict], formula: str,
                         elements: list[str], source: str = "oqmd") -> go.Figure:
    """
    Build an interactive ternary figure for *source* ('oqmd' or 'mp').

    Trace indices used by filter buttons:
      0 – non-ICSD, off-hull
      1 – non-ICSD, on-hull
      2 – ICSD, off-hull
      3 – ICSD, on-hull
      4 – target composition
    """
    if source == "oqmd":
        cmin, cmax = CMIN_OQMD, CMAX_OQMD
        source_label = "OQMD"
    elif source == "mp":
        cmin, cmax = CMIN_MP, CMAX_MP
        source_label = "Materials Project"
    elif source == "alex_pbe":
        cmin, cmax = CMIN_ALEX, CMAX_ALEX
        source_label = "Alexandria PBE"
    else:
        cmin, cmax = CMIN_ALEX, CMAX_ALEX
        source_label = "Alexandria PBEsol"
    no_entry_msg = f"(no {source_label} entry at this composition)"

    entries = _lowest_per_composition(all_entries)
    n_total = len(entries)
    # OQMD: on-hull when stability ≤ 0.  MP/Alexandria: on-hull when stability ≈ 0.
    hull_tol = 1e-6 if source == "oqmd" else 0.001
    n_hull = sum(1 for e in entries
                 if e.get("stability") is not None and e["stability"] <= hull_tol)
    n_icsd = sum(1 for e in entries if e.get("icsd"))

    target_fracs = _target_fracs_from_formula(formula, elements)

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
            on = e.get("stability") is not None and e["stability"] <= hull_tol
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
            de.append(e.get("delta_e"))
            tips.append(_make_tooltip(e, is_target=is_target_group, source=source))
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
                            colorscale=COLORSCALE, cmin=cmin, cmax=cmax,
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
                color=de, colorscale=COLORSCALE, cmin=cmin, cmax=cmax,
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

    # Attach colorbar
    cb_title = ("Formation Energy<br>ΔE (eV/atom)" if source == "oqmd"
                else "Hull Distance<br>(eV/atom)")
    for tr in fig.data:
        if tr.marker.color is not None and len(tr.marker.color) > 0:
            tr.marker.colorbar = dict(
                title=dict(text=cb_title, side="right"),
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
                color=de, colorscale=COLORSCALE, cmin=cmin, cmax=cmax,
                line=dict(color="black", width=1.5),
            ),
            text=[_make_tooltip(target_entry, is_target=True, source=source)],
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
            text=[f"★ Target: {formula}<br>{no_entry_msg}"],
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
                f"<b>{system_name} — {source_label} phases</b>"
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
# Plugin interface helpers
# ---------------------------------------------------------------------------

def _build_phase_table_html(entries: list[dict], source: str, formula: str,
                             table_id: str) -> str:
    """Build the phase table HTML (without toggleTPTable definition) for one source."""
    elements = parse_elements_from_formula(formula)
    space_label = "-".join(elements)
    n_total = len(entries)
    hull_tol = 1e-6 if source == "oqmd" else 0.001
    n_stable = sum(1 for e in entries
                   if e.get("stability") is not None and e["stability"] <= hull_tol)
    n_icsd = sum(1 for e in entries if e.get("icsd"))
    if source == "oqmd":
        source_label = "OQMD"
        stab_col  = "Stability (meV/atom)"
        entry_col = "Entry"
        note      = "OQMD: negative = below hull (stable)."
    elif source == "mp":
        source_label = "Materials Project"
        stab_col  = "Hull Dist. (meV/atom)"
        entry_col = "MP ID"
        note      = "MP: 0 = on hull, positive = above."
    elif source == "alex_pbe":
        source_label = "Alexandria PBE"
        stab_col  = "Hull Dist. (meV/atom)"
        entry_col = "Entry ID"
        note      = "Alexandria: 0 = on hull, positive = above."
    else:
        source_label = "Alexandria PBEsol"
        stab_col  = "Hull Dist. (meV/atom)"
        entry_col = "Entry ID"
        note      = "Alexandria: 0 = on hull, positive = above."

    rows_html = ""
    for e in entries:
        comp = e["composition_id"].replace(" ", "")
        stability = e.get("stability")
        delta_e = e.get("delta_e")
        space = e["space"]

        if stability is None:
            stab_str, stab_style = "—", ""
        else:
            stab_meV = round(stability * 1000)
            stab_str = f"{stab_meV:+d}"
            if stability <= hull_tol:
                stab_style = ' style="color:#28a745; font-weight:bold;"'
            elif stability <= 0.05:
                stab_style = ' style="color:#856404;"'
            else:
                stab_style = ' style="color:#721c24;"'

        delta_e_str = f"{delta_e * 1000:.1f}" if delta_e is not None else "—"
        icsd_str = "✓" if e.get("icsd") else ""

        if source == "oqmd":
            entry_id = e.get("entry_id")
            cif_rel = _oqmd_cif_rel_path(space, e["composition_id"], entry_id, stability)
            entry_cell = (
                f'<a href="https://oqmd.org/materials/entry/{entry_id}" '
                f'class="external-link" target="_blank">{entry_id}</a>'
                if entry_id else "—"
            )
        elif source == "mp":
            mp_id = e.get("mp_id")
            cif_rel = _mp_cif_rel_path(space, e["composition_id"], mp_id, stability)
            entry_cell = (
                f'<a href="https://next-gen.materialsproject.org/materials/{mp_id}" '
                f'class="external-link" target="_blank">{mp_id}</a>'
                if mp_id else "—"
            )
        else:
            alex_id = e.get("entry_id", "")
            cif_rel = _alex_cif_rel_path(source, space, e["composition_id"],
                                          alex_id, stability)
            entry_cell = f"<code>{alex_id}</code>" if alex_id else "—"

        cif_html = (f'<a href="{cif_rel}" class="cif-link" download>CIF</a>'
                    if cif_rel else "—")

        rows_html += f"""<tr>
    <td><code>{comp}</code></td>
    <td{stab_style}>{stab_str}</td>
    <td>{delta_e_str}</td>
    <td style="text-align:center; color:#28a745;">{icsd_str}</td>
    <td>{space}</td>
    <td>{entry_cell}</td>
    <td>{cif_html}</td>
</tr>\n"""

    return f"""
<p style="color:#555; font-size:0.9em; margin:12px 0 8px 0;">
    {n_total} {source_label} entries in the <strong>{space_label}</strong> element system
    ({n_stable} on hull, {n_icsd} ICSD-tagged). {note} Sorted most stable first.
</p>
<table id="{table_id}" class="fit-table">
    <thead>
        <tr>
            <th>Composition</th>
            <th>{stab_col}</th>
            <th>ΔE (meV/atom)</th>
            <th>ICSD</th>
            <th>Space</th>
            <th>{entry_col}</th>
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
    var t = document.getElementById('{table_id}');
    if (!t) return;
    var rows = t.tBodies[0].rows;
    for (var i = {DEFAULT_SHOW}; i < rows.length; i++) rows[i].style.display = 'none';
}})();
</script>
"""


# ---------------------------------------------------------------------------
# Plugin interface
# ---------------------------------------------------------------------------

def get_summary_cards(df) -> list[dict]:
    n_oqmd        = sum(1 for _ in OQMD_DIR.rglob("*.cif"))        if OQMD_DIR.exists()        else 0
    n_mp          = sum(1 for _ in MP_DIR.rglob("*.cif"))           if MP_DIR.exists()           else 0
    n_alex_pbe    = sum(1 for _ in ALEX_PBE_DIR.rglob("*.cif"))    if ALEX_PBE_DIR.exists()    else 0
    n_alex_pbesol = sum(1 for _ in ALEX_PBESOL_DIR.rglob("*.cif")) if ALEX_PBESOL_DIR.exists() else 0
    cards = [{"label": "OQMD Phase CIFs", "value": str(n_oqmd)}]
    if n_mp          > 0: cards.append({"label": "MP Phase CIFs",          "value": str(n_mp)})
    if n_alex_pbe    > 0: cards.append({"label": "Alexandria PBE CIFs",    "value": str(n_alex_pbe)})
    if n_alex_pbesol > 0: cards.append({"label": "Alexandria PBEsol CIFs", "value": str(n_alex_pbesol)})
    return cards


def get_table_columns(df) -> list[str]:
    return []


def generate(row, plots_dir: Path, results_dir: Path) -> None:
    pass


def get_detail_section(row, plots_dir: Path, results_dir: Path) -> dict | None:
    formula = row["formula"]
    safe_id = re.sub(r"[^A-Za-z0-9]", "_", formula)
    elements = parse_elements_from_formula(formula)

    entries_oqmd      = _load_all_entries(formula, OQMD_DIR)
    entries_mp        = _load_all_entries(formula, MP_DIR)
    entries_alex_pbe  = _load_all_entries(formula, ALEX_PBE_DIR)
    entries_alex_pbes = _load_all_entries(formula, ALEX_PBESOL_DIR)

    if not any([entries_oqmd, entries_mp, entries_alex_pbe, entries_alex_pbes]):
        return None

    # Determine which sources have data
    sources_available = [
        ("oqmd",      entries_oqmd,      "OQMD"),
        ("mp",        entries_mp,        "Materials Project"),
        ("alex_pbe",  entries_alex_pbe,  "Alexandria PBE"),
        ("alex_pbesol", entries_alex_pbes, "Alexandria PBEsol"),
    ]
    active_sources = [(k, e, lbl) for k, e, lbl in sources_available if e]
    use_toggle = len(active_sources) >= 2
    first_source = active_sources[0][0] if active_sources else "oqmd"

    html = ""

    # --- Source toggle bar (when ≥2 sources have data) --------------------
    if use_toggle:
        n = len(active_sources)
        btn_html = ""
        for i, (src_key, _, lbl) in enumerate(active_sources):
            is_first_btn = (i == 0)
            if n == 1:
                radius = "4px"
            elif i == 0:
                radius = "4px 0 0 4px"
            elif i == n - 1:
                radius = "0 4px 4px 0"
            else:
                radius = "0"
            if is_first_btn:
                style = (f"padding:6px 16px; cursor:pointer; border:1px solid #2c7be5; "
                         f"background:#2c7be5; color:white; font-weight:bold; "
                         f"border-radius:{radius};")
            else:
                style = (f"padding:6px 16px; cursor:pointer; border:1px solid #aaa; "
                         f"background:#f8f9fa; color:#333; font-weight:normal; "
                         f"border-radius:{radius};")
            # Remove gap between adjacent buttons
            margin = "" if i == 0 else "margin-left:-1px;"
            btn_html += (
                f'<button id="tp_src_btn_{src_key}_{safe_id}" '
                f'onclick="toggleTPSource(\'{safe_id}\', \'{src_key}\')" '
                f'style="{style}{margin}">{lbl}</button>'
            )
        html += f'<div style="margin-bottom:14px;">{btn_html}</div>\n'

    # --- Build per-source HTML blocks -------------------------------------
    source_html_map = {}
    for src_key, entries, _ in active_sources:
        src_html = ""
        if len(elements) == 3:
            fig = _make_ternary_figure(entries, formula, elements, source=src_key)
            src_html += fig.to_html(
                full_html=False,
                include_plotlyjs="cdn",
                config={"displayModeBar": False},
            )
        src_html += _build_phase_table_html(
            entries, src_key, formula, f"tp_table_{src_key}_{safe_id}"
        )
        source_html_map[src_key] = src_html

    # --- Assemble with toggle wrappers ------------------------------------
    if use_toggle:
        all_keys = [k for k, _, _ in active_sources]
        for src_key, src_html in source_html_map.items():
            display = "" if src_key == first_source else ' style="display:none"'
            html += f'<div id="tp_src_{src_key}_{safe_id}"{display}>{src_html}</div>\n'

        all_keys_js = json.dumps(all_keys)
        html += f"""
<script>
function toggleTPSource(safeId, source) {{
    {all_keys_js}.forEach(function(s) {{
        var active = (s === source);
        var div = document.getElementById('tp_src_' + s + '_' + safeId);
        var btn = document.getElementById('tp_src_btn_' + s + '_' + safeId);
        if (div) div.style.display = active ? '' : 'none';
        if (btn) {{
            if (active) {{
                btn.style.background = '#2c7be5'; btn.style.color = 'white';
                btn.style.fontWeight = 'bold'; btn.style.borderColor = '#2c7be5';
            }} else {{
                btn.style.background = '#f8f9fa'; btn.style.color = '#333';
                btn.style.fontWeight = 'normal'; btn.style.borderColor = '#aaa';
            }}
        }}
    }});
}}
function toggleTPTable(tableId, nTotal, defaultShow) {{
    var table = document.getElementById(tableId);
    var rows = table.tBodies[0].rows;
    var btn = document.getElementById(tableId + '_btn');
    var showing = document.getElementById(tableId + '_showing');
    var allVisible = rows[defaultShow] && rows[defaultShow].style.display !== 'none';
    for (var i = defaultShow; i < rows.length; i++) {{
        rows[i].style.display = allVisible ? 'none' : '';
    }}
    btn.textContent = allVisible ? 'Show all' : 'Show less';
    showing.textContent = allVisible ? Math.min(defaultShow, nTotal) : nTotal;
}}
</script>
"""
    else:
        # Single source — no toggle needed
        html += list(source_html_map.values())[0]
        html += """
<script>
function toggleTPTable(tableId, nTotal, defaultShow) {
    var table = document.getElementById(tableId);
    var rows = table.tBodies[0].rows;
    var btn = document.getElementById(tableId + '_btn');
    var showing = document.getElementById(tableId + '_showing');
    var allVisible = rows[defaultShow] && rows[defaultShow].style.display !== 'none';
    for (var i = defaultShow; i < rows.length; i++) {
        rows[i].style.display = allVisible ? 'none' : '';
    }
    btn.textContent = allVisible ? 'Show all' : 'Show less';
    showing.textContent = allVisible ? Math.min(defaultShow, nTotal) : nTotal;
}
</script>
"""

    return {"title": "Phases in Element System", "html": html}
