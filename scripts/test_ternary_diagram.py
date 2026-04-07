"""
Test script: improved ternary phase diagram visualization.

Always shows lowest-ΔE entry per unique composition (no overlap).
Four preset buttons cycle through common filter combinations.
ICSD entries use diamond markers; non-ICSD use circles.
On-hull entries (stability ≤ 0) have a black outline.
The target composition, if it matches an OQMD entry, reuses that entry's
tooltip rather than adding a silent star on top.

Usage:
    poetry run python scripts/test_ternary_diagram.py
    poetry run python scripts/test_ternary_diagram.py Mo Ta Ti 1 2 2
    poetry run python scripts/test_ternary_diagram.py Al Co Fe 4 3 1
    poetry run python scripts/test_ternary_diagram.py Gd Co Sn 3 1 6
"""
import json
import re
import sys
from itertools import combinations
from pathlib import Path

import plotly.graph_objects as go

TERNARY_DIR = Path("data/external/oqmd_ternary_phases")
OUT_DIR = Path("results/test_ternary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- CLI args --------------------------------------------------------
if len(sys.argv) >= 4:
    ELEMENTS = sys.argv[1:4]
    COUNTS = [float(x) for x in sys.argv[4:7]] if len(sys.argv) >= 7 else [1.0, 1.0, 1.0]
else:
    ELEMENTS = ["Mo", "Ta", "Ti"]
    COUNTS = [1.0, 2.0, 2.0]

SYSTEM_NAME = "-".join(ELEMENTS)
_total = sum(COUNTS)
TARGET_FRACS = {el: cnt / _total for el, cnt in zip(ELEMENTS, COUNTS)}
TARGET_FORMULA = "".join(
    f"{el}{int(cnt) if cnt == int(cnt) else cnt}"
    for el, cnt in zip(ELEMENTS, COUNTS)
)

# Color scale: blue = very negative (stable), white = 0, red = positive (unstable)
COLORSCALE = "RdBu_r"
CMIN, CMAX = -0.5, 0.3


# ---- data helpers ----------------------------------------------------

def load_all_entries(elements: list[str]) -> list[dict]:
    entries = []
    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = "-".join(sorted(combo))
            p = TERNARY_DIR / f"{space}.json"
            if not p.exists():
                continue
            for e in json.loads(p.read_text())["entries"]:
                entries.append({**e, "space": space, "order": r})
    return entries


def parse_composition(composition_id: str) -> dict[str, float]:
    """'Mo2 Ta1 Ti2' → {'Mo': 2.0, 'Ta': 1.0, 'Ti': 2.0}"""
    out = {}
    for token in composition_id.split():
        m = re.match(r"([A-Z][a-z]*)(\d+(?:\.\d+)?)", token)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def to_ternary(counts: dict[str, float], elements: list[str]) -> tuple[float, float, float]:
    vals = [counts.get(el, 0.0) for el in elements]
    total = sum(vals)
    if total == 0:
        return (0.0, 0.0, 0.0)
    return tuple(v / total for v in vals)


def fracs_match(counts: dict[str, float], target_fracs: dict[str, float],
                elements: list[str], tol: float = 0.005) -> bool:
    """Return True if a composition dict matches the target fractions."""
    fa, fb, fc = to_ternary(counts, elements)
    ta = target_fracs.get(elements[0], 0.0)
    tb = target_fracs.get(elements[1], 0.0)
    tc = target_fracs.get(elements[2], 0.0)
    return abs(fa - ta) < tol and abs(fb - tb) < tol and abs(fc - tc) < tol


def lowest_per_composition(entries: list[dict]) -> list[dict]:
    """Keep only the lowest-ΔE entry per unique composition_id."""
    by_comp: dict[str, dict] = {}
    for e in entries:
        cid = e["composition_id"]
        if cid not in by_comp or e["delta_e"] < by_comp[cid]["delta_e"]:
            by_comp[cid] = e
    return list(by_comp.values())


def make_tooltip(e: dict, is_target: bool = False) -> str:
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


# ---- figure builder --------------------------------------------------

def make_figure(all_entries: list[dict], elements: list[str]) -> go.Figure:
    """
    Build the ternary figure with four traces (ICSD × hull) plus target.
    Trace order (indices used by buttons):
      0 – non-ICSD, off-hull
      1 – non-ICSD, on-hull
      2 – ICSD, off-hull
      3 – ICSD, on-hull
      4 – target
    """
    entries = lowest_per_composition(all_entries)
    n_total = len(entries)
    n_hull = sum(1 for e in entries if e.get("stability") is not None and e["stability"] <= 1e-6)
    n_icsd = sum(1 for e in entries if e.get("icsd"))

    # Find if target matches any OQMD entry
    target_entry = None
    remaining = []
    for e in entries:
        counts = parse_composition(e["composition_id"])
        if target_entry is None and fracs_match(counts, TARGET_FRACS, elements):
            target_entry = e
        else:
            remaining.append(e)

    def split(subset):
        """Return (non_icsd_off, non_icsd_on, icsd_off, icsd_on) sub-lists."""
        non_icsd_off, non_icsd_on, icsd_off, icsd_on = [], [], [], []
        for e in subset:
            on = e.get("stability") is not None and e["stability"] <= 1e-6
            icsd = bool(e.get("icsd"))
            if icsd:
                (icsd_on if on else icsd_off).append(e)
            else:
                (non_icsd_on if on else non_icsd_off).append(e)
        return non_icsd_off, non_icsd_on, icsd_off, icsd_on

    ni_off, ni_on, ic_off, ic_on = split(remaining)

    def build_coords(subset, is_target_group=False):
        a, b, c, de, tips = [], [], [], [], []
        for e in subset:
            counts = parse_composition(e["composition_id"])
            fa, fb, fc = to_ternary(counts, elements)
            a.append(fa); b.append(fb); c.append(fc)
            de.append(e["delta_e"])
            tips.append(make_tooltip(e, is_target=is_target_group))
        return a, b, c, de, tips

    fig = go.Figure()

    def add_trace(subset, name, symbol, size, opacity, outline_color, is_target_group=False):
        outline_width = 1.5 if outline_color != "rgba(0,0,0,0)" else 0
        if not subset:
            # Empty trace — still register it for button indexing but hide from legend
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
        a, b, c, de, tips = build_coords(subset, is_target_group)
        fig.add_trace(go.Scatterternary(
            a=a, b=b, c=c,
            mode="markers",
            name=name,
            marker=dict(
                symbol=symbol,
                size=size,
                color=de,
                colorscale=COLORSCALE,
                cmin=CMIN, cmax=CMAX,
                opacity=opacity,
                line=dict(color=outline_color, width=outline_width),
            ),
            text=tips,
            hoverinfo="text",
            cliponaxis=False,
        ))

    # Trace 0: non-ICSD, off-hull
    add_trace(ni_off,  "Non-ICSD, off-hull", "circle",  7, 0.55, "rgba(0,0,0,0)")
    # Trace 1: non-ICSD, on-hull
    add_trace(ni_on,   "Non-ICSD, on-hull",  "circle", 10, 0.90, "black")
    # Trace 2: ICSD, off-hull
    add_trace(ic_off,  "ICSD, off-hull",     "diamond", 8, 0.55, "rgba(0,0,0,0)")
    # Trace 3: ICSD, on-hull
    add_trace(ic_on,   "ICSD, on-hull",      "diamond", 11, 0.90, "black")

    # Attach the shared colorbar to the first trace that has data.
    # Place it on the right side; legend will go below it.
    for tr in fig.data:
        if tr.marker.color is not None and len(tr.marker.color) > 0:
            tr.marker.colorbar = dict(
                title=dict(text="Formation Energy<br>ΔE (eV/atom)", side="right"),
                thickness=14,
                len=0.45,
                x=1.02, xanchor="left",
                y=0.98, yanchor="top",
                tickformat=".2f",
            )
            tr.marker.showscale = True
            break

    # Trace 4: target
    if target_entry is not None:
        a, b, c, de, _ = build_coords([target_entry], is_target_group=True)
        symbol = "star" if not target_entry.get("icsd") else "star-diamond"
        fig.add_trace(go.Scatterternary(
            a=a, b=b, c=c,
            mode="markers",
            name=f"Target: {TARGET_FORMULA}",
            marker=dict(
                symbol=symbol,
                size=18,
                color=de,
                colorscale=COLORSCALE,
                cmin=CMIN, cmax=CMAX,
                line=dict(color="black", width=1.5),
            ),
            text=[make_tooltip(target_entry, is_target=True)],
            hoverinfo="text",
            cliponaxis=False,
        ))
    else:
        ta = TARGET_FRACS.get(elements[0], 0.0)
        tb = TARGET_FRACS.get(elements[1], 0.0)
        tc = TARGET_FRACS.get(elements[2], 0.0)
        fig.add_trace(go.Scatterternary(
            a=[ta], b=[tb], c=[tc],
            mode="markers",
            name=f"Target: {TARGET_FORMULA}",
            marker=dict(
                symbol="star", size=18,
                color="gold",
                line=dict(color="black", width=1.5),
            ),
            text=[f"Target: {TARGET_FORMULA}<br>(no OQMD entry at this composition)"],
            hoverinfo="text",
            cliponaxis=False,
        ))

    # ---- preset filter buttons ----------------------------------------
    # Visibility arrays: [T0, T1, T2, T3, T4]
    #   T0 = non-ICSD off-hull   T2 = ICSD off-hull
    #   T1 = non-ICSD on-hull    T3 = ICSD on-hull    T4 = target
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

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{SYSTEM_NAME} — OQMD phases</b>"
                f"<br><sup>{n_total} unique compositions · "
                f"{n_hull} on hull · "
                f"{n_icsd} ICSD entries · "
                f"target: {TARGET_FORMULA}</sup>"
            ),
            x=0.5,
            y=0.97,
            yanchor="top",
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
        # Buttons sit below the ternary plot, well clear of the title
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
        # Legend sits below the colorbar in the right margin
        legend=dict(
            x=1.02, xanchor="left",
            y=0.48, yanchor="top",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#ccc", borderwidth=1,
            font=dict(size=11),
        ),
        width=820,
        height=680,
        margin=dict(l=60, r=175, t=80, b=80),
    )

    return fig


# ---- main ------------------------------------------------------------

def main():
    print(f"System: {SYSTEM_NAME}  Target: {TARGET_FORMULA}")
    entries = load_all_entries(ELEMENTS)
    print(f"  {len(entries)} raw entries loaded")
    low = lowest_per_composition(entries)
    print(f"  {len(low)} unique compositions (lowest ΔE each)")
    on_hull = sum(1 for e in low if e.get("stability") is not None and e["stability"] <= 1e-6)
    icsd = sum(1 for e in low if e.get("icsd"))
    print(f"  {on_hull} on hull,  {icsd} from ICSD")

    fig = make_figure(entries, ELEMENTS)
    out = OUT_DIR / f"{SYSTEM_NAME}_v2.html"
    fig.write_html(str(out), include_plotlyjs="cdn", full_html=True)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
