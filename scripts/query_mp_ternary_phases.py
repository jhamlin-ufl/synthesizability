# scripts/query_mp_ternary_phases.py
"""
Query the Materials Project for all phases in each unique chemical space
(unary, binary, ternary) across all samples.

Mirrors scripts/query_ternary_phases.py (OQMD) but uses the MP REST API.
Results are stored as JSON at data/external/mp_ternary_phases/<space>.json.

Schema matches the OQMD JSON schema for interoperability:
  {
    "space":    "Al-Co-Fe",
    "elements": ["Al", "Co", "Fe"],
    "entries":  [
      {
        "mp_id":          "mp-19017",
        "composition_id": "Al2 Co1",   # OQMD-style: elements sorted, space-separated
        "delta_e":        -0.424,      # formation_energy_per_atom (eV/atom)
        "stability":      0.0,         # energy_above_hull (eV/atom), 0 = on hull
        "icsd":           true         # True when theoretical=False in MP
      },
      ...
    ]
  }
"""
import json
import os
import re
import time
from itertools import combinations
from pathlib import Path
import sys

import pandas as pd
from mp_api.client import MPRester

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from synthesizability.oqmd import parse_elements_from_formula

OUTPUT_DIR = Path("data/external/mp_ternary_phases")

MP_FIELDS = [
    "material_id",
    "formula_pretty",
    "composition_reduced",
    "formation_energy_per_atom",
    "energy_above_hull",
    "is_stable",
    "theoretical",
]


def _get_api_key() -> str:
    key = os.environ.get("MP_API_KEY", "")
    if key:
        return key
    # Fall back to extracting from ~/.bashrc (for non-interactive shells)
    bashrc = Path.home() / ".bashrc"
    if bashrc.exists():
        for line in bashrc.read_text().splitlines():
            m = re.match(r"export\s+MP_API_KEY=(.+)", line.strip())
            if m:
                return m.group(1).strip().strip('"').strip("'")
    raise RuntimeError(
        "MP_API_KEY not found in environment or ~/.bashrc.\n"
        "Set it with: export MP_API_KEY=<your_key>"
    )


def _composition_to_id(composition_reduced: dict) -> str:
    """
    Convert MP reduced composition dict to OQMD-style composition_id string.
    e.g. {'Al': 2.0, 'Co': 1.0} → 'Al2 Co1'
    Elements are sorted alphabetically.
    """
    parts = []
    for el in sorted(composition_reduced.keys()):
        count = composition_reduced[el]
        # Format as integer if whole number
        if count == int(count):
            parts.append(f"{el}{int(count)}")
        else:
            parts.append(f"{el}{count}")
    return " ".join(parts)


def get_all_spaces(formulas: list[str]) -> list[list[str]]:
    """All unique chemical spaces (unary, binary, ternary) across all formulas."""
    space_set = set()
    for formula in formulas:
        elements = parse_elements_from_formula(formula)
        for r in range(1, len(elements) + 1):
            for combo in combinations(elements, r):
                space_set.add(tuple(sorted(combo)))
    return [list(s) for s in sorted(space_set)]


def query_space(mpr: MPRester, elements: list[str]) -> list[dict]:
    """
    Query MP for all phases in the exact chemical space defined by elements.
    Returns list of entry dicts in the shared JSON schema.
    """
    chemsys = "-".join(sorted(elements))
    results = mpr.materials.summary.search(
        chemsys=chemsys,
        fields=MP_FIELDS,
    )

    entries = []
    for r in results:
        try:
            comp_reduced = r.composition_reduced
            # comp_reduced is a pymatgen Composition object; get element→count dict
            comp_dict = {str(el): float(amt)
                         for el, amt in comp_reduced.items()}
            composition_id = _composition_to_id(comp_dict)

            delta_e = float(r.formation_energy_per_atom) if r.formation_energy_per_atom is not None else None
            stability = float(r.energy_above_hull) if r.energy_above_hull is not None else None
            icsd = not bool(r.theoretical)

            entries.append({
                "mp_id":          str(r.material_id),
                "composition_id": composition_id,
                "delta_e":        delta_e,
                "stability":      stability,
                "icsd":           icsd,
            })
        except Exception as e:
            print(f"    Warning: skipping entry {getattr(r, 'material_id', '?')}: {e}")

    # Sort by stability ascending (most stable = lowest energy_above_hull first)
    entries.sort(key=lambda e: (
        e["stability"] is None,
        e["stability"] if e["stability"] is not None else 0,
    ))
    return entries


def main():
    api_key = _get_api_key()
    print(f"MP API key: {api_key[:8]}...")

    print("Loading synthesis data...")
    df = pd.read_csv("data/processed/synthesis_data_no_disorder.csv")
    all_spaces = get_all_spaces(df["formula"].dropna().tolist())
    print(f"Found {len(all_spaces)} unique chemical spaces")
    by_order = {n: sum(1 for s in all_spaces if len(s) == n) for n in [1, 2, 3]}
    for n, count in by_order.items():
        print(f"  {n}-element spaces: {count}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    total_start = time.time()
    n_success = n_empty = n_skipped = 0
    n_total = len(all_spaces)

    with MPRester(api_key) as mpr:
        version = mpr.get_database_version()
        print(f"Connected to MP database version: {version}\n")

        for idx, elements in enumerate(all_spaces):
            key = "-".join(elements)
            out_path = OUTPUT_DIR / f"{key}.json"

            if out_path.exists():
                print(f"[{idx+1}/{n_total}] {key}: skipping (already exists)")
                n_skipped += 1
                continue

            print(f"[{idx+1}/{n_total}] {key}...", end=" ", flush=True)
            t0 = time.time()

            entries = query_space(mpr, elements)
            elapsed = time.time() - t0

            if entries:
                print(f"{len(entries)} entries ({elapsed:.1f}s)")
                n_success += 1
            else:
                print(f"0 entries ({elapsed:.1f}s)")
                n_empty += 1

            payload = {
                "space":    key,
                "elements": elements,
                "entries":  entries,
            }
            out_path.write_text(json.dumps(payload, indent=2))

    total_elapsed = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"Query complete in {total_elapsed:.1f}s")
    print(f"  With entries: {n_success}")
    print(f"  Empty:        {n_empty}")
    print(f"  Skipped:      {n_skipped} (already existed)")
    print(f"  Output:       {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
