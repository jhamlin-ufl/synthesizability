# scripts/query_alexandria_phases.py
"""
Query the Alexandria database for all phases in each unique chemical space
(unary, binary, ternary) across all samples.  Queries both PBE and PBEsol
functionals in one run.

Uses the public OPTIMADE API — no API key required.

Results are stored as JSON at:
  data/external/alexandria_pbe_ternary_phases/<space>.json
  data/external/alexandria_pbesol_ternary_phases/<space>.json

JSON schema matches the existing OQMD/MP schema for interoperability:
  {
    "space":    "Al-Co-Fe",
    "elements": ["Al", "Co", "Fe"],
    "entries":  [
      {
        "entry_id":      "alex-12345",
        "composition_id": "Al8 Co7 Fe1",
        "delta_e":       -0.312,
        "stability":     0.0,
        "icsd":          false
      },
      ...
    ]
  }

Field mapping from OPTIMADE response attributes:
  entry_id       ← id
  composition_id ← parsed from chemical_formula_reduced, sorted alphabetically
  delta_e        ← _alexandria_formation_energy_per_atom
  stability      ← _alexandria_hull_distance  (0 = on hull, positive = above)
  icsd           ← always false (not exposed by Alexandria OPTIMADE API)
"""
import json
import re
import time
from itertools import combinations
from pathlib import Path
import sys

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from synthesizability.oqmd import parse_elements_from_formula

FUNCTIONALS = {
    "pbe":    "data/external/alexandria_pbe_ternary_phases",
    "pbesol": "data/external/alexandria_pbesol_ternary_phases",
}

BASE_URLS = {
    "pbe":    "https://alexandria.icams.rub.de/pbe/v1",
    "pbesol": "https://alexandria.icams.rub.de/pbesol/v1",
}

RESPONSE_FIELDS = (
    "id,chemical_formula_reduced,nelements,elements,"
    "_alexandria_formation_energy_per_atom,_alexandria_hull_distance"
)
PAGE_LIMIT = 500
SLEEP_BETWEEN = 1.0   # seconds between requests
MAX_RETRIES = 1


def _composition_to_id(formula_reduced: str) -> str:
    """
    Parse OPTIMADE chemical_formula_reduced string into OQMD-style composition_id.
    e.g. "Al8Co7Fe1" → "Al8 Co7 Fe1"  (elements sorted alphabetically)
    """
    tokens = re.findall(r"([A-Z][a-z]*)(\d+)", formula_reduced)
    parts = []
    for el, count in sorted(tokens, key=lambda x: x[0]):
        parts.append(f"{el}{count}")
    return " ".join(parts)


def _fetch_page(url: str, params: dict, retry: int = 0) -> dict | None:
    """GET a single OPTIMADE page; retry once on failure."""
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        if retry < MAX_RETRIES:
            time.sleep(2)
            return _fetch_page(url, params, retry + 1)
        print(f"    ERROR fetching {url}: {e}")
        return None


def query_space(base_url: str, elements: list[str]) -> list[dict]:
    """
    Query Alexandria for all phases in the chemical space defined by elements.
    Uses exact element count filter so binary/unary subspaces are queried
    separately.  Returns list of entry dicts in the shared JSON schema.
    """
    n = len(elements)
    el_list = ", ".join(f'"{e}"' for e in sorted(elements))
    filter_str = f'elements HAS ALL {el_list} AND nelements={n}'

    url = f"{base_url}/structures"
    params = {
        "filter": filter_str,
        "page_limit": PAGE_LIMIT,
        "response_fields": RESPONSE_FIELDS,
    }

    entries = []
    page = 0
    while True:
        page += 1
        data = _fetch_page(url, params)
        if data is None:
            break

        for item in data.get("data", []):
            attrs = item.get("attributes", {})
            try:
                formula_reduced = attrs.get("chemical_formula_reduced", "")
                composition_id = _composition_to_id(formula_reduced)
                delta_e   = attrs.get("_alexandria_formation_energy_per_atom")
                stability = attrs.get("_alexandria_hull_distance")

                if delta_e   is not None: delta_e   = float(delta_e)
                if stability is not None: stability = float(stability)

                entries.append({
                    "entry_id":      str(item["id"]),
                    "composition_id": composition_id,
                    "delta_e":        delta_e,
                    "stability":      stability,
                    "icsd":           False,
                })
            except Exception as e:
                print(f"    Warning: skipping entry {item.get('id', '?')}: {e}")

        links = data.get("links", {})
        meta  = data.get("meta",  {})
        if not meta.get("more_data_available", False):
            break

        next_url = links.get("next")
        if not next_url:
            break
        # next_url is fully formed — use it as-is, clear params to avoid duplication
        url    = next_url
        params = {}
        time.sleep(SLEEP_BETWEEN)

    entries.sort(key=lambda e: (
        e["stability"] is None,
        e["stability"] if e["stability"] is not None else 0,
    ))
    return entries


def get_all_spaces(formulas: list[str]) -> list[list[str]]:
    """All unique chemical spaces (unary, binary, ternary) across all formulas."""
    space_set = set()
    for formula in formulas:
        elements = parse_elements_from_formula(formula)
        for r in range(1, len(elements) + 1):
            for combo in combinations(elements, r):
                space_set.add(tuple(sorted(combo)))
    return [list(s) for s in sorted(space_set)]


def main():
    print("Loading synthesis data...")
    df = pd.read_csv("data/processed/synthesis_data_no_disorder.csv")
    all_spaces = get_all_spaces(df["formula"].dropna().tolist())
    print(f"Found {len(all_spaces)} unique chemical spaces")
    by_order = {n: sum(1 for s in all_spaces if len(s) == n) for n in [1, 2, 3]}
    for n, count in by_order.items():
        print(f"  {n}-element spaces: {count}")
    print()

    for functional, out_dir_str in FUNCTIONALS.items():
        out_dir = Path(out_dir_str)
        out_dir.mkdir(parents=True, exist_ok=True)
        base_url = BASE_URLS[functional]

        print(f"{'='*60}")
        print(f"Querying Alexandria {functional.upper()}: {base_url}")
        print(f"Output: {out_dir}")
        print(f"{'='*60}")

        n_success = n_empty = n_skipped = 0
        n_total = len(all_spaces)

        for idx, elements in enumerate(all_spaces):
            key = "-".join(elements)
            out_path = out_dir / f"{key}.json"

            if out_path.exists():
                print(f"[{idx+1}/{n_total}] {key}: skipping (already exists)")
                n_skipped += 1
                continue

            print(f"[{idx+1}/{n_total}] {key}...", end=" ", flush=True)
            t0 = time.time()

            entries = query_space(base_url, elements)
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
            time.sleep(SLEEP_BETWEEN)

        print(f"\nDone with {functional.upper()}")
        print(f"  With entries: {n_success}")
        print(f"  Empty:        {n_empty}")
        print(f"  Skipped:      {n_skipped} (already existed)")
        print()


if __name__ == "__main__":
    main()
